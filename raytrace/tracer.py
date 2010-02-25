#    Copyright 2009, Teraview Ltd.
#
#    This file is part of Raytrace.
#
#    Raytrace is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.


from enthought.traits.api import HasTraits, Array, Float, Complex,\
            Property, List, Instance, Range, Any,\
            Tuple, Event, cached_property, Set, Int, Trait, Button,\
            self, Str, Bool, PythonValue, Enum
from enthought.traits.ui.api import View, Item, ListEditor, VSplit,\
            RangeEditor, ScrubberEditor, HSplit, VGroup, TextEditor,\
            TupleEditor, VGroup, HGroup, TreeEditor, TreeNode, TitleEditor,\
            ShellEditor, Controller
            
from enthought.traits.ui.file_dialog import save_file
            
from enthought.tvtk.api import tvtk
from enthought.tvtk.pyface.scene_model import SceneModel
from enthought.tvtk.pyface.scene_editor import SceneEditor
import numpy
import threading, os, itertools
import wx
from itertools import chain, izip, islice, count
from raytrace.sources import BaseRaySource
from raytrace.rays import RayCollection, collectRays
from raytrace.constraints import BaseConstraint
from raytrace.has_queue import HasQueue, on_trait_change
from raytrace.faces import Face
from raytrace.bases import Traceable, Probe, Result
from raytrace.utils import normaliseVector, transformNormals, transformPoints,\
        transformVectors, dotprod
from raytrace import ctracer

counter = count()
    
    
class RayTraceModelHandler(Controller):
    pass
    
    
class RayTraceModel(HasQueue):
    scene = Instance(SceneModel, (), {'background':(1,1,0.8)}, transient=True)
    
    optics = List(Traceable)
    sources = List(BaseRaySource)
    probes = List(Probe)
    constraints = List(BaseConstraint)
    results = List(Result)
    
    optical_path = Float(0.0, transient=True)
    
    update = Event()
    _updating = Bool(False)
    update_complete = Event()
    
    Self = self
    ShellObj = PythonValue({}, transient=True)
        
    recursion_limit = Int(10, desc="maximum number of refractions or reflections")
    
    save_btn = Button("Save scene")
    
    def _optics_changed(self, opticList):
        scene = self.scene
        #del scene.actor_list[:]    
        for o in opticList:
            scene.add_actors(o.get_actors(scene))
        
        for optic in opticList:
            optic.on_trait_change(self.trace_all, "update")
            optic.on_trait_change(self.render_vtk, "render")
        self.trace_all()
    
    def _rays_changed(self, rayList):
        scene = self.scene
        sources = [o.pipeline for o in rayList]
        mappers = [tvtk.PolyDataMapper(input=s.output) for s in sources]
        actors = [tvtk.Actor(mapper=m) for m in mappers]
        for actor in actors:
            property = actor.property
            property.color = (1,0.5,0)
        scene.add_actors(actors)
        self.trace_all()
        
    def _probes_changed(self, probeList):
        scene = self.scene
        #del scene.actor_list[:]    
        for p in probeList:
            scene.add_actors(p.get_actors(scene))
        for probe in probeList:
            probe.on_trait_change(self.update_probes, "update")
            probe.on_trait_change(self.render_vtk, "render")
        self.trace_all()
        
    def _constraints_changed(self, constraintsList):
        for constraint in constraintsList:
            constraint.on_trait_change(self.trace_all, "update")
        self.trace_all()
        
    def _results_changed(self, resultsList):
        pass #not yet sure what we need to do here
        
    def update_probes(self):
        if self.scene is not None:
            self.render_vtk()
        
    def trace_all(self):
        if not self._updating:
            self._updating = True
            self.update = True
        
    @on_trait_change("update", dispatch="queued")
    def do_update(self):
        optics = self.optics
        #print "trace", 
        counter.next()
        if optics is not None:
            for o in optics:
                o.intersections = []
            for ray_source in self.sources:
                self.trace_ray_source(ray_source, optics)
            for o in optics:
                o.update_complete()
            for r in self.results:
                r.calc_result(self)
        self.render_vtk()
        self._updating = False
        
    def trace_detail(self, async=False):
        optics = [o.clone_traits() for o in self.optics]
        for child, parent in izip(optics, self.optics):
            child.shadow_parent = parent
        sources = [s.clone_traits() for s in self.sources]
        for child, parent in izip(sources, self.sources):
            child.shadow_parent = parent
        probes = [p.clone_traits() for p in self.probes]
        for child, parent in izip(probes, self.probes):
            child.shadow_parent = parent
        if async:
            self.thd = threading.Thread(target=self.async_trace, 
                                args=(optics, sources, probes))
            self.thd.start()
        else:
            self.async_trace(optics, sources, probes)
        
    def async_trace(self, optics, sources, probes):
        """called in a thread to do background tracing"""
        for o in optics:
            o.intersections = []
        for ray_source in sources:
            self.trace_ray_source_detail(ray_source, optics)
        for probe in probes:
            probe.find_intersections(ray_source)
        for o in optics:
            o.update_complete()
            
        wx.CallAfter(self.on_trace_complete, optics, sources)
        
    def on_trace_complete(self, optics, sources):
        for s in sources:
            s.shadow_parent.copy_traits(s)
        for o in optics:
            o.shadow_parent.copy_traits(o)
        print "async trace complete"
        
    def render_vtk(self):
        if self.scene is not None:
            self.scene.render()
        
    def trace_ray_source(self, ray_source, optics):
        """trace a ray source"""
        rays = ray_source.InputRays #FIXME
        traced_rays = []
        limit = self.recursion_limit
        count = 0
        face_sets = [o.faces for o in optics]
        all_faces = list(itertools.chain(*(fs.faces for fs in face_sets)))
        for i, f in enumerate(all_faces):
            f.idx = i
            f.update()
        for fs in face_sets:
            fs.sync_transforms()
        
        while rays.n_rays>0 and count<limit:
            #print "count", count
            traced_rays.append(rays)
            rays = ctracer.trace_segment(rays, face_sets, all_faces)
            count += 1
        ray_source.TracedRays = traced_rays
        ray_source.data_source.modified()
        
    def trace_sequence(self, input_rays, faces_sequence):
        """
        Perform a sequential ray-trace.
        
        @param input_rays: a RayCollection instance
        @param optics_sequence: a list of Face instances or lists of Faces
        
        returns - the traced rays, as a list of RayCollections including
                the initial input rays
        """
        traced_rays = [rays]
        rays = input_rays
        for faces in faces_sequence:
            if isinstance(faces, Face):
                intersections = face.trace_rays(rays)
                mask = intersections['length']!=numpy.Infinity
                intersections = intersections[mask]
                points = intersections['point']
                children = face.eval_children(rays, points)
            else:
                intersections = numpy.column_stack([f.trace_rays(rays) for f in faces])
                shortest = numpy.argmin(intersections['length'], axis=1)
                ar = numpy.arange(size)
                lengths = intersections['length'][ar,shortest]
                
                #now remove infinite rays
                mask = lengths!=numpy.Infinity
                shortest = shortest[mask]
                ar = ar[mask]
                
            rays = children
            traces_rays.append(rays)
        return traced_rays
        
    def trace_segment(self, rays, optics):
        """trace a RayCollection"""
        if not optics:
            return None
        size = rays.origin.shape[0]
        #optic.trace_rays should return a 1d recarray with 
        #dtype=(('length','f8'),('face','O'),('point','f8',[3,]))
        intersections = numpy.column_stack([o.trace_rays(rays) for o in optics])
        
        shortest = numpy.argmin(intersections['length'], axis=1)
        ar = numpy.arange(size)
        
        intersections = intersections[ar,shortest] #reduce 2D to 1D
        lengths = intersections['length']
        
        #now find the infinite rays, to be filtered out later
        mask = lengths!=numpy.Infinity
        and_finite = lambda x: numpy.logical_and(mask, x)
        
        faces = intersections['face']
        
        ###why?
        if intersections.size==1:
            points = intersections['point']
        else:
            points = intersections['point']
            
        lengthsT = lengths.reshape(-1,1)
        #print "shape", lengthsT.shape, lengthsT.dtype, lengthsT[:5]
        rays.length = lengthsT
        
        face_mask = ((f, and_finite(faces==f)) for f in set(faces))
        face_mask = (a for a in face_mask if a[1].any())
        
        #tell the input rays which faces terminate each ray
        #If a ray has no intersection, the end_face is None
        faces[numpy.logical_not(mask)] = None
        rays.end_face = faces
        
        children = filter(None,[f.eval_children(rays, points, m) for f,m in face_mask])
        if len(children)==0:
            return None
        new_rays = collectRays(*children)
        new_rays.parent = rays
        return new_rays
    
    def _save_btn_changed(self):
        filename = save_file()
        if not filename: return
        fmap = {".stp": self.write_to_STEP,
                ".step": self.write_to_STEP,
                ".wrl": self.write_to_VRML,
                ".vrml": self.write_to_VRML}
        ext = os.path.splitext(filename)[-1].lower()
        try:
            fmap[ext](filename)
        except KeyError:
            self.write_to_STEP(filename)
    
    def write_to_VRML(self, fname):
        scene = self.scene
        if scene is not None:
            renwin = scene._renwin
            if filename:
                writer = tvtk.VRMLExporter(file_name=fname,
                                           render_window=renwin)
                writer.update()
                writer.write()
                
    def write_to_STEP(self, fname):
        from raytrace.step_export import export_shapes2 as export_shapes
        optics = self.optics
        sources = self.sources
        shapes_colors = filter(None, (o.make_step_shape() for o in optics))
        shapes_colors.extend(filter(None,[s.make_step_shape() for s in sources]))
        
        shapes = [s for s,c in shapes_colors]
        colors = [c for s,c in shapes_colors]
        export_shapes(shapes, fname, colorList=colors)
        
    def _sources_changed(self, source_list):
        scene = self.scene
        for source in source_list:
            for actor in source.actors:
                scene.add_actor(actor)
            source.on_trait_change(self.trace_all, "update")
            source.on_trait_change(self.render_vtk, "render")
        self.trace_all()
    
    
#use a singleton handler
controller = RayTraceModelHandler()
    
        
def on_dclick(*obj):
    print "objects", obj
    obj[0].edit_traits(kind="live", parent=controller.info.ui.control)
    
    
        
tree_editor = TreeEditor(
                nodes=[
                       TreeNode(
                        node_for=[RayTraceModel],
                        children='',
                        auto_open=True,
                        label="=My Model",
                        view = View()
                        ),
                       TreeNode(
                        node_for=[RayTraceModel],
                        children='optics',
                        auto_open=True,
                        label="=Components",
                        view = View(),
                        ),
                       TreeNode(
                        node_for=[RayTraceModel],
                        children='sources',
                        auto_open=True,
                        label="=Ray Sources",
                        view = View(),
                        ),
                        TreeNode(
                        node_for=[RayTraceModel],
                        children='probes',
                        auto_open=True,
                        label="=Probes",
                        view = View(),
                        ),
                       TreeNode(
                        node_for=[Traceable],
                        children='',
                        auto_open=False,
                        label="name",
                        ),
                       TreeNode(
                        node_for=[BaseRaySource],
                        children='',
                        auto_open=True,
                        label="name",
                        ),
                       TreeNode(
                        node_for=[Probe],
                        children='',
                        auto_open=True,
                        label="name",
                        ),
                       TreeNode(
                        node_for=[RayTraceModel],
                        children='constraints',
                        auto_open=True,
                        label="=Constraints",
                        view = View()
                        ),
                       TreeNode(
                        node_for=[RayTraceModel],
                        children='results',
                        auto_open=True,
                        label="=Results",
                        view = View()
                        ),
                       TreeNode(
                        node_for=[BaseConstraint],
                        children='',
                        auto_open=True,
                        label="name",
                        ),
                       TreeNode(
                        node_for=[Face],
                        children='',
                        auto_open=False,
                        label="name",
                        ),
                       TreeNode(
                        node_for=[Result],
                        children='',
                        auto_open=False,
                        label="name",
                        ),
                       ],
                orientation='vertical',
                hide_root=True,
                on_dclick=on_dclick,
                )
    
ray_tracer_view = View(
                   HSplit(
                    VSplit(
                       Item('scene',editor=SceneEditor(),
                            height=600),
                       #Item('optics@', editor=ListEditor(use_notebook=True),
                       #     width=200),
                       Item('ShellObj', editor=ShellEditor(share=False)),
                       show_labels=False,
                       dock="vertical"
                       ),
                       VGroup(Item('Self', 
                                id="TracerModelID",
                                editor=tree_editor, width=200),
                            Item('save_btn'),
                            show_labels=False
                            ),
                    show_labels=False,
                    dock="horizontal",
                    id="raytrace.model"
                   ),
                   resizable=True,
                   #height=500,
                   width=800,
                   id="raytrace.view",
                   handler=controller
                   )
    
RayTraceModel.class_trait_view("traits_view", ray_tracer_view)
    
        

