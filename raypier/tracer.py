#    Copyright 2009, Teraview Ltd.
#
#    This file is part of Raypier.
#
#    Raypier is free software: you can redistribute it and/or modify
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



from traits.api import HasTraits, Array, Float, Complex,\
            Property, List, Instance, Range, Any,\
            Tuple, Event, cached_property, Set, Int, Trait, Button,\
            self, Str, Bool, PythonValue, Enum, File
from traitsui.api import View, Item, ListEditor, VSplit,\
            RangeEditor, ScrubberEditor, HSplit, VGroup, TextEditor,\
            TupleEditor, VGroup, HGroup, TreeEditor, TreeNode, TitleEditor,\
            ShellEditor, Controller, Tabbed
            
from traitsui.menu import Menu, MenuBar, Action, Separator
            
from pyface.api import FileDialog, OK as FD_OK
            
from tvtk.api import tvtk
from tvtk.pyface.scene_model import SceneModel
from tvtk.pyface.scene_editor import SceneEditor
import numpy
import threading, os, itertools
import os
import traceback
import yaml
from contextlib import contextmanager
from itertools import chain, islice, count
from raypier.sources import BaseRaySource
from raypier.core.ctracer import Face, RayCollection
from raypier.core.tracer import trace_rays
from raypier.constraints import BaseConstraint
from raypier.has_queue import HasQueue, on_trait_change
from raypier.bases import Traceable, Probe, Result
from raypier.utils import normaliseVector, transformNormals, transformPoints,\
        transformVectors, dotprod
        
from raypier.core import ctracer
from raypier.qt_future_call import FutureCall

counter = count()

from raypier import mirrors, prisms, corner_cubes, ellipsoids, sources,\
    results, beamstop, lenses, beamsplitters, waveplates, apertures, fields, gausslet_sources

optics_classes = sorted(Traceable.subclasses, key=lambda c: c.__name__)

optics_menu = Menu(name = "Components...",
                   *[Action(name=cls.__name__, action="insert_"+cls.__name__)\
                    for cls in optics_classes]
                    )

source_classes = sorted(BaseRaySource.subclasses, key=lambda c:c.__name__)
sources_menu = Menu(name = "Sources...",
                   *[Action(name=cls.__name__, action="insert_"+cls.__name__)\
                    for cls in source_classes]
                    )

results_classes = sorted(Result.subclasses, key=lambda c:c.__name__)
results_menu = Menu(name = "Results...",
                   *[Action(name=cls.__name__, action="insert_"+cls.__name__)\
                    for cls in results_classes]
                    )

menubar = MenuBar(Menu(Action(name="Open...", action="open_file_action"),
                       Action(name="Save...", action="save_file_action"),
                       Action(name="Save As...", action="save_as_action"),
                       name="File..."),
                  Menu(optics_menu,
                       sources_menu,
                       results_menu,
                       name="Insert..."))
    
    
class RayTraceModelHandler(Controller):
    def save_file_action(self, info):
        model = info.ui.context['object']
        try:
            model.save_as_yaml()
        except IOError:
            self.save_as_action(info)
        
    def save_as_action(self, info):
        dlg = FileDialog(parent=info.ui.control,
                         action='save as', 
                         wildcard=FileDialog.create_wildcard("Raypier model", ["*.raymod"]))
        ret = dlg.open()
        if ret == FD_OK:
            fname = dlg.path
            model = info.ui.context['object']
            model.save_as_yaml(filename=fname)
        info.ui.title = fname
    
    def open_file_action(self, info):
        dlg = FileDialog(parent=info.ui.control,
                         action='open', 
                         wildcard=FileDialog.create_wildcard("Raypier model", ["*.raymod"]))
        ret = dlg.open()
        if ret == FD_OK:
            fname = dlg.path
            if os.path.exists(fname):
                model = info.ui.context['object']
                model.load_from_yaml(fname)
            info.ui.title = fname
        
    def insert_component(self, info, cls):
        tracer = info.ui.context['object']
        tracer.optics.append(cls())
        
    def insert_source(self, info, cls):
        tracer = info.ui.context['object']
        tracer.sources.append(cls())
        
    def insert_result(self, info, cls):
        tracer = info.ui.context['object']
        tracer.results.append(cls())
        

for cls in optics_classes:
    def make_action(_cls):
        def action(self, info):
            self.insert_component(info, _cls)
        return action
    setattr(RayTraceModelHandler, "insert_"+cls.__name__, make_action(cls))
    
for cls in source_classes:
    def make_action(_cls):
        def action(self, info):
            self.insert_source(info, _cls)
        return action
    setattr(RayTraceModelHandler, "insert_"+cls.__name__, make_action(cls))
    
for cls in results_classes:
    def make_action(_cls):
        def action(self, info):
            self.insert_result(info, _cls)
        return action
    setattr(RayTraceModelHandler, "insert_"+cls.__name__, make_action(cls))
    
    
class RayTraceModel(HasQueue):
    scene = Instance(SceneModel, (), {'background':(1,1,0.8)}, transient=True)
    
    optics = List(Traceable)
    sources = List(BaseRaySource)
    probes = List(Probe)
    constraints = List(BaseConstraint)
    results = List(Result)
    
    optical_path = Float(0.0, transient=True)
    
    all_faces = List(ctracer.Face, desc="global list of all faces, created automatically "
                     " when a tracing operation is initiated. Ray end_face_idx "
                     "can be used to index this list")
    face_sets = List(ctracer.FaceList, desc="list of FaceLists extracted from all "
                     "optics when a tracing operation is initiated")
    
    update = Event() #triggers a tracing operation
    _updating = Bool(False) #indicating that tracing is in progress
    update_complete = Event()
    
    _hold_off = Bool(False) #Blocks tracing (while model parameters are altered)
    _update_requested = Bool(False)
    
    Self = self
    ShellObj = PythonValue({}, transient=True)
        
    recursion_limit = Int(200, desc="maximum number of refractions or reflections")
    
    save_btn = Button("Save scene")
    
    filename = File()
    
    def load_from_yaml(self, filename):
        with open(filename, 'r') as fobj:
            model = yaml.load(fobj)
        print(model)
        self.optics = model['components']
        self.sources = model['sources']
        self.results = model['results']
        self.update = True
    
    def save_as_yaml(self, filename=None):
        if filename is None:
            filename = self.filename
            if self.filename is None:
                raise IOError("no preset filename")
        model = {"components":list(self.optics),
                 "sources": list(self.sources),
                 "results": list(self.results)}
        with open(filename, 'w') as fobj:
            yaml.dump(model, fobj)
        self.filename = filename
    
    @on_trait_change("optics[]")
    def on_optics_change(self, obj, name, removed, opticList):
        #print("adding", opticList, removed, name)
        scene = self.scene
        #del scene.actor_list[:]    
        for o in opticList:
            scene.add_actors(o.get_actors(scene))
        for o in removed:
            try:
                scene.remove_actors(o.get_actors(scene))
            except:
                pass
        
        for optic in opticList:
            optic.on_trait_change(self.trace_all, "update")
            optic.on_trait_change(self.render_vtk, "render")
        for optic in removed:
            optic.on_trait_change(self.trace_all, "update", remove=True)
            optic.on_trait_change(self.render_vtk, "render", remove=True)
        self.trace_all()
    
    def _rays_changed(self, rayList):
        scene = self.scene
        sources = [o.pipeline for o in rayList]
        mappers = [tvtk.PolyDataMapper(input_connection=s.output_port) for s in sources]
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
        if self._hold_off:
            self._update_requested=True
            return
        if not self._updating:
            self._updating = True
            self.update = True
            
    @contextmanager
    def hold_off(self):
        """A provides a context where the tracing operations are blocked so
        the user can edit parameters of the model without triggering multiple traces.
        A final trace occurs on exiting the context if required."""
        self._hold_off = True
        try:
            yield
        except:
            self._update_requested = False
            raise
        finally:
            self._hold_off = False
        if self._update_requested:
            self._update_requested = False
            self.trace_all()
        
    @on_trait_change("update", dispatch="queued")
    def do_update(self):
        optics = self.optics
        print( "trace") 
        next(counter)
        try:
            if optics is not None:
                self.prepare_to_trace()
                for o in optics:
                    o.intersections = []
                for ray_source in self.sources:
                    self.trace_ray_source(ray_source, optics)
                for probe in self.probes:
                    try:
                        probe.evaluate(self.sources)
                    except:
                        traceback.print_exc()
                for o in optics:
                    o.update_complete()
                for r in self.results:
                    try:
                        r.calc_result(self)
                    except:
                        traceback.print_exc()
            self.render_vtk()
            self._updating = False
        except:
            traceback.print_exc()
        
    def trace_detail(self, _async=False):
        optics = [o.clone_traits() for o in self.optics]
        for child, parent in zip(optics, self.optics):
            child.shadow_parent = parent
        sources = [s.clone_traits() for s in self.sources]
        for child, parent in zip(sources, self.sources):
            child.shadow_parent = parent
        probes = [p.clone_traits() for p in self.probes]
        for child, parent in zip(probes, self.probes):
            child.shadow_parent = parent
        if _async:
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
            
        #wx.CallAfter(self.on_trace_complete, optics, sources)
        FutureCall(0, self.on_trace_complete, optics, sources)
        
    def on_trace_complete(self, optics, sources):
        for s in sources:
            s.shadow_parent.copy_traits(s)
        for o in optics:
            o.shadow_parent.copy_traits(o)
        print("async trace complete")
        
    def render_vtk(self):
        if self.scene is not None:
            self.scene.render()
            
    def prepare_to_trace(self):
        """Called before a tracing operation is performed, to do
        all synchronisation between optics and their faces
        """
        face_sets = [o.faces for o in self.optics]
        for fs in face_sets:
            fs.sync_transforms()
            
        self.face_sets = face_sets
        
    def trace_ray_source(self, ray_source, optics):
        """trace a ray source asequentially, using the ctracer framework"""
        max_length = ray_source.max_ray_len
        rays = ray_source.input_rays #FIXME
        face_lists = self.face_sets
        rays.wavelengths = numpy.ascontiguousarray(ray_source.wavelength_list, numpy.double)
        try:
            traced_rays, all_faces = trace_rays(rays, face_lists, 
                                                recursion_limit=self.recursion_limit, 
                                                max_length=max_length)
            self.all_faces = all_faces
            ray_source.traced_rays = traced_rays
        finally:
            ray_source.data_source.modified()
        
    def trace_sequence(self, input_rays, faces_sequence):
        """
        Perform a sequential ray-trace. *** THIS IS BROKEN ***
        
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
    
    def _save_btn_changed(self):
        dlg = FileDialog(action="save as",
                         wildcard=FileDialog.create_wildcard("Scene model", 
                                                             ["*.stp",
                                                              "*.step",
                                                              "*.wrl",
                                                              "*.vrml"]))
        ret = dlg.open()
        if ret != FD_OK: return
        filename = dlg.path
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
        from raypier.step_export import export_shapes2 as export_shapes
        optics = self.optics
        sources = self.sources
        shapes_colors = [_f for _f in (o.make_step_shape() for o in optics) if _f]
        shapes_colors.extend([_f for _f in [s.make_step_shape() for s in sources] if _f])
        
        shapes = [s for s,c in shapes_colors]
        colors = [c for s,c in shapes_colors]
        export_shapes(shapes, fname, colorList=colors)
        
    def ipython_view(self, width, height, view={}):
        import ipywidgets as widgets
        from IPython.display import Image, display, clear_output
        
        renderer = tvtk.Renderer()
        for actor in self.scene.actor_list:
            renderer.add_actor(actor)
        renderer.background = (1,1,0.8)
#         
        renderer.reset_camera()
        camera = renderer.active_camera
        if "position" in view:
            camera.position = view['position']
        if "focal_point" in view:
            camera.focal_point = view['focal_point']
        if "view_up" in view:
            camera.view_up = view['view_up']
#         
        renderWindow = tvtk.RenderWindow()
        renderWindow.off_screen_rendering = True
        renderWindow.add_renderer(renderer)
        renderWindow.size = (width, height)
        renderWindow.render()
        
        windowToImageFilter = tvtk.WindowToImageFilter()
        windowToImageFilter.input = renderWindow
        windowToImageFilter.update()
#          
        filename = "/dev/shm/temp_vtk_put.png"
        writer = tvtk.PNGWriter()
        writer.file_name = filename
        writer.write_to_memory = False
        writer.input_connection = windowToImageFilter.output_port
        writer.write()
        
        view_out = {}
        
        def show():       
            clear_output(wait=True)
            renderer.reset_camera_clipping_range()
            renderer.modified()
            renderWindow.render()
            windowToImageFilter.input = renderWindow
            windowToImageFilter.modified()
            windowToImageFilter.update()
            writer.write()
            view_out.update({"position": tuple(camera.position), 
                             "view_up": tuple(camera.view_up),
                             "focal_point": tuple(camera.focal_point)})
            return display(Image(filename), grp)
        
        def r_up(arg):
            camera.orthogonalize_view_up()
            camera.elevation(10)
            return show()
        
        def r_down(arg):
            camera.orthogonalize_view_up()
            camera.elevation(-10)
            return show()
        
        def r_left(arg):
            camera.orthogonalize_view_up()
            camera.azimuth(10)
            return show()
        
        def r_right(arg):
            camera.orthogonalize_view_up()
            camera.azimuth(-10)
            return show()
        
        def pan_left(arg):
            camera.orthogonalize_view_up()
            camera.yaw(-2)
            return show()
        
        def pan_right(arg):
            camera.orthogonalize_view_up()
            camera.yaw(2)
            return show()
        
        def pan_up(arg):
            camera.orthogonalize_view_up()
            camera.pitch(-2)
            return show()
        
        def pan_down(arg):
            camera.orthogonalize_view_up()
            camera.pitch(2)
            return show()
        
        def roll_left(arg):
            camera.roll(10)
            return show()
        
        def roll_right(arg):
            camera.roll(-10)
            return show()
        
        def zoom_in(arg):
            camera.dolly(1.2)
            return show()
            
        def zoom_out(arg):
            camera.dolly(0.8)
            return show()
        
        
        b1 = widgets.Button(description = '\u2191')
        b1.on_click(r_up)
        b2 = widgets.Button(description = '\u2193')
        b2.on_click(r_down)
        b3 = widgets.Button(description = '\u2190')
        b3.on_click(r_left)
        b4 = widgets.Button(description = '\u2192')
        b4.on_click(r_right)
        b5 = widgets.Button(description = '\u21ba')
        b5.on_click(roll_left)
        b6 = widgets.Button(description = '\u21bb')
        b6.on_click(roll_right)
        
        b7 = widgets.Button(description = '+')
        b7.on_click(zoom_in)
        b8 = widgets.Button(description = '-')
        b8.on_click(zoom_out)
        b9 = widgets.Button(description = '\u2190')
        b9.on_click(pan_left)
        b10 = widgets.Button(description = '\u2192')
        b10.on_click(pan_right)
        b11 = widgets.Button(description = '\u2191')
        b11.on_click(pan_up)
        b12 = widgets.Button(description = '\u2193')
        b12.on_click(pan_down)
        
        grp1 = widgets.HBox(border_style="solid",
                            border_width=1,
                            border_radius=5,
                            padding=2,
                            margin=1)
        grp1.children=[b1,b2,b3,b4,b5,b6]
        
        grp2 = widgets.HBox(border_style="solid",
                            border_width=1,
                            border_radius=5,
                            padding=2,
                            margin=1)
        grp2.children=[b7,b8,b9,b10,b11,b12]
        
        
        grp = widgets.HBox(children=[grp1,grp2])
        
        show()
        return view_out
        
        
    def render_bitmap(self, width, height, filename=None,
                      azimuth=15.0, elevation=30.0, roll=0.0,
                      zoom=1.0, pan_h=0.0, pan_v=0.0
                      ):
        renderer = tvtk.Renderer()
        for actor in self.scene.actor_list:
            renderer.add_actor(actor)
        renderer.background = (1,1,0.8)
        
        renderer.reset_camera()
        camera = renderer.active_camera
        camera.roll(roll)
        camera.elevation(elevation)
        camera.azimuth(azimuth)
        camera.dolly(zoom)
        camera.yaw(pan_h)
        camera.pitch(pan_v)
        
        renderWindow = tvtk.RenderWindow()
        renderWindow.off_screen_rendering = True
        renderWindow.add_renderer(renderer)
        renderWindow.size = (width, height)
        renderWindow.render()
         
        windowToImageFilter = tvtk.WindowToImageFilter()
        windowToImageFilter.input = renderWindow
        windowToImageFilter.update()
         
        writer = tvtk.PNGWriter()
        if filename is not None:
            writer.file_name = filename
            writer.write_to_memory = False
        else:
            writer.write_to_memory = True
        writer.input_connection = windowToImageFilter.output_port
        writer.write()
        #data = numpy.asarray(writer.result).tostring()
        return writer.result
        
    @on_trait_change("sources[]")
    def on_sources_changed(self, obj, name, removed, source_list):
        scene = self.scene
        for source in source_list:
            for actor in source.actors:
                scene.add_actor(actor)
            source.on_trait_change(self.trace_all, "update")
            source.on_trait_change(self.render_vtk, "render")
            
        for source in removed:
            scene.remove_actors(source.actors)
            source.on_trait_change(self.trace_all, "update", remove=True)
            source.on_trait_change(self.render_vtk, "render", remove=True)
        self.trace_all()
    
    
#use a singleton handler
controller = RayTraceModelHandler()
    
        
def on_dclick(*obj):
    print("objects", obj)
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
                       Tabbed(
                           Item('ShellObj', editor=ShellEditor(share=False)),
                           Item('constraints', style="custom", editor=ListEditor(use_notebook=True)),
                           show_labels=False
                           ),
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
                    id="raypier.model"
                   ),
                   resizable=True,
                   #height=500,
                   width=800,
                   id="raypier.view",
                   handler=controller,
                   menubar=menubar,
                   )
    
RayTraceModel.class_trait_view("traits_view", ray_tracer_view)
    
        

