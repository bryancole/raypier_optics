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
            Property, List, Instance, on_trait_change, Range, Any,\
            Tuple, Event, cached_property, Set, Int, Trait, Button,\
            self, Str, Bool, PythonValue, Enum
from enthought.traits.ui.api import View, Item, ListEditor, VSplit,\
            RangeEditor, ScrubberEditor, HSplit, VGroup, TextEditor,\
            TupleEditor, VGroup, HGroup, TreeEditor, TreeNode, TitleEditor,\
            ShellEditor
from enthought.tvtk.api import tvtk
from enthought.tvtk.pyface.scene_model import SceneModel
from enthought.tvtk.pyface.scene_editor import SceneEditor
import numpy
import threading
import wx
from itertools import chain, izip, islice
from raytrace.sources import BaseRaySource, collectRays, RayCollection

Vector = Array(shape=(3,))

NumEditor = TextEditor(auto_set=False, enter_set=True, evaluate=float)
ComplexEditor = TextEditor(auto_set=False, enter_set=True, 
                           evaluate=float)

ROField = TextEditor(auto_set=False, enter_set=True, evaluate=float)
VectorEditor = TupleEditor(labels=['x','y','z'], auto_set=False, enter_set=True)

dotprod = lambda a,b: (a*b).sum(axis=-1)[...,numpy.newaxis]

def transformPoints(t, pts):
    """Apply a vtkTransform to a numpy array of points"""
    out = tvtk.Points()
    t.transform_points(pts, out)
    return numpy.asarray(out)

def transformVectors(t, pts):
    """Apply a vtkTransform to a numpy array of vectors.
    I.e. ignore the translation component"""
    out = tvtk.DoubleArray(number_of_components=3)
    t.transform_vectors(pts, out)
    return numpy.asarray(out)

def transformNormals(t, pts):
    """Apply a vtkTransform to a numpy array of normals.
    I.e. ignore the translation component and normalise
    the resulting magnitudes"""
    out = tvtk.DoubleArray(number_of_components=3)
    t.transform_normals(pts, out)
    return numpy.asarray(out)

def normaliseVector(a):
    """normalise a (3,) vector or a (n,3) array of vectors"""
    a= numpy.asarray(a)
    mag = numpy.sqrt((a**2).sum(axis=-1))
    return (a.T/mag).T


def Convert_to_SP(input_v, normal_v, E1_vector, E1_amp, E2_amp):
    """
    All inputs are 2D arrays
    """
    cross = numpy.cross
    
    E2_vector = cross(input_v, E1_vector)
    
    v = cross(input_v, normal_v)
    S_vector = numpy.where(numpy.all(v==0, axis=1).reshape(-1,1)*numpy.ones(3), 
                           normaliseVector(E1_vector),
                           normaliseVector(v) )
    
    v = cross(input_v, S_vector)
    P_vector = normaliseVector(v)
    
    S_amp = E1_amp*dotprod(E1_vector,S_vector) + E2_amp*dotprod(E2_vector, S_vector)
                
    P_amp = E1_amp*dotprod(E1_vector,P_vector) + E2_amp*dotprod(E2_vector, P_vector)
    
    return S_amp, P_amp, S_vector, P_vector

    
    
class Direction(HasTraits):
    x = Float
    y = Float
    z = Float


class ModelObject(HasTraits):
    name = Str("A traceable component")
    
    centre = Tuple(0.,0.,0.) #position
    
    _orientation = Tuple(float, float)
    
    orientation = Property(Range(-180.0,180.0), depends_on="_orientation")
    elevation = Property(Range(-180.,180.), depends_on="_orientation")
    
    rotation = Range(-180.0,180.0, value=0.0) #rotation around orientation axis
    
    direction_btn = Button("set")
    direction = Property(Tuple(float,float,float),depends_on="_orientation")
    
    dir_x = Property(depends_on="direction")
    dir_y = Property(depends_on="direction")
    dir_z = Property(depends_on="direction")
    
    transform = Instance(tvtk.Transform, (), transient=True)
    
    display = Enum("shaded", "wireframe", "hidden")
    
    actors = Instance(tvtk.ActorCollection, transient=True)
    render = Event() #request rerendering (but not necessarily re-tracing)
    
    def _display_changed(self, vnew):
        if vnew=="shaded":
            for actor in self.actors:
                actor.visibility = True
                actor.property.representation = "surface"
        elif vnew=="wireframe":
            for actor in self.actors:
                actor.visibility = True
                actor.property.representation = "wireframe"
        else:
            for actor in self.actors:
                actor.visibility = False
        self.render = True
    
    def get_actors(self, scene):
        return self.actors
    
    def _direction_btn_changed(self):
        d = self.direction
        D = Direction(x=d[0],y=d[1],z=d[2])
        D.edit_traits(kind='modal')
        self.direction = D.x, D.y, D.z
    
    def _get_orientation(self): return self._orientation[0]
    
    def _set_orientation(self, v): self._orientation = (v, self._orientation[1])
    
    def _get_elevation(self): return self._orientation[1]
    
    def _set_elevation(self, v): self._orientation = self._orientation[0], v
    
    def _get_dir_x(self): return self.direction[0]
    
    def _get_dir_y(self): return self.direction[1]
    
    def _get_dir_z(self): return self.direction[2]
    
    @on_trait_change("_orientation, centre, rotation")
    def on_position(self):
        trans = self.transform
        trans.identity()
        trans.translate(*self.centre)
        o,e = self._orientation
        trans.rotate_z(o)
        trans.rotate_x(e)
        trans.rotate_z(self.rotation)
        #print "set transform", self._orientation
        self.update = True
        
    @cached_property
    def _get_direction(self):
        temp = tvtk.Transform()
        o,e = self._orientation
        temp.rotate_z(o)
        temp.rotate_x(e)
        temp.rotate_z(self.rotation)
        direct = temp.transform_point(0,0,1)
        #print "get direction", direct, o, e
        return direct
    
    def _set_direction(self, d):
        x,y,z = normaliseVector(d)
        Theta = numpy.arccos(z)
        theta = 180*Theta/numpy.pi
        phi = 180*numpy.arctan2(x,y)/numpy.pi
        #print "set direction", -phi, -theta
        self._orientation = -phi, -theta
        
        
class Probe(ModelObject):
    pass
        
    
class Traceable(ModelObject):
    vtkproperty = Instance(tvtk.Property, transient=True)
    
    tolerance = Float(0.0001) #excludes ray-segments below this length

    update = Event() #request re-tracing
    
    intersections = List([])
    
    #all Traceables have a pipeline to generate a VTK visualisation of themselves
    pipeline = Any(transient=True) #a tvtk.DataSetAlgorithm ?
    
    polydata = Property(depends_on=['update', ])
    
    def _actors_default(self):
        pipeline = self.pipeline
        
        map = tvtk.PolyDataMapper(input=pipeline.output)
        act = tvtk.Actor(mapper=map)
        act.property = self.vtkproperty
        actors = tvtk.ActorCollection()
        actors.append(act)
        return actors
    
    def _get_polydata(self):
        self.pipeline.update()
        pd = self.pipeline.output
        #print pd
        return pd
    
    def trace_rays(self, rays):
        """traces a RayCollection.
        returns - a recarray of intersetions with the same size as rays
                  dtype=(('length','d'),('cell','i'),('point','d',[3,]))
                  
                  'length' is the physical length from the start of the ray
                  to the intersection
                  
                  'cell' is the ID of the intersecting face
                  
                  'point' is the position coord of the intersection
        """
        
        raise NotImplementedError
    
    def update_complete(self):
        pass
    
    def eval_children(self, rays, points, cells, mask):
        """
        actually calculates the new ray-segments. Physics here
        for Fresnel reflections.
        
        rays - a RayCollection object
        points - a (Nx3) array of intersection coordinates
        cells - a length N array of cell ids
        mask - a bool array selecting items for this Optic
        """
        raise NotImplementedError
    
    def compute_normal(self, point, cell_ids):
        """
        Computes the surface normal vectors for a given array of
        intersection points
        
        @param point: (n,3) array of points
        @param cells: (n,) array of ints giving the cell IDs for each intersection
        
        @return: (n,3) array, the normal vectors (should have unit magnitude)
        """
        raise NotImplementedError


Traceable.uigroup = VGroup(
                   Item('name', editor=TitleEditor(), springy=False,
                        show_label=False),
                   Item('display'),
                   VGroup(
                   Item('orientation', editor=ScrubberEditor()),
                   Item('elevation', editor=ScrubberEditor()),
                   Item('rotation', editor=ScrubberEditor()),
                   ),
                   HGroup(Item('centre', 
                               show_label=False, 
                               editor=VectorEditor,
                               springy=True), 
                          show_border=True,
                          label="Centre"),
                   HGroup(
                          VGroup(
                          Item("dir_x", style="readonly", label="x"),
                          Item("dir_y", style="readonly", label="y"),
                          Item("dir_z", style="readonly", label="z"),
                          springy=True
                          ),
                          Item('direction_btn', show_label=False, width=-60),
                          show_border=True,
                          label="Direction"
                          ),
                    )

    
class Optic(Traceable):
    n_inside = Complex(1.0+0.0j) #refractive
    n_outside = Complex(1.0+0.0j)
    
    all_rays = Bool(False, desc="trace all reflected rays")
    
    vtkproperty = tvtk.Property(opacity = 0.4,
                             color = (0.8,0.8,1.0))
    
    @on_trait_change("n_inside, n_outside")
    def n_changed(self):
        self.update = True
    
    def compute_normal(self, points, cell_ids):
        """
        Evaluate the surface normal in the world frame-of-reference
        @param points: flaot64 ndarray of shape (n,3) giving intersection points
        @param cell_ids: Int ndaray of shape (n,) with cell ids
        """
        raise NotImplementedError
    
    def eval_children(self, rays, points, cells, mask):
        """
        actually calculates the new ray-segments. Physics here
        for Fresnel reflections.
        
        rays - a RayCollection object
        points - a (Nx3) array of intersection coordinates
        cells - a length N array of cell ids
        mask - a bool array selecting items for this Optic
        """
        points = points[mask]
        cells = cells[mask] ###reshape not necessary
        normal = self.compute_normal(points, cells)
        input_v = rays.direction[mask]
        
        parent_ids = numpy.arange(mask.shape[0])[mask]
        optic = numpy.repeat([self,], points.shape[0] )
        
        S_amp, P_amp, S_vec, P_vec = Convert_to_SP(input_v, 
                                                   normal, 
                                                   rays.E_vector[mask], 
                                                   rays.E1_amp[mask], 
                                                   rays.E2_amp[mask])

        #this is cos(theta), where theta is the angle between the
        #normal and the incident ray
        cosTheta = dotprod(normal, input_v)
        
        origin = points
            
        fromoutside = cosTheta < 0
        n1 = numpy.where(fromoutside, self.n_outside.real, self.n_inside.real)
        n2 = numpy.where(fromoutside, self.n_inside.real, self.n_outside.real)
        flip = numpy.where(fromoutside, 1, -1)
            
        abscosTheta = numpy.abs(cosTheta)
        
        N2 = (n2/n1)**2
        N2cosTheta = N2*abscosTheta
        
        #if this is less than zero, we have Total Internal Reflection
        N2_sin2 = abscosTheta**2 + (N2 - 1)
        
        TIR = N2_sin2 < 0.0
        sqrt = numpy.sqrt
        
        cosThetaNormal = cosTheta*normal
        reflected = input_v - 2*cosThetaNormal
        sqrtN2sin2 = numpy.where(TIR, 1.0j*sqrt(-N2_sin2), sqrt(N2_sin2))
        #print "root n2.sin2", sqrtN2sin2
        
        #Fresnel equations for reflection
        R_p = (N2cosTheta - sqrtN2sin2) / (N2cosTheta + sqrtN2sin2)
        R_s = (abscosTheta - sqrtN2sin2) / (abscosTheta + sqrtN2sin2)
        #print "R_s", R_s, "R_p", R_p
        
        ###Now calculate transmitted rays
        d1 = input_v
        tangent = d1 - cosThetaNormal
        
        tan_mag_sq = ((n1*tangent/n2)**2).sum(axis=1).reshape(-1,1)        
        
        c2 = numpy.sqrt(1 - tan_mag_sq)
        transmitted = tangent*(n1/n2) - c2*normal*flip 
        #print d1, normal, tangent, transmitted, "T"
        
        cos1 = abscosTheta
        #cos of angle between outgoing ray and normal
        cos2 = abs(dotprod(transmitted, normal))
        
        Two_n1_cos1 = (2*n1)*cos1
        
        aspect = sqrt(cos2/cos1) * Two_n1_cos1
        
        #Fresnel equations for transmission
        T_p = aspect / ( n2*cos1 + n1*cos2 )
        T_s = aspect / ( n2*cos2 + n1*cos1 )
        #print "T_s", T_s, "T_p", T_p
        
        if self.all_rays:
            refl_rays = RayCollection(origin=origin,
                                       direction = reflected,
                                       max_length = rays.max_length,
                                       E_vector = S_vec,
                                       E1_amp = S_amp*R_s,
                                       E2_amp = P_amp*R_p,
                                       parent_ids = parent_ids,
                                       optic = optic,
                                       face_id = cells,
                                       refractive_index=n1)
            
            trans_rays = RayCollection(origin=origin,
                                       direction = transmitted,
                                       max_length = rays.max_length,
                                       E_vector = S_vec,
                                       E1_amp = S_amp*T_s,
                                       E2_amp = P_amp*T_p,
                                       parent_ids = parent_ids,
                                       optic = optic,
                                       face_id = cells,
                                       refractive_index=n2)
            
            allrays = collectRays(refl_rays, trans_rays)
            allrays.parent = rays
            return allrays
        else:
            TIR.shape=-1,1
            tir = TIR*numpy.ones(3)
            direction = numpy.where(tir, reflected,transmitted)
            E1_amp = S_amp*numpy.where(TIR, R_s, T_s)
            E2_amp = P_amp*numpy.where(TIR, R_p, T_p)
            refractive_index = numpy.where(TIR, n1, n2)
            
            return RayCollection(origin=origin,
                               direction = direction,
                               max_length = rays.max_length,
                               E_vector = S_vec,
                               E1_amp = E1_amp,
                               E2_amp = E2_amp,
                               parent_ids = parent_ids,
                               optic = optic,
                               face_id = cells,
                               refractive_index=refractive_index) 
    
    
class VTKOptic(Optic):
    """Polygonal optics using a vtkOBBTree to do the ray intersections"""
    data_source = Instance(tvtk.ProgrammableSource, (), transient=True)
    
    obb = Instance(tvtk.OBBTree, (), transient=True)
    
    def __getstate__(self):
        d = super(VTKOptic, self).__getstate__()
        bad = ['polydata','pipeline','data_source','obb', 'transform']
        for b in bad:
            d.pop(b, None)
        return d
    
    def __init__(self, *args, **kwds):
        super(VTKOptic, self).__init__(*args, **kwds)
        self.on_trait_change(self.on__polydata_changed, "_polydata")
        self.on__polydata_changed()
        
    def on__polydata_changed(self):
        self.data_source.modified()
        self.update = True
    
    def _polydata_changed(self):
        obb = self.obb
        obb.free_search_structure()
        obb.data_set = self.polydata
        obb.build_locator()
        self.data_source.modified()
        
    def _pipeline_default(self):
        source = self.data_source
        def execute():
            polydata = self._polydata
            output = source.poly_data_output
            output.shallow_copy(polydata)
        source.set_execute_method(execute)
        t = self.transform
        transf = tvtk.TransformFilter(input=source.output, transform=t)
        tri = tvtk.TriangleFilter(input=transf.output)
        return tri
        
    def trace_segment(self, seg, last_optic=None, last_cell=None):
        """
        Finds the intersection of the given ray-segment with the 
        object geometry data
        """
        p1 = seg.origin
        p2 = p1 + seg.MAX_RAY_LENGTH*seg.direction
        pts = tvtk.Points()
        ids = tvtk.IdList()
        ret = self.obb.intersect_with_line(p1, p2, pts, ids)
        sqrt = numpy.sqrt
        array = numpy.array
        if ret==0:
            return None
        if self is not last_optic:
            last_cell = None
        data = [ (sqrt(((array(p) - p1)**2).sum()), p, Id)
                    for p, Id in izip(pts, ids) 
                    if Id is not last_cell ]
        if not data:
            return None
        short = min(data, key=lambda a: a[0])
        return short[0], short[1], short[2], self
        
    
class BeamStop(VTKOptic):
    name = "Beam Stop"
    height = Float #distance between parallel faces
    width = Float #width of parallel faces
    
    _polydata = Property(depends_on=['height', 'width'])
    
    ray_count = Property(Int, depends_on=['intersections_items'])
    ave_ray_length = Float(0.0)
    ave_power = Float(0.0)
    
    traits_view = View(VGroup(
                       Traceable.uigroup,
                       Item('ray_count', style="readonly",height=-25),
                       Item('ave_ray_length', label="path length",
                            style="readonly",height=-25),
                       Item('ave_power', label="incident power",
                            style="readonly",height=-25)
                       )
                       )
    
    @cached_property
    def _get__polydata(self):
        h = self.height/2
        w = self.width/2
        points = [(-w,-h,0),
                  (w,-h,0),
                  (w,h,0),
                  (-w,h,0)]
        cells = [(0,1,2),
                 (2,3,0)
                 ]
        pd = tvtk.PolyData(points=numpy.array(points), 
                           polys=numpy.array(cells))
        return pd
    
    def _get_ray_count(self):
        return len(self.intersections)
    
    def eval_children(self, seg, point, cell_id):
        return []
    
    def update_complete(self):
        try:
            ave = sum(seg.cum_length for seg in self.intersections)/len(self.intersections)
        except ZeroDivisionError:
            ave = 0.0
        self.ave_ray_length = ave
        
#        for seg in self.intersections:
#            print seg.E1_amp, seg.E2_amp
        
        try:
            pwr = sum(self.eval_pwr(seg) for seg in self.intersections)/len(self.intersections)
        except ZeroDivisionError:
            pwr = 0.0
        self.ave_power = pwr
        
    @staticmethod
    def eval_pwr(seg):
        E1 = seg.E1_amp
        E2 = seg.E2_amp
        P1 = (abs(E1)**2)
        P2 = (abs(E2)**2)
        return P1 + P2
    
    
class RayTraceModel(HasTraits):
    scene = Instance(SceneModel, (), {'background':(1,1,0.8)}, transient=True)
    
    optics = List(Traceable)
    sources = List(BaseRaySource)
    probes = List(Probe)
    
    optical_path = Float(0.0, transient=True)
    
    update = Event()
    
    Self = self
    ShellObj = PythonValue(transient=True)
        
    recursion_limit = Int(10, desc="maximum number of refractions or reflections")
    
    def _optics_changed(self, opticList):
        scene = self.scene
        #del scene.actor_list[:]    
        for o in opticList:
            scene.add_actors(o.get_actors(scene))
        
        for optic in opticList:
            optic.on_trait_change(self.trace_all, "update")
            optic.on_trait_change(self.render_vtk, "render")
        self.update = True
    
    def _rays_changed(self, rayList):
        scene = self.scene
        sources = [o.pipeline for o in rayList]
        mappers = [tvtk.PolyDataMapper(input=s.output) for s in sources]
        actors = [tvtk.Actor(mapper=m) for m in mappers]
        for actor in actors:
            property = actor.property
            property.color = (1,0.5,0)
        scene.add_actors(actors)
        self.update = True
        
    def _probes_changed(self, probeList):
        scene = self.scene
        #del scene.actor_list[:]    
        for p in probeList:
            scene.add_actors(p.get_actors(scene))
        for probe in probeList:
            probe.on_trait_change(self.update_probes, "update")
            probe.on_trait_change(self.render_vtk, "render")
        self.update = True
        
    def update_probes(self):
        if self.scene is not None:
            self.render_vtk()
        
    @on_trait_change("update")
    def trace_all(self):
        optics = self.optics
        if optics:
            for o in optics:
                o.intersections = []
            for ray_source in self.sources:
                self.trace_ray_source(ray_source, optics)
            for o in optics:
                o.update_complete()
        self.render_vtk()
        
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
        rays = ray_source.InputRays
        traced_rays = []
        limit = self.recursion_limit
        count = 0
        while (rays is not None) and count<limit:
            #print "count", count
            traced_rays.append(rays)
            rays = self.trace_segment(rays, optics)
            count += 1
        ray_source.TracedRays = traced_rays
        ray_source.data_source.modified()
            
    def trace_ray_source_detail(self, ray_source, optics):
        """trace a ray source, using it's detailed rays"""
        rays = ray_source.InputDetailRays
        traced_rays = []
        limit = self.recursion_limit
        count = 0
        while (rays is not None) and count<limit:
            #print "count", count
            traced_rays.append(rays)
            rays = self.trace_segment(rays, optics)
            count += 1
        ray_source.TracedDetailRays = traced_rays
        
    def trace_segment(self, rays, optics):
        """trace a RayCollection"""
        size = rays.origin.shape[0]
        #optic.trace_rays should return a 1d recarray with 
        #dtype=(('length','f8'),('cell','i2'),('point','f8',[3,]))
        intersections = numpy.column_stack([o.trace_rays(rays) for o in optics])
        shortest = numpy.argmin(intersections['length'], axis=1)
        ar = numpy.arange(size)
        lengths = intersections['length'][ar,shortest]
        
        #now remove infinite rays
        mask = lengths!=numpy.Infinity
        shortest = shortest[mask]
        ar = ar[mask]
        
        optic_ixt = numpy.choose(shortest, optics)
        cell_ids = intersections['cell'][ar,shortest]
        if intersections.size==1:
            points = intersections['point'][ar,shortest,:]
        else:
            points = intersections['point'][:,ar,shortest].T
        lengthsT = lengths.reshape(-1,1)
        #print "shape", lengthsT.shape, lengthsT.dtype, lengthsT[:5]
        rays.length = lengthsT
        pair_mask = ((o, o==optic_ixt) for o in numpy.unique(optic_ixt))
        children = [o.eval_children(rays, points, cell_ids, m) 
                    for o,m in pair_mask]
        if len(children)==0:
            return None
        new_rays = collectRays(*children)
        new_rays.parent = rays
        return new_rays
        
    def _sources_changed(self, source_list):
        scene = self.scene
        for source in source_list:
            for actor in source.actors:
                scene.add_actor(actor)
            source.on_trait_change(self.trace_all, "update")
        self.update = True
        
        
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
                        auto_open=True,
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
                       
                       ],
                orientation='vertical',
                hide_root=True
                )
        
    
ray_tracer_view = View(
                   HSplit(
                    VSplit(
                       Item('scene',editor=SceneEditor(),
                            height=600),
                       #Item('optics@', editor=ListEditor(use_notebook=True),
                       #     width=200),
                       Item('ShellObj', editor=ShellEditor()),
                       show_labels=False,
                       dock="vertical"
                       ),
                       Item('Self', 
                            id="TracerModelID",
                            editor=tree_editor, width=200),
                    show_labels=False,
                    dock="horizontal",
                    id="raytrace.model"
                   ),
                   resizable=True,
                   #height=500,
                   width=800,
                   id="raytrace.view"
                   )
    
RayTraceModel.class_trait_view("traits_view", ray_tracer_view)
    
        

