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


from traits.api import HasTraits, Array, BaseFloat, Complex,\
            Property, List, Instance, Range, Any,\
            Tuple, Event, cached_property, Set, Int, Trait, Button,\
            self, Str, Bool, PythonValue, Enum, MetaHasTraits
from traitsui.api import View, Item, ListEditor, VSplit,\
            RangeEditor, ScrubberEditor, HSplit, VGroup, TextEditor,\
            TupleEditor, VGroup, HGroup, TreeEditor, TreeNode, TitleEditor,\
            ShellEditor
            
from traitsui.file_dialog import save_file
            
from tvtk.api import tvtk
import numpy
import threading, os, itertools
from itertools import chain, islice, count
import yaml
from abc import ABCMeta

from raypier.constraints import BaseConstraint
from raypier.has_queue import HasQueue, on_trait_change
from raypier.utils import normaliseVector, transformNormals, transformPoints,\
        transformVectors, dotprod
from raypier.core import ctracer, cmaterials
from raypier.core.ctracer import RayCollection, GaussletCollection
from .shapes import BaseShape
from raypier.vtk_algorithms import EmptyGridSource
from .editors import NumEditor, IntEditor, ComplexEditor, ROField, VectorEditor

Vector = Array(shape=(3,))

counter = count()


class BaseRayCollection(metaclass=ABCMeta):
    pass
BaseRayCollection.register(RayCollection)
BaseRayCollection.register(GaussletCollection)


class Float(BaseFloat):
    def validate(self, obj, name, value):
        return float(value)


class Direction(HasTraits):
    x = Float
    y = Float
    z = Float


class RaypierObjectMetaclass(MetaHasTraits):
    """
    The metaclass for YAMLObject.
    """
    def __init__(cls, name, bases, kwds):
        super(RaypierObjectMetaclass, cls).__init__(name, bases, kwds)
        if 'yaml_tag' in kwds and kwds['yaml_tag'] is not None:
            pass
        else:
            cls.yaml_tag = "!"+name
        cls.yaml_loader.add_constructor(cls.yaml_tag, cls.from_yaml)
        cls.yaml_dumper.add_representer(cls, cls.to_yaml)
        
        if not cls.abstract:
            cls.subclasses.add(cls)


class RaypierObject(object, metaclass=RaypierObjectMetaclass):
    """
    An object that can dump itself to a YAML stream
    and load itself from a YAML stream.
    """
    __slots__ = ()  # no direct instantiation, so allow immutable subclasses

    yaml_loader = yaml.Loader
    yaml_dumper = yaml.Dumper

    yaml_tag = None
    yaml_flow_style = None
    
    abstract = True
    subclasses = set()

    def from_yaml(cls, loader, node):
        """
        Convert a representation node to a Python object.
        """
        return loader.construct_yaml_object(node, cls)
    from_yaml = classmethod(from_yaml)

    def to_yaml(cls, dumper, data):
        """
        Convert a Python object to a representation node.
        """
        return dumper.represent_yaml_object(cls.yaml_tag, data, cls,
                flow_style=cls.yaml_flow_style)
    to_yaml = classmethod(to_yaml)
    
    @classmethod
    def class_trait_view_elements(cls):
        return None


class Renderable(HasQueue, RaypierObject, metaclass=RaypierObjectMetaclass):
    display = Enum("shaded", "wireframe", "hidden")
    
    actors = Instance(tvtk.ActorCollection, (), transient=True)
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


class ModelObject(Renderable):
    abstract=True
    name = Str("A traceable component")
    
    update = Event() #request re-tracing
    
    centre = Tuple(0.,0.,0.) #position
    
    _orientation = Tuple(Float, Float)
    
    orientation = Property(Range(-180.0,180.0), transient=True,
                           depends_on="_orientation")
    elevation = Property(Range(-180.,180.), transient=True,
                         depends_on="_orientation")
    
    rotation = Range(-180.0,180.0, value=0.0) #rotation around orientation axis
    
    direction_btn = Button("set")
    direction = Property(Tuple(float,float,float),depends_on="_orientation")
    
    x_axis = Property(Tuple(float,float,float),depends_on="_orientation, rotation")
    
    dir_x = Property(depends_on="direction")
    dir_y = Property(depends_on="direction")
    dir_z = Property(depends_on="direction")
    
    transform = Instance(tvtk.Transform, (), transient=True)
    
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
        
    def translate_local(self, du, dv, dw):
        """Translate the object along it's local axes, where
        dw is along the direction-vector, du and dv are along the orthogonal
        axes to this"""
        centre = numpy.array(self.centre)
        w = numpy.array(self.direction)
        u = numpy.array(self.x_axis)
        v = numpy.cross(w,u)
        new_centre = centre + (du*u + dv*v + dw*w)
        self.centre = tuple(new_centre)
        
    def _get_x_axis(self):
        temp = tvtk.Transform()
        o,e = self._orientation
        temp.rotate_z(o)
        temp.rotate_x(e)
        temp.rotate_z(self.rotation)
        direct = temp.transform_point(1,0,0)
        return direct
        
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
        
    def make_step_shape(self):
        """Creates an OpenCascade BRep Shape
        representation of the object, which can be
        exported to STEP format"""
        return False, None
        
        
class Probe(ModelObject):
    """An object which calculates some result from the raypier output, based on the
    probe geometry.
    
    Although probes inheret the 'update' trait from its superclass, this event doesn't trigger
    a new tracing operation.
    """ 
    
    def evaluate(self, src_list):
        """called once after all sources have been traced"""
        
    
class Traceable(ModelObject):
    vtkproperty = Instance(tvtk.Property, transient=True)
    
    intersections = List([], transient=True)
    
    material = Instance(ctracer.InterfaceMaterial)
    
    faces = Instance(ctracer.FaceList, 
                desc="Container of traceable faces (Face instances)",
                 transient = True)
    
    #all Traceables have a pipeline to generate a VTK visualisation of themselves
    pipeline = Any(transient=True) #a tvtk.DataSetAlgorithm ?
    
    polydata = Property(depends_on=['update', ])
    
    abstract=True
    subclasses = set()
    
    def _actors_default(self):
        pipeline = self.pipeline
        
        map = tvtk.PolyDataMapper(input_connection=pipeline.output_port)
        map.scalar_visibility = False
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
        
        @param rays: a RayCollection instance
        @param face_id: id of the face to trace. If None, trace all faces
        
        returns - a recarray of intersetions with the same size as rays
                  dtype=(('length','d'),('face','O'),('point','d',[3,]))
                  
                  'length' is the physical length from the start of the ray
                  to the intersection
                  
                  'cell' is the ID of the intersecting face
                  
                  'point' is the position coord of the intersection (in world coords)
        """
        max_length = rays.max_length
        p1 = rays.origin
        p2 = p1 + max_length*rays.direction
        return self.intersect(p1, p2, max_length)
        
    def intersect(self, p1, p2, max_length):
        t = self.transform
        inv_t = t.linear_inverse
        P1 = transformPoints(inv_t, p1)
        P2 = transformPoints(inv_t, p2)
        
        faces = self.faces
        if len(faces)>1:
            traces = numpy.column_stack([f.intersect(P1, P2, max_length) for f in faces])   
            nearest = numpy.argmin(traces['length'], axis=1)
            ar = numpy.arange(traces.shape[0])
            shortest = traces[ar,nearest]
        else:
            shortest = faces[0].intersect(P1, P2, max_length)
            
        t_points = shortest['point']
        points = transformPoints(t, t_points)
        shortest['point'] = points
        return shortest
    
    def intersect_line(self, p1, p2):
        """Find the nearest intersection of a single line defined by two points
        p1 and p2
        
        returns a tuple (L, F, P), where L is the scalar length from p1 to the
        intersection point, F is the intersecting face and P is the intersection
        point (a 3-vector).
        """
        p1 = numpy.asarray(p1).reshape(-1,1)
        p2 = numpy.asarray(p2).reshape(-1,1)
        max_length = ((p1-p2)**2).sum(axis=0)[0]
        nearest = self.intersect(p1, p2, max_length)[0]
        return nearest['length'], nearest['cell'], nearest['point']
    
    def update_complete(self):
        pass

OrientationEditor = RangeEditor(low=-180.0, high=180.0, mode='slider')


Traceable.uigroup = VGroup(
                   Item('name', editor=TitleEditor(), springy=False,
                        show_label=False),
                   Item('display'),
                   VGroup(
                   Item('orientation', editor=OrientationEditor),
                   Item('elevation', editor=OrientationEditor),
                   Item('rotation', editor=OrientationEditor),
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
                   label="Geometry"
                    )


class ShapedTraceable(Traceable):
    shape = Instance(BaseShape)
    
    grid_extent = Tuple((-50.,50.,-50.,50.0,-5.,5.))
    grid_resolution = Tuple((50,50,30))
    grid_in = Instance(EmptyGridSource, ())    
    
    clip = Instance(tvtk.ClipVolume, ())
    
    def eval_sag_top(self, points):
        return None
    
    def eval_sag_bottom(self, points):
        return None
    
    def eval_grid_extent(self):
        return (-50.,50.,-50.,50.0,-5.,5.)
    
    @on_trait_change("grid_extent", "grid_resolution")
    def _update_grid(self):
        grid = self.grid_in
        e = self.grid_extent
        d = self.grid_resolution
        grid.dimensions = d
        grid.spacing = (e[1]-e[0])/(d[0]-1), (e[3]-e[2])/(d[1]-1), (e[5]-e[4])/(d[2]-1)
        grid.origin = (e[0],e[2],e[4])
        grid.modified()
        self.update=True
        
    def _grid_extent_default(self):
        return self.eval_grid_extent()
    
    def build_pipeline(self, input_node):
        return input_node
    
    def _pipeline_default(self):
        grid = self.grid_in
        self._update_grid()
        
        shape = self.shape
        func = shape.impl_func
        clip = self.clip 
        clip.input_connection=grid.output_port
        clip.clip_function=func
        clip.inside_out=1
        
        attrb = tvtk.ProgrammableAttributeDataFilter(input_connection=clip.output_port)
        def update( *args ):
            in_data = attrb.get_input_data_object(0,0)
            points = in_data.points.to_array()
            top = self.eval_sag_top(points)
            bottom = -self.eval_sag_bottom(points)
            if top is None:
                top = bottom
            if bottom is None:
                bottom = top
            if top is not bottom:
                result = numpy.maximum(top, bottom) #Or maybe I need maximum
            else:
                result = top
            out = attrb.get_output_data_object(0)
            out.point_data.scalars=result
        attrb.set_execute_method(update)
        
        clip2 = tvtk.ClipDataSet(input_connection=attrb.output_port,
                                 value=0.0,
                                 inside_out=1)
        
        #topoly = tvtk.GeometryFilter(input_connection=clip2.output_port)
        topoly = tvtk.DataSetSurfaceFilter(input_connection=clip2.output_port)
        norms = tvtk.PolyDataNormals(input_connection=topoly.output_port)
        
        node = self.build_pipeline(norms)
        
        transF = tvtk.TransformFilter(input_connection=node.output_port, 
                                      transform=self.transform)
        #self.grid_extent = extent
        return transF
    
    
class Optic(Traceable):
    n_inside = Complex(1.0+0.0j) #refractive
    n_outside = Complex(1.0+0.0j)
    
    all_rays = Bool(False, desc="trace all reflected rays")
    
    vtkproperty = Instance(tvtk.Property, (), {'opacity': 0.4, 'color': (0.8,0.8,1.0)})
                             
    def _material_default(self):
        m = cmaterials.FullDielectricMaterial(n_inside = self.n_inside,
                                    n_outside = self.n_outside)
        return m
    
    def calc_refractive_index(self, wavelengths):
        """
        Evaluates an array of (complex) refractive indices.
        
        :param wavelengths: a shape=(N,1) array of wavelengths
        
        :returns: a 2-tuple representing the inside and outside
        
        refractive indices respectively. The items in the tuple can be
        either:
         
          1) an arrays with the same shape as wavelengths and with
             dtype=numpy.complex128
          2) a complex scalar
        """
        return self.n_inside, self.n_outside
    
    @on_trait_change("n_inside, n_outside")
    def n_changed(self):
        self.material.n_inside = self.n_inside
        self.material.n_outside = self.n_outside
        self.update = True
        
        
class ShapedOptic(ShapedTraceable, Optic):
    vtkproperty = Instance(tvtk.Property, (), {'opacity': 0.4, 'color': (0.8,0.8,1.0)})
    
    def _material_default(self):
        m = cmaterials.FullDielectricMaterial(n_inside = self.n_inside,
                                    n_outside = self.n_outside)
        return m
    
    def calc_refractive_index(self, wavelengths):
        """
        Evaluates an array of (complex) refractive indices.
        @param wavelengths: a shape=(N,1) array of wavelengths
        @returns: a 2-tuple representing the inside and outside
        refractive indices respectively. The items in the tuple can be
        either 1) an arrays with the same shape as wavelengths and with
                    dtype=numpy.complex128
            or 2) a complex scalar
        """
        return self.n_inside, self.n_outside
    
    @on_trait_change("n_inside, n_outside")
    def n_changed(self):
        self.material.n_inside = self.n_inside
        self.material.n_outside = self.n_outside
        self.update = True
    
    
class VTKOptic(Optic):
    """Polygonal optics using a vtkOBBTree to do the ray intersections"""
    abstract=True
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
        transf = tvtk.TransformFilter(input_connection=source.output_port, 
                                      transform=t)
        tri = tvtk.TriangleFilter(input_connection=transf.output_port)
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
                    for p, Id in zip(pts, ids) 
                    if Id is not last_cell ]
        if not data:
            return None
        short = min(data, key=lambda a: a[0])
        return short[0], short[1], short[2], self
    
    
class Result(HasTraits, RaypierObject):
    abstract=True
    subclasses = set()
    name = Str("a result")
    
    def calc_result(self, tracer):
        """
        Called at the end of a tracing operation, so the result can 
        be evaluated
        
        @param tracer: the RayTraceModel instance
        """
        raise NotImplementedError
