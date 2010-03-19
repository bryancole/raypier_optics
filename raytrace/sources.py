#    Copyright 2009, Teraview Ltd., Bryan Cole
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


import numpy
import itertools

from enthought.traits.api import HasTraits, Int, Float, \
     Bool, Property, Array, Event, List, cached_property, Str,\
     Instance, on_trait_change, Trait, Enum, Title

from enthought.traits.ui.api import View, Item, Tabbed, VGroup, Include, \
    Group

from enthought.tvtk.api import tvtk

from raytrace.ctracer import RayCollection, Ray, ray_dtype
from raytrace.utils import normaliseVector, Range, TupleVector, Tuple
from raytrace.bases import Renderable

Vector = Array(shape=(3,))



class BaseRaySource(HasTraits):
    name = Title("Ray Source")
    update = Event()
    display = Enum("pipes", "wires", "hidden")
    render = Event()
    mapper = Instance(tvtk.PolyDataMapper, ())
    
    InputRays = Property(Instance(RayCollection), depends_on="max_ray_len")
    TracedRays = List(RayCollection)
    
    InputDetailRays = Property(Instance(RayCollection), depends_on="InputRays")
    TracedDetailRays = List(RayCollection)

    detail_resolution = Int(32)

    max_ray_len = Float(200.0)
    
    scale_factor = Float(0.2)
    
    show_start = Bool(True)
    show_polarisation = Bool(False)
    show_direction = Bool(False)
    show_normals = Bool(False)
    
    #idx selector for the input ways which should be visualised
    view_ray_ids = Trait(None, Array(dtype=numpy.int))
    
    tube = Instance(tvtk.TubeFilter, (), transient=True)
    sphere = Instance(tvtk.SphereSource, (), transient=True)
    data_source = Instance(tvtk.ProgrammableSource, transient=True)
    
    normals_source = Instance(tvtk.ProgrammableSource, transient=True)
    normals_actor = Instance(tvtk.Actor, transient=True)
    
    ray_actor = Instance(tvtk.Actor, transient=True)
    start_actor = Instance(tvtk.Actor, transient=True)
    #field_actor = Instance(tvtk.Actor)
    
    export_pipes = Bool(False, desc="if set, rays are exported to STEP files as pipes")
    
    actors = Property(List, transient=True)
    
    vtkproperty = Instance(tvtk.Property, (), {'color':(1,0.5,0)}, transient=True)
    
    display_grp = VGroup(Item('display'),
                       Item('show_start'),
                       Item('show_normals'),
                       Item('scale_factor'),
                       Item('export_pipes',label="export as pipes"),
                       label="Display")
    
    traits_view = View(Item('name', show_label=False),
                       Tabbed(Include('geom_grp'),
                              display_grp)
                       )
    
    def __repr__(self):
        return self.name
    
    def get_sequence_to_face(self, face):
        """returns a list of list of Face objects, those encountered
        on the route to the target face"""
        #find the first RayCollection which contains the target face
        for rays in self.TracedRays[1:]:
            faces = list(rays.face)
            if face in faces:
                break
        else:
            raise ValueError("no rays trace to this face")
        
        seq = []
        
        #now iterate back up the ray-tree collecting only the faces
        #on the path to the target face
        parent = rays.parent
        ids = rays.parent_ids[face==rays.face]
        while True:
            if parent == self.TracedRays[0]:
                break
            faces = list(numpy.unique1d(parent.face[ids]))
            seq.append(faces)
            ids = parent.parent_ids[ids]
            rays = parent
            parent = rays.parent
        seq.reverse()
        return seq
    
    def eval_angular_spread(self, idx):
        """A helper method to evaluate the angular spread of a ray-segment.
        @param idx: the index of the RayCollection in the TracedRay list 
                    to analyse
                    
        @returns: the average angle from the mean ray, in degrees"""
        rays = self.TracedRays[idx]
        ave_dir = normaliseVector(rays.direction.mean(axis=0))
        dotprod = (ave_dir[numpy.newaxis,:] * rays.direction).sum(axis=1)
        angles = numpy.arccos(dotpod)*180 / numpy.pi
        return angles.mean()
    
    def make_step_shape(self):
        if self.export_pipes:
            from raytrace.step_export import make_rays_both as make_rays
        else:
            from raytrace.step_export import make_rays_wires as make_rays
        return make_rays(self.TracedRays, self.scale_factor), "red"
    
    def _TracedRays_changed(self):
        self.data_source.modified()
        self.normals_source.modified()
        
    def _get_InputRays(self):
        return None
    
    def _InputRays_changed(self):
        self.update=True
        
    def _scale_factor_changed(self, scale):
        self.tube.radius = scale
        self.render = True
    
    def _get_actors(self):
        actors = [self.ray_actor, self.start_actor, self.normals_actor]
        #return actors
        return actors[:2]  #makeshift turning off of normal glyphs
    
    def _normals_source_default(self):
        source = tvtk.ProgrammableSource()
        def execute():
            output = source.poly_data_output
            points = numpy.vstack([p.origin for p in self.TracedRays])
            normals = numpy.vstack([p.normals for p in self.TracedRays])
            output.points = points
            output.point_data.normals = normals
            print "calc normals GLYPH"
        source.set_execute_method(execute)
        return source
        
        
    def _normals_actor_default(self):
        source = self.normals_source
        glyph = tvtk.ArrowSource()
        glyph_filter = tvtk.Glyph3D(source=glyph.output,
                                    input = source.output,
                                    scale_factor=10.0,
                                    vector_mode='use_normal'
                                    )
        map = tvtk.PolyDataMapper(input=glyph_filter.output)
        act = tvtk.Actor(mapper=map)
        act.property.color = (0.0,1.0,0.0)
        return act
    
    def _data_source_default(self):
        source = tvtk.ProgrammableSource()
        def execute():
            output = source.poly_data_output
            pointArrayList = []
            for rays in self.TracedRays:                
                start_pos = [r.origin for r in rays]
                end_pos = [r.termination for r in rays]
                #print "start", start_pos
                #print "end", end_pos
                interleaved = numpy.array([start_pos, end_pos]).swapaxes(0,1).reshape(-1,3)
                pointArrayList.append(interleaved)
            if pointArrayList:
                points = numpy.vstack(pointArrayList)
                output.points=points
                lines = list([i,i+1] for i in xrange(0,points.shape[0],2))
                output.lines = lines
        source.set_execute_method(execute)
        return source
    
    def _ray_actor_default(self):
        tube = self.tube
        tube.input=self.data_source.output
        tube.number_of_sides = 20
        
        map = self.mapper
        act = tvtk.Actor(mapper=map)
        act.visibility = True
        display = self.display
        if display=="pipes":
            map.input=tube.output
        elif display=="wires":
            map.input = self.data_source.output
        else:
            act.visibility = False
        prop = self.vtkproperty
        if prop:
            act.property = prop
        return act
    
    def _display_changed(self, vnew):
        actors = self.actors
        for act in actors:
            act.visibility = True
        if vnew=="pipes":
            self.mapper.input = self.tube.output
        elif vnew=="wires":
            self.mapper.input = self.data_source.output
        elif vnew=="hidden":
            for act in actors:
                act.visibility = False
        self.render = True
    
    def _start_actor_default(self):
        map = tvtk.PolyDataMapper(input=self.sphere.output)
        act = tvtk.Actor(mapper=map)
        return act
    
    
class ParallelRaySource(BaseRaySource):
    origin = Tuple((0.,0.,0.))
    direction = Tuple((0.,0.,1.))
    number = Int(20)
    radius = Float(10.)
    
    view_ray_ids = numpy.arange(20)
    
    InputRays = Property(Instance(RayCollection), 
                         depends_on="origin, direction, number, radius, max_ray_len")
    
    geom_grp = VGroup(Group(Item('origin', show_label=False,resizable=True), 
                            show_border=True,
                            label="Origin position",
                            padding=0),
                       Group(Item('direction', show_label=False, resizable=True),
                            show_border=True,
                            label="Direction"),
                       Item('number'),
                       #Item('rings'),
                       Item('radius'),
                       label="Geometry")
    
    
    @on_trait_change("focus, direction, number, radius, max_ray_len")
    def on_update(self):
        self.data_source.modified()
        self.update=True
    
    @cached_property
    def _get_InputRays(self):
        print "making rays"
        origin = numpy.array(self.origin)
        direction = numpy.array(self.direction)
        count = self.number
        radius = self.radius
        
        max_axis = numpy.abs(direction).argmax()
        if max_axis==0:
            v = numpy.array([0.,1.,0.])
        else:
            v = numpy.array([1.,0.,0.])
        d1 = numpy.cross(direction, v)
        d1 = normaliseVector(d1)
        d2 = numpy.cross(direction, d1)
        d2 = normaliseVector(d2)
        angles = (i*2*numpy.pi/count for i in xrange(count))
        offsets = [radius*(d1*numpy.sin(a) + d2*numpy.cos(a)) for a in angles]
        origins = [origin + offset for offset in offsets]
#        directions = numpy.ones_like(origins) * direction
#        ray_data = numpy.zeros(count, dtype=ray_dtype)
#        ray_data['E_vector'] = [[1,0,0]]
#        ray_data['E1_amp'] = 1.0 + 0.0j
#        ray_data['E2_amp'] = 0.0
#        ray_data['refractive_index'] = 1.0+0.0j
#        ray_data['normal'] = [[0,1,0]]
#        rays = RayCollection.from_array(ray_data)
        print "made rays"
        rays = RayCollection(len(origins))
        for i in xrange(len(origins)):
            ray = Ray(origin=origins[i], direction=direction,
                        length=numpy.Inf)
            rays.add_ray(ray)
        return rays
    
    

class ConfocalRaySource(BaseRaySource):
    focus = TupleVector
    direction = TupleVector
    number = Range(1,50,20, editor_traits={'mode':'spinner'})
    theta = Range(0.0,90.0,value=30)
    working_dist = Float(100.0)
    rings = Range(1,50,3, editor_traits={'mode':'spinner'})
    
    principle_axes = Property(Tuple(Array,Array), depends_on="direction")
    
    #view_ray_ids = numpy.arange(20)
    
    InputRays = Property(Instance(RayCollection), 
                         depends_on="focus, direction, number, rings, theta, working_dist, max_ray_len")
    
    geom_grp = VGroup(Group(Item('focus', show_label=False,resizable=True), 
                            show_border=True,
                            label="Focus position",
                            padding=0),
                       Group(Item('direction', show_label=False, resizable=True),
                            show_border=True,
                            label="Direction"),
                       Item('number'),
                       Item('rings'),
                       Item('theta', label="angle"),
                       Item('working_dist', label="working dist"),
                       label="Geometry")
    
    @cached_property
    def _get_principle_axes(self):
        direction = self.direction
        max_axis = numpy.abs(direction).argmax()
        if max_axis==0:
            v = numpy.array([0.,1.,0.])
        else:
            v = numpy.array([1.,0.,0.])
        d1 = numpy.cross(direction, v)
        d1 = normaliseVector(d1)
        d2 = numpy.cross(direction, d1)
        d2 = normaliseVector(d2)
        return (d1, d2)
    
    @on_trait_change("focus, direction, number, theta, working_dist, max_ray_len, rings")
    def on_update(self):
        self.data_source.modified()
        self.update=True
    
    @cached_property
    def _get_InputRays(self):
        focus = numpy.array(self.focus)
        direction = numpy.array(self.direction)
        working_dist = self.working_dist
        origin = focus - direction*working_dist
        count = self.number
        rings = self.rings
        theta = self.theta
        radius = numpy.tan(numpy.pi * theta/180) * working_dist
        
        d1, d2 = self.principle_axes
        
        radii = [(i+1)*radius/rings for i in xrange(rings)]
        angles = [i*2*numpy.pi/count for i in xrange(count)]
        offsets = [r*(d1*numpy.sin(a) + d2*numpy.cos(a)) 
                   for a in angles
                   for r in radii]
        
        origins = [origin] + [origin + offset for offset in offsets]
        directions = normaliseVector([direction] + 
                                     [direction*working_dist - offset 
                                      for offset in offsets])
        rays = RayCollection(len(origins))
        for i in xrange(len(origins)):
            ray = Ray(origin=origins[i], direction=directions[i],
                        length=numpy.Inf)
            rays.add_ray(ray)

#        rays.set_polarisation(1, 0, 0)
#        size = origins.shape[0],1
#        rays.E1_amp = numpy.ones(size, dtype=numpy.complex128)
#        rays.E2_amp = numpy.zeros(size, dtype=numpy.complex128)
#        rays.refractive_index = numpy.ones(size, dtype=numpy.complex128)
#        rays.offset_length = numpy.sqrt(((origins - focus)**2).sum(axis=-1)).reshape(-1,1)
#        rays.normals = numpy.zeros_like(origins)
#        rays.normals[:,1]=1.0
        return rays
    
    @cached_property
    def _get_InputDetailRays(self):
        focus = numpy.array(self.focus)
        direction = numpy.array(self.direction)
        working_dist = self.working_dist
        origin = focus - direction*working_dist
        count = self.number
        theta = self.theta
        radius = numpy.tan(numpy.pi * theta/180)
        
        d1, d2 = self.principle_axes
        
        det_res = self.detail_resolution
        slen = det_res * 1j
        newaxis = numpy.newaxis
        
        d1o, d2o = numpy.ogrid[-radius:radius:slen,-radius:radius:slen]
        
        index = numpy.arange(det_res*det_res).reshape(det_res,det_res)
        
        cells = numpy.dstack([index[:-1,:-1],
                                    index[1:,:-1],
                                    index[1:,1:],
                                    index[:-1,1:]]).reshape(-1,4)
        
        offsets = (d1o[:,:,newaxis]*d1[newaxis,newaxis,:] \
                    + d2o[:,:,newaxis]*d2[newaxis,newaxis,:]).reshape(-1,3)
                    
        in_rad = (offsets**2).sum(axis=-1) < radius**2
        
        #cull points outside the radius
        offsets = offsets[in_rad,:]
        
        #now cull cells with points outside the radius
        map = numpy.ones(in_rad.shape, numpy.int) * -1
        map[in_rad] = numpy.arange(in_rad.sum())
        
        mask = in_rad[cells].all(axis=-1)
        cells = map[cells[mask,:]]

        origins = origin[newaxis,:] + offsets*working_dist
        directions = normaliseVector(direction[newaxis,:] - offsets)
        
        rays = RayCollection(origin=origins, direction=directions,
                             max_length=self.max_ray_len)
        rays.set_polarisation(1, 0, 0)
        size = origins.shape[0],1
        rays.E1_amp = numpy.ones(size, dtype=numpy.complex128)
        rays.E2_amp = numpy.zeros(size, dtype=numpy.complex128)
        rays.refractive_index = numpy.ones(size, dtype=numpy.complex128)
        rays.cells = cells
        rays.offset_length = numpy.sqrt(((origins - focus)**2).sum(axis=-1)).reshape(-1,1)
        #print cells.max(), rays.number, "check"
        return rays