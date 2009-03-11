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
from enthought.traits.api import HasTraits, Int, Float, Range,\
     Bool, Property, Array, Event, List, cached_property, Str,\
     Instance, Tuple, on_trait_change, Trait

from enthought.traits.ui.api import View, Item

from enthought.tvtk.api import tvtk


Vector = Array(shape=(3,))

def normaliseVector(a):
    a= numpy.asarray(a)
    mag = numpy.sqrt((a**2).sum(axis=-1)).reshape(-1,1)
    return a/mag

def normaliseVector1d(a):
    a= numpy.asarray(a)
    mag = numpy.sqrt((a**2).sum())
    return a/mag


VectorArray = Array(shape=(None,3), dtype=numpy.double)

def collectRays(*rayList):
    """Combines a sequence of RayCollections into a single larger one"""
    origin = numpy.vstack([r.origin for r in rayList])
    direction = numpy.vstack([r.direction for r in rayList])
    length = numpy.vstack([r.length for r in rayList])
    optic = numpy.concatenate([r.optic for r in rayList])
    face_id = numpy.concatenate([r.face_id for r in rayList])
    refractive_index = numpy.vstack([r.refractive_index for r in rayList])
    E_vector = numpy.vstack([r.E_vector for r in rayList])
    E1_amp = numpy.vstack([r.E1_amp for r in rayList])
    E2_amp = numpy.vstack([r.E2_amp for r in rayList])
    parent_ids = numpy.concatenate([r.parent_ids for r in rayList])
    
    #parent = rayList[0].parent
    max_length = max(r.max_length for r in rayList)
    
    newRays = RayCollection(origin=origin, direction=direction,
                            length=length, optic=optic,
                            face_id=face_id, refractive_index=refractive_index,
                            E_vector=E_vector, E1_amp=E1_amp,
                            E2_amp = E2_amp, parent_ids=parent_ids,
                            max_length=max_length)
    return newRays


class RayCollection(HasTraits):
    origin = VectorArray
    direction = Property(VectorArray)
    length = Array(shape=(None,1), dtype=numpy.double)
    termination = Property(depends_on="origin, direction, length")
    
    optic = Array(shape=(None,), dtype=numpy.object)
    face_id = Array(shape=(None,), dtype=numpy.uint32)
    
    refractive_index = Array(shape=(None,1), dtype=numpy.complex128)
    
    E_vector = VectorArray
    E1_amp = Array(shape=(None,1), dtype=numpy.complex128)
    E2_amp = Array(shape=(None,1), dtype=numpy.complex128)
    
    parent = Instance("RayCollection")
    parent_ids = Array(shape=(None,), dtype=numpy.uint32)
    
    max_length = Float
    
    def set_polarisation(self, Ex, Ey, Ez):
        E = numpy.array([Ex, Ey, Ez])
        
        orth = numpy.cross(E, self.direction)
        para = numpy.cross(orth, self.direction)
        para /= numpy.sqrt((para**2).sum(axis=1)).reshape(-1,1)
        
        self.E_vector = para
    
    def _length_default(self):
        return numpy.ones((self.origin.shape[0],1)) * numpy.Infinity
    
    def _get_termination(self):
        length = self.length.clip(0, self.max_length).reshape(-1,1)
        try:
            return self.origin + length*self.direction
        except:
            pass
    
    def _get_direction(self):
        return self._direction
    
    def _set_direction(self, d):
        self._direction = normaliseVector(d)


class BaseRaySource(HasTraits):
    name = Str("Ray Source")
    update = Event()
    
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
    
    #idx selector for the input ways which should be visualised
    view_ray_ids = Trait(None, Array(dtype=numpy.int))
    
    tube = Instance(tvtk.TubeFilter, (), transient=True)
    sphere = Instance(tvtk.SphereSource, (), transient=True)
    data_source = Instance(tvtk.ProgrammableSource, transient=True)
    
    ray_actor = Instance(tvtk.Actor, transient=True)
    start_actor = Instance(tvtk.Actor, transient=True)
    #field_actor = Instance(tvtk.Actor)
    
    actors = Property(List, transient=True)
    
    vtkproperty = Instance(tvtk.Property, (), {'color':(1,0.5,0)}, transient=True)
    
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
    
    def _TracedRays_changed(self):
        self.data_source.modified()
        
    def _get_InputRays(self):
        return None
    
    def _InputRays_changed(self):
        self.update=True
        
    def _scale_factor_changed(self, scale):
        self.tube.radius = scale
        self.update = True
    
    def _get_actors(self):
        actors = [self.ray_actor, self.start_actor]
        prop = self.vtkproperty
        if prop:
            for actor in actors:
                actor.property = prop
        return actors
    
    def _data_source_default(self):
        source = tvtk.ProgrammableSource()
        def execute():
            output = source.poly_data_output
            
            pointArrayList = []
            for rays in self.TracedRays:
                try:
                    id_mask = numpy.setmember1d(rays.parent_ids, ids)
                except (NameError, AttributeError):
                    id_mask = numpy.zeros(rays.length.shape[0], dtype=numpy.bool)
                    view_ray_ids = self.view_ray_ids
                    if view_ray_ids is not None:
                        id_mask[self.view_ray_ids] = True
                    else:
                        id_mask[:] = True
                
                start_pos = rays.origin[id_mask]
                end_pos = rays.termination[id_mask]
                ids = numpy.argwhere(id_mask)[:,0]
                interleaved = numpy.array([start_pos, end_pos]).swapaxes(0,1).reshape(-1,3)
                pointArrayList.append(interleaved)
                
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
        
        map = tvtk.PolyDataMapper(input=tube.output)
        act = tvtk.Actor(mapper=map)
        return act
    
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
    
    traits_view = View(Item('name'),
                       Item('focus'),
                       Item('direction'),
                       Item('number'),
                       Item('radius'),
                       Item('show_start'),
                       Item('scale_factor'),
                       )
    
    @on_trait_change("focus, direction, number, radius, max_ray_len")
    def on_update(self):
        self.data_source.modified()
        self.update=True
    
    @cached_property
    def _get_InputRays(self):
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
        d1 = normaliseVector1d(d1)
        d2 = numpy.cross(direction, d1)
        d2 = normaliseVector1d(d2)
        angles = (i*2*numpy.pi/count for i in xrange(count))
        offsets = [radius*(d1*numpy.sin(a) + d2*numpy.cos(a)) for a in angles]
        
        origins = numpy.array([origin + offset for offset in offsets])
        directions = numpy.ones_like(origins) * direction
        rays = RayCollection(origin=origins, direction=directions,
                             max_length=self.max_ray_len)
        rays.set_polarisation(1, 0, 0)
        size = origins.shape[0],1
        rays.E1_amp = numpy.ones(size, dtype=numpy.complex128)
        rays.E2_amp = numpy.zeros(size, dtype=numpy.complex128)
        rays.refractive_index = numpy.ones(size, dtype=numpy.complex128)
        return rays
    
    

class ConfocalRaySource(BaseRaySource):
    focus = Tuple((0.,0.,0.))
    direction = Tuple((0.,0.,1.))
    number = Int(20)
    theta = Range(0.0,90.0,value=30)
    working_dist = Float(100.0)
    
    principle_axes = Property(Tuple(Array,Array), depends_on="direction")
    
    #view_ray_ids = numpy.arange(20)
    
    InputRays = Property(Instance(RayCollection), 
                         depends_on="focus, direction, number, theta, working_dist, max_ray_len")
    
    traits_view = View(Item('name'),
                       Item('focus'),
                       Item('direction'),
                       Item('number'),
                       Item('theta', label="angle"),
                       Item('working_dist', label="working dist"),
                       Item('show_start'),
                       Item('scale_factor'),
                       )
    
    @cached_property
    def _get_principle_axes(self):
        direction = self.direction
        max_axis = numpy.abs(direction).argmax()
        if max_axis==0:
            v = numpy.array([0.,1.,0.])
        else:
            v = numpy.array([1.,0.,0.])
        d1 = numpy.cross(direction, v)
        d1 = normaliseVector1d(d1)
        d2 = numpy.cross(direction, d1)
        d2 = normaliseVector1d(d2)
        return (d1, d2)
    
    @on_trait_change("focus, direction, number, theta, working_dist, max_ray_len")
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
        theta = self.theta
        radius = numpy.tan(numpy.pi * theta/180) * working_dist
        
        d1, d2 = self.principle_axes
        
        angles = (i*2*numpy.pi/count for i in xrange(count))
        offsets = [radius*(d1*numpy.sin(a) + d2*numpy.cos(a)) for a in angles]
        
        origins = numpy.array([origin + offset for offset in offsets])
        directions = normaliseVector([direction*working_dist - offset 
                                      for offset in offsets])
        rays = RayCollection(origin=origins, direction=directions,
                             max_length=self.max_ray_len)
        rays.set_polarisation(1, 0, 0)
        size = origins.shape[0],1
        rays.E1_amp = numpy.ones(size, dtype=numpy.complex128)
        rays.E2_amp = numpy.zeros(size, dtype=numpy.complex128)
        rays.refractive_index = numpy.ones(size, dtype=numpy.complex128)
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
        
        slen = self.detail_resolution * 1j
        newaxis = numpy.newaxis
        
        d1o, d2o = numpy.ogrid[-radius:radius:slen,-radius:radius:slen]
        
        offsets = (d1o[:,:,newaxis]*d1[newaxis,newaxis,:] \
                    + d2o[:,:,newaxis]*d2[newaxis,newaxis,:]).reshape(-1,3)
                    
        in_rad = (offsets**2).sum(axis=-1) < radius**2
        
        offsets = offsets[in_rad,:]

        origins = origin[newaxis,:] + offsets*working_dist
        directions = normaliseVector(direction[newaxis,:] - origins)
        
        rays = RayCollection(origin=origins, direction=directions,
                             max_length=self.max_ray_len)
        rays.set_polarisation(1, 0, 0)
        size = origins.shape[0],1
        rays.E1_amp = numpy.ones(size, dtype=numpy.complex128)
        rays.E2_amp = numpy.zeros(size, dtype=numpy.complex128)
        rays.refractive_index = numpy.ones(size, dtype=numpy.complex128)
        return rays