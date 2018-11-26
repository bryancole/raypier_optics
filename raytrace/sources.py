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

from traits.api import HasTraits, Int, Float, \
     Bool, Property, Array, Event, List, cached_property, Str,\
     Instance, on_trait_change, Trait, Enum, Title, Complex

from traitsui.api import View, Item, Tabbed, VGroup, Include, \
    Group

from tvtk.api import tvtk

from raytrace.ctracer import RayCollection, Ray, ray_dtype
from raytrace.utils import normaliseVector, Range, TupleVector, Tuple, \
            UnitTupleVector, UnitVectorTrait
from raytrace.bases import Renderable, RaytraceObject, NumEditor

Vector = Array(shape=(3,))


class BaseBase(HasTraits, RaytraceObject):
    pass


class BaseRaySource(BaseBase):
    abstract=True
    subclasses = set()
    name = Title("Ray Source", allow_selection=False)
    update = Event()
    display = Enum("pipes", "wires", "hidden")
    render = Event()
    mapper = Instance(tvtk.PolyDataMapper, (), transient=True)
    
    wavelength_list = List() #List of wavelengths (floats) given in microns
    
    InputRays = Property(Instance(RayCollection), depends_on="max_ray_len")
    TracedRays = List(RayCollection, transient=True)
    
    InputDetailRays = Property(Instance(RayCollection), depends_on="InputRays")
    TracedDetailRays = List(RayCollection, transient=True)

    detail_resolution = Int(32)

    max_ray_len = Float(200.0)
    
    scale_factor = Float(0.2)
    
    show_start = Bool(True)
    show_polarisation = Bool(False)
    show_direction = Bool(False)
    show_normals = Bool(False)
    
    #idx selector for the input ways which should be visualised
    #making this transient because pyYAML fails to serialise arrays
    view_ray_ids = Trait(None, Array(dtype=numpy.int), transient=True)
    
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


    def get_ray_list_by_id(self):
        """return a list of lists of dictionaries where the index of the outer list is the ray id (defined by order encountered)
        , and the inner index is the recursion # backwards (0 being the ray that dies, 1 is its parent, and so on) 
        and the dictionary keys are the attributes of the ray object"""
        result = []
        #keys is copy and pasted from the dtype of a ray object
        keys = ['origin','direction','normal','E_vector','refractive_index','E1_amp','E2_amp','length','wavelength','parent_idx','end_face_idx']
        #start at the last ray collection and work backwards
        temp = list(self.TracedRays)    #otherwise its a TraitListObject
        raycollections = []             #make data a list of lists of lists
        for raycollection in temp:
            tmp = []
            for x in raycollection.copy_as_array():
                tmp.append(list(x))
            raycollections.append(tmp)

        length = len(raycollections)-1

        for i,collection in enumerate(reversed(raycollections)):
            for ray in collection:
                if ray:                  #rays already accounted for are set to none
                    lst = []             #list to contain all relevant dictionaries for each ray
                    r = {}               #the first dictionary
                    for x,att in enumerate(ray):
                        r[keys[x]] = att 
                    lst.append(r)
                    cntr = length-i
                    cntr -= 1
                    parent_idx = r['parent_idx']
                    while cntr >= 0:     #iterate back through the ray collections
                        r = {}          #and repeat what we did for the first ray
                        arr = raycollections[cntr][parent_idx]
                        for x,att in enumerate(arr):
                            r[keys[x]] = att 
                        lst.append(r.copy())
                        raycollections[cntr][parent_idx] = None     #mark this ray done
                        cntr -= 1
                        parent_idx = r['parent_idx']
                    result.append(lst)

        
        return result

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
        glyph_filter = tvtk.Glyph3D(input_connection = source.output_port,
                                    scale_factor=10.0,
                                    vector_mode='use_normal'
                                    )
        glyph_filter.set_source_connection(glyph.output_port)
        map = tvtk.PolyDataMapper(input_connection=glyph_filter.output_port)
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
        tube.input_connection=self.data_source.output_port
        tube.number_of_sides = 20
        
        map = self.mapper
        act = tvtk.Actor(mapper=map)
        act.visibility = True
        display = self.display
        if display=="pipes":
            map.input_connection=tube.output_port
        elif display=="wires":
            map.input_connection = self.data_source.output_port
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
            self.mapper.input_connection = self.tube.output_port
        elif vnew=="wires":
            self.mapper.input_connection = self.data_source.output_port
        elif vnew=="hidden":
            for act in actors:
                act.visibility = False
        self.render = True
    
    def _start_actor_default(self):
        map = tvtk.PolyDataMapper(input_connection=self.sphere.output_port)
        act = tvtk.Actor(mapper=map)
        return act
    
    
class SingleRaySource(BaseRaySource):
    abstract = False
    origin = Tuple((0.,0.,0.))
    direction = UnitTupleVector
    
    E_vector = UnitVectorTrait((1.,0.,0.), editor_traits={'cols':3,
                                'labels':['x','y','z']})
    E1_amp = Complex(1.0+0.0j)
    E2_amp = Complex(0.0+0.0j)
    wavelength = Float(0.78)
    
    view_ray_ids = numpy.arange(1)
    
    InputRays = Property(Instance(RayCollection), 
                         depends_on="origin, direction, max_ray_len, E_vector, E1_amp, E2_amp")
    
    geom_grp = VGroup(Group(Item('origin', show_label=False,resizable=True), 
                            show_border=True,
                            label="Origin position",
                            padding=0),
                       Group(Item('direction', show_label=False, resizable=True),
                            show_border=True,
                            label="Direction"),
                       label="Geometry")
    
    @on_trait_change("wavelength")
    def _do_wavelength_changed(self):
        self.wavelength_list = [self.wavelength]
    
    @on_trait_change("direction, max_ray_len")
    def on_update(self):
        self.data_source.modified()
        self.update=True
    
    @cached_property
    def _get_InputRays(self):
        origin = numpy.array(self.origin)
        direction = numpy.array(self.direction)
        max_axis = numpy.abs(direction).argmax()
        if max_axis==0:
            v = numpy.array([0.,1.,0.])
        else:
            v = numpy.array([1.,0.,0.])
        d1 = numpy.cross(direction, v)
        d1 = normaliseVector(d1)
        d2 = numpy.cross(direction, d1)
        d2 = normaliseVector(d2)

        E_vector = numpy.cross(direction, self.E_vector)
        E_vector = numpy.cross(E_vector, direction)

        rays = RayCollection(1)
        ray = Ray()
        ray.origin = origin
        ray.direction = direction
        ray.wavelength_idx = 0
        ray.E_vector = E_vector
        ray.E1_amp = self.E1_amp
        ray.E2_amp = self.E2_amp
        ray.refractive_index = 1.0
        ray.normals = (0.,1.,0.) #What is this for?
        ray.parent_idx = -1
        ray.end_face_idx = -1
        
        rays.add_ray(ray)
        return rays
    
    
class BroadbandRaySource(SingleRaySource):
    """Creates a number of coincident rays covering a range of wavelengths.
    The wavelength spacings may be uniform in wavelength, or uniform in frequency.
    Wavelength is given in microns and is the in-vacuo wavelength.
    """
    number = Int(200, auto_set=False, enter_set=True)
    wavelength_start = Float(1.50, editor=NumEditor)
    wavelength_end = Float(1.60, editor=NumEditor)
    uniform_deltaf = Bool(True)
    
    view_ray_ids = numpy.arange(200)
    
    InputRays = Property(Instance(RayCollection), 
                         depends_on="origin, direction, number, wavelength_start, "
                         "wavelength_end, uniform_deltaf, max_ray_len, E_vector")
    
    geom_grp = VGroup(Group(Item('origin', show_label=False,resizable=True), 
                            show_border=True,
                            label="Origin position",
                            padding=0),
                       Group(Item('direction', show_label=False, resizable=True),
                            show_border=True,
                            label="Direction"),
                       Item('number'),
                       Item('wavelength_start'),
                       Item('wavelength_end'),
                       Item("uniform_deltaf"),
                       label="Geometry")
    
    @on_trait_change("direction, number, wavelength_start, wavelength_end, uniform_deltaf, max_ray_len")
    def on_update(self):
        self.data_source.modified()
        self.update=True
        
    @on_trait_change("number, wavelength_start, wavelength_end, uniform_deltaf")
    def _do_wavelength_changed(self):
        if self.uniform_deltaf:
            f_start = 1./self.wavelength_start
            f_end = 1./self.wavelength_end
            wavelengths = 1./numpy.linspace(f_start, f_end, self.number)
        else:
            wavelengths = numpy.linspace(self.wavelength_start, self.wavelength_end,
                                         self.number)
        self.wavelength_list = list(wavelengths)
        
    @cached_property
    def _get_InputRays(self):
        origin = numpy.array(self.origin)
        direction = numpy.array(self.direction)
        count = self.number

        E_vector = numpy.cross(self.E_vector, direction)
        E_vector = numpy.cross(E_vector, direction)
        
        wavelen_idx = numpy.arange(count)

        ray_data = numpy.zeros(count, dtype=ray_dtype)            
        ray_data['origin'] = origin.reshape(-1,3)
        ray_data['direction'] = direction.reshape(-1,3)
        ray_data['wavelength_idx'] = wavelen_idx
        ray_data['E_vector'] = [normaliseVector(E_vector)]
        ray_data['E1_amp'] = self.E1_amp
        ray_data['E2_amp'] = self.E2_amp
        ray_data['refractive_index'] = 1.0+0.0j
        ray_data['normal'] = [[0,1,0]]
        rays = RayCollection.from_array(ray_data)
        return rays
    
    
class ParallelRaySource(SingleRaySource):
    
    number = Int(20, auto_set=False, enter_set=True)
    radius = Float(10.,editor=NumEditor)
    rings = Range(0,50,3, editor_traits={'mode':'spinner'})
    
    view_ray_ids = numpy.arange(20)
    
    InputRays = Property(Instance(RayCollection), 
                         depends_on="origin, direction, number, rings, radius, max_ray_len, E_vector")
    
    geom_grp = VGroup(Group(Item('origin', show_label=False,resizable=True), 
                            show_border=True,
                            label="Origin position",
                            padding=0),
                       Group(Item('direction', show_label=False, resizable=True),
                            show_border=True,
                            label="Direction"),
                       Item('number'),
                       Item('rings'),
                       Item('radius'),
                       label="Geometry")
    
    
    @on_trait_change("direction, number, rings, radius, max_ray_len")
    def on_update(self):
        self.data_source.modified()
        self.update=True
    
    @cached_property
    def _get_InputRays(self):
        origin = numpy.array(self.origin)
        direction = numpy.array(self.direction)
        count = self.number
        rings = self.rings
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

        E_vector = numpy.cross(self.E_vector, direction)
        E_vector = numpy.cross(E_vector, direction)

        ray_data = numpy.zeros((rings*count)+1, dtype=ray_dtype)
        if rings:
            radii = ((numpy.arange(rings)+1)*(radius/rings))[:,None,None]
            angles = (numpy.arange(count)*(2*numpy.pi/count))[None,:,None]
            offsets = radii*(d1*numpy.sin(angles) + d2*numpy.cos(angles)) 
            offsets.shape = (-1,3)
            
            ray_data['origin'][1:] = offsets
        ray_data['origin'] += origin
        ray_data['direction'] = direction
        ray_data['wavelength_idx'] = 0
        ray_data['E_vector'] = [normaliseVector(E_vector)]
        ray_data['E1_amp'] = self.E1_amp
        ray_data['E2_amp'] = self.E2_amp
        ray_data['refractive_index'] = 1.0+0.0j
        ray_data['normal'] = [[0,1,0]]
        rays = RayCollection.from_array(ray_data)
        return rays
    
    
class RectRaySource(BaseRaySource):
    """ rays from a rectangular aperture """ 
    origin = Tuple((0.,0.,0.))
    abstract = False
    direction = UnitTupleVector
    number = Int(20,auto_set=False,enter_set=True)		#total rays is n^2
    length = Float(10.,editor=NumEditor)
    width = Float(10,editor=NumEditor)
    rings = Range(1,50,3, editor_traits={'mode':'spinner'})
    theta = Float(0.,editor=NumEditor)
    randomness = Bool(False)
    
    view_ray_ids = numpy.arange(20)
    
    InputRays = Property(Instance(RayCollection), 
                         depends_on="origin, direction, number, length, width, theta, max_ray_len")
    
    geom_grp = VGroup(Group(Item('origin', show_label=False,resizable=True), 
                            show_border=True,
                            label="Origin position",
                            padding=0),
                       Group(Item('direction', show_label=False, resizable=True),
                            show_border=True,
                            label="Direction"),
                       Item('number'),
                       Item('length'),
                       Item('width'),
                       Item('theta'),
		       Item('randomness'),
                       label="Geometry")
    
    
    @on_trait_change("focus, direction, number, length, width, theta, max_ray_len")
    def on_update(self):
        self.data_source.modified()
        self.update=True
    
    @cached_property
    def _get_InputRays(self):
        from utils import z_rotation
        origin = numpy.array(self.origin)
        direction = numpy.array(self.direction)
        count = self.number
        length = self.length
        width = self.width
        randomness = self.randomness
        max_axis = numpy.abs(direction).argmax()
        if max_axis==0:
            v = numpy.array([0.,1.,0.])
        else:
            v = numpy.array([1.,0.,0.])
        d1 = numpy.cross(direction, v)
        d1 = normaliseVector(d1)
        d2 = numpy.cross(direction, d1)
        d2 = normaliseVector(d2)
        a = length/2
        b = width/2
        X_range = numpy.linspace(-a, a, count)
        Y_range = numpy.linspace(-b, b, count)

        ray_data = numpy.zeros(count**2, dtype=ray_dtype)
        offsets = numpy.zeros([X_range.size*Y_range.size,3])
        
        if randomness:
            from random import uniform
        
        for i,x in enumerate(X_range):
            for j,y in enumerate(Y_range):
                if randomness:
                    x = x + uniform(-a/count, a/count)
                    y = y + uniform(-b/count, b/count)
                point = d1 * x + d2 * y
                block = i * Y_range.size
                offsets[block + j] = point
        
        origins = numpy.array([origin + offset for offset in offsets])
        
        dirmatrix = numpy.matrix(direction)
        #print "dirmatrix: ",z_rotation(numpy.radians(self.theta))
        raydir = z_rotation(numpy.radians(self.theta))*(dirmatrix.T)
        directions = numpy.ones_like(origins) * numpy.array(raydir.T)

        ray_data['origin'] = origins
        ray_data['direction'] = directions
        ray_data['E_vector'] = [[1,0,0]]
        ray_data['E1_amp'] = 1.0 + 0.0j
        ray_data['E2_amp'] = 0.0
        ray_data['refractive_index'] = 1.0+0.0j
        ray_data['normal'] = [[0,1,0]]
        rays = RayCollection.from_array(ray_data)
        return rays
    

class ConfocalRaySource(SingleRaySource):
    abstract = False
    focus = TupleVector
    direction = UnitTupleVector
    number = Range(1,50,20, editor_traits={'mode':'spinner'})
    theta = Range(0.0,90.0,value=30)
    working_dist = Float(100.0,editor=NumEditor)
    rings = Range(0,50,3, editor_traits={'mode':'spinner'})
    reverse = 1    #-1 is converging, +1 is diverging
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
        rev = self.reverse #to determine converging or diverging  
        focus = numpy.array(self.focus)
        direction = numpy.array(self.direction)
        working_dist = self.working_dist
        origin = focus - rev*direction*working_dist
        count = self.number
        rings = self.rings
        theta = self.theta
        radius = numpy.tan(numpy.pi * theta/180) * working_dist
        
        d1, d2 = self.principle_axes
        ray_data = numpy.zeros((rings*count)+1, dtype=ray_dtype)
        
        radii = ((numpy.arange(rings)+1)*(radius/rings))[:,None,None]
        angles = (numpy.arange(count)*(2*numpy.pi/count))[None,:,None]
        offsets = radii*(d1*numpy.sin(angles) + d2*numpy.cos(angles)) 
        offsets.shape = (-1,3)

        ray_data['origin'][1:] = offsets
        ray_data['origin'] += origin
        ray_data['direction'][0] = normaliseVector(rev*direction)
        ray_data['direction'][1:] = normaliseVector((rev*direction*working_dist) - offsets)
        ray_data['length'] = self.max_ray_len
        ray_data['E_vector'] = self.E_vector
        ray_data['E1_amp'] = self.E1_amp
        ray_data['E2_amp'] = self.E2_amp
        ray_data['wavelength_idx'] = 0
        ray_data['refractive_index'] = 1.0
        ray_data['normal'] = [0,1,0]
        rays = RayCollection.from_array(ray_data)
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

class AdHocSource(BaseRaySource):
    '''create a source by specifying the input rays yourself''' 
    
    InputRays = Property(Instance(RayCollection), 
                         depends_on="origin, direction, number, rings, radius, max_ray_len")

    def __init__(self,*args,**kwargs):
        self.rays = kwargs.pop('rays')
        super(AdHocSource,self).__init__(*args,**kwargs)

    @cached_property
    def _get_InputRays(self):
        #i actually don't understand traits well at all.  this whole class is a hack.
        return self.rays

