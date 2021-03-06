#    Copyright 2009, Teraview Ltd., Bryan Cole
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


import numpy
import time
import itertools
import traceback

from traits.api import HasTraits, Int, Float, \
     Bool, Property, Array, Event, List, cached_property, Str,\
     Instance, on_trait_change, Trait, Enum, Title, Complex, Dict

from traitsui.api import View, Item, Tabbed, VGroup, Include, \
    Group

from tvtk.api import tvtk

from raypier.core.ctracer import RayCollection, GaussletCollection, Ray, ray_dtype, GAUSSLET_, PARABASAL_
from raypier.utils import normaliseVector, Range, TupleVector, Tuple, \
            UnitTupleVector, UnitVectorTrait
from raypier.bases import RaypierObject, NumEditor, BaseRayCollection

Vector = Array(shape=(3,))


class BaseBase(HasTraits, RaypierObject):
    pass


class BaseRaySource(BaseBase):
    abstract=True
    subclasses = set()
    name = Title("Ray Source", allow_selection=False)
    
    ### Updated to the current timestamp whenever the traced rays are changed (i.e. after a tracing operation)
    _mtime = Float(0.0)
    
    ### The update event trigger a new tracing operation
    update = Event()
    display = Enum("pipes", "wires", "hidden")
    
    ### The render event signals that the VTK view needs to be re-rendered
    render = Event()
    mapper = Instance(tvtk.PolyDataMapper, (), transient=True)
    
    wavelength_list = List() #List of wavelengths (floats) given in microns
    
    input_rays = Property(Instance(BaseRayCollection), depends_on="max_ray_len")
    traced_rays = List(BaseRayCollection, transient=True)
    
    InputDetailRays = Property(Instance(BaseRayCollection), depends_on="input_rays")
    TracedDetailRays = List(BaseRayCollection, transient=True)

    detail_resolution = Int(32)

    max_ray_len = Float(200.0)
    
    scale_factor = Float(0.1)
    
    show_start = Bool(True)
    show_polarisation = Bool(False)
    show_direction = Bool(False)
    show_normals = Bool(False)
    opacity = Range(0.0,1.0,1.0)
    
    ### A dict which maps each RayCollection instance in the traced_rays
    ### list to a list/array of bools of the same length as the RayCollection
    ### This bool array indicates if a ray should be hidden in the display
    ray_mask = Dict()
    
    #idx selector for the input ways which should be visualised
    #making this transient because pyYAML fails to serialise arrays
    view_ray_ids = Trait(None, Array(dtype=numpy.int), transient=True)
    
    tube = Instance(tvtk.TubeFilter, (), transient=True)
    sphere = Instance(tvtk.SphereSource, (), transient=True)
    data_source = Instance(tvtk.ProgrammableSource, transient=True)
    
    normals_source = Instance(tvtk.ProgrammableSource, transient=True)
    normals_actor = Instance(tvtk.Actor, transient=True)
    normals_glyph = Instance(tvtk.Glyph3D, transient=True)
    
    ray_actor = Instance(tvtk.Actor, transient=True)
    start_actor = Instance(tvtk.Actor, transient=True)
    #field_actor = Instance(tvtk.Actor)
    
    export_pipes = Bool(False, desc="if set, rays are exported to STEP files as pipes")
    
    actors = Property(List, transient=True)
    
    vtkproperty = Instance(tvtk.Property, (), {'color':(1,0.5,0), 'opacity': 1.0}, transient=True)
    
    display_grp = VGroup(Item('display'),
                       Item('show_start'),
                       Item('show_normals'),
                       Item('scale_factor'),
                       Item('export_pipes',label="export as pipes"),
                       Item('opacity', editor=NumEditor),
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
        for rays in self.traced_rays[1:]:
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
            if parent == self.traced_rays[0]:
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
        temp = list(self.traced_rays)    #otherwise its a TraitListObject
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
        
        :param idx: the index of the RayCollection in the TracedRay list 
                    to analyse
                    
        :returns: the average angle from the mean ray, in degrees
        
        """
        rays = self.traced_rays[idx]
        ave_dir = normaliseVector(rays.direction.mean(axis=0))
        dotprod = (ave_dir[numpy.newaxis,:] * rays.direction).sum(axis=1)
        angles = numpy.arccos(dotpod)*180 / numpy.pi
        return angles.mean()
    
    def make_step_shape(self):
        if self.export_pipes:
            from raypier.step_export import make_rays_both as make_rays
        else:
            from raypier.step_export import make_rays_wires as make_rays
        return make_rays(self.traced_rays, self.scale_factor), "red"
    
    def _traced_rays_changed(self):
        self.data_source.modified()
        self.normals_source.modified()
        self._mtime = time.monotonic()
        print("SET MTIME")
        
    def _get_input_rays(self):
        return None
    
    def _input_rays_changed(self):
        self.update=True
        
    def _scale_factor_changed(self, scale):
        self.tube.radius = scale
        self.normals_glyph.scale_factor=self.scale_factor*10
        self.sphere.radius = self.scale_factor*3
        self.render = True
        
    def _show_normals_changed(self):
        self.normals_source.modified()
        self.render=True
        
    def _ray_mask_changed(self):
        self.update = True
        
    def _opacity_changed(self):
        self.vtkproperty.opacity = self.opacity
        self.render = True
    
    def _get_actors(self):
        actors = [self.ray_actor, self.start_actor, self.normals_actor]
        return actors
        #return actors[:2]  #makeshift turning off of normal glyphs
        
    def _normals_glyph_default(self):
        glyph =  tvtk.Glyph3D(scale_factor=self.scale_factor*10,
                                    vector_mode='use_normal'
                                    )
        return glyph
    
    def _normals_source_default(self):
        source = tvtk.ProgrammableSource()
        def execute():
            if not self.show_normals:
                return
            output = source.poly_data_output
            points = numpy.vstack([p.origin for p in self.traced_rays[1:]])
            normals = numpy.vstack([p.normal for p in self.traced_rays[1:]])
            output.points = points
            output.point_data.normals = normals
            #print("calc normals GLYPH")
        source.set_execute_method(execute)
        return source
        
    def _normals_actor_default(self):
        source = self.normals_source
        glyph = tvtk.ArrowSource(tip_radius=0.3, shaft_radius=0.1)
        glyph_filter = self.normals_glyph
        glyph_filter.set_input_connection(source.output_port)
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
            mask = self.ray_mask
            gausslet = isinstance(self.input_rays, GaussletCollection)
            for rays in self.traced_rays:
                vis = mask.get(rays, [])
                if len(vis) != len(rays):
                    vis = itertools.repeat(True)
                if gausslet:
                    start_pos = [r.base_ray.origin for r,v in zip(rays, vis) if v]
                    end_pos = [r.base_ray.termination for r,v in zip(rays, vis) if v]
                else:
                    start_pos = [r.origin for r,v in zip(rays, vis) if v]
                    end_pos = [r.termination for r,v in zip(rays, vis) if v]
                #print "start", start_pos
                #print "end", end_pos
                interleaved = numpy.array([start_pos, end_pos]).swapaxes(0,1).reshape(-1,3)
                pointArrayList.append(interleaved)
            if pointArrayList:
                points = numpy.vstack(pointArrayList)
                output.points=points
                lines = list([i,i+1] for i in range(0,points.shape[0],2))
                output.lines = lines
        source.set_execute_method(execute)
        return source
    
    def _ray_actor_default(self):
        tube = self.tube
        tube.radius = self.scale_factor
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
        prop.opacity = self.opacity
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
    wavelength = Float(0.78) # in microns
    
    view_ray_ids = numpy.arange(1)
    
    input_rays = Property(Instance(RayCollection), 
                         depends_on="origin, direction, max_ray_len, E_vector, E1_amp, E2_amp")
    
    geom_grp = VGroup(Group(Item('origin', show_label=False,resizable=True), 
                            show_border=True,
                            label="Origin position",
                            padding=0),
                       Group(Item('direction', show_label=False, resizable=True),
                            show_border=True,
                            label="Direction"),
                       label="Geometry")
    
    def _wavelength_list_default(self):
        return [self.wavelength,]
    
    @on_trait_change("wavelength")
    def _do_wavelength_changed(self):
        self.wavelength_list = [self.wavelength]
    
    @on_trait_change("origin, direction, max_ray_len")
    def on_update(self):
        self.data_source.modified()
        self.update=True
    
    @cached_property
    def _get_input_rays(self):
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
    
    input_rays = Property(Instance(RayCollection), 
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
    def _get_input_rays(self):
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
    
    input_rays = Property(Instance(RayCollection), 
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
    def _get_input_rays(self):
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
    
    input_rays = Property(Instance(RayCollection), 
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
    def _get_input_rays(self):
        from .utils import z_rotation
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
    
    
class GaussianBeamRaySource(SingleRaySource):
    """
    Launches a set of skew rays simulating the propagation of a Gaussian beam.
    See Paul D Colbourne - Proc. SPIE 9293, International Optical Design Conference 2014, 92931S
    and refs theirin.
    """
    beam_waist = Float(100.0) #in microns
    number = Range(2,64,16)
    working_distance = Float(0.0) #the distance from the beam waist to the source location
    
    input_rays = Property(Instance(RayCollection), 
                         depends_on="origin, direction, number, beam_waist, wavelength, "
                         "max_ray_len, E_vector, working_distance")
    
    geom_grp = VGroup(Group(Item('origin', show_label=False,resizable=True), 
                            show_border=True,
                            label="Origin position",
                            padding=0),
                       Group(Item('direction', show_label=False, resizable=True),
                            show_border=True,
                            label="Direction"),
                       Item('number'),
                       Item("wavelength"),
                       Item('beam_waist', tooltip="1/e^2 beam intensity radius, in microns"),
                       Item('working_distance', tooltip="The distance from the source to the beam waist"),
                       label="Geometry")
    
    
    @on_trait_change("direction, number, beam_waist, max_ray_len, working_distance")
    def on_update(self):
        self.data_source.modified()
        self.update=True
    
    @cached_property
    def _get_input_rays(self):
        try:
            origin = numpy.array(self.origin)
            direction = numpy.array(self.direction)
            count = self.number
            radius = self.beam_waist/2000.0 #convert to mm
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
    
            ray_data = numpy.zeros((2*count)+1, dtype=ray_dtype)
            
            theta0 = self.wavelength/(numpy.pi*radius*1000.0)
            lr = numpy.array([1,-1])[:,None,None]
            eye = numpy.array([1,1])[:,None,None]
            alpha = (numpy.arange(count)*(2*numpy.pi/count))[None,:,None]
            origins = eye*radius*(d1*numpy.sin(alpha) + d2*numpy.cos(alpha))
            directions = theta0*lr*(d1*numpy.cos(alpha) - d2*numpy.sin(alpha))
            
            origins.shape = (-1,3)
            directions.shape = (-1,3)
            directions += direction
            directions = normaliseVector(directions)
            origins = origins - (self.working_distance*directions)
                
            ray_data['origin'][1:] = origins
            ray_data['origin'] += origin + (self.working_distance*direction)
            ray_data['direction'][1:] = directions
            ray_data['direction'][0] = direction
            ray_data['wavelength_idx'] = 0
            ray_data['E_vector'] = [normaliseVector(E_vector)]
            ray_data['E1_amp'] = self.E1_amp
            ray_data['E2_amp'] = self.E2_amp
            ray_data['refractive_index'] = 1.0+0.0j
            ray_data['normal'] = [[0,1,0]]
            rays = RayCollection.from_array(ray_data)
            rays.wavelengths = numpy.asarray(self.wavelength_list)
        except:
            import traceback
            traceback.print_exc()
            return RayCollection(5)
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
    
    input_rays = Property(Instance(RayCollection), 
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
    def _get_input_rays(self):
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
    
    input_rays = Property(Instance(RayCollection), 
                         depends_on="origin, direction, number, rings, radius, max_ray_len")

    def __init__(self,*args,**kwargs):
        self.rays = kwargs.pop('rays')
        super(AdHocSource,self).__init__(*args,**kwargs)

    @cached_property
    def _get_input_rays(self):
        #i actually don't understand traits well at all.  this whole class is a hack.
        return self.rays
    
    
class RayFieldSource(SingleRaySource):    
    show_mesh = Bool(True)
    mesh_source = Instance(tvtk.ProgrammableSource, transient=True)
    mesh_mapper = Instance(tvtk.PolyDataMapper, (), transient=True)
    mesh_actor = Instance(tvtk.Actor, transient=True)
    
    def _get_actors(self):
        actors = [self.ray_actor, self.start_actor, self.normals_actor, self.mesh_actor]
        return actors
        
    def _mesh_source_default(self):
        source = tvtk.ProgrammableSource()
        def execute():
            if not self.show_mesh:
                return
            output = source.poly_data_output
            rays = self.input_rays
            points = rays.origin
            neighbours = rays.neighbours
            fz = frozenset
            tris = set()
            for i in range(points.shape[0]):
                nb = neighbours[i]
                if any(a<=0 for a in nb):
                    continue
                nb1, nb2, nb3, nb4, nb5, nb6 = nb
                these = ((i, nb1, nb2),
                             (i, nb2, nb3),
                             (i, nb3, nb4),
                             (i, nb4, nb5),
                             (i, nb5, nb6),
                             (i, nb6, nb1))
                for tri in these:
                    tris.add( fz(tri) )
            output.points = points
            output.polys = [list(a) for a in tris]
            
            E1 = rays.E1_amp
            E2 = rays.E2_amp
            I = E1.real**2 + E1.imag**2 + E2.real**2 + E2.imag**2
            output.point_data.scalars = I
            self.mesh_mapper.scalar_range = (I.min(), I.max())

        source.set_execute_method(execute)
        return source
    
    def _mesh_actor_default(self):
        source = self.mesh_source
        map = self.mesh_mapper
        map.input_connection = source.output_port
        act = tvtk.Actor(mapper=map)
        #act.property.color = (0.0,1.0,0.0)
        act.property.representation = "surface"
        return act
    
    @on_trait_change("origin, direction")
    def _on_update_mesh(self):
        self.mesh_source.modified()
    

class HexagonalRayFieldSource(RayFieldSource):
    abstract = False
    radius = Float(10.,editor=NumEditor)
    resolution = Float(10.0, editor=NumEditor)
    gauss_width = Float(2.0, editor=NumEditor)
    
    input_rays = Property(Instance(RayCollection), 
                         depends_on="origin, direction, resolution, wavelength, gauss_width, radius, max_ray_len, E_vector")
    
    geom_grp = VGroup(Group(Item('origin', show_label=False,resizable=True), 
                            show_border=True,
                            label="Origin position",
                            padding=0),
                       Group(Item('direction', show_label=False, resizable=True),
                            show_border=True,
                            label="Direction"),
                       Item('radius'),
                       Item('resolution'),
                       Item('gauss_width'),
                       label="Geometry")

    @on_trait_change("direction, spacing, radius, max_ray_len")
    def on_update(self):
        self.data_source.modified()
        self.mesh_source.modified()
        self.update=True
    
    @cached_property
    def _get_input_rays(self):
        origin = numpy.array(self.origin)
        direction = numpy.array(self.direction)
        radius = self.radius
        spacing = radius/self.resolution
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
        
        cosV = numpy.cos(numpy.pi*30/180.)
        sinV = numpy.sin(numpy.pi*30/180.)
        d2 = cosV*d2 + sinV*d1
        
        nsteps = int(radius/(spacing*numpy.cos(numpy.pi*30/180.)))
        i = numpy.arange(-nsteps-1, nsteps+2)
        j = numpy.arange(-nsteps-1, nsteps+2)
        
        vi, vj = numpy.meshgrid(i,j)
        label = numpy.arange(vi.size).reshape(*vi.shape) #give each point a unique index
        neighbours = numpy.full((label.shape[0], label.shape[1], 6), -1, 'i')
        neighbours[1:,:,3] = label[:-1,:]
        neighbours[:-1,:,0] = label[1:,:]
        neighbours[:,1:,4] = label[:,:-1]
        neighbours[:,:-1,1] = label[:,1:]
        neighbours[1:,:-1,2] = label[:-1,1:]
        neighbours[:-1,1:,5] = label[1:,:-1]
        
        neighbours.shape = -1,6
        
        vi.shape = -1
        vj.shape = -1
        
        ri = ((vi + sinV*vj)**2 + (cosV*vj**2))
        select = ri < (radius/spacing)**2
        
        xi = vi[select]
        yj = vj[select]
        
        backmap = numpy.full(vi.shape[0]+1, -1, 'i')
        backmap[:-1][select] = numpy.arange(len(xi))
        selnb = backmap[neighbours]
        neighbours = selnb[select,:]
        
        offsets = xi[:,None]*d1 + yj[:,None]*d2
        offsets *= spacing
        
        gauss = numpy.exp(-(ri[select]*(spacing**2)/((self.gauss_width)**2))) 
        
        ray_data = numpy.zeros(offsets.shape[0], dtype=ray_dtype)
            
        ray_data['origin'] = offsets
        ray_data['origin'] += origin
        ray_data['direction'] = direction
        ray_data['wavelength_idx'] = 0
        ray_data['E_vector'] = [normaliseVector(E_vector)]
        ray_data['E1_amp'] = self.E1_amp * gauss
        ray_data['E2_amp'] = self.E2_amp * gauss
        ray_data['refractive_index'] = 1.0+0.0j
        ray_data['normal'] = [[0,1,0]]
        rays = RayCollection.from_array(ray_data)
        rays.neighbours = neighbours
        return rays
    
    
class ConfocalRayFieldSource(HexagonalRayFieldSource):
    abstract = False
    angle = Float(10.0)
    
    ### Number of angle steps along the radius
    resolution = Range(0.0, None, 10.0)
    
    working_dist = Float(100.0)
    
    numerical_aperture = Float(0.12)
    
    input_rays = Property(Instance(RayCollection), 
                         depends_on="origin, direction, angle, resolution, "
                                    "working_dist, max_ray_len, E_vector, wavelength, "
                                    "E1_amp, E2_amp, numerical_aperture")
    
    geom_grp = VGroup(Group(Item('origin', show_label=False,resizable=True), 
                            show_border=True,
                            label="Origin position",
                            padding=0),
                       Group(Item('direction', show_label=False, resizable=True),
                            show_border=True,
                            label="Direction"),
                       Item('wavelength', editor=NumEditor),
                       Item('numerical_aperture', editor=NumEditor),
                       Item('angle', editor=NumEditor),
                       Item('resolution', editor=NumEditor),
                       Item('working_dist', editor=NumEditor),
                       label="Geometry")

    @on_trait_change("angle, resolution, numerical_aperture, working_dist, max_ray_len")
    def on_update(self):
        self.data_source.modified()
        self.mesh_source.modified()
    
    @cached_property
    def _get_input_rays(self):
        try:
            origin = numpy.array(self.origin)
            direction = numpy.array(self.direction)
            
            radius = numpy.tan(self.angle*numpy.pi/180)
            spacing = radius /( self.resolution + 0.5 )
            
            max_axis = numpy.abs(direction).argmax()
            if max_axis==0:
                v = numpy.array([0.,1.,0.])
            else:
                v = numpy.array([1.,0.,0.])
            d1 = numpy.cross(direction, v)
            d1 = normaliseVector(d1)
            d2 = numpy.cross(direction, d1)
            d2 = normaliseVector(d2)
            
            cosV = numpy.cos(numpy.pi*30/180.)
            sinV = numpy.sin(numpy.pi*30/180.)
            d2 = cosV*d2 + sinV*d1
            
            nsteps = int(radius/(spacing*numpy.cos(numpy.pi*30/180.)))
            i = numpy.arange(-nsteps-1, nsteps+2)
            j = numpy.arange(-nsteps-1, nsteps+2)
            
            vi, vj = numpy.meshgrid(i,j)
            label = numpy.arange(vi.size).reshape(*vi.shape) #give each point a unique index
            neighbours = numpy.full((label.shape[0], label.shape[1], 6), -1, 'i')
            neighbours[1:,:,3] = label[:-1,:]
            neighbours[:-1,:,0] = label[1:,:]
            neighbours[:,1:,4] = label[:,:-1]
            neighbours[:,:-1,1] = label[:,1:]
            neighbours[1:,:-1,2] = label[:-1,1:]
            neighbours[:-1,1:,5] = label[1:,:-1]
            
            neighbours.shape = -1,6
            
            vi.shape = -1
            vj.shape = -1
            
            ri = ((vi + sinV*vj)**2 + (cosV*vj**2))
            select = ri < (radius/spacing)**2
            
            xi = vi[select]
            yj = vj[select]
            
            backmap = numpy.full(vi.shape[0]+1, -1, 'i')
            backmap[:-1][select] = numpy.arange(len(xi))
            selnb = backmap[neighbours]
            neighbours = selnb[select,:]
            
            offsets = xi[:,None]*d1 + yj[:,None]*d2
            offsets *= spacing
            
            directions = direction[None,:] - offsets
            directions = normaliseVector(directions)
            
            E_vector = numpy.cross(self.E_vector, directions)
            E_vector = numpy.cross(E_vector, directions)
            E_vector = normaliseVector(E_vector)
            
            r_na = numpy.tan(numpy.arcsin(self.numerical_aperture))
            gauss = numpy.exp(-(ri[select]*(spacing**2))/(r_na**2))*spacing
            
            ray_data = numpy.zeros(offsets.shape[0], dtype=ray_dtype)
                
            wl = self.wavelength
            path_length = numpy.sqrt(((offsets + direction[None,:])**2).sum(axis=-1))*self.working_dist
            phase = -2*numpy.pi*((1000*path_length/wl))#%1.0)
            ray_data['origin'] = offsets*self.working_dist
            ray_data['origin'] += origin
            ray_data['direction'] = directions
            ray_data['wavelength_idx'] = 0
            ray_data['E_vector'] = E_vector
            ray_data['E1_amp'] = self.E1_amp * gauss
            ray_data['E2_amp'] = self.E2_amp * gauss
            ray_data['refractive_index'] = 1.0+0.0j
            ray_data['normal'] = [[0,1,0]]
            ray_data['phase'] = phase
            rays = RayCollection.from_array(ray_data)
            rays.neighbours = neighbours
        except:
            traceback.print_exc()
            return RayCollection(0)
        print("Made input rays.")
        return rays
    
    
class SingleGaussletSource(SingleRaySource):
    beam_waist = Float(100.0) #in microns
    working_distance = Float(0.0) #the distance from the beam waist to the source location
    
    input_rays = Property(Instance(RayCollection), 
                         depends_on="origin, direction, beam_waist, wavelength, "
                         "max_ray_len, E_vector, working_distance")
    
    geom_grp = VGroup(Group(Item('origin', show_label=False,resizable=True), 
                            show_border=True,
                            label="Origin position",
                            padding=0),
                       Group(Item('direction', show_label=False, resizable=True),
                            show_border=True,
                            label="Direction"),
                       Item("wavelength"),
                       Item('beam_waist', tooltip="1/e^2 beam intensity radius, in microns"),
                       Item('working_distance', tooltip="The distance from the source to the beam waist"),
                       label="Geometry")
    
    
    @on_trait_change("direction, beam_waist, max_ray_len, working_distance")
    def on_update(self):
        self.data_source.modified()
        self.update=True
    
    @cached_property
    def _get_input_rays(self):
        try:
            origin = numpy.array(self.origin)
            direction = numpy.array(self.direction)
            radius = self.beam_waist/2000.0 #convert to mm
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
    
            count = 3
            ray_data = numpy.zeros((2*count)+1, dtype=ray_dtype)
            
            theta0 = self.wavelength/(numpy.pi*radius*1000.0)
            lr = numpy.array([1,-1])[:,None,None]
            eye = numpy.array([1,1])[:,None,None]
            alpha = (numpy.arange(count)*(2*numpy.pi/count))[None,:,None]
            origins = eye*radius*(d1*numpy.sin(alpha) + d2*numpy.cos(alpha))
            directions = theta0*lr*(d1*numpy.cos(alpha) - d2*numpy.sin(alpha))
            
            origins.shape = (-1,3)
            directions.shape = (-1,3)
            directions += direction
            directions = normaliseVector(directions)
            origins = origins - (self.working_distance*directions)
                
            ray_data['origin'] = origin
            ray_data['origin'][1:] += origins + (self.working_distance*direction)
            ray_data['direction'][1:] = directions
            ray_data['direction'][0] = direction
            ray_data['wavelength_idx'] = 0
            ray_data['E_vector'] = [normaliseVector(E_vector)]
            ray_data['E1_amp'] = self.E1_amp
            ray_data['E2_amp'] = self.E2_amp
            ray_data['refractive_index'] = 1.0+0.0j
            ray_data['normal'] = [[0,1,0]]
            ray_data['ray_type_id'] = GAUSSLET_
            ray_data['ray_type_id'][1:] |= PARABASAL_ 
            rays = RayCollection.from_array(ray_data)
            
            neighbours = numpy.full((7, 6), -1, 'i')
            neighbours[0,:] = [1,2,3,4,5,6]
            rays.neighbours = neighbours
        except:
            import traceback
            traceback.print_exc()
            return RayCollection(5)
        return rays
    
    
