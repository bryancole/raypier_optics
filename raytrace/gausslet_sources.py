
import numpy
import itertools
from traits.api import Instance, Title, Float, Tuple, Complex, Property, cached_property

from tvtk.api import tvtk

from raytrace.sources import BaseRaySource, UnitTupleVector, UnitVectorTrait
from raytrace.ctracer import GaussletCollection, gausslet_dtype
from raytrace.tracer import normaliseVector


class BaseGaussletSource(BaseRaySource):
    name = Title("Gausslet Source", allow_selection=False)
    para_data_source = Instance(tvtk.ProgrammableSource, transient=True)
    para_tube = Instance(tvtk.TubeFilter, (), transient=True)
    para_ray_actor = Instance(tvtk.Actor, transient=True)
    para_property = Instance(tvtk.Property, (), {'color':(0.0,0.5,0.5), 'opacity': 1.0}, transient=True)
    para_mapper = Instance(tvtk.PolyDataMapper, (), transient=True)
    
    def _TracedRays_changed(self):
        self.data_source.modified()
        self.para_data_source.modified()
        self.normals_source.modified()
        
    def _opacity_changed(self):
        self.vtkproperty.opacity = self.opacity
        self.para_property.opacity = self.opacity
        self.render = True
        
    def _get_actors(self):
        actors = [self.ray_actor, self.para_ray_actor, self.start_actor, self.normals_actor]
        return actors
    
    def _normals_source_default(self):
        source = tvtk.ProgrammableSource()
        
        def get_origins(gc_list):
            for gc in gc_list:
                yield gc.origin.reshape(-1,3)
                yield gc.para_origin.reshape(-1,3)
                
        def get_normals(gc_list):
            for gc in gc_list:
                yield gc.normal.reshape(-1,3)
                yield gc.para_normal.reshape(-1,3)
                
        def execute():
            if not self.show_normals:
                return
            output = source.poly_data_output
            points = numpy.vstack( list( get_origins(self.TracedRays[1:]) ) )
            normals = numpy.vstack( list( get_normals(self.TracedRays[1:]) ) )
            output.points = points
            output.point_data.normals = normals
            #print("calc normals GLYPH")
        source.set_execute_method(execute)
        return source
    
    def _para_data_source_default(self):
        source = tvtk.ProgrammableSource()
        def execute():
            output = source.poly_data_output
            pointArrayList = []
            mask = self.ray_mask
            for rays in self.TracedRays:
                vis = mask.get(rays, [])
                if len(vis) != len(rays):
                    vis = itertools.repeat(True)
                start_pos_ = (r.parabasal_rays for r,v in zip(rays, vis) if v)
                start_pos = [p.origin for p in itertools.chain.from_iterable(start_pos_)]
                end_pos_ = (r.parabasal_rays for r,v in zip(rays, vis) if v)
                end_pos = [p.termination for p in itertools.chain.from_iterable(end_pos_)]
                interleaved = numpy.array([start_pos, end_pos]).swapaxes(0,1).reshape(-1,3)
                pointArrayList.append(interleaved)
            if pointArrayList:
                points = numpy.vstack(pointArrayList)
                output.points=points
                lines = list([i,i+1] for i in range(0,points.shape[0],2))
                output.lines = lines
        source.set_execute_method(execute)
        return source
    
    def _para_ray_actor_default(self):
        tube = self.para_tube
        tube.radius = self.scale_factor
        tube.input_connection=self.para_data_source.output_port
        tube.number_of_sides = 20
        
        map = self.para_mapper
        act = tvtk.Actor(mapper=map)
        act.visibility = True
        display = self.display
        if display=="pipes":
            map.input_connection=tube.output_port
        elif display=="wires":
            map.input_connection = self.para_data_source.output_port
        else:
            act.visibility = False
        prop = self.para_property
        prop.opacity = self.opacity
        if prop:
            act.property = prop
        return act
    
    def _scale_factor_changed(self, scale):
        self.tube.radius = scale
        self.para_tube.radius = scale
        self.normals_glyph.scale_factor=self.scale_factor*10
        self.sphere.radius = self.scale_factor*3
        self.render = True
        
    def _display_changed(self, vnew):
        actors = self.actors
        for act in actors:
            act.visibility = True
        if vnew=="pipes":
            self.mapper.input_connection = self.tube.output_port
            self.para_mapper.input_connection = self.para_tube.output_port
        elif vnew=="wires":
            self.mapper.input_connection = self.data_source.output_port
            self.para_mapper.input_connection = self.para_data_source.output_port
        elif vnew=="hidden":
            for act in actors:
                act.visibility = False
        self.render = True
        
        
class SingleGaussletSource(BaseGaussletSource):
    abstract=False
    
    origin = Tuple((0.,0.,0.))
    direction = UnitTupleVector
    
    E_vector = UnitVectorTrait((1.,0.,0.), editor_traits={'cols':3,
                                'labels':['x','y','z']})
    E1_amp = Complex(1.0+0.0j)
    E2_amp = Complex(0.0+0.0j)
    wavelength = Float(0.78) # in microns
    
    beam_waist = Float(1.0)
    
    InputRays = Property(Instance(GaussletCollection), 
                         depends_on="origin, direction, max_ray_len, E_vector, E1_amp, E2_amp, beam_waist")
    
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
        
        angles = (numpy.arange(6)*numpy.pi/3.)[:,None]
        offsets = numpy.cos(angles)*d1[None,:] + numpy.sin(angles)*d2[None,:]
        offsets *= self.beam_waist
        
        data = numpy.zeros((1,), dtype=gausslet_dtype)
        base_ray = data['base_ray']
        base_ray['origin'] = origin[None,:]
        base_ray['direction'] = direction[None,:]
        base_ray['E_vector'] = E_vector[None,:]
        base_ray['E1_amp'] = 1.0
        base_ray['length'] = self.max_ray_len
        
        para = data['para_rays']
        para['direction'] = direction[None,None,:]
        para['origin'] = offsets[None,:,:]
        para['length'] = self.max_ray_len
        para['normal'] = -direction[None,:]
        
        gc = GaussletCollection.from_array(data)
        #print("A", [g.length for g in gc[0].parabasal_rays] , [self.max_ray_len]*6 )
        assert [g.length for g in gc[0].parabasal_rays] == [self.max_ray_len]*6
        return gc
        