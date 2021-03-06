
import numpy
import itertools
import time
from traits.api import Instance, Title, Float, Tuple, Complex, Property, cached_property, observe, Bool, Range
from traitsui.api import View, Item, VGroup, Tabbed, Include, Group
from tvtk.api import tvtk

from raypier.sources import BaseRaySource, UnitTupleVector, UnitVectorTrait
from raypier.core.ctracer import GaussletCollection, gausslet_dtype, ray_dtype
from raypier.tracer import normaliseVector
from raypier.editors import NumEditor
from raypier.decompositions import decompose_angle, make_hexagonal_grid


class BaseGaussletSource(BaseRaySource):
    """
    Abstract base class for Gausslet Sources.
    """
    #:object title
    name = Title("Gausslet Source", allow_selection=False)
    
    para_data_source = Instance(tvtk.ProgrammableSource, transient=True)
    para_tube = Instance(tvtk.TubeFilter, (), transient=True)
    para_ray_actor = Instance(tvtk.Actor, transient=True)
    para_property = Instance(tvtk.Property, (), {'color':(0.0,0.5,0.5), 'opacity': 1.0}, transient=True)
    para_mapper = Instance(tvtk.PolyDataMapper, (), transient=True)
    
    #: set True to render the parabasal rays
    show_paras = Bool(False)
    
    #: wavelength, in microns.
    wavelength = Float(0.78) # in microns
    
    #: Total power in the source beam, in Watts, calculated as the sum of the power in each individual
    #: gausslet. This is, in fact,  *not* the true beam power. Due to the overlap between the modes,
    #: the actual beam power due to the coherent sum of the modes is greater than the sum of the 
    #: single-mode powers.
    beam_power = Float(1.0)
    
    display_grp = VGroup(Item('display'),
                       Item('show_start'),
                       Item('show_normals'),
                       Item('show_paras'),
                       Item('scale_factor'),
                       Item('export_pipes',label="export as pipes"),
                       Item('opacity', editor=NumEditor),
                       label="Display")
    
    params_grp = VGroup(
                       Item("max_angle"),
                       Item('beam_waist', tooltip="1/e^2 beam intensity radius, in microns"),
                       Item('sample_spacing', tooltip="Sample spacing for the field in the aperture"),
                       label="Parameters")
    
    traits_view = View(Item('name', show_label=False),
                       Tabbed(VGroup(Group(Item('origin', show_label=False,resizable=True), 
                                           show_border=True,
                                    label="Origin position",
                                    padding=0),
                               Group(Item('direction', show_label=False, resizable=True),
                                    show_border=True,
                                    label="Direction"),
                               Item("wavelength"),
                               Include('params_grp'),
                            ),
                              display_grp)
                       )
    
    def _traced_rays_changed(self):
        self.data_source.modified()
        self.para_data_source.modified()
        self.normals_source.modified()
        self._mtime = time.monotonic()
        
    def _input_rays_changed(self):
        self.update=True
        
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
                yield gc.base_rays.origin.reshape(-1,3)
                if self.show_paras:
                    yield gc.para_origin.reshape(-1,3)
                
        def get_normals(gc_list):
            for gc in gc_list:
                yield gc.base_rays.normal.reshape(-1,3)
                if self.show_paras:
                    yield gc.para_normal.reshape(-1,3)
                
        def execute():
            if not self.show_normals:
                return
            output = source.poly_data_output
            points = numpy.vstack( list( get_origins(self.traced_rays[1:]) ) )
            normals = numpy.vstack( list( get_normals(self.traced_rays[1:]) ) )
            output.points = points
            output.point_data.normals = normals
            #print("calc normals GLYPH")
        source.set_execute_method(execute)
        return source
    
    @observe("show_paras")
    def _on_show_paras_changed(self, evt):
        self.para_data_source.modified()
        self.render=True
    
    def _para_data_source_default(self):
        source = tvtk.ProgrammableSource()
        def execute():
            if not self.show_paras:
                return
            output = source.poly_data_output
            pointArrayList = []
            mask = self.ray_mask
            for rays in self.traced_rays:
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
    """
    A source which launches a single Gausslet. Mainly used for testing.
    """
    abstract=False
    
    #: The origin of the Gausslet base-ray
    origin = Tuple((0.,0.,0.))
    
    #: Base-ray direction unit-vector.
    direction = UnitTupleVector
    
    #: Polarisation unit-vector for the E1 amplitude. Orthogonal to the direction vector.
    E_vector = UnitVectorTrait((1.,0.,0.), editor_traits={'cols':3,
                                'labels':['x','y','z']})
    
    #: Complex amplitude for the E1_vector polarisation direction. 
    E1_amp = Complex(1.0+0.0j)
    
    #: Complex amplitude for the E2_vector, implicitly orthogonal to both `direction` and `E1_vector`.
    E2_amp = Complex(0.0+0.0j)
    
    #: Gaussian 1/e**2 beam radius, in microns, at the beam waist.
    beam_waist = Float(1.0) #radius, in microns
    
    #: The GaussletCollection launched by this object.
    input_rays = Property(Instance(GaussletCollection), 
                         depends_on="origin, direction, max_ray_len, E_vector, E1_amp, E2_amp, beam_waist")
    
    params_grp = VGroup(
                       Item("wavelength"),
                       Item("max_angle"),
                       Item('beam_waist', tooltip="1/e^2 beam intensity radius, in microns"),
                       Item('sample_spacing', tooltip="Sample spacing for the field in the aperture"),
                       label="Parameters")
    
    def _wavelength_list_default(self):
        return [self.wavelength,]
    
    @observe("wavelength")
    def _do_wavelength_changed(self, evt):
        self.wavelength_list = [self.wavelength]
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
    
    
class TopHatGaussletSource(SingleGaussletSource):
    """
    Launches Gausslets forming a top-hat intensity profile at the beam-waist.
    """
    
    #: The diameter of the beam at the beam-waist, in microns.
    beam_waist = Float(100.0) #diameter, in microns
    
    #: The sample-spacing across the beam waist. 
    sample_spacing = Float(10.0) #in microns
    
    #: The maximum divergence angle for the Gausslet base-rays.
    max_angle = Float(10.0)
    
    input_rays = Property(Instance(GaussletCollection), 
                         depends_on="origin, direction, beam_waist, wavelength, "
                         "max_ray_len, E_vector, sample_spacing, max_angle")
    
    geom_grp = VGroup(Group(Item('origin', show_label=False,resizable=True), 
                            show_border=True,
                            label="Origin position",
                            padding=0),
                       Group(Item('direction', show_label=False, resizable=True),
                            show_border=True,
                            label="Direction"),
                       Item("wavelength"),
                       Item("max_angle"),
                       Item('beam_waist', tooltip="1/e^2 beam intensity radius, in microns"),
                       Item('sample_spacing', tooltip="Sample spacing for the field in the aperture"),
                       label="Geometry")
    
    
    @observe("direction, beam_waist, max_ray_len, sample_spacing, max_angle")
    def on_update(self, evt):
        self.data_source.modified()
        self.update=True
    
    @cached_property
    def _get_input_rays(self):
        try:
            origin = numpy.array(self.origin)
            direction = numpy.array(self.direction)
            radius = self.beam_waist/2.0 #convert to mm
            
            max_axis = numpy.abs(direction).argmax()
            if max_axis==0:
                v = numpy.array([0.,1.,0.])
            else:
                v = numpy.array([1.,0.,0.])
                
            max_angle=self.max_angle
            input_spacing = self.sample_spacing
            wavelength = self.wavelength
            oversample=4
            
            x = numpy.linspace(-radius,radius, int(radius/input_spacing)+1)
            X,Y = numpy.meshgrid(x,x)
            select = (X**2 + Y**2) < radius
            E_field = numpy.zeros((2,len(x),len(x)), dtype=numpy.complex128)
            E1 = E_field[0,:,:]
            E1[select] = 1.0 
                
            rays = decompose_angle(origin, direction, v, E_field, input_spacing, max_angle, wavelength, oversample)
                
        except:
            import traceback
            traceback.print_exc()
            return GaussletCollection(5)
        return rays
    
    
class CollimatedGaussletSource(SingleGaussletSource):
    """
    Launches a quasi-collimated beam of Gausslets.
    """
    
    #: The max radius of the Gausslets at the objects origin.
    radius = Float(10.,editor=NumEditor)
    
    #: Sets the number of Gausslets found along the radius of the beam, 
    #: I.e. Controls the spacing of the Gausslets, being radius/resolution.  
    resolution = Float(10.0, editor=NumEditor)
    
    #: Controls the overlap between adjacent Gausslets. If blending==1.0, the 
    #: Gaussian 1/e^2 intensity width is equal to the Gausslet spacing.
    #: The Gausslet width = spacing/blending. I.e. a higher value reduces the Gausslet overlap.
    blending = Float(1.5)
    
    input_rays = Property(Instance(GaussletCollection), 
                         depends_on="origin, direction, max_ray_len, E_vector, E1_amp, "
                                    "E2_amp, radius, resolution, wavelength, beam_waist, blending")
    
    params_grp = VGroup(
                        Item("radius", editor=NumEditor),
                        Item("resolution", editor=NumEditor),
                       Item('beam_waist', tooltip="1/e^2 beam intensity radius, in microns", editor=NumEditor),
                       Item("blending", editor=NumEditor),
                       label="Parameters")
    
    @cached_property
    def _get_input_rays(self):
        self.wavelength_list = [self.wavelength]
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
        
        vi.shape = -1
        vj.shape = -1
        
        ri = ((vi + sinV*vj)**2 + (cosV*vj**2))
        select = ri < (radius/spacing)**2
        
        xi = vi[select]
        yj = vj[select]
        
        offsets = xi[:,None]*d1 + yj[:,None]*d2
        offsets *= spacing
        gauss = numpy.exp(-(ri[select]*(spacing**2)/((self.beam_waist)**2))) 
        
        ray_data = numpy.zeros(offsets.shape[0], dtype=ray_dtype)
        
        E1_amp = self.E1_amp
        E2_amp = self.E2_amp
        #print("E_vals", E1_amp, E2_amp)
        P1 = (E1_amp * E1_amp.conjugate()).real
        P2 = (E2_amp * E1_amp.conjugate()).real
        
        ### The true E_field is E1_amp (or E2_amp) * sqrt(area of beam)###
        beamlet_area = (spacing**2) * numpy.pi
        total_power = ((P1 + P2) * (gauss**2)).sum() #* beamlet_area
        scaling = numpy.sqrt(self.beam_power/total_power)
            
        ray_data['origin'] = offsets
        ray_data['origin'] += origin
        ray_data['direction'] = direction
        ray_data['wavelength_idx'] = 0
        ray_data['E_vector'] = [normaliseVector(E_vector)]
        ray_data['E1_amp'] = scaling * self.E1_amp * gauss
        ray_data['E2_amp'] = scaling * self.E2_amp * gauss
        ray_data['refractive_index'] = 1.0+0.0j
        ray_data['normal'] = [[0,1,0]]
        ray_data['ray_type_id'] = (1<<1) #indicates Gausslets
        
        rays = GaussletCollection.from_rays(ray_data)
        wl = numpy.array(self.wavelength_list)
        rays.wavelengths = wl
        working_dist = 0.0
        rays.config_parabasal_rays(wl, spacing/self.blending, working_dist)
        return rays
        
        
class GaussianPointSource(SingleGaussletSource):
    """
    Creates a low-F# point source with symmetric Gaussian profile.
    
        resolution - the number of rays along the angular radius.
        
        working_dist - the distance in front of the origin plane at which the focus is found
        
        numerical_aperture - equal of the sin of the maximum ray angle
        
        blending - determines the mode overlap
    """
    abstract = False
    
    #: Number of angle steps along the radius
    resolution = Range(0.0, None, 10.0)
    
    #: The distance from the object origin to the beam waist. A negative value
    #: means that the focus is behind the origin. 
    working_dist = Float(100.0)
    
    #: The maximum numerical aperture of the base-rays in the beam.
    numerical_aperture = Float(0.2)
    
    #: Sets the overlap between adjacent Gausslets *in the far field*.
    #: I.e. a higher value increases the Gausslet beam waist, which
    #: reduces the divergence of the Gaussian modes such that their
    #: overlap in the far field is reduced.
    blending = Float(1.0)
    
    ### Traits Properties seems to be broken. No good way to flush the cache so I'm doing it myself.
    _input_rays = Instance(GaussletCollection)
    input_rays = Property(Instance(GaussletCollection))
    
    params_grp = VGroup(
                        Item("numerical_aperture", editor=NumEditor),
                        Item("resolution", editor=NumEditor),
                       Item('beam_waist', tooltip="1/e^2 beam intensity radius, in microns", editor=NumEditor),
                       Item('working_dist', editor=NumEditor),
                       Item("blending", editor=NumEditor),
                       label="Parameters")
    
    @observe("origin, direction, max_ray_len, E_vector, E1_amp, "\
            "E2_amp, numerical_aperture, resolution, wavelength, "\
            "beam_waist, blending, working_dist")
    def _clear_input_rays(self, evt):
        del self._input_rays
        self.update = True
    
    def _get_input_rays(self):
        rays = self._input_rays
        if rays is None:
            rays = self._eval_input_rays()
            self._input_rays = rays
        return rays
    
    def _eval_input_rays(self):
        self.wavelength_list = [self.wavelength]
        origin = numpy.array(self.origin)
        direction = normaliseVector(numpy.array(self.direction))
        
        d2 = normaliseVector(numpy.cross(direction, self.E_vector))
        d1 = numpy.cross(direction, d2)
        
        focal_radius = self.beam_waist/1000.0 #convert to mm
        
        n = 1.0 + 0.0j #refractive index of medium
        k_abs = (2000.0*numpy.pi)/self.wavelength #k now in inverse-mm
        k_radius = 2./focal_radius
        k_max = self.numerical_aperture * k_abs
        k_grid_spacing = k_max / self.resolution
        
        kx,ky = make_hexagonal_grid(k_max, spacing=k_grid_spacing)
        
        krsq = kx**2 + ky**2
        kz = numpy.sqrt(k_abs**2 - krsq)
        
        waist_location = origin + self.working_dist*direction 
        directions = (kx[:,None]*d1 + ky[:,None]*d2 + kz[:,None]*direction)/k_abs
        directions = normaliseVector(directions)
        
        gauss = numpy.exp(-(krsq/(k_radius**2)))
        E1_amp = gauss * self.E1_amp
        E2_amp = gauss * self.E2_amp
        
        ray_data = numpy.zeros(gauss.shape[0], dtype=ray_dtype)
            
        E_vectors = numpy.cross(directions, d2[None,:])
            
        ray_data['origin'] = waist_location[None,:]
        ray_data['direction'] = directions
        ray_data['wavelength_idx'] = 0
        ray_data['E_vector'] = E_vectors
        ray_data['E1_amp'] = E1_amp
        ray_data['E2_amp'] = E2_amp
        ray_data['refractive_index'] = 1.0+0.0j
        ray_data['normal'] = [[0,1,0]]
        ray_data['ray_type_id'] = (1<<1) #indicates Gausslets
        
        rays = GaussletCollection.from_rays(ray_data)
        wl = numpy.array(self.wavelength_list)
        rays.wavelengths = wl
        working_dist = 0.0
        gauss_radius = 2.*self.blending/(k_grid_spacing)
        rays.config_parabasal_rays(wl, gauss_radius, working_dist)
        rays.project_to_plane(origin, direction)
        return rays
    