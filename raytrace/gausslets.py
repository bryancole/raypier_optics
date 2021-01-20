"""
Functions relating to the creation and manipulations of gausslets
"""

import numpy

from numpy import fft
from scipy.interpolate import RectBivariateSpline

from raytrace.utils import normaliseVector
from raytrace.ctracer import ray_dtype, GaussletCollection, FaceList
from raytrace.bases import Traceable
from raytrace.cmaterials import ResampleGaussletMaterial
from raytrace.cfaces import CircularFace
from raytrace.fields import eval_Efield_from_gausslets

from traits.api import Range, Float, Array, Property, Instance, observe, Str
from traitsui.api import View, Group, Item

from tvtk.api import tvtk
from raytrace.editors import NumEditor


def next_power_of_two(n):
    return int(2**numpy.ceil(numpy.log2(n)))


def make_hexagonal_grid(radius, spacing=1.0):
    cosV = numpy.cos(numpy.pi*30/180.)
    sinV = numpy.sin(numpy.pi*30/180.)
    
    nsteps = int(radius/(spacing*numpy.cos(numpy.pi*30/180.)))
    i = numpy.arange(-nsteps-1, nsteps+2)
    j = numpy.arange(-nsteps-1, nsteps+2)
    
    vi, vj = numpy.meshgrid(i,j)
    
    vi.shape = -1
    vj.shape = -1
    
    x = (vi + sinV*vj)*spacing
    y = cosV*vj*spacing
    
    r2 = x**2 + y**2
    select = r2 < (radius**2)
    return x[select], y[select]


def decompose_angle(origin, direction, axis1, E_field, input_spacing, max_angle, wavelength, oversample=4):
    """
    Compute a set of Gausslets over a range of directions, where the gausslet amplitude is
    obtained by FFT of the input E-field profile. The Gausslets all have an origin at the given
    origin point.
    
    params:
        origin - a (x,y,z) position vector for the origin of the Gausslets and the centre of the E-field
                distribution
        direction - a direction vector giving the output optical axis
        axis1 - a vector orthogonal to the direction giving the E1 polarisation axis and the direction of the E-field
                array 1st axis.
        E_field - a complex array of shape (N,M,3). The vector E-field.
        input_spacing - a scalar giving the sample spacing for the input E-field, in microns
        max_angle - sets the maximum output angle for the outgoing gausslets. Rays outside this angle are omitted. Units=degrees
        wavelength - sets the wavelength for the source, in microns
        oversample - sets the amount of padding added to the E_field data to increase the oversampling of the output.
    """
    input_spacing /= 1000.0
    wavelength /= 1000.0
    direction = normaliseVector(direction)
    d2 = normaliseVector(numpy.cross(direction, axis1))
    d1 = numpy.cross(direction, d2)
    
    N,M = E_field.shape[:2]
    target_size = next_power_of_two(max(N,M)) * oversample
    
    n_start = int((target_size/2) - (N/2))
    m_start = int((target_size/2) - (M/2))
    
    data_in = numpy.zeros((target_size,target_size,3), dtype=numpy.complex128)
    
    data_in[n_start:n_start+N, m_start:m_start+M, :] = E_field
    data_in = fft.ifftshift(data_in, axes=(0,1))
    
    data_out = fft.fft2(data_in, axes=(0,1)) / (target_size**2)
    data_out = fft.fftshift(data_out, axes=(0,1))
        
    ###The individual gausslet beam-waist radii at the origin are determined from the
    ### angular spacing of the output gausslets, and the wavelength.
    kmax = 2.0*numpy.pi/(2*input_spacing) #
    
    k_abs = 2.0*numpy.pi/wavelength #where wavelength is in microns
    
    kr = numpy.linspace(-kmax,kmax,data_out.shape[0])
    
    k_limit = k_abs*numpy.sin(numpy.pi*max_angle/180.0)
    k_grid_spacing = (kr[1]-kr[0]) #*oversample #The spacing for the hexagonal grid
    
    kx,ky = make_hexagonal_grid(k_limit, spacing=k_grid_spacing)
    _x = kx
    _y = ky
    
    kz = numpy.sqrt(k_abs**2 - kx**2 - ky**2)
        
    directions = (kx[:,None]*d1 + ky[:,None]*d2 + kz[:,None]*direction)/k_abs
    #directions = normaliseVector(directions)
    #print("interp:", x.min(),x.max(),y.min(),y.max(), kmax, kz, k_limit)
    
    Ex_real = RectBivariateSpline(kr,kr,data_out[:,:,0].real)(_x, _y, grid=False)
    Ex_imag = RectBivariateSpline(kr,kr,data_out[:,:,0].imag)(_x, _y, grid=False)
    Ey_real = RectBivariateSpline(kr,kr,data_out[:,:,1].real)(_x, _y, grid=False)
    Ey_imag = RectBivariateSpline(kr,kr,data_out[:,:,1].imag)(_x, _y, grid=False)
    Ez_real = RectBivariateSpline(kr,kr,data_out[:,:,2].real)(_x, _y, grid=False)
    Ez_imag = RectBivariateSpline(kr,kr,data_out[:,:,2].imag)(_x, _y, grid=False)
    
    E_real = numpy.column_stack((Ex_real, Ey_real, Ez_real))
    E_imag = numpy.column_stack((Ex_imag, Ey_imag, Ez_imag))
    
    E_full = E_real - 1.0j*E_imag
    
    E_vectors = normaliseVector(numpy.cross(directions, d2))
    E2_vectors = normaliseVector(numpy.cross(E_vectors, directions))
    
    E1_amp = (E_vectors * E_full).sum(axis=1)
    E2_amp = (E2_vectors * E_full).sum(axis=1)
    
    ray_data = numpy.zeros(directions.shape[0], dtype=ray_dtype)
    
    ray_data['origin'] = origin
    ray_data['direction'] = directions
    ray_data['wavelength_idx'] = 0
    ray_data['E_vector'] = E_vectors
    ray_data['E1_amp'] = E1_amp
    ray_data['E2_amp'] = E2_amp
    ray_data['refractive_index'] = 1.0+0.0j
    ray_data['normal'] = [[0,1,0]]
    ray_data['phase'] = 0.0
    rays = GaussletCollection.from_rays(ray_data)
    rays.wavelengths = wl = numpy.array([wavelength])
    working_dist=0.0
    #calculate beam-waists in microns
    #gausslet_radius = 2.0*numpy.pi/(k_grid_spacing * oversample)
    gausslet_radius = wavelength*k_abs/(k_grid_spacing*numpy.pi)
    print("Gausslet radius:", gausslet_radius)
    rays.config_parabasal_rays(wl, gausslet_radius, working_dist)
    rays.wavelengths = wl
    return rays, data_out


def decompose_position():
    pass


class AngleDecompositionPlane(Traceable):
    name = Str("Angle Decomposition")
    sample_spacing = Float(1.0) #in microns
    
    width = Range(0,512,64)
    height = Range(0,512,64)
    
    _mask = Instance(numpy.ndarray)
    mask = Property()
    
    ### Sets the geometry of the intersection face
    diameter = Float(25.0)
    offset = Float(0.0)
    
    ### Maximum ray ange in degrees
    max_angle = Float(30.0)
    
    material = Instance(ResampleGaussletMaterial)
    
    _E_field = Array() #store the E_field at the aperture for debugging purposes
    _k_field = Array()
    
    traits_view = View(Group(
                    Traceable.uigroup,
                    Item("sample_spacing", editor=NumEditor),
                    Item("width", editor=NumEditor),
                    Item("height", editor=NumEditor),
                    Item("max_angle", editor=NumEditor)
                    ))
    
    def evaluate_decomposed_rays(self, input_rays):
        origin = numpy.asarray(self.centre)
        direction = numpy.asarray(self.direction)
        axis1 = numpy.asarray(self.x_axis)
        axis2 = numpy.cross(axis1, direction)
        spacing = self.sample_spacing
        
        wavelengths = input_rays.wavelengths
        
        x = numpy.arange(self.width)*spacing/1000.0
        x -= x.mean()
        y = numpy.arange(self.height)*spacing/1000.0
        y -= y.mean()
        
        points = x[:,None,None]*axis1[None,None,:] + y[None,:,None]*axis2[None,None,:]
        points += origin[None,None,:]
        
        E_field = eval_Efield_from_gausslets(input_rays, points.reshape(-1,3), wavelengths).reshape(*points.shape)
        self._E_field = E_field
        mask = self._mask
        if mask is not None:
            E_field *= mask[:,:,None]
        print("Efield:", E_field.shape)
        print("spacing:", spacing)
        max_angle = self.max_angle
        new_rays, debug_k = decompose_angle(origin, direction, axis1, E_field, spacing, 
                                   max_angle, wavelengths[0], oversample=4)
        print("Decomposed rays count:", len(new_rays))
        self._k_field = debug_k
        return new_rays
    
    def _get_mask(self):
        return self._mask
    
    def _set_mask(self, obj):
        w,h = obj.shape
        self._mask = numpy.clip(obj.astype('d'), 0.0, 1.0)
        self.width = w
        self.height = h
        
    @observe("width, height")
    def _chk_dims(self, evt):
        mask = self._mask
        if mask is None:
            return
        if (self.width, self.height) != mask.shape:
            raise ValueError("Width and height must match shape of mask.")
    
    def _material_default(self):
        m = ResampleGaussletMaterial(eval_func=self.evaluate_decomposed_rays)
        return m
    
    def _faces_default(self):
        fl = FaceList(owner=self)
        fl.faces =  [CircularFace(owner=self, z_plane=0.0, material=self.material)]
        return fl
        
    def _pipeline_default(self):
        radius = self.diameter/2.
        circle = tvtk.ArcSource(use_normal_and_angle=True,
                                center=[0.0,0.0,0.0],
                                polar_vector=[0.0, radius, 0.0],
                                normal=[0.0,0.0,1.0],
                                angle=360.0
                                )
        
        trans = tvtk.TransformFilter(input_connection=circle.output_port,
                                     transform=self.transform)
        return trans
