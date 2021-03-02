

import numpy

from numpy import fft
from scipy.interpolate import RectBivariateSpline
from scipy.sparse.linalg import lsqr, lsmr
from scipy.sparse import coo_matrix

from .ctracer import ray_dtype, GaussletCollection
from .utils import normaliseVector
from .fields import eval_Efield_from_gausslets, EFieldSummation
from .unwrap2d import unwrap2d
from .cfields import calc_mode_curvature, build_interaction_matrix, apply_mode_curvature


root_pi = numpy.sqrt(numpy.pi)

def next_power_of_two(n):
    return int(2**numpy.ceil(numpy.log2(n)))

def lagrange_invariant(gausslet):
    base_ray = gausslet.base_ray
    origin = numpy.asarray(base_ray.origin)
    direction = numpy.asarray(base_ray.direction)
    axis1 = numpy.asarray(base_ray.E_vector)
    axis2 = numpy.cross(axis1, direction)
    
    paras = gausslet.parabasal_rays
    
    def get_origin(porigin):
        return numpy.asarray(porigin) - origin
    
    def get_direction(pdir):
        return numpy.asarray(pdir)# - direction
    
    h1 = [(get_origin(p.origin)*axis1).sum() for p in paras]
    h2 = [(get_origin(p.origin)*axis2).sum() for p in paras]
    u1 = [(get_direction(p.direction)*axis1).sum() for p in paras]
    u2 = [(get_direction(p.direction)*axis2).sum() for p in paras]
    
    def HU(i,j):
        return (h1[i]*u1[j] + h2[i]*u2[j]) - ( h1[j]*u1[i] + h2[j]*u2[i] )
    
    #out = [HU(0,3), HU(1,2), HU(0,5), HU(1,4), HU(3,4), HU(2,5)]
    out = [HU(0,5), HU(1,2), HU(3,4), HU(1,4), HU(0,3), HU(2,5)]
#     out = numpy.zeros((6,6),'d')
#     for i in range(6):
#         for j in range(6):
#             out[i,j] = 
            
    return out


def make_hexagonal_grid(radius, spacing=1.0, connectivity=False):
    """Creates a 2d hexagonal grid.
    Radius and spacing have the same units.
    Returns - x,y: 1d arrays of the respective coordinate.
    if connectivity=True:
        Returns - x, y, neighbours: where neighbours is an (N,6) array of indices giving the adjacent neighbours to each coordinate
    """
    
    cosV = numpy.cos(numpy.pi*30/180.)
    sinV = numpy.sin(numpy.pi*30/180.)
    
    nsteps = int(radius/(spacing*numpy.cos(numpy.pi*30/180.)))
    i = numpy.arange(-nsteps-1, nsteps+2)
    j = numpy.arange(-nsteps-1, nsteps+2)
    
    vi, vj = numpy.meshgrid(i,j)
    
    ni, nj = vi.shape
    vi.shape = -1
    vj.shape = -1
    
    x = (vi + sinV*vj)*spacing
    y = cosV*vj*spacing
    
    r2 = x**2 + y**2
    select = r2 < (radius**2)
    
    if connectivity:
        labels = -numpy.ones((ni+2, nj+2), numpy.int)
        labels[1:-1, 1:-1] = numpy.arange(ni*nj).reshape(ni,nj)
        
        centres = -numpy.ones((vi.size,), numpy.int)
        neighbours = numpy.dstack([labels[2:,1:-1], labels[:-2,1:-1],
                                   labels[1:-1,2:], labels[1:-1,:-2],
                                   labels[2:,:-2], labels[:-2,2:]])
        neighbours.shape = (-1,6)
        
        sel = centres
        centres[select] = numpy.arange(select.sum())
        nb = centres[numpy.compress(select, neighbours, axis=0)]
        
        return x[select], y[select], nb
    else:
        return x[select], y[select]


def decompose_angle(origin, direction, axis1, E_field, input_spacing, max_angle, wavelength, 
                    oversample=4, E_max=None, pos_max=None):
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
    
    kr = numpy.linspace(-kmax,kmax,data_out.shape[0]+1)[:-1]
    
    k_limit = k_abs*numpy.sin(numpy.pi*max_angle/180.0)
    k_grid_spacing = (kr[1]-kr[0])*oversample #The spacing for the hexagonal grid
    
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
    
    E_full = -(-1.0j*E_real + E_imag)
    
    E_vectors = normaliseVector(numpy.cross(directions, d2))
    E2_vectors = normaliseVector(numpy.cross(directions, E_vectors))
    
    E1_amp = (E_vectors * E_full).sum(axis=1)
    E2_amp = (E2_vectors * E_full).sum(axis=1)
    
    ray_data = numpy.zeros(directions.shape[0], dtype=ray_dtype)
    
    gausslet_radius = wavelength*k_abs/(k_grid_spacing*numpy.pi)
    
    ### Some sort of guess to get the output power roughly right.
    scaling = numpy.sqrt(2)*input_spacing*target_size
    
    ray_data['origin'] = origin + ((d1+d2)*(input_spacing*0.5))
    ray_data['direction'] = directions
    ray_data['wavelength_idx'] = 0
    ray_data['E_vector'] = E_vectors
    ray_data['E1_amp'] = E1_amp*scaling
    ray_data['E2_amp'] = E2_amp*scaling
    ray_data['refractive_index'] = 1.0+0.0j
    ray_data['normal'] = [[0,1,0]]
    ray_data['phase'] = 0.0
    rays = GaussletCollection.from_rays(ray_data)
    rays.wavelengths = wl = numpy.array([wavelength])
    working_dist=0.0
    
    print("Gausslet radius:", gausslet_radius)
    rays.config_parabasal_rays(wl, gausslet_radius, working_dist)
    rays.wavelengths = wl
    
    if E_max is not None and pos_max is not None:
        print("E_max:", E_max, "pos_max:", pos_max)
        e_test = eval_Efield_from_gausslets(rays, pos_max[None,:], blending=1.0)
        print("e_test:", e_test)
        scaling = (E_max/e_test).mean()
        #ray_data['E1_amp'] *= scaling
        #ray_data['E2_amp'] *= scaling
        
        rays = GaussletCollection.from_rays(ray_data)
        rays.wavelengths = wl
        rays.config_parabasal_rays(wl, gausslet_radius, working_dist)
    
    return rays, data_out


def decompose_position(input_rays, origin, direction, axis1, radius, resolution, curvature=None, blending=1.5):
    """
    Spatially decompose the input Gausslets into a new set of Gausslets defined over a circular area of the 
    given radius. The output ray density is given by the resolution parameter.
    
    Parameters
    ----------
    
    input_rays : GaussletCollection
                 The input rays for the decomposition
                 
    origin : tuple, (float, float, float)
             The centre of the decomposition plane
             
    direction : tuple (float, float, float)
                The normal-vector for the decomposition plane
                
    axis1 : tuple (float, float, float)
            A vector defining the orientation of the decomposition plane (orthogonal to direction)
                
    radius : float
             The maximum radius giving the limit of the output rays
             
    resolution : float
                 The number of rays generated between the decomposition plane centre and radius edge
                 
    curvature : float, optional
                The approximate radius-of-curvature of the wavefront. I.e. distance to focus. None=Inf.
                
    Returns
    -------
    
    rays : GaussletCollection
            The out-going Gausslets
    """
    spacing = radius / resolution
    wavelengths = input_rays.wavelengths
    
    if len(numpy.unique(wavelengths)) > 1:
        raise ValueError("Can't decompose wavefront with more than one wavelength present.")
    
    wavelength = wavelengths[0]
    
    origin = numpy.asarray(origin)
    direction = normaliseVector(numpy.asarray(direction))
    axis1 = normaliseVector(numpy.asarray(axis1))
    axis2 = normaliseVector(numpy.cross(axis1, direction))
    axis1 = numpy.cross(axis2, direction)
    
    ###In this version, we will evaluate the field on a cartesian grid
    _radius = radius * (1 + 2./resolution)
    x_ = y_ = numpy.linspace(-_radius,_radius, int((1.1*resolution)*2))
    
    x,y = numpy.meshgrid(x_, y_)
    
    origins_in = origin[None,:] + x.reshape(-1,)[:,None]*axis1[None,:] + y.reshape(-1)[:,None]*axis2[None,:]
    
    efield = EFieldSummation(input_rays)
    E_field = efield.evaluate(origins_in)
    
    ### Project onto local axes to get orthogonal polarisations and choose the one with the most power
    E1_amp = (E_field*axis1[None,:]).sum(axis=1)
    E2_amp = (E_field*axis2[None,:]).sum(axis=1)
    
    P1 = (E1_amp.real**2 + E1_amp.imag**2).sum()
    P2 = (E2_amp.real**2 + E2_amp.imag**2).sum()
    
    if P1 > P2:
        ### Always in the range -pi to +pi
        phase = numpy.arctan2(E1_amp.imag, E1_amp.real)
    else:
        phase = numpy.arctan2(E2_amp.imag, E2_amp.real)
    
    nx = len(x_)
    ny = len(y_)
    phase.shape = (nx, ny)
    
    if curvature is not None:
        sphz = -(curvature - numpy.sqrt(curvature*curvature - x*x - y*y))
        sph = sphz*2000.0*numpy.pi/wavelength - numpy.pi
        phase = phase - sph
        phase = phase%(2*numpy.pi)
        phase = phase - numpy.pi
    else:
        sph = 0
    
    uphase, residuals = unwrap2d(phase, anchor=(nx//2, ny//2))
    uphase2 = uphase + sph
    
    ### Setting s=0 results in oscillatory behaviour near the edge along the x-idrection.
    wavefront = RectBivariateSpline(x_, y_, -uphase2*wavelength/(2000*numpy.pi), kx=3, ky=3, s=0.001)
    
    ### Now make a hex-grid of new ray start-points
    rx,ry, nb = make_hexagonal_grid(radius, spacing=spacing, connectivity=True)
    
    N = len(rx)
    
    origins = origin[None,:] + rx[:,None]*axis1[None,:] + ry[:,None]*axis2[None,:]
    
    ### First derivatives give us the ray directions
    dx = wavefront(rx,ry,dx=1, grid=False)
    dy = wavefront(rx,ry,dy=1, grid=False)
    
    ### Get second derivatives, to get wavefront curvature
    dx2 = wavefront(rx,ry,dx=2, grid=False)
    dy2 = wavefront(rx,ry,dy=2, grid=False)
    dxdy = wavefront(rx,ry,dx=1,dy=1, grid=False)

    ###Convert to A,B,C coefficients
    A,B,C,xl,yl,zl = calc_mode_curvature(rx, ry, dx, dy, dx2, dy2, dxdy)
    
    E_in = efield.evaluate(origins) #should be a (N,3) array of complex values
    
    directions = axis1[None,:]*zl[:,0,None] + axis2[None,:]*zl[:,1,None] + direction[None,:]*zl[:,2,None]
    E_vectors = axis1[None,:]*xl[:,0,None] + axis2[None,:]*xl[:,1,None] + direction[None,:]*xl[:,2,None]
    H_vectors = axis1[None,:]*yl[:,0,None] + axis2[None,:]*yl[:,1,None] + direction[None,:]*yl[:,2,None]
        
    E1_amp = (E_in*E_vectors).sum(axis=-1)
    E2_amp = (E_in*H_vectors).sum(axis=-1)
    
    ray_data = numpy.zeros(N, dtype=ray_dtype)
    
    ray_data['origin'] = origins
    ray_data['direction'] = directions
    ray_data['E_vector'] = E_vectors
    ray_data['refractive_index'] = 1.0
    ray_data['E1_amp'] = E1_amp
    ray_data['E2_amp'] = E2_amp
    
    gausslets = GaussletCollection.from_rays(ray_data)
    gausslets.wavelengths = wavelengths
    gausslets.config_parabasal_rays(wavelengths, spacing/blending, 0.0)
    apply_mode_curvature(gausslets, -A, -B, -C)
    
    E_test = eval_Efield_from_gausslets(gausslets, origins)
    
    power_scaling = (E_in.real**2 + E_in.imag**2).sum() / (E_test.real**2 + E_test.imag**2).sum()
    
    gausslets.scale_amplitude(numpy.sqrt(power_scaling))
    
    return gausslets, E_field, uphase
    