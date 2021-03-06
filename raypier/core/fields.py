
"""
Functions relating to the evaluation of the optical E-field by
summation of General Astigmatic Gaussian Beams
"""


from .cfields import sum_gaussian_modes, evaluate_modes as evaluate_modes_c
from .utils import normaliseVector, dotprod
from .ctracer import RayCollection

import numpy
scipy = None



def find_ray_gen(probe_centre, traced_rays):
    probe = numpy.asarray(probe_centre)
    waypoints = [rays.origin.mean(axis=0) for rays in traced_rays]
    waypoints.append(traced_rays[-1].termination.mean(axis=0))
    
    dist = [ ((wpt-probe)**2).sum() for wpt in waypoints ]
    
    imin = numpy.argmin(dist)
    if imin + 1 >= len(waypoints):
        return -1
    
    closest = waypoints[imin]
    
    dir1 = closest - waypoints[imin-1]
    dir1 = normaliseVector(dir1)
    
    dir2 = waypoints[imin+1] - closest
    dir2 = normaliseVector(dir2)
    
    normal = normaliseVector(dir1 + dir2)
    if dotprod(probe-closest, normal) > 0:
        idx = imin
    else:
        idx = imin-1
    if idx < 0:
        idx = 0
    if idx > len(waypoints)-2:
        idx = -1
    return idx
    
    
    ### unfinished ###


def project_to_sphere(rays, centre=(0,0,0), radius=10.0):
    """project the given set of rays back to their intercept 
    with a sphere at the given centre and radius.
    
    rays - an array of ray_t dtype
    """
    centre = numpy.asarray(centre).reshape(1,3)
    origin = rays['origin'] - centre
    direction = rays['direction'] #assume this is already normalised
    c = (origin**2).sum(axis=1) - (radius**2)
    b = 2*(direction * origin).sum(axis=1)
    a = 1 
    
    d = b**2 - 4*c
    selector = (d>=0.0)
    d = numpy.sqrt(d[selector])
    b = b[selector]
    root1 = (-b + d)/2
    root2 = (-b - d)/2
    #choose most negative path to intersection
    alpha = numpy.column_stack((root1, root2)).min(axis=1)
    
    rays = rays[selector]
    rays['origin'] += alpha[:,None]*direction[selector]
    
    rays['accumulated_path'] += alpha*rays['refractive_index'].real
    return rays 
    

def evaluate_neighbours(rays, neighbours_idx):
    """For each of a rays neighbours, we need to project the neighbour back
    onto the plane containing the main rays origin, to obtain the (x,y)
    cordinate for the neighbour ray, relative to the main ray origin
    
    rays - a ray_t array of length N
    neighbours_idx - a array on ints of shape (N,6)
    
    returns - a 5-tuple (rays, x,y,dx,dy) where x,y are N*6 arrays for the coordinate of each neighbour.
              dx and dy represent the change in direction of the neighbouring rays (i.e. 
              curvature of the wavefront). The returns rays are the subset of the input rays with
              6 neighbours (i.e. edge-rays are dropped).
              """
    
    mask = (neighbours_idx>=0).all(axis=1)
    nb = neighbours_idx[mask,:]
    n_origin = rays['origin'][nb,:]
    n_direction = rays['direction'][nb,:]
    origin = rays['origin'][mask,None,:]
    direction = rays['direction'][mask,None,:]
    E = rays['E_vector'][mask,None,:] #Assume normalised
    H = numpy.cross(E, direction) #Should also be unit length
    offsets = n_origin - origin
    alpha = -(offsets * direction).sum(axis=-1)/(n_direction * direction).sum(axis=-1)
    projected = (offsets + alpha[:,:,None]*n_direction)
    x = (projected * E).sum(axis=-1)
    y = (projected * H).sum(axis=-1)
    
    dz = (n_direction * direction).sum(axis=-1)
    dx = (n_direction * E).sum(axis=-1)/dz 
    dy = (n_direction * H).sum(axis=-1)/dz
    return (rays[mask], x,y,dx,dy)


def evaluate_neighbours_gc(gausslets):
    ga=gausslets
    origin = ga['base_ray']['origin'][:,None,:]
    direction = ga['base_ray']['direction'][:,None,:]
    n_direction = ga['para_rays']['direction']
    E = ga['base_ray']['E_vector'][:,None,:]
    H = numpy.cross(E, direction)
    
    offsets = ga['para_rays']['origin'] - origin
    #alpha = -(offsets * direction).sum(axis=-1)/(n_direction * direction).sum(axis=-1)
    ### Don't need to project onto origin normal plane as all para-rays start on this plane.
    x = (offsets * E).sum(axis=-1)
    y = (offsets * H).sum(axis=-1)
    dz = (n_direction * direction).sum(axis=-1)
    dx = (n_direction * E).sum(axis=-1)/dz 
    dy = (n_direction * H).sum(axis=-1)/dz
    
    ###This doesn't help
    #     r = x**2 + y**2
    #     r2 = (dx*x + dy*y)/r
    #     dx2 = x*r2
    #     dy2 = y*r2
    return (ga['base_ray'], x, y, dx, dy)
    
    
def evaluate_modes(neighbour_x, neighbour_y, dx, dy, blending=1.0):
    """For each ray in rays, use its nearest neighbours
    to compute the best fit for the Astigmatic Gaussian Beam
    parameters.
    
    rays - an array of ray_t with shape (N,)
    x,y,dx,dy - array of shape (N,6)
    
    For N rays, return a Nx3 complex array of coeffs"""
    global scipy
    if scipy is None:
        global lsqr
        global block_diag
        import scipy
        from scipy.sparse.linalg import lsqr
        from scipy.sparse import block_diag
    
    ### Do linear least squares on each ray and neighbours
    
    x = neighbour_x
    y = neighbour_y
    
    M = block_diag(numpy.dstack((x**2, 2*x*y, y**2)))
    b = numpy.ones(M.shape[0])
    fit = lsqr(M,b)
    im_coefs = fit[0].reshape(x.shape[0],3)

    x = neighbour_x
    y = neighbour_y
    O = numpy.zeros_like(x)
    M_ = numpy.dstack((x,y,O, O,x,y)).reshape(-1,12,3)
    b = numpy.dstack( (dx,dy) ).reshape(-1)
    M = block_diag(M_)
    fit = lsqr(M,b)
    re_coefs = fit[0].reshape(x.shape[0],3)
    
    return ((blending*1j)*im_coefs) + re_coefs


def Gamma(z, A, B, C):
    """Used in testing"""
    detG0 = (A*C) - (B*B)
    denom = (1 + (z*(A+C)) + (z*z)*detG0) #*2
    AA = (A + z*detG0)/denom
    BB = B/denom
    CC = (C + z*detG0)/denom
    return AA, BB, CC


def ExtractGamma(gausslet_collection, blending=1.0):
    """Used in Testing"""
    gc = gausslet_collection.copy_as_array() 
    rays, x, y, dx, dy = evaluate_neighbours_gc(gc)
    modes = evaluate_modes_c(x, y, dx, dy, blending=blending)
    return modes

    
def eval_Efield_from_rays(ray_collection, points, wavelengths, 
                          blending=1.0,
                          exit_pupil_offset=0.0, 
                          exit_pupil_centre=(0.0,0.0,0.0)):
    rays = ray_collection.copy_as_array() 
    radius = exit_pupil_offset
        
    if radius:
        projected = project_to_sphere(rays, exit_pupil_centre, radius)
    else:
        projected = rays
    
    neighbours_idx = ray_collection.neighbours
    rays, x, y, dx, dy = evaluate_neighbours(projected, neighbours_idx)

    modes = evaluate_modes_c(x, y, dx, dy, blending=blending)
    
    _rays = RayCollection.from_array(rays)
    
    E = sum_gaussian_modes(_rays, modes, wavelengths, points)
    
    return E


class EFieldSummation(object):
    """
    For situations where you wish to evaluate the E-field from a set of Gausslets with different sets of evaluation points,
    this class provides a small optimisation by performing the maths to convert ray-intercepts to Gaussian mode parameters
    up front.
    """
    def __init__(self, gausslet_collection, wavelengths=None, blending=1.0 ):
        if wavelengths is None:
            wavelengths = numpy.asarray(gausslet_collection.wavelengths)
        if wavelengths is None:
            raise ValueError("No wavelengths supplied")
        self.wavelengths = wavelengths
        self.gc = gc = gausslet_collection.copy_as_array() 
        rays, x, y, dx, dy = evaluate_neighbours_gc(gc)
        self.modes = evaluate_modes_c(x, y, dx, dy, blending=blending)
        self.base_rays = RayCollection.from_array(rays)
        
    def evaluate(self, points):
        """
        Called to calculate the E-field for the given points.
        """
        points = numpy.ascontiguousarray(points)
        shape = points.shape
        points.shape=(-1,3)
        E = sum_gaussian_modes(self.base_rays, 
                              self.modes, 
                              self.wavelengths, points)
        E.shape = shape
        return E


def eval_Efield_from_gausslets(gausslet_collection, points, 
                               wavelengths = None,
                               blending=1.0, **kwds):
    """
    Calculates the vector E-field is each of the points given. The returned 
    array of field-vectors will have the same length as `points` and 
    has `numpy.complex128` dtype.
    
    :param GaussletCollection gc: The set of Gausslets for which the field should be calculated
    :param ndarray[N,3] points: An array of shape (N,3) giving the points at which the field will be evaluated.
    :param ndarray[] wavelengths: A 1d array containing the wavelengths to be used for the field calculation, 
                                    overriding the wavelengths data contained by the GaussletCollection object.
    :param float blending: The 1/width of each Gaussian mode at the evaluation points. A value of unity (the default),
                            means the parabasal rays are determined to be the 1/e point in the field amplitude.
    """
    gc = gausslet_collection.copy_as_array() 
    if wavelengths is None:
        wavelengths = numpy.asarray(gausslet_collection.wavelengths)
    rays, x, y, dx, dy = evaluate_neighbours_gc(gc)
    modes = evaluate_modes_c(x, y, dx, dy, blending=blending)
    _rays = RayCollection.from_array(rays)
    E = sum_gaussian_modes(_rays, modes, wavelengths, points)
    return E

