"""
Functions relating to the evaluation of the optical E-field by
summation of General Astigmatic Gaussian Beams
"""

from traits.api import on_trait_change, Float, Instance,Event, Int,\
        Property, Str, Array, cached_property, List, Bool, observe, Button

from traitsui.api import View, Item, VGroup, Tabbed

from tvtk.api import tvtk

from raypier.bases import Probe, Traceable, NumEditor
from raypier.probes import BaseCapturePlane
from raypier.sources import RayCollection, BaseRaySource
from raypier.core.cfields import sum_gaussian_modes, evaluate_modes as evaluate_modes_c
from raypier.find_focus import find_ray_focus
from raypier.utils import normaliseVector, dotprod
from raypier.core.ctracer import GaussletCollection, RayCollection

import numpy
from raypier.editors import IntEditor

scipy = None
    
import time
import traceback


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


def eval_Efield_from_gausslets(gausslet_collection, points, wavelengths,
                               blending=1.0, **kwds):
    gc = gausslet_collection.copy_as_array() 
    rays, x, y, dx, dy = evaluate_neighbours_gc(gc)
    modes = evaluate_modes_c(x, y, dx, dy, blending=blending)
    _rays = RayCollection.from_array(rays)
    E = sum_gaussian_modes(_rays, modes, wavelengths, points)
    return E
    
    
class EFieldPlane(Probe):
    name = Str("E-Field Probe")
    detector = Instance(BaseCapturePlane)
    align_detector = Bool(True)
    
    gen_idx = Int(-1)
    width = Float(0.5) #in mm
    height = Float(0.5)
    size = Int(30)
            
    exit_pupil_offset = Float(10.0) #in mm
    blending = Float(1.0)
    
    centre_on_focus_btn = Button()
    
    ###The output of the probe
    E_field = Array()
    intensity = Property(Array, depends_on="E_field")
    total_power = Property(depends_on="intensity, width, height, size")
    
    ### The mean real part of the refractive index for the material where the probe plane lies
    ### Get get this from the selected/captured rays.
    refractive_index = Float(1.0)
        
    _mtime = Float(0.0)
    _src_list = List()
    _plane_src = Instance(tvtk.PlaneSource, (), 
                    {"x_resolution":1, "y_resolution":1},
                    transient=True)
    
    _attrib = Instance(tvtk.ProgrammableAttributeDataFilter,(), transient=True)
                    
    traits_view = View(Tabbed(
                        VGroup(
                       Traceable.uigroup,
                       Item("total_power", style="readonly"),
                       Item('size', editor=IntEditor),
                       Item('width', editor=NumEditor),
                       Item('height', editor=NumEditor),
                       Item('exit_pupil_offset', editor=NumEditor),
                       Item('blending', editor=NumEditor),
                       Item('gen_idx', editor=IntEditor),
                       Item('centre_on_focus_btn', show_label=False, label="Centre on focus")
                   )))
    
    @observe("centre")
    def on_move(self, evt):
        detector = self.detector
        if detector is not None and self.align_detector:
            detector.centre = evt.new
        self._mtime = 0.0
        self.on_change()
    
    @on_trait_change("orientation, size, width, height, exit_pupil_offset, blending, gen_idx")
    def config_pipeline(self):
        src = self._plane_src
        size = self.size
        side = self.width/2.
        yside = self.height/2
        
        src.origin = (-side,-yside,0)
        src.point1 = (-side,yside,0)
        src.point2 = (side,-yside,0)
        
        src.x_resolution = size
        src.y_resolution = size
        
        self._mtime = 0.0
        self.on_change()
        
    @on_trait_change("update")
    def on_change(self):
        try:
            self.evaluate(None)
        except:
            traceback.print_exc()
    
    def _actors_default(self):
        source = self._plane_src
        attr = self._attrib
        attr.input_connection = source.output_port
        
        trans_f = tvtk.TransformFilter(input_connection=attr.output_port,
                        transform = self.transform)
        map = tvtk.PolyDataMapper(input_connection=trans_f.output_port)
        
        def execute():
            dataobj = attr.poly_data_output
            data = self.intensity.T
            dataobj.cell_data.scalars = data.ravel()
            map.scalar_range = (data.min(), data.max()) 
        attr.set_execute_method(execute)
        
        act = tvtk.Actor(mapper=map)
        actors = tvtk.ActorCollection()
        actors.append(act)
        return actors
    
    @cached_property
    def _get_intensity(self):        
        E = self.E_field
        return (E.real**2).sum(axis=-1) + (E.imag**2).sum(axis=-1)
    
    def _get_total_power(self):
        power = self.intensity.sum()*(self.width*self.height/(self.size**2))
        power *= self.refractive_index
        return power
        
    def evaluate(self, src_list):
        mtime = self._mtime
        detector = self.detector
        if src_list is None:
            src_list = self._src_list
        else:
            self._src_list = src_list
        
        start = time.monotonic()
            
        if detector is not None:
            detector.evaluate(src_list)
            if detector.captured is None:
                return
            if mtime > detector._mtime:
                return
            ray_list = [detector.captured,]
        else:
            idx = self.gen_idx
            ray_list = [src.TracedRays[idx] for src in src_list if src.TracedRays]
        
        n_list = []
        for ct, rays in enumerate(ray_list):
            wavelengths = rays.wavelengths
            
            size = self.size
            side = self.width/2.
            yside = self.height/2
            px = numpy.linspace(-side,side,size)
            py = numpy.linspace(-yside, yside, size)
            points = numpy.dstack(numpy.meshgrid(px,py,0)).reshape(-1,3)
            ###What a chore
            trns = self.transform
            pts_in = tvtk.Points()
            pts_out = tvtk.Points()
            pts_in.from_array(points)
            trns.transform_points(pts_in, pts_out)
            points2 = pts_out.to_array().astype('d')
            
            if isinstance(rays, GaussletCollection):
                n_list.append(rays.base_rays.refractive_index.real)
                E = eval_Efield_from_gausslets(rays, points2, wavelengths, 
                                               blending=self.blending)
            else:
                n_list.append(rays.refractive_index.real)
                E = eval_Efield_from_rays(rays, points2, wavelengths, 
                                          blending=self.blending,
                                          exit_pupil_offset=self.exit_pupil_offset,
                                          exit_pupil_centre=self.centre)
            
            if ct==0:
                self.E_field = E.reshape(self.size, self.size, 3)
            else:
                self.E_field += E.reshape(self.size, self.size, 3)
                
        self.refractive_index = numpy.concatenate(n_list).mean()
        
        self._attrib.modified()
        end = time.monotonic()
        self._mtime = end
        print(f"Field calculation took: {end-start} s")
        
            
    def intersect_plane(self, rays):
        """
        @param rays: a numpy array of ray_t dtype
        
        intersect rays with plane, returning the indices of the intersecting rays
        """
        #raise NotImplementedError
        pass
    
    @observe("centre_on_focus_btn")
    def do_centre_on_focus(self, evt):
        self.centre_on_focus()
    
    def centre_on_focus(self, idx=-1, offset=(0,0,0)):
        """
        Locate the centre of the field-probe on the point of closest approach 
        for the indicated RayCollection (by default, the last one).
        """
        detector = self.detector
        if detector is None:
            src = self.source
            rays = src.TracedRays[idx]
        else:
            rays = detector.captured
        self.centre = tuple(a+b for a,b in zip(find_ray_focus(rays), offset))
    
    
    