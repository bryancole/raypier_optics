"""
Functions relating to the evaluation of the optical E-field by
summation of General Astigmatic Gaussian Beams
"""

from traits.api import on_trait_change, Float, Instance,Event, Int,\
        Property, Button, Str, Array

from traitsui.api import View, Item, VGroup, DropEditor

from tvtk.api import tvtk

from raytrace.bases import Probe, Traceable, NumEditor, Vector
from raytrace.sources import RayCollection, BaseRaySource

import numpy


def project_to_sphere(rays, centre=(0,0,0), radius=10.0):
    """project the given set of rays back to their intercept 
    with a sphere at the given centre and radius.
    Returns a new RayCollection"""
    
    
def evaluate_modes(rays, neighbours):
    """For each ray in rays, use its nearest neighbours
    to compute the best fit for the Astigmatic Gaussian Beam
    parameters.
    For N rays, return a Nx3 complex array of coeffs"""
    ### Do linear least squares on each ray and neighbours
    
    
def evaluate__E(rays, mode_coeffs, points):
    """Evaluate the E-field at each position in the
    _points_ array.
    """
    
    
class EFieldPlane(Probe):    
    width = Float(0.5) #in mm
    size = Int(30)
            
    ray_source = Instance(klass="raytrace.sources.BaseRaySource")
    
    exit_pupil_offset = Float(10.0) #in mm
    
    ###The output of the probe
    E_field = Array()
    
    _eval_btn = Button("calculate")
    
    update = Event() #request re-tracing
    
    _plane_src = Instance(tvtk.PlaneSource, (), 
                    {"x_resolution":1, "y_resolution":1},
                    transient=True)
                    
    traits_view = View(VGroup(
                       Traceable.uigroup,
                       Item('size', editor=NumEditor),
                       Item('width', editor=NumEditor),
                        ),
                   )
    
        
    def __eval_btn_fired(self):
        source = self.ray_source
        
        pass
    
    @on_trait_change("size, width")
    def config_pipeline(self):
        src = self._plane_src
        size = self.size
        side = self.width/2.
        
        src.origin = (-side,-side,0)
        src.point1 = (-side,side,0)
        src.point2 = (side,-side,0)
        
        src.x_resolution = size
        src.y_resolution = size
        
        self.update = True
        
    @on_trait_change("centre, direction, orientation, exit_pupil_offset")
    def on_change(self):
        self.update = True
    
    def _actors_default(self):
        source = self._plane_src
        trans_f = tvtk.TransformFilter(input_connection=source.output_port,
                        transform = self.transform)
        map = tvtk.PolyDataMapper(input_connection=trans_f.output_port)
        act = tvtk.Actor(mapper=map)
        actors = tvtk.ActorCollection()
        actors.append(act)
        return actors
        
    def evaluate(self, ray_src):
        c = 2.99792458e8 * 1e-9 #convert to mm/ps
        traced_rays = ray_src.TracedRays
        all_wavelengths = numpy.asarray(ray_src.wavelength_list)
        all_rays = [r.copy_as_array() for r in traced_rays]
        neighbours = list(ray_src.iter_neighbours())
        
        intersections = [self.intersect_plane(rays) for rays in all_rays]
            
    def intersect_plane(self, rays):
        """
        @param rays: a numpy array of ray_t dtype
        
        intersect rays with plane, returning the indices of the intersecting rays
        """
        #raise NotImplementedError
        pass
    
    
    