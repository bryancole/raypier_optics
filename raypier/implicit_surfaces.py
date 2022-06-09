

from raypier.core import cimplicit_surfs as csurfs
from raypier.core.ctracer import ImplicitSurface
from raypier.core.cimplicit_surfs import Intersection, Union, Difference

from tvtk.api import tvtk

from traits.api import HasStrictTraits, Instance, Tuple, observe, Float


class BaseImplicitSurface(HasStrictTraits):
    c_surface = Instance(ImplicitSurface)
    vtk_impl_func = Instance(tvtk.ImplicitFunction)
    
    
class Plane(BaseImplicitSurface):
    origin = Tuple(0.,0.,0.)
    normal = Tuple(0.,0.,1.)
    
    def _c_surface_default(self):
        return csurfs.Plane(origin=self.origin, normal=self.normal)
    
    def _vtk_impl_func_default(self):
        return tvtk.Plane(origin=self.origin, normal=self.normal)
    
    @observe(["origin","normal"])
    def observe_params(self, ch):
        setattr(self.c_surface, ch.name, ch.new)
        setattr(self.vtk_impl_func, ch.name, ch.new)
        
        
class Sphere(BaseImplicitSurface):
    centre = Tuple(0.,0.,0.)
    radius = Float(1.0)
    
    def _c_surface_default(self):
        return csurfs.Sphere(centre=self.centre, radius=self.radius)
    
    def _vtk_impl_func_default(self):
        return tvtk.Sphere(center=self.centre, radius=self.radius)
    
    @observe("centre")
    def observe_centre(self, ch):
        self.c_surface.centre = ch.new
        self.vtk_impl_func.center = ch.new
        
    @observe("radius")
    def observe_radius(self, ch):
        self.c_surface.radius = ch.new
        self.vtk_impl_func.radius = ch.new
        
        
class Cylinder(BaseImplicitSurface):
    origin = Tuple(0.,0.,0.)
    axis = Tuple(0.,0.,1.)
    radius = Float(10.)
    
    def _c_surface_default(self):
        return csurfs.Cylinder(origin=self.origin, axis=self.axis, radius=self.radius)
    
    def _vtk_impl_func_default(self):
        return tvtk.Cylinder(center=self.origin, axis=self.axis, radius=self.radius)
    
    @observe("origin")
    def observe_centre(self, ch):
        self.c_surface.origin = ch.new
        self.vtk_impl_func.center = ch.new
    
    @observe(["axis", "radius"])
    def observe_params(self, ch):
        setattr(self.c_surface, ch.name, ch.new)
        setattr(self.vtk_impl_func, ch.name, ch.new)
        