

from raypier.core import cimplicit_surfs as csurfs

from tvtk.api import tvtk

from traits.api import HasStrictTraits, Instance, Tuple, observe, Float


class BaseImplicitSurface(HasStrictTraits):
    c_surface = Instance(csurfs.ImplicitSurface)
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
        setattr(self.c_surface, ch.name, ch.value)
        setattr(self.vtk_impl_func, ch.name, ch.value)
        
        
class Sphere(BaseImplicitSurface):
    centre = Tuple(0.,0.,0.)
    radius = Float(1.0)
    
    def _c_surface_default(self):
        return csurfs.Sphere(centre=self.centre, radius=self.radius)
    
    def _vtk_impl_func_default(self):
        return tvtk.Sphere(center=self.centre, radius=self.radius)
    
    @observe("centre")
    def observe_centre(self, ch):
        self.c_surface.centre = ch.value
        self.vtk_impl_func.center = ch.value
        
    @observe("radius")
    def observe_radius(self, ch):
        self.c_surface.radius = ch.value
        self.vtk_impl_func.radius = ch.value
        
        
class Cylinder(BaseImplicitSurface):
    origin = Tuple(0.,0.,0.)
    axis = Tuple(0.,0.,1.)
    radius = Float(10.)
    
    def _c_surface_default(self):
        return csurfs.Cylinder(origin=self.origin, axis=self.axis, radius=self.radius)
    
    def _vtk_impl_func_default(self):
        return tvtk.Cylinder(origin=self.origin, axis=self.axis, radius=self.radius)
    
    @observe(["origin", "axis", "radius"])
    def observe_params(self, ch):
        setattr(self.c_surface, ch.name, ch.value)
        setattr(self.vtk_impl_func, ch.name, ch.value)
        