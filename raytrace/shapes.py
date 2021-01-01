

from traits.api import Instance, HasTraits, Property, on_trait_change, cached_property,\
        Tuple, Float, Enum

from tvtk.api import tvtk

from .core import cshapes
from .core.ctracer import Shape


class BaseShape(HasTraits):
    cshape = Instance(Shape)
    
    def __and__(self, other):
        return LogicalAND(self, other)
    
    def __or__(self, other):
        return LogicalOR(self, other)
    
    def __xor__(self, other):
        return LogicalXOR(self, other)
    
    def __invert__(self):
        return InvertedShape(self)
    
    
class InvertedShape(BaseShape):
    shape = Instance(BaseShape)
    cshape = Property(depends_on="shape")
    impl_func = Property(tvtk.ImplicitFunction, depends_on="shape")
    
    @cached_property
    def _get_cshape(self):
        return cshapes.InvertShape(self.shape)
        
    @cached_property
    def _get_impl_func(self):
        ### Can't figure out how to invert a VTK implicit function! crap.
        raise NotImplementedError
    

class BooleanShape(BaseShape):
    shape1 = Instance(BaseShape)
    shape2 = Instance(BaseShape)
    
    cshape = Property(depends_on="shape1, shape2")
    impl_func = Property(depends_on="shape1, shape2, operation_type")
    
    operation_type = Enum("union", "difference", "intersection", "union_of_magnitudes")
    cshape_type = Enum(cshapes.BooleanAND, cshapes.BooleanOR, cshapes.BooleanXOR, cshapes.InvertShape)
    
    def __init__(self, shape1, shape2):
        self.shape1=shape1
        self.shape2=shape2
        
    @cached_property
    def _get_cshape(self):
        return self.cshape_type(self.shape1.cshape,
                                  self.shape2.cshape)
        
    @cached_property
    def _get_impl_func(self):
        func = tvtk.ImplicitBoolean()
        func.add_function(self.shape1.impl_func)
        func.add_function(self.shape2.impl_func)
        func.operation_type = self.operation_type
        return func


class LogicalAND(BooleanShape):
    """Meaning inside A AND inside B"""
    operation_type = "intersection"
    cshape_type = cshapes.BooleanAND
        
        
class LogicalOR(BooleanShape):
    """Inside A OR inside B"""
    operation_type = "union"
    cshape_type = cshapes.BooleanOR
        
        
class LogicalXOR(BooleanShape):
    """(Inside A and outside B) or (Inside B and outside A)"""
    operation_type = "difference"
    cshape_type = cshapes.BooleanXOR
        
        
class BasicShape(BaseShape):
    centre = Tuple(Float(0.0),Float(0.0))
    cshape = Property(depends_on="centre")
    
    
class CircleShape(BasicShape):
    radius = Float(1.0)
    cshape = Property(depends_on="centre, radius")
    impl_func = Instance(tvtk.ImplicitFunction)
    
    @cached_property
    def _get_cshape(self):
        return cshapes.CircleShape(centre=self.centre, radius=self.radius)
    
    def _impl_func_default(self):
        x,y = self.centre
        func = tvtk.Cylinder(axis=(0,0,1), center=(x, y, 0),
                             radius=self.radius)
        return func
    
    @on_trait_change("radius, centre")
    def _update_imp_func(self):
        x,y = self.centre
        func = self.impl_func
        func.radius = self.radius
        func.center = (x,y,0)
        
        

    