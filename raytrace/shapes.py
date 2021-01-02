

from traits.api import Instance, HasTraits, Property, on_trait_change, cached_property,\
        Tuple, Float, Enum, observe
        
from traitsui.api import View, Item, VGroup, Label, Group

from tvtk.api import tvtk

from .core import cshapes
from .core.ctracer import Shape
from .editors import NumEditor


class BaseShape(HasTraits):
    cshape = Instance(Shape)
    
    def __and__(self, other):
        return BooleanShape(self, other, operation="intersection")
    
    def __or__(self, other):
        return BooleanShape(self, other, operation="union")
    
    def __xor__(self, other):
        return BooleanShape(self, other, operation="difference")
    
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
    
    cshape = Property(depends_on="shape1, shape2, cshape_type")
    impl_func = Property(depends_on="shape1, shape2, operation_type")
    
    operation_type = Enum("difference", "union", "intersection", "union_of_magnitudes")
    cshape_type = Enum(cshapes.BooleanAND, cshapes.BooleanOR, cshapes.BooleanXOR, cshapes.InvertShape)
    
    traits_view = View(VGroup(
        Group(Item("shape1", style="custom", show_label=False, padding=-5),
              show_border=True,
              label="Shape 1")
              ,
        Item("operation_type", style="simple", show_label=False),
        Group(Item("shape2", style="custom", show_label=False, padding=-5),
              show_border=True,
              label="Shape 2"),
        padding = 0,
        show_border=True,
        label="Boolean Combination"
        ),
    )
    
    cshape_map = {"union": cshapes.BooleanOR, 
                      "difference": cshapes.BooleanXOR,
                      "union_of_magnitudes": cshapes.BooleanOR, 
                      "intersection": cshapes.BooleanAND}
    
    def __init__(self, shape1, shape2, operation="difference"):
        self.operation_type = operation
        self.shape1=shape1
        self.shape2=shape2
        
    @cached_property
    def _get_cshape(self):
        return self.cshape_type(self.shape1.cshape,
                                  self.shape2.cshape)
        
    @observe("operation_type")
    def on_operation_changed(self, evt):
        self.cshape_type = self.cshape_map[evt.new]
        
    @cached_property
    def _get_impl_func(self):
        func = tvtk.ImplicitBoolean()
        func.add_function(self.shape1.impl_func)
        func.add_function(self.shape2.impl_func)
        func.operation_type = self.operation_type
        return func
        
        
class BasicShape(BaseShape):
    centre = Tuple(Float(0.0),Float(0.0))
    cshape = Instance(Shape)
    impl_func = Instance(tvtk.ImplicitFunction)
    
    traits_view = View(Item("centre"))
    
    
class CircleShape(BasicShape):
    radius = Float(1.0)
    
    traits_view = View(
            Item("centre"),
            Item("radius", editor=NumEditor))
    
    def _cshape_default(self):
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
        self.cshape.radius = self.radius
        self.cshape.centre = (x,y)
        
        
class RectangleShape(BasicShape):
    width = Float(5.0)
    height = Float(7.0)
    
    _zextent = Float(10000.)
    
    traits_view = View(VGroup(Item("centre"),
                            Item("width", editor=NumEditor),
                              Item("height", editor=NumEditor)))

    def _cshape_default(self):
        return cshapes.RectangleShape(centre=self.centre,
                                      width=self.width,
                                      height=self.height)
        
    def _impl_func_default(self):
        x,y = self.centre
        w=self.width/2.
        h = self.height/2.
        z = self._zextent
        func = tvtk.Box(bounds=(x-w,x+w,y-h,y+y,-z,z))
        return func
    
    @on_trait_change("width, height, centre")
    def _update_imp_func(self):
        func = self.impl_func
        x,y = self.centre
        w=self.width/2.
        h = self.height/2.
        z = self._zextent
        func.bounds = (x-w,x+w,y-h,y+h,-z,z)
        cshape = self.cshape
        cshape.width = self.width
        cshape.height = self.height
        cshape.centre = (x,y)
        
    