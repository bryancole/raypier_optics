

from traits.api import Instance, HasStrictTraits, Property, cached_property,\
        Tuple, Float, Enum, Event, Str, observe
        
from traitsui.api import View, Item, VGroup, Label, Group

from tvtk.api import tvtk

from .core import cshapes
from .core.ctracer import Shape
from .editors import NumEditor


class BaseShape(HasStrictTraits):
    updated = Event()
    
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
    ###Inputs
    shape1 = Instance(BaseShape)
    shape2 = Instance(BaseShape)
    operation_type = Enum("difference", "union", "intersection", "union_of_magnitudes")
    
    ###Outputs
    bounds = Property(depends_on="shape1.bounds, shape2.bounds")
    cshape = Property(depends_on="shape1, shape2, operation_type")
    impl_func = Property(depends_on="shape1, shape2, operation_type")
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
    
    def __init__(self, shape1, shape2, **kwds):
        operation = kwds.pop("operation", "difference")
        super().__init__(**kwds)
        self.operation_type = operation
        self.shape1=shape1
        self.shape2=shape2
        
    @observe("shape1.updated, shape2.updated, operation_type")
    def _on_changed(self, evt):
        self.updated = True
        
    @cached_property
    def _get_bounds(self):
        b1 = self.shape1.bounds
        b2 = self.shape2.bounds
        return ( min(b1[0],b2[0]), max(b1[1],b2[1]), min(b1[2],b2[2]), max(b1[3],b2[3]) )
        
    @cached_property
    def _get_cshape(self):
        return self.cshape_map[self.operation_type](self.shape1.cshape,
                                  self.shape2.cshape)
        
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
    
    bounds = Property()
    
    traits_view = View(Item("centre"))
    
    
class CircleShape(BasicShape):
    radius = Float(1.0)
    
    bounds = Property(depends_on="radius, centre")
    
    traits_view = View(
            Item("centre"),
            Item("radius", editor=NumEditor))
    
    @cached_property
    def _get_bounds(self):
        x,y = self.centre
        r = self.radius
        return (x-r,x+r,y-r,y+r)
    
    def _cshape_default(self):
        return cshapes.CircleShape(centre=self.centre, radius=self.radius)
    
    def _impl_func_default(self):
        x,y = self.centre
        func = tvtk.Cylinder(axis=(0,0,1), center=(x, y, 0),
                             radius=self.radius)
        return func
    
    @observe("radius, centre")
    def _update_imp_func(self, evt):
        x,y = self.centre
        func = self.impl_func
        func.radius = self.radius
        func.center = (x,y,0)
        self.cshape.radius = self.radius
        self.cshape.centre = (x,y)
        self.updated = True
        
        
class RectangleShape(BasicShape):
    width = Float(5.0)
    height = Float(7.0)
    
    bounds = Property(depends_on="width, height, centre")
    
    _zextent = Float(10000.)
    
    traits_view = View(VGroup(Item("centre"),
                            Item("width", editor=NumEditor),
                              Item("height", editor=NumEditor)))
    
    @cached_property
    def _get_bounds(self):
        x,y = self.centre
        w=self.width/2.
        h = self.height/2.
        return (x-w,x+w,y-h,y+h)

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
    
    @observe("width, height, centre")
    def _update_imp_func(self, evt):
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
        self.updated = True
    