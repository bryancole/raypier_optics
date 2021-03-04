

from traits.api import Instance, HasStrictTraits, Property, cached_property,\
        Tuple, Float, Enum, Event, Str, observe, List, ComparisonMode
        
from traitsui.api import View, Item, VGroup, Label, Group, TupleEditor, ListEditor

from tvtk.api import tvtk

import numpy as np

from .core import cshapes
from .core.ctracer import Shape
from .editors import NumEditor


class BaseShape(HasStrictTraits):
    """Abstract base class for all Shape objects.
    """
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
    """Performs a unary NOT operation on the given Shape
    """
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
    """
    Abstract base class for shape boolean operations.
    
    A BooleanShape owns two subshapes, `shape1` and `shape2`.
    These can be combined using any of the logical operations: difference,
    union, intersection and union-of-magnitudes.
    """
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
    """
    Abstract base class for Shape primitives.
    """
    cshape = Instance(Shape)
    impl_func = Instance(tvtk.ImplicitFunction)
    
    bounds = Property()
    
    
    
class CircleShape(BasicShape):
    """
    A Shape representing a circle.
    
    Params
    ------
    
    centre - a 2-tuple giving the local centre coordinate
    
    radius - a float giving the circle radius
    
    """
    centre = Tuple(Float(0.0),Float(0.0))
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
    """
    A Shape representing a rectangle.
    
    Params
    ------
    
    centre - a 2-tuple of floats giving the centre of the rectangle in local coordinates.
    
    width - the width of the rectangle
    
    height - the height of the rectangle
    """
    centre = Tuple(Float(0.0),Float(0.0))
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
        
        
coord_editor = TupleEditor(cols=2, labels=["x", "y"])
    
    
class PolygonShape(BasicShape):
    """
    A general polygon shape outline, defined in terms of a list of 2-tuples for the 
    (x,y) corner points in the local coordinate system of the optic.
    """
    
    #: A list of (x,y) coordinates. There should be at least 3 points in the list.
    #: The shape boundary should be non-self-intersecting.
    coordinates = List(trait=Tuple(Float,Float), comparison_mode=ComparisonMode.identity) #List()
    
    _coords_evt = Event()
    
    #: A property giving the bounding box for the shape.
    bounds = Property(depends_on="_coords_evt")
    
    traits_view = View(Item("coordinates", editor=ListEditor(editor=coord_editor)))
    
    
    def _cshape_default(self):
        return cshapes.PolygonShape(coordinates=np.asarray(self.coordinates))
    
    def _coordinates_default(self):
        return [(-1.,-1.),(0.,1.),(1.,-1.)]
    
    @cached_property
    def _get_bounds(self):
        coords = self.coordinates
        xmin = min(a[0] for a in coords)
        xmax = max(a[0] for a in coords)
        ymin = min(a[1] for a in coords)
        ymax = max(a[1] for a in coords)
        return (xmin,xmax,ymin,ymax)
    
    def _impl_func_default(self):
        points = [(a[0],a[1],0.0) for a in self.coordinates]
        func = tvtk.ImplicitSelectionLoop(loop=points,
                                          normal=(0.,0.,1.))
        return func
    
    @observe("coordinates.items")
    def _update_coords(self, evt):
        self._coords_evt = True 
        coords = self.coordinates
        points = [(a[0],a[1],0.0) for a in coords]
        self.impl_func.loop = points
        self.cshape.coordinates = np.array(coords)
        self.updated=True
        
        
class HexagonShape(PolygonShape):
    """A Hexagon shape, provided as a convenience. Subclasses PolygonShape
    """
    
    #: The centre of the hexagon
    centre = Tuple(Float(0.0),Float(0.0))
    
    #: The distance from the centre to the corners of the hexagon
    radius = Float(25.0)
    
    #: A rotation angle to be applied to the hexagon, in degrees
    rotation = Float(0.0)
    
    traits_view = View(VGroup(
            Item("centre", editor=TupleEditor(cols=2,labels=["x","y"])),
            Item("radius", editor=NumEditor),
            Item("rotation", editor=NumEditor)
            ))
    
    def make_coords(self):
        r = self.radius
        x,y = self.centre
        rot = self.rotation*np.pi/180.0
        angle = (rot + i*np.pi/3 for i in range(6))
        coords = [(r*np.cos(th)+x, r*np.sin(th)+y) for th in angle]
        return coords
    
    def _coordinates_default(self):
        return self.make_coords()
    
    @observe("radius, rotation, centre")
    def _radius_updated(self, evt):
        self.coordinates = self.make_coords()