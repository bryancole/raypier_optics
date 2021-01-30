
"""
Traited wrappers for various dispersion-curve definitions
"""


from .dispersion import NondispersiveCurve, NamedDispersionCurve

from traits.api import HasStrictTraits, Instance, Float, Str, Enum, Property, cached_property, Bool

from traitsui.api import View, VGroup, Item

from .editors import NumEditor


GlassNames = sorted(NamedDispersionCurve.get_glass_names())



class OpticalMaterial(HasStrictTraits):
    name = Str("Part A")
    from_database = Bool(True)
    refractive_index = Float(1.0)
    glass_name=Enum(GlassNames)
    absorption=Float(0.0)
    
    dispersion_curve = Property(depends_on="glass_name, absorption, refractive_index, from_database")
    
    traits_view = View(VGroup(
            Item("from_database"),
            Item("glass_name", style="simple", enabled_when="from_database"),
            Item("refractive_index", editor=NumEditor, enabled_when="not from_database"),
            Item("absorption", editor=NumEditor),
        ))
    
    def __init__(self, *args, **kwds):
        if "from_database" not in kwds:
            if ("refractive_index" in kwds) and ("glass_name" not in kwds):
                kwds['from_database'] = False
        return super().__init__(*args, **kwds)
    
    @cached_property
    def _get_dispersion_curve(self):
        if self.from_database:
            return NamedDispersionCurve(self.glass_name, absorption=self.absorption)
        else:
            return NondispersiveCurve(refractive_index=self.refractive_index,
                                  absorption=self.absorption)
    
    
class AirMaterial(OpticalMaterial):
    from_database = False
    name = "Air"
    
    
air = AirMaterial()