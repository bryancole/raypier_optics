
from traits.api import Float, Str

from .lenses import AsphericLens
from .cmaterials import CoatedDispersiveMaterial
from .dispersion import NamedDispersionCurve, NondispersiveCurve


class DispersiveAsphericLens(AsphericLens):
    glass_type = Str("D-ZLAF52LA")
    coating_thickness = Float(0.283) #microns
    coating_refractive_index = Float(1.37)
    
    def _material_default(self):
        glass = NamedDispersionCurve(self.glass_type)
        cri = self.coating_refractive_index
        coating = NondispersiveCurve(refractive_index=cri)
        air = NondispersiveCurve(1.0)
        
        m = CoatedDispersiveMaterial(dispersion_inside=glass,
                                     dispersion_outside=air,
                                     dispersion_coating=coating)
        return m


class Thorlabs355440(DispersiveAsphericLens):
    name = "Thorlabs#355440"
    diameter = 4.7
    CT = 3.8
    A_curvature = -3.200655
    A_conic = -4.321649
    A4 = -5.521153E-3
    A6 = 1.981378E-3
    A8 = -4.782553E-4
    A10 = 7.328134E-5
    B_curvature = 3.200655
    B_conic = -4.321649
    B4 = 5.521153E-3
    B6 = -1.981378E-3
    B8 = 4.782553E-4
    B10 = -7.328134E-5
    glass_type = "D-ZLAF52LA"
    
    
    
    