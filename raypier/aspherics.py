
from traits.api import Float, Str

from .lenses import AsphericLens
from .core.cmaterials import CoatedDispersiveMaterial
from .dispersion import NamedDispersionCurve, NondispersiveCurve


class DispersiveAsphericLens(AsphericLens):
    glass_type = Str("N-BK7")
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
    A12 = -5.920460E-6
    A14 = 9.104334E-8
    A16 = 1.291935E-8
    
    B_curvature = 3.200655
    B_conic = -4.321649
    B4 = 5.521153E-3
    B6 = -1.981378E-3
    B8 = 4.782553E-4
    B10 = -7.328134E-5
    B12 = 5.920460E-6
    B14 = -9.104334E-8
    B16 = -1.291935E-8
    glass_type = "D-ZLAF52LA"
    
    
class Thorlabs355440_B(AsphericLens):
    name = "Thorlabs#355440"
    diameter = 4.7
    CT = 3.8
    A_curvature = -3.200655
    A_conic = -4.321649
    A4 = -5.521153E-3
    A6 = 1.981378E-3
    A8 = -4.782553E-4
    A10 = 7.328134E-5
    A12 = -5.920460E-6
    A14 = 9.104334E-8
    A16 = 1.291935E-8
    
    B_curvature = 3.200655
    B_conic = -4.321649
    B4 = 5.521153E-3
    B6 = -1.981378E-3
    B8 = 4.782553E-4
    B10 = -7.328134E-5
    B12 = 5.920460E-6
    B14 = -9.104334E-8
    B16 = -1.291935E-8
    
    coating_thickness = Float(0.283) #microns
    coating_refractive_index = Float(1.37)
    
    def _material_default(self):
        #glass = NondispersiveCurve(1.786) #980nm
        glass = NondispersiveCurve(1.7935) #780nm
        cri = self.coating_refractive_index
        coating = NondispersiveCurve(refractive_index=cri)
        air = NondispersiveCurve(1.0)
        
        m = CoatedDispersiveMaterial(dispersion_inside=glass,
                                     dispersion_outside=air,
                                     dispersion_coating=coating)
        return m
    
    
class Thorlabs352440(AsphericLens):
    """This is an obsolete lens but one we have a lot of experience with which makes it a good test-piece"""
    name = "Thorlabs#352440"
    diameter = 4.7
    CT = 4.06
    A_curvature = -2.40
    A_conic = -4.456395
    A4 = -5.2675614E-03
    A6 = 1.0753165E-03
    A8 = -3.2335319E-04
    A10 = 2.5517886E-05
    B_curvature = 2.39
    B_conic = -0.8371770
    B4 = -7.3100129E-03
    B6 = -3.9100155E-04
    B8 = 1.2597638E-04
    B10 = -2.6451531E-05
    
    coating_thickness = Float(0.283) #microns
    coating_refractive_index = Float(1.37)
    
    ###glass_type = "ECO-550"
    
    def _material_default(self):
        ### ECO-550 is not in the database :(
        glass = NondispersiveCurve(1.594) #Value for 780nm
        cri = self.coating_refractive_index
        coating = NondispersiveCurve(refractive_index=cri)
        air = NondispersiveCurve(1.0)
        
        m = CoatedDispersiveMaterial(dispersion_inside=glass,
                                     dispersion_outside=air,
                                     dispersion_coating=coating)
        return m
    
    
    
    