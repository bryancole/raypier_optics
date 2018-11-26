'''
Support for achromat lenses 

Created on 23 Nov 2018

@author: bryan
'''
from traits.api import Float, Instance, on_trait_change

from raytrace.lenses import BaseLens
from raytrace.dispersion import BaseDispersionCurve, NamedDispersionCurve, \
        NondispersiveCurve
from raytrace.cmaterials import CoatedDispersiveMaterial
from raytrace.cfaces import SphericalFace
from raytrace.ctracer import FaceList
from raytrace.custom_sources import EmptyGridSource

from tvtk.api import tvtk
from numpy import sqrt


class Doublet(BaseLens):
    abstract = False
    name = "Achromatic Doublet"
    
    CT1 = Float(6.0, desc="first centre thickness")
    CT2 = Float(4.0, desc="second centre thickness")
    diameter = Float(25.0)
    curvature1 = Float(43.96, desc="radius of curvature for first outer face")
    curvature2 = Float(-42.90, desc="radius of curvature for the inner face")
    curvature3 = Float(-392.21, desc="radius of curvature for the second outer face")
    dispersion1 = Instance(BaseDispersionCurve, desc="Dispersion curve for first material")
    dispersion2 = Instance(BaseDispersionCurve, desc="Dispersion curve for second material")
    dispersion_coating = Instance(BaseDispersionCurve, desc="Dispersion curve for the coating")
    coating_thickness = Float(0.25, desc="Thickness of the AR coating, in microns")
    
    _material1 = Instance(CoatedDispersiveMaterial, ())
    _material2 = Instance(CoatedDispersiveMaterial, ())
    _material3 = Instance(CoatedDispersiveMaterial, ())
    
    vtk_grid = Instance(EmptyGridSource, ())
    vtk_cylinder = Instance(tvtk.Cylinder, ())
    vtk_sphere1 = Instance(tvtk.Sphere, ())
    vtk_sphere2 = Instance(tvtk.Sphere, ())
    vtk_sphere3 = Instance(tvtk.Sphere, ())
    
    clip2 = Instance(tvtk.ClipDataSet, ())
    
    @on_trait_change("dispersion1, dispersion2, dispersion_coating, coating_thickness")
    def on_materials_changed(self):
        self.apply_materials()
        self.update = True
    
    def apply_materials(self):
        air = NondispersiveCurve(1.0)
        self._material1.dispersion_inside = self.dispersion1
        self._material1.dispersion_outside = air
        self._material1.dispersion_coating = self.dispersion_coating
        self._material1.coating_thickness = self.coating_thickness
        
        self._material2.dispersion_inside = self.dispersion1
        self._material2.dispersion_outside = self.dispersion2
        self._material2.dispersion_coating = self.dispersion_coating
        self._material2.coating_thickness = 0.0
        
        self._material3.dispersion_inside = self.dispersion2
        self._material3.dispersion_outside = air
        self._material3.dispersion_coating = self.dispersion_coating
        self._material3.coating_thickness = self.coating_thickness
        
    @on_trait_change("diameter, CT1, CT2, curvature1, curvature2, curvature3")
    def on_geometry_changed(self):
        face1, face2, face3 = self.faces.faces
        face1.diameter = self.diameter
        face1.z_height = self.CT1
        face1.curvature = self.curvature1
        
        face2.diameter = self.diameter
        face2.curvature = self.curvature2
        
        face3.diameter = self.diameter
        face3.z_height = self.CT2
        face3.curvature = self.curvature3
        self.update = True
    
    def _faces_default(self):
        self.on_apply_materials()
        fl = FaceList(owner=self)
        fl.faces = [ 
                  SphericalFace(owner=self, 
                                diameter=self.diameter,
                                material=self._material1,
                                z_height=self.CT1, 
                                curvature=self.curvature1),
                  SphericalFace(owner=self, 
                                diameter=self.diameter,
                                material=self._material2,
                                z_height=0.0, 
                                curvature=self.curvature2),
                  SphericalFace(owner=self, 
                                diameter=self.diameter,
                                material=self._material3,
                                z_height=-self.CT2, 
                                curvature=self.curvature3),
                ]
        return fl
    
    @on_trait_change("CT1, CT2, diameter, curvature1, curvature3, curvature2")
    def config_pipeline(self):
        ct1 = self.CT1
        ct2 = self.CT2
        rad = self.diameter/2
        curve1 = self.curvature1
        curve2 = self.curvature2
        curve3 = self.curvature3
        
        size = 41
        spacing = 2*rad / (size-1)
        if curve1 >= 0.0:
            extra1=0.0
        else:
            extra1 = -curve1 - sqrt(curve1**2 - rad**2)
            
        if curve3 <= 0.0:
            extra2=0.0
        else:
            extra2 = curve2 - sqrt(curve2**2 - rad**2)
            
        lsize1 = int((ct1+extra1)/spacing) + 2
        lsize2 = int((ct2+extra2)/spacing) + 2
        
        grid = self.vtk_grid
        grid.dimensions = (size,size,lsize1+lsize2)
        grid.origin = (-rad, -rad, -lsize2*spacing)
        grid.spacing = (spacing, spacing, spacing)
                      
        cyl = self.vtk_cylinder
        cyl.center = (0,0,0)
        cyl.radius = rad
        
        s1 = self.vtk_sphere1
        s1.center = (0,0,ct1 - curve1)
        s1.radius = abs(curve1)
        
        s3 = self.vtk_sphere3
        s3.center = (0,0,ct2 - curve2)
        
        self.clip2.inside_out = bool(curve >= 0.0)
        
        self.vtk_grid.modified()
        self.update=True
        
                                         
    def _pipeline_default(self):
        grid = self.vtk_grid
        #grid.set_execute_method(self.create_grid)
        grid.modified()
        
        trans = tvtk.Transform()
        trans.rotate_x(90.)
        cyl = self.vtk_cylinder
        cyl.transform = trans
        
        clip1 = tvtk.ClipVolume(input_connection=grid.output_port,
                                 clip_function=self.vtk_cylinder,
                                 inside_out=1)
        
        self.clip2.set(input_connection = clip1.output_port,
                      clip_function=self.vtk_sphere,
                      inside_out=1)
        
        topoly = tvtk.GeometryFilter(input_connection=self.clip2.output_port)
        norms = tvtk.PolyDataNormals(input_connection=topoly.output_port)
        
        transF = tvtk.TransformFilter(input_connection=norms.output_port, 
                                      transform=self.transform)
        self.config_pipeline()
        grid.modified()
        return transF
    
    
class EdmundOptic45805(Doublet):
    name = "Edumnd Optics Achromat#45805"
    CT1 = 6.0
    CT2 = 4.0
    diameter = 25.4
    offset = 0.0
    curvature1 = 43.96
    curvature2 = -42.90
    curvature3 = -392.21
    dispersion1 = NamedDispersionCurve("N-LAK22")
    dispersion2 = NamedDispersionCurve("N-SF6")
    dispersion_coating = NondispersiveCurve(refractive_index=1.37)
    coating_thickness = 0.283 #microns