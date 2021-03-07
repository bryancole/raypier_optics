'''
Support for achromat lenses 

Created on 23 Nov 2018

:author: bryan
'''
from traits.api import Float, Instance, on_trait_change
from traitsui.api import View, Item, VGroup

from raypier.bases import NumEditor, Traceable
from raypier.lenses import BaseLens
from raypier.dispersion import BaseDispersionCurve, NamedDispersionCurve, \
        NondispersiveCurve
from raypier.core.cmaterials import CoatedDispersiveMaterial
from raypier.core.cfaces import SphericalFace
from raypier.core.ctracer import FaceList
from raypier.vtk_algorithms import EmptyGridSource

from tvtk.api import tvtk
from numpy import sqrt


class Singlet(BaseLens):
    """
    A dispersive singlet lens with spherical surfaces. Not really achromatic at all.
    """
    abstract = False
    name = "sperical singlet"
    
    #: Lens diameter in mm
    diameter = Float(25.0)
    
    #: Center thickness in mm
    CT = Float(6.0, desc="first centre thickness")
    
    #: Radius of curvature for the first surface, in mm
    curvature1 = Float(43.96, desc="radius of curvature for first outer face")
    
    #: Radius of curvature for the second surface, in mm
    curvature2 = Float(-42.90, desc="radius of curvature for the inner face")
    
    #: Dispersion curve for the lens bulk :type raypier.core.cmaterials.BaseDispersionCurve:
    dispersion = Instance(BaseDispersionCurve, desc="Dispersion curve for the glass")
    
    #: Dispersion curve for the lens coating :type raypier.core.cmaterials.BaseDispersionCurve:
    dispersion_coating = Instance(BaseDispersionCurve, desc="Dispersion curve for the coating")
    
    #: Coating thickness, in microns
    coating_thickness = Float(0.25, desc="Thickness of the AR coating, in microns")
    
    _material1 = Instance(CoatedDispersiveMaterial, ())
    _material2 = Instance(CoatedDispersiveMaterial, ())
    
    vtk_grid = Instance(EmptyGridSource, ())
    vtk_cylinder = Instance(tvtk.Cylinder, ())
    vtk_sphere1 = Instance(tvtk.Sphere, ())
    vtk_sphere2 = Instance(tvtk.Sphere, ())
    
    clip2 = Instance(tvtk.ClipDataSet, ())
    clip3 = Instance(tvtk.ClipDataSet, ())
    
    traits_view = View(VGroup(
                       Traceable.uigroup,  
                       Item('CT', editor=NumEditor),
                       Item('diameter', editor=NumEditor),
                       Item('curvature1', editor=NumEditor),
                       Item('curvature2', editor=NumEditor),
                       )
                    )
    
    @on_trait_change("dispersion, dispersion_coating, coating_thickness")
    def on_materials_changed(self):
        self.apply_materials()
        self.update = True
    
    def apply_materials(self):
        air = NondispersiveCurve(1.0)
        self._material1.dispersion_inside = self.dispersion
        self._material1.dispersion_outside = air
        self._material1.dispersion_coating = self.dispersion_coating
        self._material1.coating_thickness = self.coating_thickness
        
        self._material2.dispersion_inside = air
        self._material2.dispersion_outside = self.dispersion
        self._material2.dispersion_coating = self.dispersion_coating
        self._material2.coating_thickness = self.coating_thickness
        
    @on_trait_change("diameter, CT, curvature1, curvature2")
    def on_geometry_changed(self):
        face1, face2 = self.faces.faces
        face1.diameter = self.diameter
        face1.z_height = self.CT/2
        face1.curvature = self.curvature1
        
        face2.diameter = self.diameter
        face2.z_height = -self.CT/2
        face2.curvature = self.curvature2
        
        self.update = True
    
    def _faces_default(self):
        self.apply_materials()
        fl = FaceList(owner=self)
        fl.faces = [ 
                  SphericalFace(owner=self, 
                                diameter=self.diameter,
                                material=self._material1,
                                z_height=self.CT/2, 
                                curvature=self.curvature1),
                  SphericalFace(owner=self, 
                                diameter=self.diameter,
                                material=self._material2,
                                z_height=-self.CT/2, 
                                curvature=self.curvature2),
                ]
        return fl
    
    @on_trait_change("CT, diameter, curvature1, curvature2")
    def config_pipeline(self):
        ct = self.CT
        rad = self.diameter/2
        curve1 = self.curvature1
        curve2 = self.curvature2
        
        size = 41
        spacing = 2*rad / (size-1)
        if curve1 >= 0.0:
            extra1=0.0
        else:
            extra1 = -curve1 - sqrt(curve1**2 - rad**2)
            
        if curve2 <= 0.0:
            extra2=0.0
        else:
            extra2 = curve2 - sqrt(curve2**2 - rad**2)
            
        lsize1 = int((ct/2+extra1)/spacing) + 2
        lsize2 = int((ct/2+extra2)/spacing) + 2
        
        grid = self.vtk_grid
        grid.dimensions = (size,size,lsize1+lsize2)
        grid.origin = (-rad, -rad, -lsize2*spacing)
        grid.spacing = (spacing, spacing, spacing)
                      
        cyl = self.vtk_cylinder
        cyl.center = (0,0,0)
        cyl.radius = rad
        
        s1 = self.vtk_sphere1
        s1.center = (0,0,ct/2 - curve1)
        s1.radius = abs(curve1)
        
        s2 = self.vtk_sphere2
        s2.center = (0,0,-ct/2-curve2)
        s2.radius = abs(curve2)
        
        self.clip2.inside_out = bool(curve1 >= 0.0)
        self.clip3.inside_out = bool(curve2 <= 0.0)
        
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
                      clip_function=self.vtk_sphere1,
                      inside_out=1)
        
        self.clip3.set(input_connection = self.clip2.output_port,
                       clip_function=self.vtk_sphere2,
                       inside_out=1)
        
        topoly = tvtk.GeometryFilter(input_connection=self.clip3.output_port)
        
        norms = tvtk.PolyDataNormals(input_connection=topoly.output_port)
        
        transF = tvtk.TransformFilter(input_connection=norms.output_port, 
                                      transform=self.transform)
        
        self.on_geometry_changed()
        self.config_pipeline()
        grid.modified()
        return transF
    

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
    clip3 = Instance(tvtk.ClipDataSet, ())
    
    traits_view = View(VGroup(
                       Traceable.uigroup,  
                       Item('CT1', editor=NumEditor),
                       Item('CT2', editor=NumEditor),
                       Item('diameter', editor=NumEditor),
                       Item('curvature1', editor=NumEditor),
                       Item('curvature2', editor=NumEditor),
                       Item('curvature3', editor=NumEditor),
                       )
                    )
    
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
        
        self._material2.dispersion_inside = self.dispersion2
        self._material2.dispersion_outside = self.dispersion1
        self._material2.dispersion_coating = self.dispersion_coating
        self._material2.coating_thickness = 0.0
        
        ### We reverse inside and outside because really, outside means
        ### "the direction the surface normals point".
        self._material3.dispersion_inside = air
        self._material3.dispersion_outside = self.dispersion2
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
        face3.z_height = -self.CT2
        face3.curvature = self.curvature3
        self.update = True
    
    def _faces_default(self):
        self.apply_materials()
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
            extra3=0.0
        else:
            extra3 = curve3 - sqrt(curve3**2 - rad**2)
            
        lsize1 = int((ct1+extra1)/spacing) + 2
        lsize2 = int((ct2+extra3)/spacing) + 2
        
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
        
        s2 = self.vtk_sphere2
        s2.center = (0,0,-curve2)
        s2.radius = abs(curve2)
        
        s3 = self.vtk_sphere3
        s3.center = (0,0,-ct2 - curve3)
        s3.radius = abs(curve3)
        
        self.clip2.inside_out = bool(curve1 >= 0.0)
        self.clip3.inside_out = bool(curve3 <= 0.0)
        
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
                      clip_function=self.vtk_sphere1,
                      inside_out=1)
        
        self.clip3.set(input_connection = self.clip2.output_port,
                       clip_function=self.vtk_sphere3,
                       inside_out=1)
        
        clip4 = tvtk.ClipDataSet(input_connection = self.clip3.output_port,
                       clip_function=self.vtk_sphere2,
                       inside_out=1)
        
        clip5 = tvtk.ClipDataSet(input_connection = self.clip3.output_port,
                       clip_function=self.vtk_sphere2,
                       inside_out=0)
        
        topoly = tvtk.GeometryFilter(input_connection=clip4.output_port)
        topoly2 = tvtk.GeometryFilter(input_connection=clip5.output_port)
        
        append = tvtk.AppendPolyData()
        append.add_input_connection(topoly.output_port)
        append.add_input_connection(topoly2.output_port)
        
        norms = tvtk.PolyDataNormals(input_connection=append.output_port)
        
        transF = tvtk.TransformFilter(input_connection=norms.output_port, 
                                      transform=self.transform)
        
        self.on_geometry_changed()
        self.config_pipeline()
        grid.modified()
        return transF
    
    def make_step_shape(self):
        """Creates an OpenCascade BRep Shape
        representation of the object, which can be
        exported to STEP format"""
        from raypier.step_export import make_spherical_lens2, make_compound
        centre = self.centre
        direction = self.direction
        x_axis = self.x_axis
        diameter = self.diameter
        CT1 = self.CT1
        CT2 = 0.0
        CT3 = -self.CT2
        curvature1 = self.curvature1
        curvature2 = self.curvature2
        curvature3 = self.curvature3
        
        lens1 = make_spherical_lens2(CT1, CT2, diameter, curvature1, 
                                     curvature2, centre, direction, x_axis)
        lens2 = make_spherical_lens2(CT2, CT3, diameter, curvature2, 
                                     curvature3, centre, direction, x_axis)
        doublet = make_compound([lens1, lens2])
        return doublet, "yellow"
    
    
class EdmundOptic45805(Doublet):
    name = "Edmund Optics Achromat#45805"
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
    
    
class EdmundOptic45806(Doublet):
    name = "Edmund Optics Achromat#45-806"
    CT1 = 6.0
    CT2 = 4.0
    diameter = 25.4
    offset = 0.0
    curvature1 = 57.01
    curvature2 = -56.72
    curvature3 = -656.51
    dispersion1 = NamedDispersionCurve("N-LAK22")
    dispersion2 = NamedDispersionCurve("N-SF6")
    dispersion_coating = NondispersiveCurve(refractive_index=1.37)
    coating_thickness = 0.283 #microns
    