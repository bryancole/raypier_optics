#    Copyright 2009, Teraview Ltd., Bryan Cole
#
#    This file is part of Raypier.
#
#    Raypier is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.

from traits.api import Float, Instance, on_trait_change, Property, cached_property

from traitsui.api import View, Item, VGroup

from tvtk.api import tvtk


from raypier.bases import Traceable, normaliseVector, NumEditor,\
     ComplexEditor, VectorEditor, Optic, ShapedTraceable
     
from raypier.utils import transformPoints, dotprod
from raypier.sources import RayCollection
from raypier.core.cfaces import CircularFace, RectangularFace, ShapedSphericalFace
from raypier.core.ctracer import FaceList
from raypier.core.cmaterials import DielectricMaterial, \
        CoatedDispersiveMaterial, PECMaterial
from raypier.dispersion import BaseDispersionCurve, NondispersiveCurve
from raypier.shapes import CircleShape, RectangleShape

import math, numpy



class BaseMirror(Traceable):
    abstract=True
    def _vtkproperty_default(self):
        return tvtk.Property(opacity=1.0,
                             color=(0.8,0.8,0.8),
                                representation="surface")


class PECMirror(BaseMirror):
    name = "PEC Mirror"
    abstract = False
    diameter = Float(25.4)
    thickness = Float(5.0, desc="purely for visualisation purposes")
    offset = Float(0.0)
    
    vtk_cylinder = Instance(tvtk.CylinderSource, (),
                            dict(resolution=32), transient=True)
    
    cyl_trans = Instance(tvtk.Transform, (), transient=True)
    
    traits_view = View(VGroup(
                       Traceable.uigroup,
                       Item('diameter', editor=NumEditor),
                       Item('thickness', editor=NumEditor),
                       Item('offset', editor=NumEditor)
                        ),
                   )
    
    def _material_default(self):
        return PECMaterial()
    
    def _faces_default(self):
        m = self.material
        fl = FaceList(owner=self)
        fl.faces = [CircularFace(owner=self, material=m)]
        return fl
    
    def make_step_shape(self):
        from raypier.step_export import make_cylinder
        cyl = make_cylinder(self.centre, 
                             self.direction, 
                             self.diameter/2, 
                             self.thickness,
                             self.offset,
                             self.x_axis)
        return cyl, "green"
    
    @on_trait_change("diameter, thickness, offset")
    def config_pipeline(self):
        cyl = self.vtk_cylinder
        cyl.radius = self.diameter/2.
        thick = self.thickness
        cyl.height = thick
        
        t = self.cyl_trans
        t.identity()
        t.translate(self.offset,0,thick/2.)
        t.rotate_x(90.)
        
        self.update = True
    
    def _pipeline_default(self):
        cyl = self.vtk_cylinder
        norms = tvtk.PolyDataNormals(input_connection=cyl.output_port)
        transF1 = tvtk.TransformFilter(input_connection=norms.output_port, 
                                       transform=self.cyl_trans)
        transF2 = tvtk.TransformFilter(input_connection=transF1.output_port, 
                                       transform=self.transform)
        self.config_pipeline()
        return transF2
    
    
class SphericalMirrorWithHole(ShapedTraceable):
    name = "Spherical Mirror + Hole"
    abstract = False
    
    diameter = Float(25.4)
    hole_diameter = Float(5.0)
    hole_offset = Float(0.0)
    
    curvature = Float(100.0)
    
    thickness = Float(5.0, desc="centre thickness, purely for visualisation purposes")
    
    
    traits_view = View(VGroup(
                       Traceable.uigroup,
                       Item('diameter', editor=NumEditor),
                       Item('hole_diameter', editor=NumEditor),
                       Item('thickness', editor=NumEditor),
                       Item('curvature', editor=NumEditor),
                       Item('hole_offset', editor=NumEditor)
                        ),
                   )
        
    
    def _material_default(self):
        return PECMaterial()
    
    def _shape_default(self):
        return CircleShape(radius=self.diameter/2.) ^ CircleShape(radius=self.hole_diameter/2.,
                                                                  centre=(self.hole_offset,0))
    
    def _faces_default(self):
        m = self.material
        fl = FaceList(owner=self)
        fl.faces = [ShapedSphericalFace(owner=self, material=m, shape=self.shape.cshape)]
        return fl
    
    @on_trait_change("diameter, hole_diameter, thickness, curvature, hole_offset")
    def config_params(self):
        self.grid_extent = self.eval_grid_extent()
        self.shape.shape1.radius = self.diameter/2.
        self.shape.shape2.radius = self.hole_diameter/2.
        self.shape.shape2.centre = (self.hole_offset,0)
        face = self.faces.faces[0]
        face.z_height = self.thickness
        face.curvature = -self.curvature
        self.update = True
        
    def eval_grid_extent(self):
        r = self.diameter/2
        c = self.curvature
        h = c - numpy.sqrt(c**2 - r**2)
        return (-r,r,-r,r,0,self.thickness + h)
    
    def eval_sag_top(self, points):
        centre = numpy.array([0,0,(self.curvature+self.thickness)])
        func = ((points - centre)**2).sum(axis=1) - (self.curvature**2)
        return -func
    
    def eval_sag_bottom(self, points):
        return None
    
    

class PlanarWindow(PECMirror, Optic):
    n_inside = 1.5
    name = "Planar window"
    abstract = False
    
    traits_view = View(VGroup(
                       Traceable.uigroup,
                       Item('diameter', editor=NumEditor),
                       Item('thickness', editor=NumEditor),
                       Item('offset', editor=NumEditor),
                       Item('n_inside', editor=NumEditor),
                        ),
                   )
                   
    vtkproperty = tvtk.Property(opacity = 0.4,
                             color = (0.8,0.8,1.0))
                             
    def _material_default(self):
        return Optic._material_default(self)
                   
    def _thickness_changed(self, new_t):
        self.faces[1].z_plane = new_t
    
    def _faces_default(self):
        m = self.material
        fl = FaceList(owner=self)
        fl.faces = [CircularFace(owner=self, z_plane=0, material=m),
                    CircularFace(owner=self, z_plane=self.thickness, 
                        invert_normal=True, material=m)]
        return fl
    

class PlanarDispersiveWindow(PECMirror, Optic):
    name = "Planar dispersive window"
    abstract = False
    
    material = Instance(CoatedDispersiveMaterial, ())
    
    dispersion = Instance(BaseDispersionCurve, desc="Dispersion curve for the glass")
    dispersion_coating = Instance(BaseDispersionCurve, desc="Dispersion curve for the coating")
    coating_thickness = Float(0.25, desc="Thickness of the AR coating, in microns")
    
    traits_view = View(VGroup(
                       Traceable.uigroup,
                       Item('diameter', editor=NumEditor),
                       Item('thickness', editor=NumEditor),
                       Item('offset', editor=NumEditor),
                       Item('n_inside', editor=NumEditor),
                        ),
                   )
                   
    vtkproperty = tvtk.Property(opacity = 0.4,
                             color = (0.6,0.9,0.5))
                             
    @on_trait_change("dispersion, dispersion_coating, coating_thickness")
    def on_materials_changed(self):
        self.apply_materials()
        self.update = True
    
    def apply_materials(self):
        air = NondispersiveCurve(1.0)
        self.material.dispersion_inside = self.dispersion
        self.material.dispersion_outside = air
        self.material.dispersion_coating = self.dispersion_coating
        self.material.coating_thickness = self.coating_thickness
                   
    def _thickness_changed(self, new_t):
        self.faces[1].z_plane = new_t
        self.update = True
    
    def _faces_default(self):
        self.apply_materials()
        m = self.material
        fl = FaceList(owner=self)
        fl.faces = [CircularFace(owner=self, z_plane=0, material=m),
                    CircularFace(owner=self, z_plane=self.thickness, 
                        invert_normal=True, material=m)]
        return fl


class RectMirror(BaseMirror):
    name = "Rectangular Mirror"
    abstract = False
    length = Float(25.4)
    width = Float(25.4)
    thickness = Float(5.0, desc="purely for visualisation purposes")
    offset = Float(0.0)
    
    vtk_cube = Instance(tvtk.CubeSource, (), transient=True)
    
    cube_trans = Instance(tvtk.Transform, (), transient=True)
    
    traits_view = View(VGroup(
                       Traceable.uigroup,
                       Item('length', editor=NumEditor),
                       Item('width', editor=NumEditor),
                       Item('thickness', editor=NumEditor),
                       Item('offset', editor=NumEditor)
                        ),
                   )
    
    def _faces_default(self):
        fl = FaceList(owner=self)
        fl.faces = [RectangularFace(owner=self, material = self.material)]
        return fl
    
    ''' #copied from cylinder, probably easy to modify for a cube
    def make_step_shape(self):
        from raypier.step_export import make_cylinder
        cyl = make_cylinder(self.centre, 
                             self.direction, 
                             self.diameter/2, 
                             self.thickness,
                             self.offset,
                             self.x_axis)
        return cyl, "green"
    '''
    @on_trait_change("length, width, offset")
    def config_pipeline(self):
        cube = self.vtk_cube
        cube.x_length = self.length
        cube.y_length = self.width
        thick = self.thickness
        cube.z_length = thick
        
        t = self.cube_trans
        t.identity()
        t.translate(self.offset,0,-thick/2) #not sure about self.offset
        t.rotate_x(180.)
        
        self.update = True
    
    def _pipeline_default(self):
        cube = self.vtk_cube
        norms = tvtk.PolyDataNormals(input_connection=cube.output_port)
        transF1 = tvtk.TransformFilter(input_connection=norms.output_port, 
                                       transform=self.cube_trans)
        transF2 = tvtk.TransformFilter(input_connection=transF1.output_port, 
                                       transform=self.transform)
        self.config_pipeline()
        return transF2
    
    def make_step_shape(self):
        """Creates an OpenCascade BRep Shape
        representation of the object, which can be
        exported to STEP format"""
        from raypier.step_export import make_box
        position = self.centre
        direction = self.direction
        x_axis = self.x_axis
        dz = self.thickness
        dx = self.length
        dy = self.width
        
        box = make_box(position, direction, x_axis, dx, dy, dz)
        
        return box, "green"
    

class RectWindow(RectMirror, Optic):
    n_inside = 1.5
    name = "Rectangular window"
    abstract = False
    
    traits_view = View(VGroup(
                       Traceable.uigroup,
                       Item('length', editor=NumEditor),
                       Item('width',editor=NumEditor),
                       Item('thickness', editor=NumEditor),
                       Item('offset', editor=NumEditor),
                       Item('n_inside', editor=NumEditor),
                        ),
                   )
                   
    vtkproperty = tvtk.Property(opacity = 0.4,
                             color = (0.8,0.8,1.0))
                             
    def _material_default(self):
        return Optic._material_default(self)
                   
    def _thickness_changed(self, new_t):
        self.faces[1].z_plane = new_t
    
    def _faces_default(self):
        m = self.material
        fl = FaceList(owner=self)
        fl.faces = [RectangularFace(owner=self, z_plane=0, material=m),
                    RectangularFace(owner=self, z_plane=-self.thickness, 
                        invert_normal=True, material=m)]
        return fl
    
            
if __name__=="__main__":
    from ray_tracer import RayTraceModel, BuildRaySet, BuildConfocalRaySet
    
    mirror = PECMirror(centre=(0,0,0),
                       orientation=0.0,
                       elevation=45.0)
    
#    input_rays = BuildRaySet(origin = (0,0,-20),
#                         direction = (0,0,1),
#                         radius=5.0,
#                         count=20)
    input_rays = BuildConfocalRaySet(focus=(0,0,5), 
                                     direction=(0,0,1),
                                     theta=20, 
                                     working_dist=20, 
                                     count=20)
    
    model = RayTraceModel(optics=[mirror], rays=input_rays)
    model.configure_traits()
