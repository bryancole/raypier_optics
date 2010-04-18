#    Copyright 2009, Teraview Ltd., Bryan Cole
#
#    This file is part of Raytrace.
#
#    Raytrace is free software: you can redistribute it and/or modify
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

from enthought.traits.api import Float, Instance, on_trait_change

from enthought.traits.ui.api import View, Item, VGroup

from enthought.tvtk.api import tvtk


from raytrace.bases import Traceable, normaliseVector, NumEditor,\
     ComplexEditor, VectorEditor, Optic
     
from raytrace.utils import transformPoints, dotprod
from raytrace.sources import RayCollection
from raytrace.cfaces import CircularFace, RectangularFace
from raytrace.ctracer import FaceList
from raytrace.cmaterials import DielectricMaterial

import math, numpy



class BaseMirror(Traceable):
    def _vtkproperty_default(self):
        return tvtk.Property(opacity=1.0,
                             color=(0.8,0.8,0.8),
                                representation="surface")


class PECMirror(BaseMirror):
    name = "PEC Mirror"
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
    
    def _faces_default(self):
        fl = FaceList(owner=self)
        fl.faces = [CircularFace(owner=self)]
        return fl
    
    def make_step_shape(self):
        from raytrace.step_export import make_cylinder
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
        norms = tvtk.PolyDataNormals(input=cyl.output)
        transF1 = tvtk.TransformFilter(input=norms.output, transform=self.cyl_trans)
        transF2 = tvtk.TransformFilter(input=transF1.output, transform=self.transform)
        self.config_pipeline()
        return transF2
    

class PlanarWindow(PECMirror, Optic):
    n_inside = 1.5
    name = "Planar window"
    
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

class RectMirror(BaseMirror):
    name = "Rectangular Mirror"
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
        fl.faces = [RectangularFace(owner=self)]
        return fl
    
    ''' #copied from cylinder, probably easy to modify for a cube
    def make_step_shape(self):
        from raytrace.step_export import make_cylinder
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
        norms = tvtk.PolyDataNormals(input=cube.output)
        transF1 = tvtk.TransformFilter(input=norms.output, transform=self.cube_trans)
        transF2 = tvtk.TransformFilter(input=transF1.output, transform=self.transform)
        self.config_pipeline()
        return transF2

class RectWindow(RectMirror, Optic):
    n_inside = 1.5
    name = "Rectangular window"
    
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
