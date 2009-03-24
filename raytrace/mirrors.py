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


from raytrace.tracer import Traceable, normaliseVector, NumEditor,\
     ComplexEditor, VectorEditor
     
from raytrace.utils import transformPoints, dotprod
from raytrace.sources import RayCollection
from faces import CircularFace, PECFace

import math, numpy


class PECCircularFace(CircularFace, PECFace):
    pass


class BaseMirror(Traceable):
    def compute_normal(self, points, cells):
        t = self.transform
        n = numpy.asarray([t.transform_vector(0,0,-1),])
        return numpy.ones(points.shape) * n
    
    def eval_children(self, rays, points, cells, mask=slice(None,None,None)):
        """
        actually calculates the new ray-segments. Physics here
        for Fresnel reflections.
        
        rays - a RayCollection object
        points - a (Nx3) array of intersection coordinates
        cells - a length N array of cell ids
        mask - a bool array selecting items for this Optic
        """
        points = points[mask]
        n = rays.refractive_index[mask]
        cells = cells[mask]
        normal = self.compute_normal(points, cells)
        input_v = rays.direction[mask]
        
        parent_ids = numpy.arange(mask.shape[0])[mask]
        optic = numpy.repeat([self,], points.shape[0] )
        
        S_amp, P_amp, S_vec, P_vec = Convert_to_SP(input_v, 
                                                   normal, 
                                                   rays.E_vector[mask], 
                                                   rays.E1_amp[mask], 
                                                   rays.E2_amp[mask])

        #this is cos(theta), where theta is the angle between the
        #normal and the incident ray
        cosTheta = dotprod(normal, input_v)
        cosThetaNormal = cosTheta*normal
        reflected = input_v - 2*cosThetaNormal
        
        refl_rays = RayCollection(origin=points,
                                   direction = reflected,
                                   max_length = rays.max_length,
                                   E_vector = S_vec,
                                   E1_amp = -S_amp,
                                   E2_amp = -P_amp,
                                   parent = rays,
                                   parent_ids = parent_ids,
                                   optic = optic,
                                   face_id = cells,
                                   refractive_index=n)
        return refl_rays


class PECMirror(BaseMirror):
    name = "PEC Mirror"
    diameter = Float(25.4)
    thickness = Float(5.0, desc="purely for visualisation purposes")
    
    vtk_cylinder = Instance(tvtk.CylinderSource, (),
                            dict(resolution=32), transient=True)
    
    cyl_trans = Instance(tvtk.Transform, (), transient=True)
    
    vtkproperty = tvtk.Property(color=(0.8,0.8,0.8),
                                representation="surface")
    
    traits_view = View(VGroup(
                       Traceable.uigroup,
                       Item('diameter', editor=NumEditor),
                       Item('thickness', editor=NumEditor),
                        ),
                   )
    
    def _diameter_changed(self, dnew):
        self.faces[0].diameter = dnew
    
    def _faces_default(self):
        return [PECCircularFace(transform=self.transform,
                                diameter = self.diameter)]
    
    def make_step_shape(self):
        from raytrace.step_export import make_cylinder
        return make_cylinder(self.centre, 
                             self.direction, 
                             self.diameter/2, 
                             self.thickness)
    
    @on_trait_change("diameter, thickness")
    def config_pipeline(self):
        cyl = self.vtk_cylinder
        cyl.radius = self.diameter/2.
        thick = self.thickness
        cyl.height = thick
        
        t = self.cyl_trans
        t.identity()
        t.translate(0,0,thick/2.)
        t.rotate_x(90.)
        
        self.update = True
    
    def _pipeline_default(self):
        cyl = self.vtk_cylinder
        norms = tvtk.PolyDataNormals(input=cyl.output)
        transF1 = tvtk.TransformFilter(input=norms.output, transform=self.cyl_trans)
        transF2 = tvtk.TransformFilter(input=transF1.output, transform=self.transform)
        self.config_pipeline()
        return transF2
    
            
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
