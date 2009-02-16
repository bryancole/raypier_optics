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
     ComplexEditor, VectorEditor, Convert_to_SP, transformPoints, dotprod
     
from raytrace.sources import RayCollection

import math, numpy


class BaseMirror(Traceable):
    def compute_normal(self, points, cells):
        t = self.transform
        n = numpy.asarray([t.transform_vector(0,0,-1),])
        return numpy.ones(points.shape) * n
    
    def eval_children(self, rays, points, cells, mask):
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
    
#    def trace_segment(self, seg, last_optic=None, last_cell=None):
#        """
#        only one cell
#        """
#        p1 = seg.origin
#        p2 = p1 + seg.MAX_RAY_LENGTH*seg.direction
#        t = self.transform
#        inv_t = self.transform.linear_inverse
#        x1,y1,z1 = inv_t.transform_point(p1)
#        x2,y2,z2 = inv_t.transform_point(p2)
#        
#        r = self.diameter/2
#
#        if (z1 >= 0) and (z2 > 0):
#            return
#        elif (z1 <= 0) and (z2 < 0):
#            return
#        if abs(z1) < self.tolerance:
#            return
#
#        h = -z1/(z2-z1)
#        X = x1 + h*(x2-x1)
#        Y = y1 + h*(y2-y1)
#        #print X, Y, X**2 + Y**2, r**2
#        if (X**2 + Y**2) > r**2:
#            return
#
#        pt = numpy.asarray(t.transform_point(X,Y,0))
#        cell = 0
#        dist = numpy.sqrt(((pt-p1)**2).sum())
#        return dist, pt, cell, self
    
    def trace_rays(self, rays):
        """re-implemented from base class"""
        max_length = rays.max_length
        p1 = rays.origin
        p2 = p1 + max_length*rays.direction
        t = self.transform
        inv_t = t.linear_inverse
        P1 = transformPoints(inv_t, p1)
        P2 = transformPoints(inv_t, p2)
        
        r = self.diameter/2
        
        x1,y1,z1 = P1.T
        x2,y2,z2 = P2.T
        
        h = -z1/(z2-z1)
        X = x1 + h*(x2-x1)
        Y = y1 + h*(y2-y1)
        
        length = max_length*h
        
        mask = (X**2 + Y**2) > r**2
        mask = numpy.logical_or(mask, h < self.tolerance)
        mask = numpy.logical_or(mask, h > 1.0)
        
        length[mask] = numpy.Infinity
        
        t_points = numpy.column_stack((X, Y, numpy.zeros_like(X)))
        points = transformPoints(t, t_points)
    
        dtype=([('length','f8'),('cell','i2'),('point','f8',3)])
        result = numpy.zeros(p1.shape[0], dtype=dtype)
        result['length'] = length
        result['point'] = points
        
        return result
            
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
