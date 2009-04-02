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

from enthought.traits.ui.api import View, Item, ListEditor, VSplit,\
            RangeEditor, ScrubberEditor, HSplit, VGroup

from enthought.tvtk.api import tvtk
     
from raytrace.bases import Optic, normaliseVector, NumEditor,\
    ComplexEditor, Traceable, transformPoints, transformNormals
    
from raytrace.faces import CircularFace, SphericalFace, DielectricFace

import math, numpy


class PlanarFace(CircularFace, DielectricFace):
    pass

class ConvexFace(SphericalFace, DielectricFace):
    pass


class BaseLens(Optic):
    pass


class PlanoConvexLens(BaseLens):
    name = "Plano-Convex Lens"
    n_inside = 1.5
    
    CT = Float(5.0, desc="centre thickness")
    diameter = Float(15.0)
    curvature = Float(11.7, desc="radius of curvature for spherical face")
    
    #vtkproperty = tvtk.Property(representation="wireframe")
    
    vtk_grid = Instance(tvtk.ProgrammableSource, ())
    vtk_cylinder = Instance(tvtk.Cylinder, ())
    vtk_sphere = Instance(tvtk.Sphere, ())
    
    clip2 = Instance(tvtk.ClipDataSet, ())
    
    traits_view = View(VGroup(
                       Traceable.uigroup,  
                       Item('n_inside', editor=ComplexEditor),
                       Item('CT', editor=NumEditor),
                       Item('diameter', editor=NumEditor),
                       Item('curvature', editor=NumEditor)
                       )
                    )
    
    def _faces_default(self):
        return [PlanarFace(owner=self), ConvexFace(owner=self)]
    
    def create_grid(self):        
        ct = self.CT
        r = self.diameter/2
        curve = self.curvature
        
        source = self.vtk_grid
        sp = source.structured_points_output
        sp.initialize()
        size = 40
        spacing = 2*r / (size-1)
        if curve >= 0.0:
            extra=0.0
        else:
            extra = - curve - math.sqrt(curve**2 - r**2)
        lsize = int((ct+extra)/spacing) + 1
        print "lsize", lsize, lsize * spacing
        sp.dimensions = (size,size,lsize)
        sp.whole_extent=(0,size,0,size,0,lsize)
        sp.origin = (-r, -r, 0)
        sp.spacing = (spacing, spacing, spacing)
        sp.set_update_extent_to_whole_extent()
    
    @on_trait_change("CT, diameter, curvature")
    def config_pipeline(self):
        ct = self.CT
        rad = self.diameter/2
        curve = self.curvature
                      
        cyl = self.vtk_cylinder
        cyl.center = (0,0,0)
        cyl.radius = rad
        
        s = self.vtk_sphere
        s.center = (0,0,ct - curve)
        s.radius = abs(curve)
        
        self.clip2.inside_out = bool(curve >= 0.0)
        
        self.vtk_grid.modified()
        self.update=True
        
                                         
    def _pipeline_default(self):
        grid = self.vtk_grid
        grid.set_execute_method(self.create_grid)
        grid.modified()
        
        trans = tvtk.Transform()
        trans.rotate_x(90.)
        cyl = self.vtk_cylinder
        cyl.transform = trans
        
        clip1 = tvtk.ClipVolume(input=grid.structured_points_output,
                                 clip_function=self.vtk_cylinder,
                                 inside_out=1)
        
        self.clip2.set(input = clip1.output,
                      clip_function=self.vtk_sphere,
                      inside_out=1)
        
        topoly = tvtk.GeometryFilter(input=self.clip2.output)
        norms = tvtk.PolyDataNormals(input=topoly.output)
        
        transF = tvtk.TransformFilter(input=norms.output, transform=self.transform)
        self.config_pipeline()
        grid.modified()
        return transF
    
#    def trace_rays(self, rays, face_id=None):
#        p1 = rays.origin
#        p2 = p1 + rays.max_length*rays.direction
#        t = self.transform
#        inv_t = t.linear_inverse
#        P1 = transformPoints(inv_t, p1)
#        P2 = transformPoints(inv_t, p2)
#        
#        dtype=([('length','f8'),('cell','i2'),('point','f8',3)])
#        result = numpy.zeros(p1.shape[0], dtype=dtype)
#        
#        if face_id is None:
#            plane = self.intersect_with_planeT(P1, P2)
#            sphere = self.intersect_with_sphereT(P1, P2)
#            
#            choice = numpy.argmin([plane[0], sphere[0]], axis=0)
#            
#            result['length'] = numpy.choose(choice, [plane[0], sphere[0]]) * rays.max_length
#            result['cell'] = choice
#            points = numpy.choose(choice.reshape(-1,1)*numpy.ones((1,3),'i'), [plane[1], sphere[1]])
#            result['point'] = transformPoints(t, points)
#        else:
#            intersect = {0: self.intersect_with_planeT,
#                         1: self.intersect_with_sphereT}[face_id]
#            result['length'] = intersect[0] * rays.max_length
#            results['cell'] = face_id
#            result['point'] = transformPoints(t, intersect[1])
#        return result
#    
#    def intersect_with_planeT(self, p1, p2):
#        r = self.diameter/2
#        x1,y1,z1 = p1.T
#        x2,y2,z2 = p2.T
#        #check the ray crosses zero
#        mask = (z1>=0) == (z2>=0)
#        #check it's not a self-intersection
#        mask = numpy.logical_or(mask, abs(z1) < self.tolerance)
#
#        h = -z1/(z2-z1)
#        X = x1 + h*(x2-x1)
#        Y = y1 + h*(y2-y1)
#        #print X, Y, X**2 + Y**2, r**2
#        mask = numpy.logical_or(mask, (X**2 + Y**2) > r**2)
#        
#        h[mask] = numpy.Infinity
#        
#        return h, numpy.column_stack((X, Y, numpy.zeros_like(X)))
#            
#    def intersect_with_sphereT(self, p1, p2):
#        rad = self.curvature
#        dot = lambda a,b: (a*b).sum(axis=-1)
#        r = p1
#        s = p2 - r
#        c = numpy.array([[0,0,self.CT - self.curvature]])
#        d = r - c
#        
#        A = dot(s,s)
#        B = 2*(dot(s,d))
#        C = dot(d,d) - rad**2
#        
#        D = B**2 - 4*A*C
#        
#        mask = D < 0
#
#        E = numpy.sqrt(D)
#        roots = ((-B+E)/(2*A), (-B-E)/(2*A))
#        
#        for R in roots:
#            m = numpy.logical_or(R>1.0, R<0.01)
#            m = numpy.logical_or(m, mask)
#            pt = r + R.reshape(-1,1)*s
#            m = numpy.logical_or(m, (pt[:,0]**2 + pt[:,1]**2) > (self.diameter/2.)**2)
#            m = numpy.logical_or(m, ((pt[:,2]-c[:,2])/rad)<0)
#            R[m] = numpy.Infinity
#        R1,R2 = roots
#        selected = numpy.choose(R2<R1,[R1,R2])
#        
#        return selected, r + selected.reshape(-1,1)*s
#        
#    def compute_normal(self, points, cell_ids):
#        trans = self.transform
#        inv_trans = trans.linear_inverse
#        
#        if isinstance(cell_ids, int):
#            if cell_id==0:
#                #intersection with plane face
#                return normaliseVector(trans.transform_normal(0,0,-1))
#            elif cell_id==1:
#                #spherical surface
#                P = transformPoints(inv_trans, points)
#                centre = numpy.array(self.vtk_sphere.center)
#                radius = self.vtk_sphere.radius
#                curve = self.curvature
#                grad = (P - centre)*2
#                if curve < 0.0:
#                    grad *= -1
#                return transformNormals(trans, grad)
#            else:
#                raise ValueError("unknown cell ID")
#        else:            
#            P = transformPoints(inv_trans, points)
#            centre = numpy.array(self.vtk_sphere.center)
#            radius = self.vtk_sphere.radius
#            curve = self.curvature
#            grad = (P - centre)*2
#            if curve < 0.0:
#                grad *= -1
#            
#            normals = transformNormals(trans, grad)
#            
#            normals[cell_ids==0] = numpy.array(trans.transform_normal(0,0,-1))
#            return normals

        
if __name__=="__main__":
    from ray_tracer import RayTraceModel, BuildRaySet
    
    lens = PlanoConvexLens(orientation=0.0,
                           elevation=0.0,
                           CT=5.,
                           curvature=12.)
    
    input_rays = BuildRaySet(origin = (0,0,-20),
                         direction = (0,0,1),
                         radius=5.0,
                         count=20)
    
    model = RayTraceModel(optics=[lens], rays=input_rays)
    model.configure_traits()
