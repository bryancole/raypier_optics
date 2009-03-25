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

from enthought.traits.api import HasTraits, Instance, Float, Complex,\
        Tuple, Property, on_trait_change, PrototypedFrom, Any

from enthought.tvtk.api import tvtk

import numpy

from raytrace.utils import Convert_to_SP, dotprod, transformPoints,\
        transformNormals
from raytrace.rays import RayCollection

class Face(HasTraits):
    owner = Any #Instance(klass="raytrace.tracer.Traceable")
    
    transform = PrototypedFrom('owner')
    
    tolerance = Float(0.0001) #excludes ray-segments below this length
    
    def trace_rays(self, rays):
        raise NotImplementedError
    
    def intersect(self, P1, P2):
        raise NotImplementedError
        
    def compute_normal(self, points):
        """
        Computes the surface normal vectors for a given array of
        intersection points, in the global coordinate system
        
        @param point: (n,3) array of points
        
        @return: (n,3) array, the normal vectors (should have unit magnitude)
        """
        raise NotImplementedError
                       
    def eval_children(self, rays, points, mask=slice(None,None,None)):
        """
        actually calculates the new ray-segments. Physics here
        for Fresnel reflections.
        
        rays - a RayCollection object
        points - a (Nx3) array of intersection coordinates
        mask - a bool array selecting items for this Optic. (optional)
        """
        raise NotImplementedError
    

class CircularFace(Face):
    diameter = PrototypedFrom('owner')
    
    def trace_rays(self, rays):
        raise NotImplementedError
    
    def intersect(self, P1, P2, max_lenth):
        """
        Calculate intersection point for a ray segment between two points
        P1 and P2, in the local coordinate system.
        """
        max_length = numpy.sqrt(((P2 - P1)**2).sum(axis=1))
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
    
        dtype=([('length','f8'),('face', 'O'),('point','f8',3)])
        result = numpy.empty(P1.shape[0], dtype=dtype)
        result['length'] = length
        result['face'] = self
        result['point'] = t_points
        return result
    
    def compute_normal(self, points):
        """computes the surface normal in the global coordinate system"""
        t = self.transform
        n = numpy.asarray([t.transform_vector(0,0,-1),])
        return numpy.ones(points.shape) * n
    

class EllipsoidFace(Face):
    focus1 = PrototypedFrom('owner')
    focus2 = PrototypedFrom('owner')
    size = PrototypedFrom('owner')
    
    axes = PrototypedFrom('owner')
    
    ellipse_trans = PrototypedFrom('owner')
    combined_trans = PrototypedFrom('owner')
    
    X_bounds = PrototypedFrom('owner')
    Y_bounds = PrototypedFrom('owner')
    Z_bounds = PrototypedFrom('owner')
    
    def intersect(self, P1, P2, max_length):
        et = self.ellipse_trans
        inv_et = et.linear_inverse
        PP1 = transformPoints(et, P1)
        PP2 = transformPoints(et, P2)
        
        r = PP1
        s = PP2 - PP1
        
        rx, ry, rz = r.T
        sx, sy, sz = s.T
        
        xa,xb = self.axes
        
        B = xb**2
        A = xa**2
        
        a = A*(sz*sz + sy*sy) + B*sx*sx
        b = 2*( A*(rz*sz + ry*sy) + B*rx*sx )
        c = A*(rz*rz + ry*ry) + B*rx*rx - A*B
        
        d = b*b - 4*a*c
        sqd = numpy.sqrt(d)
        roots = numpy.array([(-b + sqd)/(2*a), 
                                    (-b - sqd)/(2*a)]) #shape=(2,N)
        newaxis = numpy.newaxis
        points = r[newaxis,:,:] + roots[:,:,newaxis]*s[newaxis,:,:] #shape=(2,N,3)
        
        #pts = numpy.empty_like(points)
        #pts[0,:,:] = transformPoints(inv_et, points[0,:,:]) #shape=(2,N,3)
        #pts[1,:,:] = transformPoints(inv_et, points[1,:,:])
        pts = transformPoints(inv_et, points.reshape(-1,3)).reshape(2,-1,3)
        
        xmin, xmax = min(self.X_bounds), max(self.X_bounds)
        ymin, ymax = min(self.Y_bounds), max(self.Y_bounds)
        zmin, zmax = min(self.Z_bounds), max(self.Z_bounds)
        
        mask = (d<0.0)[newaxis,:] #shape=(1,N)
        mask = numpy.logical_or(mask, pts[:,:,0]>xmax) #shape=(2,N)
        mask = numpy.logical_or(mask, pts[:,:,0]<xmin)
        mask = numpy.logical_or(mask, pts[:,:,1]>ymax)
        mask = numpy.logical_or(mask, pts[:,:,1]<ymin)
        mask = numpy.logical_or(mask, pts[:,:,2]>zmax)
        mask = numpy.logical_or(mask, pts[:,:,2]<zmin)
        
        mask = numpy.logical_or(mask, roots < self.tolerance)
        mask = numpy.logical_or(mask, roots > 1.0)
        
        roots[mask] = numpy.Infinity
        
        nearest_id = roots.argmin(axis=0) #shape = (N,)
        idx = numpy.arange(nearest_id.shape[0])
        nearest = pts[nearest_id, idx, :]
        h_nearest = roots[nearest_id, idx]
        
        dtype=([('length','f8'),('face', 'O'),('point','f8',3)])
        result = numpy.empty(P1.shape[0], dtype=dtype)
        result['length'] = h_nearest*max_length
        result['face'] = self
        result['point'] = nearest
        return result
    
    def compute_normal(self, points):
        """computes the surface normal in the global coordinate system"""
        t = self.combined_trans
        inv_t = t.linear_inverse
        
        #transform to ellipse reference frame
        P = transformPoints(inv_t, points)
        a,b = self.axes
        nx = P[:,0]/-(a**2)
        ny = P[:,1]/-(b**2)
        nz = P[:,2]/-(b**2)
        n = numpy.column_stack([nx, ny, nz]) 
        normals = transformNormals(t, n)
        return normals
    
    
    


class SphericalFace(Face):
    sphere_centre = Tuple(Float, Float, Float)
    curvature = Float
    diameter = Float

    def trace_rays(self, rays):
        raise NotImplementedError
    
    def compute_normal(self, points):
        return self.normal


class PECFace(Face):
    def eval_children(self, rays, points, mask=slice(None,None,None)):
        """
        actually calculates the new ray-segments. Physics here
        for Fresnel reflections.
        
        rays - a RayCollection object
        points - a (Nx3) array of intersection coordinates
        mask - a bool array selecting items for this Optic
        """
        points = points[mask]
        n = rays.refractive_index[mask]
        normal = self.compute_normal(points)
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
        
        faces = numpy.array([self,] * points.shape[0])
        
        refl_rays = RayCollection(origin=points,
                                   direction = reflected,
                                   max_length = rays.max_length,
                                   E_vector = S_vec,
                                   E1_amp = -S_amp,
                                   E2_amp = -P_amp,
                                   parent = rays,
                                   parent_ids = parent_ids,
                                   face = faces,
                                   refractive_index=n)
        return refl_rays
    
    
class DielectricFace(Face):
    n_inside = Complex
    n_outside = Complex
    
    def eval_children(self, rays, points, cells, mask=slice(None,None,None)):
        points = points[mask]
        if isinstance(cells, int):
            normal = self.compute_normal(points, cells)
            cells = numpy.ones(points.shape[0], numpy.int) * cells
        else:
            cells = cells[mask] ###reshape not necessary
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
        
        origin = points
            
        fromoutside = cosTheta < 0
        n1 = numpy.where(fromoutside, self.n_outside.real, self.n_inside.real)
        n2 = numpy.where(fromoutside, self.n_inside.real, self.n_outside.real)
        flip = numpy.where(fromoutside, 1, -1)
            
        abscosTheta = numpy.abs(cosTheta)
        
        N2 = (n2/n1)**2
        N2cosTheta = N2*abscosTheta
        
        #if this is less than zero, we have Total Internal Reflection
        N2_sin2 = abscosTheta**2 + (N2 - 1)
        
        TIR = N2_sin2 < 0.0
        sqrt = numpy.sqrt
        
        cosThetaNormal = cosTheta*normal
        reflected = input_v - 2*cosThetaNormal
        sqrtN2sin2 = numpy.where(TIR, 1.0j*sqrt(-N2_sin2), sqrt(N2_sin2))
        #print "root n2.sin2", sqrtN2sin2
        
        #Fresnel equations for reflection
        R_p = (N2cosTheta - sqrtN2sin2) / (N2cosTheta + sqrtN2sin2)
        R_s = (abscosTheta - sqrtN2sin2) / (abscosTheta + sqrtN2sin2)
        #print "R_s", R_s, "R_p", R_p
        
        ###Now calculate transmitted rays
        d1 = input_v
        tangent = d1 - cosThetaNormal
        
        tan_mag_sq = ((n1*tangent/n2)**2).sum(axis=1).reshape(-1,1)        
        
        c2 = numpy.sqrt(1 - tan_mag_sq)
        transmitted = tangent*(n1/n2) - c2*normal*flip 
        #print d1, normal, tangent, transmitted, "T"
        
        cos1 = abscosTheta
        #cos of angle between outgoing ray and normal
        cos2 = abs(dotprod(transmitted, normal))
        
        Two_n1_cos1 = (2*n1)*cos1
        
        aspect = sqrt(cos2/cos1) * Two_n1_cos1
        
        #Fresnel equations for transmission
        T_p = aspect / ( n2*cos1 + n1*cos2 )
        T_s = aspect / ( n2*cos2 + n1*cos1 )
        #print "T_s", T_s, "T_p", T_p
        
        if self.all_rays:
            refl_rays = RayCollection(origin=origin,
                                       direction = reflected,
                                       max_length = rays.max_length,
                                       E_vector = S_vec,
                                       E1_amp = S_amp*R_s,
                                       E2_amp = P_amp*R_p,
                                       parent_ids = parent_ids,
                                       optic = optic,
                                       face_id = cells,
                                       refractive_index=n1)
            
            trans_rays = RayCollection(origin=origin,
                                       direction = transmitted,
                                       max_length = rays.max_length,
                                       E_vector = S_vec,
                                       E1_amp = S_amp*T_s,
                                       E2_amp = P_amp*T_p,
                                       parent_ids = parent_ids,
                                       optic = optic,
                                       face_id = cells,
                                       refractive_index=n2)
            
            allrays = collectRays(refl_rays, trans_rays)
            allrays.parent = rays
            return allrays
        else:
            TIR.shape=-1,1
            tir = TIR*numpy.ones(3)
            direction = numpy.where(tir, reflected,transmitted)
            E1_amp = S_amp*numpy.where(TIR, R_s, T_s)
            E2_amp = P_amp*numpy.where(TIR, R_p, T_p)
            refractive_index = numpy.where(TIR, n1, n2)
            
            return RayCollection(origin=origin,
                               direction = direction,
                               max_length = rays.max_length,
                               E_vector = S_vec,
                               E1_amp = E1_amp,
                               E2_amp = E2_amp,
                               parent_ids = parent_ids,
                               optic = optic,
                               face_id = cells,
                               refractive_index=refractive_index) 
            

    