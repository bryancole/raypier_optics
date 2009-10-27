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
        Tuple, Property, on_trait_change, PrototypedFrom, Any, Str, Bool
        
from enthought.traits.ui.api import View, Item

from enthought.tvtk.api import tvtk

import numpy

from raytrace.utils import Convert_to_SP, dotprod, transformPoints,\
        transformNormals
        
from raytrace.rays import RayCollection

class Face(HasTraits):
    owner = Any #Instance(klass="raytrace.tracer.Traceable")
    name = Str("a traceable face")
    transform = PrototypedFrom('owner')
    
    tolerance = Float(0.0001) #excludes ray-segments below this length
    
    traits_view = View()
    
    def __repr__(self):
        cls_name = str(self.__class__.__name__)
        return "%s of %s"%(cls_name, self.owner.name)
    
    def trace_rays(self, rays):
        owner = self.owner
        max_length = rays.max_length
        p1 = rays.origin
        p2 = p1 + max_length*rays.direction
        t = owner.transform
        inv_t = t.linear_inverse
        P1 = transformPoints(inv_t, p1)
        P2 = transformPoints(inv_t, p2)
        
        intersections = self.intersect(P1, P2, max_length)
            
        rays.length = intersections['length'].reshape(-1,1)
        t_points = intersections['point']
        points = transformPoints(t, t_points)
        mask = numpy.ones(rays.number, numpy.bool) #intersections['length'] != numpy.Infinity
        print self, rays.origin.min(axis=0), rays.origin.max(axis=0)
        print points.min(axis=0), points.max(axis=0)
        children = self.eval_children(rays, points, mask)
        return children
    
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
        return None
    

class CircularFace(Face):
    name = "circular face"
    diameter = PrototypedFrom('owner')
    offset = PrototypedFrom('owner')
    
    def intersect(self, P1, P2, max_lenth):
        """
        Calculate intersection point for a ray segment between two points
        P1 and P2, in the local coordinate system.
        """
        max_length = numpy.sqrt(((P2 - P1)**2).sum(axis=1))
        r = self.diameter/2
        offset = self.offset
        
        x1,y1,z1 = P1.T
        x2,y2,z2 = P2.T
        
        h = -z1/(z2-z1)
        X = x1 + h*(x2-x1) - offset
        Y = y1 + h*(y2-y1)
        
        length = max_length*h
        
        mask = (X**2 + Y**2) > r**2
        mask = numpy.logical_or(mask, h < self.tolerance)
        mask = numpy.logical_or(mask, h > 1.0)
        
        length[mask] = numpy.Infinity
        
        t_points = numpy.column_stack((X+offset, Y, numpy.zeros_like(X)))
    
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
    
class RectangularFace(Face):
    name = "rectangular face"
    length = PrototypedFrom('owner')
    width = PrototypedFrom('owner')
    #offset = PrototypedFrom('owner')
    
    def intersect(self, P1, P2, max_lenth):
        """
        Calculate intersection point for a ray segment between two points
        P1 and P2, in the local coordinate system.
        """
        max_length = numpy.sqrt(((P2 - P1)**2).sum(axis=1))
        len = self.length/2.
        wid = self.width
        #offset = self.offset
        
        x1,y1,z1 = P1.T
        x2,y2,z2 = P2.T
        
        h = -z1/(z2-z1)
        X = x1 + h*(x2-x1) #- offset
        Y = y1 + h*(y2-y1)
        
        length = max_length*h
        
        mask = Y > len
        mask = numpy.logical_or(mask, Y < -len)
        mask = numpy.logical_or(mask, X > wid/2)
        mask = numpy.logical_or(mask, X < -wid/2)
        mask = numpy.logical_or(mask, h < self.tolerance)
        mask = numpy.logical_or(mask, h > 1.0)
        
        length[mask] = numpy.Infinity
        
        # was x + offset before i removed offsets
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
        
class PolygonFace(Face):
    name = "polygon face"
    diameter = PrototypedFrom('owner')
    offset = PrototypedFrom('owner')


class EllipsoidFace(Face):
    name = "ellipsoid face"
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
        """computes the surface normal in the global coordinate system
        
        @param points: a Nx3 array
        """
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
    
    
class ParaxialLensFace(CircularFace):
    name = "paraxial face"
    diameter = PrototypedFrom('owner')
    focal_length = PrototypedFrom('owner')
    
    def compute_normal(self, points):
        raise NotImplementedError

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
        #normal = self.compute_normal(points)
        input_v = rays.direction[mask]
        parent_ids = numpy.arange(mask.shape[0])[mask]

        t = self.transform
        inv_t = t.linear_inverse
        
        #transform to face reference frame
        P = transformPoints(inv_t, points)
        z_f = self.focal_length
        
        r_sq = P[:,0]**2 + P[:,1]**2
        
        r2z2 = numpy.sqrt(r_sq + z_f**2)
        
        input_v_t = transformNormals(inv_t, input_v)
        
        dkx = -P[:,0]/r2z2
        dky = -P[:,1]/r2z2
        
        sign = numpy.where(input_v_t[:,2]>=0, 1,-1)
        
        output_v_t = numpy.empty_like(input_v_t)
        output_v_t[:,0] = input_v_t[:,0] + dkx
        output_v_t[:,1] = input_v_t[:,1] + dky
        output_v_t[:,2] = numpy.sqrt(1 - dkx**2 - dky**2) * sign
        
        output_v = transformNormals(t, output_v_t)
        
        S_amp, P_amp, S_vec, P_vec = Convert_to_SP(input_v, 
                                                   output_v, 
                                                   rays.E_vector[mask], 
                                                   rays.E1_amp[mask], 
                                                   rays.E2_amp[mask])
        
        faces = numpy.array([self,] * points.shape[0])
        
        z_axis = numpy.array([[0,0,1]])
        cosTheta_in = dotprod(input_v_t, z_axis)
        cosTheta_out = dotprod(output_v_t, z_axis)
        
        aspect = numpy.sqrt(cosTheta_in/cosTheta_out)
        
        output_rays = RayCollection(origin=points,
                                   direction = output_v,
                                   max_length = rays.max_length,
                                   E_vector = S_vec,
                                   E1_amp = S_amp*aspect,
                                   E2_amp = P_amp*aspect,
                                   parent = rays,
                                   parent_ids = parent_ids,
                                   face = faces,
                                   refractive_index=n)
        return output_rays
    
    

class SphericalFace(Face):    
    CT = PrototypedFrom("owner")
    diameter = PrototypedFrom("owner")
    curvature = PrototypedFrom("owner")

    def trace_rays(self, rays):
        raise NotImplementedError
    
    def compute_normal(self, points):
        t = self.transform
        inv_t = t.linear_inverse
        P = transformPoints(inv_t, points)
        vtk_sphere = self.owner.vtk_sphere
        centre = numpy.array(vtk_sphere.center)
        radius = vtk_sphere.radius
        curve = self.curvature
        grad = (P - centre)*2
        if curve < 0.0:
            grad *= -1
        return transformNormals(t, grad)
    
    def intersect(self, p1, p2, max_length):
        rad = self.curvature
        dot = lambda a,b: (a*b).sum(axis=-1)
        r = p1
        s = p2 - r
        c = numpy.array([[0,0,self.CT - self.curvature]])
        d = r - c
        
        A = dot(s,s)
        B = 2*(dot(s,d))
        C = dot(d,d) - rad**2
        
        D = B**2 - 4*A*C
        
        mask = D < 0

        E = numpy.sqrt(D)
        roots = ((-B+E)/(2*A), (-B-E)/(2*A))
        
        for R in roots:
            m = numpy.logical_or(R>1.0, R<0.01)
            m = numpy.logical_or(m, mask)
            pt = r + R.reshape(-1,1)*s
            m = numpy.logical_or(m, (pt[:,0]**2 + pt[:,1]**2) > (self.diameter/2.)**2)
            m = numpy.logical_or(m, ((pt[:,2]-c[:,2])/rad)<0)
            R[m] = numpy.Infinity
        R1,R2 = roots
        selected = numpy.choose(R2<R1,[R1,R2])
        
        pts = r + selected.reshape(-1,1)*s
    
        dtype=([('length','f8'),('face', 'O'),('point','f8',3)])
        result = numpy.empty(p1.shape[0], dtype=dtype)
        result['length'] = selected*max_length
        result['face'] = self
        result['point'] = pts
        return result
    
    
class OffAxisParabolicFace(Face):
    name = "OAP Face"
    EFL = PrototypedFrom("owner")
    diameter = PrototypedFrom("owner")
    
    def compute_normal(self, points):
        """
        evaluate normalised Normal vector
        """
        n = points.shape[0]
        t = self.transform
        inv_t = t.linear_inverse
        t_points =transformPoints(inv_t, points)
        
        twoA = 1/self.EFL
        
        x,y,z = t_points.T
        r_sq = x**2 + y**2
        
        n_x = -x*twoA*r_sq
        n_y = -y*twoA*r_sq
        n_z = y**2 + x**2
        
        t_normal = numpy.column_stack((n_x, n_y, n_z))
        
#        coefs = self.owner.vtk_quadric.coefficients
#        ax = 2*coefs[0]/coefs[8]
#        ay = 2*coefs[1]/coefs[8]
#        t_normal = numpy.column_stack((ax*t_points[:,0], 
#                                       ay*t_points[:,1],
#                                       -numpy.ones(n)))
        transformed = transformNormals(t, t_normal)
        return transformed
    
    def intersect(self, P1, P2, max_length):
        """
        
        @param p1: a (n,3) array of points, start of each ray
        @param p2: a (n,3) array of point, ends of the rays
        """
        efl = self.EFL #scalar
        A = 1 / (2*efl)
        s = P2 - P1
        r = P1
        r[:,2] += efl/2.
        
        sx, sy, sz = s.T
        rx, ry, rz = r.T
        
        a = A*(sx**2 + sy**2)
        b = 2*A*(rx*sx + ry*sy) - sz
        c = A*(rx**2 + ry**2) - rz
        
        d = b**2 - 4*a*c
        
        e = numpy.sqrt(d)
        roots = [(-b+e)/(2*a), (-b-e)/(2*a)]
        
        singular = numpy.abs(a)<1e-10
        approx = (-c/b)[singular]
        
        mask = d<0
        for root in roots:
            root[mask] = numpy.Infinity
        
        root1, root2 = roots
        
        root1[singular] = approx
        root2[singular] = numpy.Infinity
        
        tol = self.tolerance
        root1[root1<tol] = numpy.Infinity
        root2[root2<tol] = numpy.Infinity
        
        P1 = r + root1[:,numpy.newaxis]*s
        P1[:,2] -= efl/2.
        
        P2 = r + root2[:,numpy.newaxis]*s
        P2[:,2] -= efl/2.
        
        rad = self.diameter/2.
        mask1 = (P1[:,0]-efl)**2 + P1[:,1]**2 > rad**2
        mask2 = (P2[:,0]-efl)**2 + P2[:,1]**2 > rad**2
        
        root1[mask1] = numpy.Infinity
        root2[mask2] = numpy.Infinity
        
        select = root1<root2
        shortest = numpy.where(select, root1, root2)

        select2 = select[:,numpy.newaxis]*numpy.ones(3,'i')
        P = numpy.choose(select2, [P2, P1])
        
        dtype=([('length','f8'),('face', 'O'),('point','f8',3)])
        result = numpy.empty(P1.shape[0], dtype=dtype)
        result['length'] = shortest*max_length
        result['face'] = self
        result['point'] = P
        return result


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
                                   refractive_index=n,
                                   normals = normal)
        return refl_rays
    
    
class DielectricFace(Face):
    n_inside = PrototypedFrom("owner")
    n_outside = PrototypedFrom("owner")
    
    all_rays = Bool(False)
    
    def eval_children(self, rays, points, mask=slice(None,None,None)):
        points = points[mask]
        normal = self.compute_normal(points)
        input_v = rays.direction[mask]
        faces = numpy.array([self,] * points.shape[0])
        parent_ids = numpy.arange(mask.shape[0])[mask]
        
        S_amp, P_amp, S_vec, P_vec = Convert_to_SP(input_v, 
                                                   normal, 
                                                   rays.E_vector[mask], 
                                                   rays.E1_amp[mask], 
                                                   rays.E2_amp[mask])

        #this is cos(theta), where theta is the angle between the
        #normal and the incident ray
        cosTheta = dotprod(normal, input_v)
        
        origin = points
        
        wavelen = rays.wavelength
        n_inside, n_outside = self.owner.calc_refractive_index(wavelen)
            
        fromoutside = cosTheta < 0
        n1 = numpy.where(fromoutside, n_outside.real, n_inside.real)
        n2 = numpy.where(fromoutside, n_inside.real, n_outside.real)
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
                                       faces=faces,
                                       refractive_index=n1,
                                       normals=normal)
            
            trans_rays = RayCollection(origin=origin,
                                       direction = transmitted,
                                       max_length = rays.max_length,
                                       E_vector = S_vec,
                                       E1_amp = S_amp*T_s,
                                       E2_amp = P_amp*T_p,
                                       parent_ids = parent_ids,
                                       faces = faces,
                                       refractive_index=n2,
                                       normals=normal)
            
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
                               parent=rays,
                               parent_ids = parent_ids,
                               faces=faces,
                               refractive_index=refractive_index) 
            

    