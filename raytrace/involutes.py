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

"""
Maths for a cylinder and also its involute.
"""
from enthought.traits.api import HasTraits, Array, Float, Complex,\
            Property, List, Instance, on_trait_change, Range, Any,\
            Tuple, Event, cached_property, Set, Int, Trait, PrototypedFrom
from enthought.traits.ui.api import View, Item, ListEditor, VSplit,\
            RangeEditor, ScrubberEditor, HSplit, VGroup
from enthought.tvtk.api import tvtk
import numpy
import math
from itertools import izip

from raytrace.bases import Traceable, normaliseVector, NumEditor,\
                transformNormals, transformPoints
from raytrace.mirrors import BaseMirror
from raytrace.faces import PECFace, Face
from raytrace.more_utils import compute_length, interpolate_z
from scipy.optimize import fsolve

#to do: take self.compute_length out of involut and use imported
# also add interpolate Z to more_utils.
      
class CylinderFace(Face):
    name = "Cylinder Face"
    length = PrototypedFrom("owner")
    radius = PrototypedFrom("owner")
    
    """angles must be between 0 and 360!"""
    begin_angle = PrototypedFrom("owner")
    end_angle = PrototypedFrom("owner")
    
    def compute_normal(self, points):
        """
        evaluate surface normal vector
        """
        n = points.shape[0]         #returns how many points there are
        t = self.transform
        inv_t = t.linear_inverse
        t_points =transformPoints(inv_t, points)

        #its just a circle, so the vector pointing to the surface from the orgin
        # already is the surface normal vector! y component = 0
        
        t_points[:,1] = numpy.zeros(t_points[:,1].shape)
        
        L = compute_length(numpy.zeros(t_points.shape),t_points)
        
        n_x = t_points[:,0]
        n_y = t_points[:,1]  # all = 0
        n_z = t_points[:,2]

        t_normal = numpy.column_stack((n_x, n_y, n_z))
        
        return transformNormals(t, t_normal)
        
    def intersect(self, P1, P2, max_length):
        """ takes lines from P1 to P2 and returns first intersection in local 
        coordinates.  z axis is normal to surface at x=0 which is also angle =0.
        However, y and z are switched for easy reuse of code.
        
        if this code ever needs optimization, one good step is to parameterize 
        circle first, and solve for intersection and theta of intersection to
        begin with and then use that theta for boundry conditions
        """
        
        n = P1.shape[0]
        R = self.radius
        
        #This was originally written with the Z axis as the long axis of the trough,
        #but inorder for the direction parameter to be useful and point from
        #vertex to focus, the y axis must be the long axis.  So, y and z are
        #here switched for all the following calculations and then switched back
        # right before the function returns its points
        
        P1[:,1:] = numpy.fliplr(P1[:,1:]).copy()
        P2[:,1:] = numpy.fliplr(P2[:,1:]).copy()
        
        #turn array of points into y = mx + q
        m = (P1[:,1]-P2[:,1])/(P1[:,0]-P2[:,0]) # m = y1-y2 / x1-x2
        q = P1[:,1]-m*P1[:,0]                   #q = y - mx
        
        #solve system of equations: y = mx+b and y^2 + x^2 = R^2 
        
        a = (m**2 + 1)
        b = 2*m*q
        c = q**2 - R**2*numpy.ones(q.shape)
        
        d = b**2 - 4*a*c
        
        miss_mask = d < 0

        E = numpy.sqrt(d)
        roots = numpy.array([(-b+E)/(2*a), (-b-E)/(2*a)])
        for root in roots:
            root[miss_mask] = numpy.inf
                        
        
        root1, root2 = roots
        
        #print "roots: ", roots
        
        #put these roots into a list of intersection points using y = mx + q
        #I make these 3d with z=0, Which i'll fix later
        inter1 = numpy.array([root1,m*root1 + q, numpy.zeros(n)]).T
        inter2 = numpy.array([root2,m*root2 + q, numpy.zeros(n)]).T
        
        #Where the slope was infinite these values are wrong:
        
        pos_result = numpy.array([P1[:,0], numpy.sqrt(R**2 - P1[:,0]**2) ,numpy.zeros(n)]).T
        neg_result = numpy.array([P1[:,0], -numpy.sqrt(R**2 - P1[:,0]**2) ,numpy.zeros(n)]).T
        perp_fix = numpy.array([numpy.abs(m) > 1000000.**2]*3).T
        
        #print "?: ",m, "\n",perp_fix
        
        inter1 = numpy.where(perp_fix,pos_result, inter1)
        inter2 = numpy.where(perp_fix,neg_result, inter2)
        
        #print "1: ",inter1,"\n",inter2

        #Where the ray was parallel to the long axis, the above fix fixes wrong
        
        parallel_result = numpy.array(numpy.ones([n,3])*numpy.inf).T
        parallel_cond = numpy.logical_and([P1[:,0] == P2[:,0]],[P1[:,1] == P2[:,1]]).T
        parallel_fix = numpy.zeros((n,3),dtype=bool)
        for i,z in enumerate(parallel_cond):
            parallel_fix[i] = z
            
        inter1[parallel_fix] = numpy.inf
        inter2[parallel_fix] = numpy.inf


        #and where there is a total miss, we want an inf, not a NaN
        miss_result = numpy.array(numpy.ones([n,3])*numpy.inf).T
        miss_fix = d<0
        
        inter1[miss_fix] = numpy.inf
        inter2[miss_fix] = numpy.inf
        
        
        # Now, are the intersections along the direction of travel?
        
        s = P2[:,:2] - P1[:,:2]     #the 2d summed vector: v1 + vs = v2
        len_s = compute_length(P1[:,:2], P2[:,:2]) #and 2d length
        dead_ray = len_s == 0       #rays of length = 0 have nonsense normals 
        s_n = s         #initialize the array
        for i,z in enumerate(s):           #normalize the vectors
            if dead_ray[i]:
                s_n[i] = numpy.zeros(s_n.shape[1])
            else:
                a = s[i,:]/len_s[i]
                s_n[i] = a

        s1 = inter1[:,:2] - P1[:,:2]
        len_s1 = compute_length(P1[:,:2], inter1[:,:2]) 
        dead_ray = len_s1 == 0       
        s1_n = s1         
        for i,z in enumerate(s1):           
            if dead_ray[i]:
                s1_n[i] = numpy.zeros(s1_n.shape[1])
            else:
                a = s1[i,:]/len_s1[i]
                s1_n[i] = a
        
        s2 = inter2[:,:2] - P1[:,:2]
        len_s2 = compute_length(P1[:,:2], inter2[:,:2]) 
        dead_ray = len_s2 == 0       
        s2_n = s2         
        for i,z in enumerate(s1):           
            if dead_ray[i]:
                s2_n[i] = numpy.zeros(s2_n.shape[1])
            else:
                a = s2[i,:]/len_s2[i]
                s2_n[i] = a

        #now use the normals to filter out intersections that are in the wrong direction
        
        backwards1 = numpy.zeros(s_n.shape,bool)
        backwards2 = backwards1.copy()
        
        #print inter1
        #print inter2
        
        # since both are vectors of length one in same dir or 180 deg apart,
        # addition should have len 2 or 0.
        
        for i,z in enumerate(s_n):
            
            temp = (s_n[i] + s1_n[i])**2
            backwards1[i] = sum(temp.T) < 1
            
            temp2 = (s_n[i] + s2_n[i])**2
            backwards2[i] = sum(temp2.T) < 1
            
        inter1[backwards1]=numpy.inf
        inter2[backwards2]=numpy.inf
        
        #print inter1
        #print inter2
        
        #now the z values can easily be interpolated:
        #change in z is proportional to total offest
        
        z1 = interpolate_z (P1, P2, inter1)
        z1[z1 -z1 != 0] = numpy.inf     #is z is a number, this will be false
        z2 = interpolate_z (P1, P2, inter2)
        z2[z2 - z2 != 0] = numpy.inf
        inter1[:,2]=z1
        inter2[:,2]=z2    
        
        print "2: ",inter1,"\n",inter2
        #now apply boundries based on begin and end angle.
        # for any given x or y value, there are two possible angles.
        # but only one will be resposible for both x and y.
        # so find the possible angles between 0 and 2pi and compare
        r = 1/R
        pi = numpy.pi
        
        y1_sin = numpy.arcsin(r*inter1[:,1])
        y1_thetas = numpy.column_stack((y1_sin, pi-y1_sin, (2*pi)+y1_sin))
        
        x1_cos = numpy.arccos(r*inter1[:,0])
        x1_thetas = numpy.column_stack((x1_cos, 2*pi-x1_cos))

        print "x und y: ",inter1[:,1],r*inter1[:,0]
        inter1_theta = numpy.ones(n)*numpy.inf
        for i,set in enumerate(y1_thetas):
            for phi in set:
                if any(abs(phi - x1_thetas[i]) < .01):
                    inter1_theta[i] = phi

        y2_sin = numpy.arcsin(r*inter2[:,1])
        y2_thetas = numpy.column_stack((y2_sin, pi-y2_sin, (2*pi)+y2_sin))
        
        x2_cos = numpy.arccos(r*inter2[:,0])
        x2_thetas = numpy.column_stack((x2_cos, 2*pi-x2_cos))
        
        inter2_theta = numpy.ones(n)*numpy.inf
        for i,set in enumerate(y2_thetas):
            for phi in set:
                if any(abs(phi - x2_thetas[i]) < .01):
                    inter2_theta[i] = phi
        
        # and use that angle to check for boundries
        #print "thetas: ",inter1_theta,inter2_theta
        #print "sections: ", inter2_theta
        small_angle = numpy.radians(self.begin_angle)
        big_angle = numpy.radians(self.end_angle)

        length = self.length
        
        boundry_mask1 = inter1_theta < small_angle
        boundry_mask1 = numpy.logical_or(boundry_mask1,inter1_theta>big_angle)
        boundry_mask1 = numpy.logical_or(boundry_mask1,inter1[:,2]>length)
        boundry_mask1 = numpy.logical_or(boundry_mask1,inter1[:,2]<0)
        
        boundry_mask1 = numpy.array([boundry_mask1]*3).T
        inter1[boundry_mask1] = numpy.inf
        
        boundry_mask2 = inter2_theta < small_angle
        boundry_mask2 = numpy.logical_or(boundry_mask2,inter2_theta>big_angle)
        boundry_mask2 = numpy.logical_or(boundry_mask2,inter2[:,2]>length)
        boundry_mask2 = numpy.logical_or(boundry_mask2,inter2[:,2]<0)
        
        boundry_mask2 = numpy.array([boundry_mask2]*3).T
        inter2[boundry_mask2] = numpy.inf

        print "3: ",inter1,"\n",inter2
        #print "boundy masks: \n", boundry_mask2, "\n", inter2  
        # next, use the distance from start to intersection to select the first 
        # intersections if there are multiple
        
        select = compute_length(P1, inter2) < compute_length(P1, inter1)

        #shortest = numpy.where(select, root1, root2)
        #mmm, numpy.where didn't like selecting vectors for some reason
        # So, I'll do it long hand
        select = compute_length(P1, inter2) < compute_length(P1, inter1)
        actual = inter1.copy()
        for i,n in enumerate(inter1):
            if select[i]:
                actual[i] = inter2[i,:]
            else:
                actual[i] = inter1[i,:]
                
        #finally, be sure the ray length to intersection is longer than the tolerance
        
        #tol_mask = self.compute_length(P1, actual) < self.tolerance
        
        #tol_mask = numpy.array([tol_mask]*3).T
        #actual[tol_mask] = numpy.inf
        
        dtype=([('length','f8'),('face', 'O'),('point','f8',3)])
        result = numpy.empty(P1.shape[0], dtype=dtype)
        result['length'] = compute_length(P1,actual)
        
        #flip y and z back to the standard order
        actual[:,1:] = numpy.fliplr(actual[:,1:]).copy()
        result['face'] = self
        result['point'] = actual
        return result

class CylinderMirrorFace(CylinderFace, PECFace):
    pass


class Cylinder(BaseMirror):
    """
    just a tube.  Can be a semi-tube, but begin and end angle must be between 0 
    and 360.  also, begin must be smaller than end angle.  This limits 
    functionality, and will be fixed in the future if thisfunctionality is ever 
    needed.
    """
    
    name = "Cylinder"
    length = Float(100.0, desc="length of tube")
    radius = Float(50.8, desc="radius of cylinder")
    begin_angle = Float(0., desc="angle at which to start making tube")
    end_angle = Float(90., desc="angle at which to stop making tube")
    resolution = Int(20, desc="domain of angle is broken into this many values")
    max_length = Float(1000.0)
    
    body = Instance(tvtk.ProgrammableSource, ())
    extruder = Instance(tvtk.LinearExtrusionFilter, ())
    
    traits_view = View(VGroup(
                        Traceable.uigroup,
                       Item('length', editor=NumEditor),
                       Item('width', editor=NumEditor),
                       Item('radius', editor=NumEditor),
                       Item('max_length', editor=NumEditor),
                       ),
                       )
    
    
    def calc_profile(self):
        output = self.body.poly_data_output
        
        start_angle = self.begin_angle
        end_angle = self.end_angle
        R = self.radius
        
        #create the 2d profile of the involute
        #r1 is vector from center of tube to outside of tube)
        #r2 is vector tangent to tube at r1 with length of a string wound around tube
        #surface is r1+r2
        #this is 2d so z = 0
        
        size = self.resolution
        theta = numpy.linspace(numpy.radians(start_angle),numpy.radians(end_angle),size).T
        
        #define involute (function of theta and starting parameters only)
        r1_x, r1_z = R*numpy.cos(theta), R*numpy.sin(theta)
        
        ys = numpy.zeros(r1_x.shape)
        
        profile = numpy.column_stack((r1_x, ys, r1_z))
        
        points = profile
        cells = [[i,i+1] for i in xrange(size-1)]
        output.points = points
        output.lines = cells
        return output
    
    def _vtkproperty_default(self):
        return tvtk.Property(opacity = 0.7,
                             color = (0.8,0.8,0))
    
    def _faces_default(self):
        return [CylinderMirrorFace(owner=self)]

# This is not necessary unless you want to export STEP files    
#    def make_step_shape(self):
#        from raytrace.step_export import make_OAP
#        return make_OAP(self.EFL, self.diameter, self.height,
#                        self.centre, self.direction, self.x_axis), "yellow"
    
                                         
    def _pipeline_default(self):

        
        self.body.set_execute_method(self.calc_profile)

        extrude = self.extruder
        extrude.input = self.body.output
        extrude.extrusion_type = "vector"
        extrude.vector = (0,1,0)
        extrude.scale_factor = self.length
        
        # cut parabolics.py here and inserted from prisms.py
        t = self.transform
        transF = tvtk.TransformFilter(input=extrude.output, transform=t)


        return transF
    
    def trace_segment(self, seg, last_optic=None, last_cell=None):
        """
        Don't care about last_optic or last_cell. We filter out
        intersections too close to the segment origin
        """
        p1 = seg.origin
        p2 = p1 + seg.MAX_RAY_LENGTH*seg.direction
        i = self.intersect_with_line(p1, p2)
        if i is None:
            return None
        i = numpy.array(i)
        dist = numpy.sqrt(((i-p1)**2).sum())
        return dist, i, 0, self
        
        
        
        
class InvoluteFace(Face):
    name = "Cylindical Involute Face"
    length = PrototypedFrom("owner")
    width = PrototypedFrom("owner")
    tube_radius = PrototypedFrom("owner")
    begin_angle = PrototypedFrom("owner")
    end_angle = PrototypedFrom("owner")
    resolution = PrototypedFrom("owner")
    max_length = PrototypedFrom("owner")
    
    intel_guess = numpy.array([0])
    #length = 100.
    #width = 30.
    #tube_radius = 6.
    #begin_angle = 0.
    #end_angle = 500.
    #resolution = 30
    #max_length = 100.
    
        
    
    def compute_length(self, start, end):
        #takes [n,3] vectors and returns [1,n] array of distances
        a = start - end
        d = a**2
        e = sum(d.T)
        distance = numpy.sqrt(e)
        
        mask = distance < self.tolerance
        distance[mask]=numpy.Infinity

        return distance      
        
    #def eval_children(self, rays, points, mask=slice(None,None,None)):
     #   return None
        
    def interpolate_z(self, first, threeD, twoD):
        # takes an 3d origin, a 3d point on the line, and intrpolates the third
        # dimesnsion of another point on that line, fow which x and y are given 
        # and z is 0
        
        # len2d1/len2d2 = LenInZ1/LenInZ2
    
        len2d1 = self.compute_length(first[:,:2], threeD[:,:2])
        len2d2 = self.compute_length(first[:,:2], twoD[:,:2])
        
        k = len2d2 / len2d1
        
        #Zf = Z1 - k*ChangeinZ         
        z = first[:,2] - (first[:,2]-threeD[:,2])*k
        
        return  z

    def involine(self, theta, args):
        """function for fsolve to solve for theta given an intersecting line"""
        X = theta
        m, b, E, R = args
        
        sin = numpy.sin
        cos = numpy.cos
        
        z = numpy.array((E-R*X + m*R)*sin(X) + (m*E - m*R*X - R)*cos(X) - b)
        return z
    
    def pts2theta(self, theta, xyz):
        """function for fsolve to find theta given an intersection point"""        
        W=self.width
        R=self.tube_radius
        
        x = xyz[0][0]
        y = xyz[0][1]
            
        r1_x, r1_y = -R*numpy.sin(theta),-R*numpy.cos(theta)
        r2_mag = numpy.abs(W - R * theta)
        r2_x, r2_y = -r2_mag*numpy.cos(theta), r2_mag*numpy.sin(theta)
        
        x0 = x - (r1_x + r2_x)
        y0 = y - (r1_y + r2_y)
        
        #distance = numpy.array([x0**2 + y0**2])
        return numpy.array([x0[0],y0[0]])
        #return numpy.array([distance])

    def compute_normal(self, points):
        """
        evaluate normalised surface Normal vector
        """
        n = points.shape[0]         #returns how many points there are
        W = self.width
        R = self.tube_radius
        t = self.transform
        
        #The intelligent guess relies on the thetas used to find these points 
        #in the intesect method.
        angles = self.intel_guess
        
        #this is a terrible guess:
        #guess = (self.begin_angle+self.end_angle)/2
        #guess = numpy.radians(guess)
        inv_t = t.linear_inverse
        t_points =transformPoints(inv_t, points)
        
        r2_x = numpy.zeros(n)
        r2_y = r2_x.copy()
        
        for j,point in enumerate(t_points):
            result = fsolve(self.pts2theta,angles[j],[point], full_output=1)
            theta = result[0]
            #print "theta:\n",theta
            #print result[2]  #error if this != 1
            r2_mag = numpy.abs(W - R * theta)
            r2_x[j], r2_y[j] = -r2_mag*numpy.cos(theta), r2_mag*numpy.sin(theta)        
        
        n_x = -r2_x
        n_y = -r2_y
        n_z = numpy.zeros(n_x.shape)
        t_normal = numpy.column_stack((n_x, n_y, n_z))
        
        
        #print t_normal
#        coefs = self.owner.vtk_quadric.coefficients
#        ax = 2*coefs[0]/coefs[8]
#        ay = 2*coefs[1]/coefs[8]
#        t_normal = numpy.column_stack((ax*t_points[:,0], 
#                                       ay*t_points[:,1],
#                                       -numpy.ones(n)))
        return transformNormals(t, t_normal)
    
        

            
    def intersect(self, P1, P2, max_length):
        """
        
        @param p1: a (n,3) array of points, start of each ray
        @param p2: a (n,3) array of point, ends of the rays
        """
        R = self.tube_radius
        W = self.width
        intel_guess = self.intel_guess
        begin_angle = numpy.radians(self.begin_angle)
        end_angle = numpy.radians(self.end_angle)
        n = P1.shape[0]         #returns how many points there are

        #turn array of points into y = mx + q
        m = (P1[:,1]-P2[:,1])/(P1[:,0]-P2[:,0]) # m = y1-y2 / x1-x2
        q = P1[:,1]-m*P1[:,0]                   #q = y - mx
        
        #define involute (function of theta and starting parameters only)
        #r1_x, r1_y = -R*numpy.sin(theta),-R*numpy.cos(theta)
        #r2_mag = numpy.abs(W - R * theta)
        #r2_x, r2_y = -r2_mag*numpy.cos(theta), r2_mag*numpy.sin(theta)
        
        #Where do these lines intersect with the involute? solve for theta when 
        #line = involute
        
        #there are only two solutions per 2Pi, so guessing every pi will give
        #every possible result.
        # Fixed: fsolve wants a better guess than that, so I used higher resolution
        guesses = numpy.arange(begin_angle,end_angle, numpy.pi/8)
        guesses = numpy.append(guesses, end_angle)
        
    
        #give output array the right shape
        intersect = numpy.ones([guesses.shape[0],n,1])*numpy.Inf
        

        #calculate all the thetas at which lines intersect involute
        for j,guess in enumerate(guesses):
            for i,z in enumerate(m):
                result = fsolve(self.involine,[guess],[m[i],q[i],W,R], full_output=1)
                #print "guess, m, q, W, R", guess,m[i],q[i],W,R
                #print "??:",result[0]
                bounds = numpy.logical_and(result[0] >= begin_angle, result[0]<= end_angle)
                #make sure fsolve found an acceptable answer
                if numpy.logical_and(result[2]==1, bounds):                
                    intersect[j,i] = result[0]
                else:
                    intersect[j,i] = numpy.Inf
                    
        #intersect now contains all the possible intersections by their theta value
        #To calculate distances and directions of travel, we need x,y and z.
        #actually, leave z out of it for now.
    
        intersection_points = numpy.zeros([guesses.shape[0],n,2])
        #print "points shape:",intersection_points.shape
        intel_guess = intersect.copy()
        #print "guess shape:",intel_guess.shape
        
        for j,theta in enumerate(intersect):
            #print "theta for ",j
            #print theta
            r1_x, r1_y = -R*numpy.sin(theta),-R*numpy.cos(theta)
            r2_mag = numpy.abs(W - R * theta)
            r2_x, r2_y = -r2_mag*numpy.cos(theta), r2_mag*numpy.sin(theta)
            x = r1_x + r2_x
            y = r1_y + r2_y
            z = 0
            for i,z in enumerate(intersection_points[j]):
                intersection_points[j,i] = [x[i][0],y[i][0]]
        #print "intersect:"
        #print intersection_points
        # Now, are the intersections along the direction of travel?
        
        s = P2[:,:2] - P1[:,:2]     #the 2d summed vector: v1 + vs = v2
        len_s = self.compute_length(P1[:,:2], P2[:,:2]) #and 2d length        
        dead_ray = len_s == 0       #rays of length = 0 have nonsense normals 
        s_n = s.copy()         #initialize the array
        for i,z in enumerate(s):           #normalize the vectors
            if dead_ray[i]:
                s_n[i] = numpy.zeros(s_n.shape[1])
            else:
                a = s[i,:]/len_s[i]
                s_n[i] = a

        s1 = numpy.ones(intersection_points.shape)*numpy.Inf
        len_s1 = numpy.ones([guesses.shape[0],n])*numpy.Inf
        s1_n = s1.copy()
        
        for j,inter in enumerate(intersection_points):
            s1[j] = inter - P1[:,:2]
            len_s1[j] = self.compute_length(P1[:,:2], inter)
            dead_ray = len_s1[j] == 0       
            for i,z in enumerate(s1[j]):
                if dead_ray[i]:
                    s1_n[j,i] = numpy.zeros(s1_n.shape[1])
                else:
                    a = s1[j,i,:]/len_s1[j,i]
                    s1_n[j,i] = a
        #print "s1_n:"
        #print s1_n
        #now use the normals to filter out intersections that are in the wrong direction
        
        backwards = numpy.zeros(s1_n.shape,bool)
        
        # since both are vectors of length one in same dir or 180 deg apart,
        # addition should have len 2 or 0.
        
        for j,v in enumerate (s1_n):
            for i,z in enumerate(s1_n[j]):
                 temp = (s_n[i] + s1_n[j,i])**2
                 backwards[j,i] = sum(temp.T) < 1
            
        intersection_points[backwards]=numpy.inf
        guess_mask = backwards[:,:,:1]
        intel_guess[guess_mask] = numpy.Inf
        
        #now the z values can easily be interpolated:
        #change in z is proportional to total offest
        third_d = numpy.zeros([guesses.shape[0],n,1])
        intersection_points = numpy.append(intersection_points,third_d,axis=2)
        
        for j,point in enumerate(intersection_points):
            z1 = self.interpolate_z (P1, P2, point)
            z1[z1 -z1 != 0] = numpy.inf     #is z is a number, this will be false
            point[:,2]=z1
        
        #Some of these supposed intersections dont actually hit the shape
        #within the given bounds.
        #x and y bounds were taken care of by limits on theta earlier.
        
        Z_bounds = numpy.array([0,self.length])
    
        zmin, zmax = min(Z_bounds), max(Z_bounds)
        
        for i,points in enumerate(intersection_points): 
            bounds_mask1 = numpy.zeros(points[:,0].shape,dtype=bool)
        
            bounds_mask1 = numpy.logical_or(bounds_mask1, points[:,2]<zmin)
            bounds_mask1 = numpy.logical_or(bounds_mask1, points[:,2]>zmax)
            bounds_mask1 = numpy.array([bounds_mask1]*3).T
            points[bounds_mask1] = numpy.inf
            guess_mask = bounds_mask1[:,:1]
            intel_guess[i][guess_mask] = numpy.inf
        # next, use the distance from start to intersection to select the first 
        # intersections if there are multiple
        
        
        actual = intersection_points[0]*numpy.inf
        best_guess = intel_guess[0].copy()
        for j,points in enumerate(intersection_points):
            select = self.compute_length(P1, actual) < self.compute_length(P1, points)
            for i,n in enumerate(actual):
                if select[i]:
                    pass
                else:
                    actual[i] = points[i,:]
                    best_guess[i] = intel_guess[j,i]
         
        
        best_guess = best_guess[numpy.where(best_guess != numpy.inf)]
        self.intel_guess = best_guess
        #finally, be sure the ray length to intersection is longer than the tolerance
        
        #tol_mask = self.compute_length(P1, actual) < self.tolerance
        
        #tol_mask = numpy.array([tol_mask]*3).T
        #actual[tol_mask] = numpy.inf
        
        
        dtype=([('length','f8'),('face', 'O'),('point','f8',3)])
        result = numpy.empty(P1.shape[0], dtype=dtype)
        result['length'] = self.compute_length(P1,actual)
        result['face'] = self
        result['point'] = actual
        return result
    
        
class CyInvoluteMirrorFace(InvoluteFace, PECFace):
    pass


class CylindricalInvolute(BaseMirror):
    """
    An involute of a cylinder.  Starts at (0,-y) (theta=0) and computes the involute
    till end_angle.  If user wraps too short a string around the cylinder, the
    magnitude of the string becomes negative and a meaningless shape results.  Though
    no legitimate use of an involute would come across this.
    """
    name = "Cylindrical Involute"
    length = Float(100.0, desc="length of trough")
    width = Float(25.4, desc="radius of aperature")
    tube_radius = Float(50.8, desc="radius of cylinder involute is formed from")
    begin_angle = Float(0., desc="angle at which to start making involute")
    end_angle = Float(90., desc="angle at which to stop taking involute")
    resolution = Int(20, desc="domain of angle is broken into this many values")
    max_length = Float(1000.0)
    
    body = Instance(tvtk.ProgrammableSource, ())
    extruder = Instance(tvtk.LinearExtrusionFilter, ())    
    
    traits_view = View(VGroup(
                        Traceable.uigroup,
                       Item('length', editor=NumEditor),
                       Item('width', editor=NumEditor),
                       Item('tube_radius', editor=NumEditor),
                       Item('max_length', editor=NumEditor),
                       ),
                       )
    
    
    def calc_profile(self):
        output = self.body.poly_data_output
        
        start_angle = self.begin_angle
        end_angle = self.end_angle
        R = self.tube_radius
        W = self.width
        
        #create the 2d profile of the involute
        #r1 is vector from center of tube to outside of tube)
        #r2 is vector tangent to tube at r1 with length of a string wound around tube
        #surface is r1+r2
        #this is 2d so z = 0
        
        size = self.resolution
        theta = numpy.linspace(numpy.radians(start_angle),numpy.radians(end_angle),size).T
        
        #define involute (function of theta and starting parameters only)
        r1_x, r1_y = -R*numpy.sin(theta),-R*numpy.cos(theta)
        r2_mag = numpy.abs(W - R * theta)
        r2_x, r2_y = -r2_mag*numpy.cos(theta), r2_mag*numpy.sin(theta)
        
        zs = numpy.zeros(r1_x.shape)
        
        r1 = numpy.column_stack((r1_x, r1_y, zs))
        r2 = numpy.column_stack((r2_x, r2_y, zs))
        
        profile = r1 + r2
        
        
    
        points = profile
        cells = [[i,i+1] for i in xrange(size-1)]
        output.points = points
        output.lines = cells
        return output
    
    def _vtkproperty_default(self):
        return tvtk.Property(opacity = 0.7,
                             color = (0.8,0.8,0))
    
    def _faces_default(self):
        return [CyInvoluteMirrorFace(owner=self)]

# This is not necessary unless you want to export STEP files    
#    def make_step_shape(self):
#        from raytrace.step_export import make_OAP
#        return make_OAP(self.EFL, self.diameter, self.height,
#                        self.centre, self.direction, self.x_axis), "yellow"
    
                                         
    def _pipeline_default(self):

        
        self.body.set_execute_method(self.calc_profile)

        extrude = self.extrude
        extrude.input=self.body.output
        extrude.extrusion_type = "vector"
        extrude.vector = (0,0,1)
        extrude.scale_factor = self.length
        
        # cut parabolics.py here and inserted from prisms.py
        t = self.transform
        transF = tvtk.TransformFilter(input=extrude.output, transform=t)


        return transF
    
    def trace_segment(self, seg, last_optic=None, last_cell=None):
        """
        Don't care about last_optic or last_cell. We filter out
        intersections too close to the segment origin
        """
        p1 = seg.origin
        p2 = p1 + seg.MAX_RAY_LENGTH*seg.direction
        i = self.intersect_with_line(p1, p2)
        if i is None:
            return None
        i = numpy.array(i)
        dist = numpy.sqrt(((i-p1)**2).sum())
        return dist, i, 0, self
        
if __name__=="__main__":

    p1 = numpy.array([[1.,1.,1.],[2.,3.,2.],[-3,0,0]])
    p2 = numpy.array([[1.,12.,12.],[10.,11.,10],[-3,4,6]])
    d = CylinderFace()
    c = d.intersect(p1,p2,100)
    