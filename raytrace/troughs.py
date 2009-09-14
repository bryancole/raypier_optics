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
Parabolic troughs and also rectangular mirrors
"""
from enthought.traits.api import HasTraits, Array, Float, Complex,\
            Property, List, Instance, on_trait_change, Range, Any,\
            Tuple, Event, cached_property, Set, Int, Trait, PrototypedFrom
from enthought.traits.ui.api import View, Item, ListEditor, VSplit,\
            RangeEditor, ScrubberEditor, HSplit, VGroup, TupleEditor
from enthought.tvtk.api import tvtk
import numpy
import math
from itertools import izip

from raytrace.bases import Traceable, normaliseVector, NumEditor,\
                transformNormals, transformPoints
from raytrace.mirrors import BaseMirror
from raytrace.faces import PECFace, RectangularFace, Face




class TroughFace(Face):
    name = "Parabolic Trough Face"
    EFL = PrototypedFrom("owner")
    #EFL = 5.
    length = PrototypedFrom("owner")
    #length = 100
    X_bounds = PrototypedFrom("owner")
    #width = 20
    
    #takes [n,3] vectors and returns [1,n] array of distances
    def compute_length(self, start, end):
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

        
    def compute_normal(self, points):
        """
        evaluate surface normal vector
        """
        n = points.shape[0]         #returns how many points there are
        t = self.transform
        inv_t = t.linear_inverse
        t_points =transformPoints(inv_t, points)
        h = 1/ (4 * self.EFL)
        
        #y and z are flipped so that y axis is the long axis
        t_points[:,1:] = numpy.fliplr(t_points[:,1:]).copy()
        
        # y = hx^2 - efl
        # dy/dx = 2hx
        
        surface_slope = 2*h*t_points[:,0]

        # perpendicular to line of slope m has slope -1/m
        
        perp_slope = -1 / surface_slope
        
        #then find a line of length 1 with this slope:
        #start with (1,m) and normalize
        
        r_sq =  1 + perp_slope**2
        L = numpy.sqrt(r_sq)
        
        #to keep normals pointing inward, points with pos x will have neg x normals
        x_sign = -1 * t_points[:,0] / numpy.abs(t_points[:,0])
        
        n_x = x_sign / L
        n_y = numpy.abs(perp_slope) / L
        n_z = numpy.zeros(n_x.shape)                 #surface doesn't point into Z at all
        
        #fix where shape was flat and normal has slope = inf
        oops = surface_slope == 0
        n_x[oops] = 0.
        n_y[oops] = 1.
        
        #notice, y and z are flipped back.
        t_normal = numpy.column_stack((n_x, n_z, n_y))
        
        
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
        n = P1.shape[0]         #returns how many points there are
        efl = self.EFL #scalar
        h = 1 / (4 * efl)
        
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

        #solve intersection of y = mx + b and y = h x^2 - EFL
        # 0 = hx^2 - mx - (b + EFL)
        # h = 1/(4*EFL)
        
        a = h
        b = -m
        c = - q - efl
        
        d = b**2 - 4*a*c
        
        e = numpy.sqrt(d)
        roots = [(-b+e)/(2*a), (-b-e)/(2*a)]
        
        root1, root2 = roots
        
        #put these roots into a list of intersection points using y = mx + q
        #I make these 3d with z=0, Which i'll fix later
        inter1 = numpy.array([root1,h*(root1**2) - efl, numpy.zeros(n)]).T
        inter2 = numpy.array([root2,h*(root2**2) - efl, numpy.zeros(n)]).T
        
        #Where the slope was infinite these values are wrong:
        
        perp_result = numpy.array([P1[:,0],h * P1[:,0]**2 - efl ,numpy.zeros(n)]).T
        perp_fix = numpy.array([P1[:,0] == P2[:,0]]*3).T
        
        inter1 = numpy.where(perp_fix,perp_result, inter1)
        inter2[perp_fix] = numpy.inf
        
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
        len_s = self.compute_length(P1[:,:2], P2[:,:2]) #and 2d length
        dead_ray = len_s == 0       #rays of length = 0 have nonsense normals 
        s_n = s         #initialize the array
        for i,z in enumerate(s):           #normalize the vectors
            if dead_ray[i]:
                s_n[i] = numpy.zeros(s_n.shape[1])
            else:
                a = s[i,:]/len_s[i]
                s_n[i] = a

        s1 = inter1[:,:2] - P1[:,:2]
        len_s1 = self.compute_length(P1[:,:2], inter1[:,:2]) 
        dead_ray = len_s1 == 0       
        s1_n = s1         
        for i,z in enumerate(s1):           
            if dead_ray[i]:
                s1_n[i] = numpy.zeros(s1_n.shape[1])
            else:
                a = s1[i,:]/len_s1[i]
                s1_n[i] = a
        
        s2 = inter2[:,:2] - P1[:,:2]
        len_s2 = self.compute_length(P1[:,:2], inter2[:,:2]) 
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
        
        z1 = self.interpolate_z (P1, P2, inter1)
        z1[z1 -z1 != 0] = numpy.inf     #is z is a number, this will be false
        z2 = self.interpolate_z (P1, P2, inter2)
        z2[z2 - z2 != 0] = numpy.inf
        inter1[:,2]=z1
        inter2[:,2]=z2
        
        #Some of these supposed intersections dont actually hit the shape
        #within the given bounds
        
        X_bounds = self.X_bounds
        # Y_bounds, determined by x bounds: y = x^2
        Z_bounds = numpy.array([0,self.length])
    
        xmin, xmax = min(X_bounds), max(X_bounds)
        zmin, zmax = min(Z_bounds), max(Z_bounds)
        
        bounds_mask1 = numpy.zeros(inter1[:,0].shape,dtype=bool)
        bounds_mask2 = numpy.zeros(inter2[:,0].shape,dtype=bool)
        
        bounds_mask1 = inter1[:,0] < xmin
        bounds_mask1 = numpy.logical_or(bounds_mask1, inter1[:,0]>xmax)
        bounds_mask1 = numpy.logical_or(bounds_mask1, inter1[:,2]<zmin)
        bounds_mask1 = numpy.logical_or(bounds_mask1, inter1[:,2]>zmax)
        bounds_mask1 = numpy.array([bounds_mask1]*3).T
        inter1[bounds_mask1] = numpy.inf
        
        bounds_mask2 = inter2[:,0] < xmin
        bounds_mask2 = numpy.logical_or(bounds_mask2, inter2[:,0]>xmax)
        bounds_mask2 = numpy.logical_or(bounds_mask2, inter2[:,2]<zmin)
        bounds_mask2 = numpy.logical_or(bounds_mask2, inter2[:,2]>zmax)
        
        
        bounds_mask2 = numpy.array([bounds_mask2]*3).T
        inter2[bounds_mask2] = numpy.inf
        
        # next, use the distance from start to intersection to select the first 
        # intersections if there are multiple
        
        select = self.compute_length(P1, inter2) < self.compute_length(P1, inter1)

        #shortest = numpy.where(select, root1, root2)
        #mmm, numpy.where didn't like selecting vectors for some reason
        # So, I'll do it long hand
        select = self.compute_length(P1, inter2) < self.compute_length(P1, inter1)
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
        result['length'] = self.compute_length(P1,actual)
        
        #flip y and z back to the standard order
        actual[:,1:] = numpy.fliplr(actual[:,1:]).copy()
        result['face'] = self
        result['point'] = actual
        return result

class TroughMirrorFace(TroughFace, PECFace):
    pass


class TroughParabloid(BaseMirror):
    """
    An trough mirror object
    """
    name = "trough"
    length = Float(100.0, desc="length of trough")
    X_bounds = Tuple((-25.5, 25.5), desc="width of parabolic profile")
    EFL = Float(50.8, desc="effective focal length")
    
    max_length = Float(1000.0)
    
    body = Instance(tvtk.ProgrammableSource, ())
    extrude = Instance(tvtk.LinearExtrusionFilter, ())
    
    traits_view = View(VGroup(
                        Traceable.uigroup,
                       Item('length', editor=NumEditor),
                       Item('width', editor=NumEditor),
                       Item('EFL', editor=NumEditor),
                       Item('X_bounds', editor=TupleEditor(cols=2, auto_set=False,
                                                           editors=[NumEditor, NumEditor],
                                                           labels=['min','max'])),
                       Item('max_length', editor=NumEditor),
                       ),
                       )
    
    @on_trait_change("X_bounds, EFL")
    def change_params(self):
        self.body.modified()
        self.update=True
        
    def _length_changed(self, l):
        self.extrude.scale_factor = l
        self.update=True
    
    def calc_profile(self):
        output = self.body.poly_data_output
        xmin, xmax = min(self.X_bounds), max(self.X_bounds)
        a = 1 / (4 * self.EFL)
        
        #create the 2d profile of parabolic trough
        size = 20
        x = numpy.linspace(xmin, xmax, size)
        z = a * (x**2) -self.EFL
        y = numpy.zeros_like(x)         #this is a 2d profile.  so, no Y
    
        points = numpy.array([x,y,z]).T 
        cells = [[i,i+1] for i in xrange(size-1)]
        output.points = points
        output.lines = cells
        return output
    
    def _vtkproperty_default(self):
        return tvtk.Property(opacity = 0.7,
                             color = (0.8,0.8,0))
    
    def _faces_default(self):
        return [TroughMirrorFace(owner=self)]

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
        
#
# ---------------------------------------------------
#

class RectMirrorFace(RectangularFace, PECFace):
    pass


class RectMirror(BaseMirror):
    """
    A rectangular mirror object
    """
    name = "rectangular mirror"
    length = Float(100.0, desc="length of trough")
    width = Float(-25.5, desc="width of parabolic profile")
    
    max_length = Float(1000.0)
    
    body = tvtk.ProgrammableSource() 
    
    traits_view = View(VGroup(
                        Traceable.uigroup,
                       Item('length', editor=NumEditor),
                       Item('width', editor=NumEditor),
                       Item('max_length', editor=NumEditor),
                       ),
                       )
    
    
    def calc_profile(self):
        output = self.body.poly_data_output
        xmin, xmax = -self.width/2, self.width/2
        size = 2
        #create the 2d profile, just a line.
        x = numpy.array([xmin, xmax])
        z = numpy.zeros_like(x)
        y = numpy.zeros_like(x)         #this is a 2d profile.  so, no Y
    
        points = numpy.array([x,y,z]).T 
        print points
        cells = [[i,i+1] for i in xrange(size-1)]
        output.points = points
        output.lines = cells
        return output
    
    def _vtkproperty_default(self):
        return tvtk.Property(opacity = 0.7,
                             color = (0.8,0.8,0))
    
    def _faces_default(self):
        return [RectMirrorFace(owner=self)]

# This is not necessary unless you want to export STEP files    
#    def make_step_shape(self):
#        from raytrace.step_export import make_OAP
#        return make_OAP(self.EFL, self.diameter, self.height,
#                        self.centre, self.direction, self.x_axis), "yellow"
    
                                         
    def _pipeline_default(self):

        
        self.body.set_execute_method(self.calc_profile)

        extrude = tvtk.LinearExtrusionFilter(input=self.body.output)
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
    
if __name__=="__main__":
    oap = TroughParabloid()
    
    mapper = tvtk.PolyDataMapper(input = oap.pipeline.output)
    actor = tvtk.Actor(mapper=mapper)
    ren = tvtk.Renderer()
    ren.add_actor(actor)
    
    ax = tvtk.Axes(origin=(0,0,0))
    axes_map = tvtk.PolyDataMapper(input=ax.output)
    axes_act = tvtk.Actor(mapper=axes_map)
    ren.add_actor(axes_act)
    
    ren.background=(0.7,0.6,0.5)
    
    renwin = tvtk.RenderWindow()
    renwin.add_renderer(ren)
    
    iren = tvtk.RenderWindowInteractor(render_window=renwin)
    iren.start()
