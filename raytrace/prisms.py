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


from enthought.traits.api import HasTraits, Array, Float, Complex,\
            Property, List, Instance, on_trait_change, Range, Any,\
            Tuple, Event, cached_property, Set, Int, Trait
from enthought.traits.ui.api import View, Item, ListEditor, VSplit,\
            RangeEditor, ScrubberEditor, HSplit, VGroup, Heading
from enthought.tvtk.api import tvtk
from enthought.tvtk.pyface.scene_model import SceneModel
from enthought.tvtk.pyface.scene_editor import SceneEditor
import numpy
from itertools import chain, izip, islice

from raytrace.tracer import Optic, VTKOptic, normaliseVector, RaySegment,\
             Traceable, NumEditor, dotprod, transformPoints, transformNormals


class Extrusion(Optic):
    """a general flat-faced optic formed by extrusion 
    of a 2D polygon"""
    profile = Array(shape=(None,2), dtype=numpy.double)
    length = Float(10.0)
    
    _face_normals = Property(Array(shape=(None,3),dtype=numpy.double),
                             depends_on="profile")
    
    data_source = Instance(tvtk.ProgrammableSource, ())
    
    extrude = Instance(tvtk.LinearExtrusionFilter, (), 
                       {'capping': True, 
                        'extrusion_type':'vector',
                        'vector': (0.,0.,1.)})
    
    @cached_property
    def _get__face_normals(self):
        profile = self.profile
        ends = numpy.roll(profile, -1, axis=0)
        t = ends - profile #edge direction
        normals = numpy.column_stack((t[:,1], -t[:,0], numpy.zeros(t.shape[0])))
        return normaliseVector(normals)
    
    @on_trait_change("length")
    def config_pipeline(self):
        self.extrude.scale_factor = self.length
        self.update=True
        
    def _profile_changed(self):
        self.data_source.modified()
        self.update=True
    
    def _pipeline_default(self):
        source = self.data_source
        def execute():
            xy = self.profile
            z = numpy.zeros(xy.shape[0])
            points = numpy.column_stack((xy,z))
            
            cells = [range(len(z)),]
            
            output = source.poly_data_output
            output.points = points
            output.polys = cells
        source.set_execute_method(execute)
        
        extrude = self.extrude
        extrude.input = source.output
        
        t = self.transform
        transf = tvtk.TransformFilter(input=extrude.output, transform=t)
        return transf
    
    def compute_normal(self, points, cells):
        """reimplemented from base class"""
        t = self.transform
        face_normals = transformNormals(t, self._face_normals)
        normals = face_normals[cells,:]
        return normals
    
    def trace_rays(self, rays):
        """re-implemented from base class:
        traces 'rays', a RayCollection instance.
        returns - a recarray of intersetions with the same size as rays
                  dtype=(('length','f8'),('cell','i2'),('point','f8',[3,]))
                  
                  'length' is the physical length from the start of the ray
                  to the intersection
                  
                  'cell' is the ID of the intersecting face
                  
                  'point' is the position coord of the intersection
        """
        max_length = rays.max_length
        p1 = rays.origin
        p2 = p1 + max_length*rays.direction
        t = self.transform
        inv_t = t.linear_inverse
        P1 = transformPoints(inv_t, p1)
        P2 = transformPoints(inv_t, p2)
        
        height = self.length
        
        ###TODO: not bothering with the end faces for the time being
        #h_top, intersect_top = self.intersect_with_end(P1, P2, height)
        #h_btm, intersect_btm = self.intersect_with_end(P1, P2, 0.0)
        
        PP1 = self.profile
        PP2 = numpy.roll(PP1, -1, axis=0)
        
        h, intersections, cid = self.intersect_with_side(P1, P2, PP1, PP2, height)
        
        length = max_length*h
        
        points = transformPoints(t, intersections)
    
        dtype=([('length','f8'),('cell','i2'),('point','f8',3)])
        result = numpy.empty(p1.shape[0], dtype=dtype)
        result['length'] = length
        result['point'] = points
        result['cell'] = cid
        
        return result
    
    def intersect_with_end(self, P1, P2, height):
        """
        Find the line intersection with a section across the extrusion,
        at the given height
        """
        x1,y1,z1 = P1.T
        x2,y2,z2 = P2.T
        
        h = (height-z1)/(z2-z1)
        X = x1 + h*(x2-x1)
        Y = y1 + h*(y2-y1)
        
        mask = h > 1.0
        mask = numpy.logical_or(mask, h < self.tolerance)
        
        edge_start = self.profile
        edge_end = numpy.roll(edge_start, -1, axis=0)
        
        ax = edge_start[:,0].reshape(-1,1)
        ay = edge_start[:,1].reshape(-1,1)
        
        dx = (edge_end[:,0]-edge_start[:,0]).reshape(-1,1)
        dy = (edge_end[:,1]-edge_start[:,1]).reshape(-1,1)
        
        y = (Y - ay)/dy
        x = ax + (y * dx)
        i_count = (x<X).sum(axis=0)
        half = i_count/2.
        outside = half==half.round()
        
        mask = numpy.logical_or(mask, outside)
        
        h[mask] = numpy.Infinity
        
        return h, numpy.column_stack((X, Y, numpy.ones_like(X)*height))
        
    def intersect_with_side(self, P1, P2, PP1, PP2, height):
        """
        
        @param P1: (n,3) array of n points, in extrusion reference frame
        @param P2: (n,3) array of points
        @param PP1: (m,2) array of m 2D start points of edge (of side)
        @param PP2: (m,2) array of m 2D end points of edge
        @param height: height of extrusion (scalar)
        """
        newaxis = numpy.newaxis
        a = PP1[newaxis,:,:] #shape = 1,m,2
        D = PP2-PP1 #shape = m,2
        d1 = normaliseVector(D) #shape = m,2
        
        P1_ = P1[:,newaxis,:]
        P2_ = P2[:,newaxis,:]
        
        r1 = P1_[:,:,:2] - a #shape = n,m,2
        h1 = numpy.cross(r1,d1[newaxis,:,:]) #shape = n,m
        
        r2 = P2_[:,:,:2] - a #shape = n,m,2
        h2 = numpy.cross(r2,d1[newaxis,:,:]) #shape = n,m
        
        h = h1/(h1-h2) #shape = n,m
        
        intersect = P1_ + h[:,:,newaxis] * (P2_-P1_) #shape = n,m,3
        
        z = intersect[:,:,2] #shape = n,m
        
        xy_intersect = intersect[:,:,:2] #shape = n,m,2
        edge_length = numpy.sqrt((D**2).sum(axis=-1)) #shape = m,
        fraction = dotprod(xy_intersect - a, d1[newaxis,:,:])[:,:,0]/edge_length[newaxis,:]
        
        mask = h > 1.0
        mask = numpy.logical_or(mask, h < self.tolerance) #shape=n,m
        mask = numpy.logical_or(mask, z < 0) #shape=n,m
        mask = numpy.logical_or(mask, z > height) #shape=n,m
        mask = numpy.logical_or(mask, fraction < 0) #shape=n,m
        mask = numpy.logical_or(mask, fraction > 1.0) #shape=n,m
        
        h[mask] = numpy.Infinity
        
        cell_id = h.argmin(axis=1) #shape = n,
        
        h_nearest = h.min(axis=1) #shape = n,
        
        n,m = intersect.shape[:2]
        nearest_intersections = intersect[numpy.arange(n),cell_id,:]
        
        return h_nearest, nearest_intersections, cell_id


class Prism(Extrusion):
    name = "prism"
    height = Float #distance from front face to apex
    width = Float #width of front face
    
    traits_view = View(VGroup(
                       Traceable.uigroup,
                       Item('all_rays'),
                       Item('n_inside'),
                       Item('length'),
                       Item('height'),
                       Item('width'),
                       )
                       )

    @on_trait_change("height, width")
    def config_profile(self):
        h = self.height
        w = self.width/2
        self.profile = [(-w,0),
                        (w,0),
                        (0,h)]
        
        
class Rhomboid(Extrusion):
    name = "rhomboid"
    
    height = Float #distance between parallel faces
    width = Float #width of parallel faces
    slant = Float #angle of the oblique faces
    
    traits_view = View(VGroup(
                       Traceable.uigroup,
                       Item('all_rays'),
                       Item('n_inside'),
                       Item('length'),
                       Item('height'),
                       Item('width'),
                       Item('slant'),
                       )
                       )

    @on_trait_change("height, width, slant")
    def config_profile(self):
        h = self.height/2
        angle = self.slant*numpy.pi/180
        s = h / numpy.tan(angle)
        w = self.width/2
        points = [(-w-s,-h),
                  (-w+s,h),
                  (w+s,h),
                  (w-s,-h)]
        points.reverse()
        self.profile = points


if __name__=="__main__":
    from ray_tracer import BuildRaySet, BeamStop, RayTraceModel
    
    input_rays = BuildRaySet(origin = (-7.42,15,0),
                             direction = (0,-1,0),
                             radius=1.0,
                             count=20)
    
    rhomboid = Rhomboid(n_inside = 1.764+0.0j,
                        orientation=0.0,
                        elevation=0.0,
                        centre=(0,0,0),
                        rotation=0,
                        length=10.0,
                        height=7.07,
                        width=14.0,
                        slant=45.0)
    
    beamstop = BeamStop(width=10,height=10,
                        centre=(7,-10,0),
                        direction=(0,1,0))
    
    #print "rhomboid", rhomboid.polydata
    print "beamstop", beamstop.polydata
    
    model = RayTraceModel(optics=[rhomboid, beamstop], rays=input_rays)
    model.configure_traits()
