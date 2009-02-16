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
a module for parabolic optics e.g. OAPs
"""
from enthought.traits.api import HasTraits, Array, Float, Complex,\
            Property, List, Instance, on_trait_change, Range, \
            Tuple, Event, cached_property, Set, Int, Trait, Bool
from enthought.traits.ui.api import View, Item, ListEditor, VSplit,\
            RangeEditor, ScrubberEditor, HSplit, VGroup
from enthought.tvtk.api import tvtk
import numpy

from raytrace.tracer import Traceable, normaliseVector, NumEditor,\
            VectorEditor, transformPoints, transformNormals
#from raytrace.sources import RayCollection
from raytrace.mirrors import BaseMirror


class Ellipsoid(BaseMirror):
    name = "Ellipsoid"
    focus1 = Tuple(-50,0,0)
    focus2 = Tuple(0, 50, 0)
    size = Float(100.0)
    X_bounds = Tuple(-25., 25.)
    Y_bounds = Tuple(-25., 25.)
    Z_bounds = Tuple(0., 50.)
    
    show_foci = Bool(True)
    
    axes = Property(Tuple, depends_on="focus1, focus2, size",
                 desc="(major, minor) axis lengths")
    
    ellipse_trans = Instance(tvtk.Transform, (), transient=True)
    
    combined_trans = Instance(tvtk.Transform, transient=True)
    
    f1_glyph = Instance(tvtk.SphereSource, (), transient=True)
    f2_glyph = Instance(tvtk.SphereSource, (), transient=True)
    
    f1_act = Instance(tvtk.Follower, (), transient=True)
    f2_act = Instance(tvtk.Follower, (), transient=True)
    
    vtk_grid = Instance(tvtk.ProgrammableSource, (), transient=True)
    vtk_quadric = Instance(tvtk.Quadric, (), transient=True)
    
    vtkproperty = tvtk.Property(opacity = 0.7,
                             color = (0.8,0.8,0))
    
    traits_view = View(VGroup(
                        Traceable.uigroup,
                        Item("focus1", editor=VectorEditor),
                        Item("focus2", editor=VectorEditor),
                        Item("show_foci"),
                        Item("size", editor=NumEditor),
                        Item("X_bounds"),
                        Item("Y_bounds"),
                        Item("Z_bounds")
                       ),
                       )
    
    @cached_property
    def _get_axes(self):
        f1 = numpy.asarray(self.focus1)
        f2 = numpy.asarray(self.focus2)
        delta = f2 - f1
        size = self.size
        h = numpy.sqrt((delta**2).sum())
        a = size/2.
        b = 0.5 * numpy.sqrt(size**2 - h**2)
        return (a,b)
    
    def _combined_trans_default(self):
        t = tvtk.Transform()
        t.concatenate(self.transform)
        t.concatenate(self.ellipse_trans.linear_inverse)
        return t
    
    def _show_foci_changed(self, val):
        self.f1_act.visibility = val
        self.f2_act.visibility = val
        self.render = True
    
    def create_grid(self):
        xmin, xmax = self.X_bounds
        ymin, ymax = self.Y_bounds
        zmin, zmax = self.Z_bounds
        
        source = self.vtk_grid
        sp = source.structured_points_output
        size = 20
        dx = (xmax-xmin) / (size-1)
        dy = (ymax-ymin) / (size-1)
        dz = (zmax-zmin) / (size-1)
        sp.dimensions = (size,size,size)
        sp.whole_extent=(0,size,0,size,0,size)
        sp.origin = (xmin, ymin, zmin)
        sp.spacing = (dx, dy, dz)
        sp.set_update_extent_to_whole_extent()
        
    @on_trait_change("X_bounds, Y_bounds, Z_bounds")
    def change_bounds(self):
        self.vtk_grid.modified()
        self.update=True
    
    @on_trait_change("focus1, focus2, size")
    def config_pipeline(self, *args):
        tp = self.transform.transform_point
        self.f1_act.position = tp(self.focus1)
        self.f2_act.position = tp(self.focus2)
        
        f1 = numpy.asarray(self.focus1)
        f2 = numpy.asarray(self.focus2)
        
        self.f1_glyph.center = f1
        self.f2_glyph.center = f2
        
        centre = (f2+f1)/2.
        ellipse_t = self.ellipse_trans
        ellipse_t.identity()
        #ellipse major axis along the X axis
        delta = f2 - f1
        ax = normaliseVector(delta)
        axy = numpy.sqrt(ax[0]**2 + ax[1]**2)
        ellipse_t.rotate_y(numpy.arctan2(ax[2], axy)*180/numpy.pi)
        
        ellipse_t.rotate_z(-numpy.arctan2(ax[1],ax[0])*180/numpy.pi)
        
        ellipse_t.translate(*-centre)
        
        a,b = self.axes
        #b = 0.5 * numpy.sqrt(size**2 - h**2)  ##not required
        
        q = self.vtk_quadric
        A1 = A2 = 1/(b**2)
        A0 = 1/(a**2)  
        #A8 = 1 
        A9 = -1
        q.coefficients = (A0, A1, A2, 
                          0, 0, 0, 
                          0, 0, 0, 
                          A9)
        self.update=True
        
                                         
    def _pipeline_default(self):
        grid = self.vtk_grid
        grid.set_execute_method(self.create_grid)
        grid.modified()
        
        quad = self.vtk_quadric
        quad.transform = self.ellipse_trans
        
        clip = tvtk.ClipVolume(input=grid.structured_points_output,
                                 clip_function=quad,
                                 inside_out=0)
        
        topoly = tvtk.GeometryFilter(input=clip.output)
        norm = tvtk.PolyDataNormals(input=topoly.output)
        transF = tvtk.TransformFilter(input=norm.output, transform=self.transform)
        self.config_pipeline()
        grid.modified()
        return transF
    
    def get_actors(self, scene):
        actors = self.actors
        
        sList = [self.f1_glyph, self.f2_glyph]
        cList = [(0,1,0),(1,0,0)]
        
        for s,c in zip(sList, cList):
            s.radius = 1.0
            map = tvtk.PolyDataMapper(input=s.output)
            act = tvtk.Actor(mapper=map, user_transform=self.transform)
            act.property.color=c
            actors.append(act)
        
        line = tvtk.LineSource(point1=(-100,0,0),
                               point2=(100,0,0))
        t_line = tvtk.TransformFilter(input=line.output,
                                      transform=self.ellipse_trans.linear_inverse)
        map = tvtk.PolyDataMapper(input=t_line.output)
        act = tvtk.Actor(mapper=map, user_transform=self.transform)
        act.property.color=(0,0,0)
        actors.append(act)
        
        l1 = tvtk.VectorText(text="F1")
        l2 = tvtk.VectorText(text="F2")
        m1 = tvtk.PolyDataMapper(input=l1.output)
        m2 = tvtk.PolyDataMapper(input=l2.output)
        act1 = self.f1_act
        act2 = self.f2_act
        
        act1.mapper = m1
        act2.mapper = m2
        
        scale = (5,5,5)
        act1.scale = scale
        act2.scale = scale

        act1.property.color=(0,0,0)
        act2.property.color=(0,0,0)
        
        act1.position = self.focus1
        act2.position = self.focus2
        
        def on_editor(new_ed):
            if new_ed is not None:
                act1.camera = new_ed._camera
                act2.camera = new_ed._camera
        
        scene.on_trait_change(on_editor, "scene_editor")
        
        actors.append(act1)
        actors.append(act2)
        return actors
    
    def trace_rays(self, rays):
        """traces a RayCollection. Reimplemented from base class
        
        returns - a recarray of intersetions with the same size as rays
                  dtype=(('length','d'),('cell','i'),('point','d',[3,]))
                  
                  'length' is the physical length from the start of the ray
                  to the intersection
                  
                  'cell' is the ID of the intersecting face
                  
                  'point' is the position coord of the intersection
        """
        p1 = rays.origin
        p2 = p1 + rays.max_length*rays.direction
        
#        t = self.transform
#        inv_t = t.linear_inverse
        
#        P1 = transformPoints(inv_t, p1)
#        P2 = transformPoints(inv_t, p2)
#        
#        #now transform into ellipse frame...
        et = self.ellipse_trans
        inv_et = et.linear_inverse
#        PP1 = transformPoints(et, P1)
#        PP2 = transformPoints(et, P2)
        
        t = self.combined_trans
        inv_t = t.linear_inverse
        PP1 = transformPoints(inv_t, p1)
        PP2 = transformPoints(inv_t, p2)
        
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
        
        xmin, xmax = self.X_bounds
        ymin, ymax = self.Y_bounds
        zmin, zmax = self.Z_bounds
        
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
        nearest = points[nearest_id, idx, :]
        h_nearest = roots[nearest_id, idx]
        
        new_pts = transformPoints(t, nearest)
        
        dtype=([('length','f8'),('cell','i2'),('point','f8',3)])
        result = numpy.empty(p1.shape[0], dtype=dtype)
        result['length'] = h_nearest*rays.max_length
        result['point'] = new_pts
        result['cell'] = 0
        return result
    
    def compute_normal(self, points, cells):
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
    
