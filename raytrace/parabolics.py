#    Copyright © 2009, Teraview Ltd., Bryan Cole
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
            Property, List, Instance, on_trait_change, Range, Any,\
            Tuple, Event, cached_property, Set, Int, Trait
from enthought.traits.ui.api import View, Item, ListEditor, VSplit,\
            RangeEditor, ScrubberEditor, HSplit, VGroup
from enthought.tvtk.api import tvtk
import numpy
import math
from itertools import izip

from raytrace.tracer import Traceable, normaliseVector, RaySegment, NumEditor,\
                transformNormals, transformPoints
from raytrace.mirrors import BaseMirror

class OffAxisParabloid(BaseMirror):
    """
    An OAP mirror object
    """
    name = "OAP"
    height = Float(25.4)
    EFL = Float(50.8)
    diameter = Float(50.8)
    
    max_length = Float(100.0)
    
    vtk_grid = Instance(tvtk.ProgrammableSource, ())
    vtk_cylinder = Instance(tvtk.Cylinder, ())
    vtk_quadric = Instance(tvtk.Quadric, ())
    
    vtkproperty = tvtk.Property(opacity = 1.0,
                             color = (0.8,0.8,0))
    
    traits_view = View(VGroup(
                        Traceable.uigroup,
                       Item('height', editor=NumEditor),
                       Item('diameter', editor=NumEditor),
                       Item('EFL', editor=NumEditor),
                       Item('max_length', editor=NumEditor),
                       ),
                       )
    
    def create_grid(self):
        EFL = self.EFL
        r = self.diameter/2
        h = self.height
        l = self.max_length
        source = self.vtk_grid
        sp = source.structured_points_output
        size = 20
        spacing = 2*r / (size-1)
        lsize = int(l/spacing)
        sp.dimensions = (size,size,lsize)
        sp.whole_extent=(0,size,0,size,0,lsize)
        sp.origin = (EFL - r, -r, -h)
        sp.spacing = (spacing, spacing, spacing)
        sp.set_update_extent_to_whole_extent()
    
    @on_trait_change("EFL, diameter, height, max_length")
    def config_pipeline(self):
        print "config"
        EFL = self.EFL
        rad = self.diameter/2
        h = self.height
        
        self.vtk_grid.modified()
                      
        cyl = self.vtk_cylinder
        cyl.center = (EFL,0,0)
        cyl.radius = rad
        
        q = self.vtk_quadric
        A0 = A1 = -1/(2*EFL)
        A8 = 1 
        A9 = EFL/2
        q.coefficients = (A0, A1, 0, 
                          0, 0, 0, 
                          0, 0, A8, 
                          A9)
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
        
        clip2 = tvtk.ClipDataSet(input = clip1.output,
                                 clip_function=self.vtk_quadric,
                                 inside_out=1)
        
        topoly = tvtk.GeometryFilter(input=clip2.output)
        norms = tvtk.PolyDataNormals(input=topoly.output)
        
        transF = tvtk.TransformFilter(input=norms.output, transform=self.transform)
        self.config_pipeline()
        #clip1.update()
        grid.modified()
        return transF
    
    def intersect_with_line(self, p1, p2):
        """
        
        @param p1: a (n,3) array of points, start of each ray
        @param p2: a (n,3) array of point, ends of the rays
        """
        trans = self.transform
        inv_t = trans.linear_inverse
        P1 = transformPoints(inv_t, p1)
        P2 = transformPoints(inv_t, p2)
        efl = self.EFL #scalar
        A = 1 / (2*efl)
        s = P2 - P1
        r = P1
        r[:,2] += self.EFL/2.
        
        sx, sy, sz = s.T
        rx, ry, rz = r.T
        
        a = A*(sx**2 + sy**2)
        b = 2*A*(rx*sx + ry*sy) - sz
        c = A*(rx**2 + ry**2) - rz
        
        d = b**2 - 4*a*c
        ###FIXME
        if d<0: #no intersection
            return None
        if a < 1e-10: #approximate to zero if we're close to the parabolic axis
            roots = [-c/b]
        else:
            e = math.sqrt(d)
            roots = [(-b+e)/(2*a), (-b-e)/(2*a)]

        pos_roots = [root for root in roots if 1.0 >= root > 0.01]
        if pos_roots:
            alpha = min(pos_roots)
        else:
            return None #no intersections in range
        
        P = [g+alpha*h for g,h in izip(P1,s)]
        
        rad = self.diameter/2.
        if (P[0]-efl)**2 + P[1]**2 > rad**2:
            return None
        
        return self.transform.transform_point(P)
        
    
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
    
    def compute_normal(self, points, cell_ids):
        """
        evaluate normalised Normal vector
        """
        n = points.shape[0]
        t = self.transform
        inv_t = t.linear_inverse
        t_points =transformPoints(inv_t, points)
        coefs = self.vtk_quadric.coefficients
        ax = 2*coefs[0]/coefs[8]
        ay = 2*coefs[1]/coefs[8]
        t_normal = numpy.column_stack((ax*t_points[:,0], 
                                       ay*t_points[:,1],
                                       -numpy.ones(n)))
        return transformNormals(t, t_normal)
    
#    def eval_children(self, seg, point, cell_id):
#        """
#        actually calculates the new ray-segments. Physics here.
#        Easy for a mirror.
#        """
#        #return []
#        normal = self.compute_normal(point)
#        
#        proj = numpy.dot(normal, seg.direction)
#        reflected = seg.direction - 2*proj*normal
#        
#        origin = numpy.asarray(point)
#            
#        refl_ray = RaySegment(origin=origin,
#                           direction = reflected,
#                           parent = seg,
#                           refractive_index=seg.refractive_index)
#        return [refl_ray,]
        
    
if __name__=="__main__":
    oap = OffAxisParabloid()
    
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
