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
from traits.api import HasTraits, Array, Float, Complex,\
            Property, List, Instance, on_trait_change, Range, Any,\
            Tuple, Event, cached_property, Set, Int, Trait
from traitsui.api import View, Item, ListEditor, VSplit,\
            RangeEditor, ScrubberEditor, HSplit, VGroup
from tvtk.api import tvtk
import numpy
import math
from itertools import izip

from raytrace.bases import Traceable, normaliseVector, NumEditor,\
                transformNormals, transformPoints
from raytrace.mirrors import BaseMirror
from raytrace.faces import OffAxisParabolicFace, PECFace


class OAPMirrorFace(OffAxisParabolicFace, PECFace):
    pass


class OffAxisParabloid(BaseMirror):
    """
    An OAP mirror object
    """
    name = "OAP"
    abstract = False
    height = Float(25.4, desc="distance between cylinder base and focus")
    EFL = Float(50.8, desc="effective focal length")
    diameter = Float(50.8, desc="outside diameter")
    
    max_length = Float(100.0)
    
    vtk_grid = Instance(tvtk.ProgrammableSource, ())
    vtk_cylinder = Instance(tvtk.Cylinder, ())
    vtk_quadric = Instance(tvtk.Quadric, ())
    
    traits_view = View(VGroup(
                        Traceable.uigroup,
                       Item('height', editor=NumEditor),
                       Item('diameter', editor=NumEditor),
                       Item('EFL', editor=NumEditor),
                       Item('max_length', editor=NumEditor),
                       ),
                       )
    
    def _vtkproperty_default(self):
        return tvtk.Property(opacity = 0.7,
                             color = (0.8,0.8,0))
    
    def _faces_default(self):
        return [OAPMirrorFace(owner=self)]
    
    def make_step_shape(self):
        from raytrace.step_export import make_OAP
        return make_OAP(self.EFL, self.diameter, self.height,
                        self.centre, self.direction, self.x_axis), "yellow"
    
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
