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


from raytrace.bases import Traceable, normaliseVector, NumEditor,\
     ComplexEditor, VectorEditor
     
from raytrace.utils import transformPoints, dotprod
from raytrace.sources import RayCollection
from faces import RectangularFace, ApertureFace

import math, numpy

class BaseAperture(Traceable):
    def _vtkproperty_default(self):
        return tvtk.Property(opacity=.35,
                             color=(0.8,0.8,0.8),
                                representation="surface")
                                
class RectApertureFace(RectangularFace, ApertureFace):
    pass



class RectAperture(BaseAperture):
    name = "Rectangular aperture"
    length = Float(100.0)
    width = Float(-25.5)
    
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
        return [RectApertureFace(owner=self)]


                                         
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