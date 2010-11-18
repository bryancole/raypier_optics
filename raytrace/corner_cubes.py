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
            Tuple, Event, cached_property, Set, Int, Trait, Bool
from enthought.traits.ui.api import View, Item, ListEditor, VSplit,\
            RangeEditor, ScrubberEditor, HSplit, VGroup, Heading
from enthought.tvtk.api import tvtk
from enthought.tvtk.pyface.scene_model import SceneModel
from enthought.tvtk.pyface.scene_editor import SceneEditor
import numpy
from itertools import chain, izip, islice, tee

#from raytrace.tracer import Optic, VTKOptic, normaliseVector, RaySegment,\
#             Traceable, NumEditor, dotprod, transformPoints, transformNormals

from raytrace.bases import Optic, Traceable, NumEditor
from raytrace.cfaces import ElipticalPlaneFace
from raytrace.ctracer import FaceList
from raytrace.mirrors import BaseMirror


class HollowRetroreflector(BaseMirror):
    name = "Hollow Retroreflector"
    diameter = Float(25.4)
    #thickness = Float(5.0, desc="purely for visualisation purposes")
    #offset = Float(0.0)
    
    imp_cyl = Instance(tvtk.Cylinder, (), transient=True)
    
    #cyl_trans = Instance(tvtk.Transform, (), transient=True)
    
    scale = Instance(tvtk.Transform, (), transient=True)
    
    traits_view = View(VGroup(
                       Traceable.uigroup,
                       Item('diameter', editor=NumEditor),
                        ),
                   )
    
    def _faces_default(self):
        fl = FaceList(owner=self)
        sqrt2 = numpy.sqrt(2)
        x = sqrt2*numpy.cos(numpy.pi*2./3)
        y = sqrt2*numpy.sin(numpy.pi*2./3)
        fl.faces = [ElipticalPlaneFace(owner=self, g_x=sqrt2, g_y=0),
                    ElipticalPlaneFace(owner=self, g_x=x, g_y=y),
                    ElipticalPlaneFace(owner=self, g_x=x, g_y=-y),
                    ]
        return fl
    
#    def make_step_shape(self):
#        from raytrace.step_export import make_cylinder
#        cyl = make_cylinder(self.centre, 
#                             self.direction, 
#                             self.diameter/2, 
#                             self.thickness,
#                             self.offset,
#                             self.x_axis)
#        return cyl, "green"
    
    @on_trait_change("diameter")
    def config_pipeline(self):
        r = self.diameter/2.
        cyl = self.imp_cyl
        cyl.radius = r
        #thick = self.thickness
        #cyl.height = thick
        
        #t = self.cyl_trans
        #t.identity()
        ##t.translate(self.offset,0,thick/2.)
        #t.rotate_x(90.)
        s=self.scale
        s.identity()
        s.scale(r,r,r)
        
        self.update = True
    
    def _pipeline_default(self):
        sqrt2 = numpy.sqrt(2)
        x = -sqrt2*numpy.cos(numpy.pi*2./3)
        y = sqrt2*numpy.sin(numpy.pi*2./3)
        p1 = tvtk.PlaneSource(origin=(0,0,0), \
                              point1=(-sqrt2,0,1), point2=(x,y,1))
        p2 = tvtk.PlaneSource(origin=(0,0,0), \
                               point1=(x,y,1), point2=(x,-y,1))
        p3 = tvtk.PlaneSource(origin=(0,0,0), \
                               point1=(x,-y,1), point2=(-sqrt2,0,1))
        for p in (p1,p2,p3):
            p.set_resolution(32,32)
        append = tvtk.AppendPolyData()
        append.add_input(p1.output)
        append.add_input(p2.output)
        append.add_input(p3.output)
        
        scale_f = tvtk.TransformFilter(input=append.output, transform=self.scale)
        
        trans = tvtk.Transform()
        trans.rotate_x(90.0)
        self.imp_cyl.transform = trans
        
        clip = tvtk.ClipPolyData(input=scale_f.output, clip_function = self.imp_cyl,
                                 inside_out=True)
                               
        #cyl = self.vtk_cylinder
        #norms = tvtk.PolyDataNormals(input=cyl.output)
        #transF1 = tvtk.TransformFilter(input=norms.output, transform=self.cyl_trans)
        transF2 = tvtk.TransformFilter(input=clip.output, transform=self.transform)
        self.config_pipeline()
        return transF2