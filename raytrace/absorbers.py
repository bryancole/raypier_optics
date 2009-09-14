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
     ComplexEditor, VectorEditor,Face
     
from raytrace.utils import Convert_to_SP, dotprod, transformPoints,\
        transformNormals
from raytrace.sources import RayCollection
from faces import CircularFace, RectangularFace

import math, numpy

class AbsorberFace(Face):
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
        
        print "called!"
        
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
                                   max_length = numpy.inf,
                                   E_vector = S_vec,
                                   E1_amp = -S_amp,
                                   E2_amp = -P_amp,
                                   parent = rays,
                                   parent_ids = parent_ids,
                                   face = faces,
                                   refractive_index=n,
                                   normals = normal)
        return 0


class CircularAbsorberFace(CircularFace, AbsorberFace):
    pass


class BaseAbsorber(Traceable):
    def _vtkproperty_default(self):
        return tvtk.Property(opacity=1.0,
                             color=(0.,0.,0.),
                                representation="surface")


class AbsorberDisk(BaseAbsorber):
    name = "PEC Mirror"
    diameter = Float(25.4)
    thickness = Float(5.0, desc="purely for visualisation purposes")
    offset = Float(0.0)
    
    vtk_cylinder = Instance(tvtk.CylinderSource, (),
                            dict(resolution=32), transient=True)
    
    cyl_trans = Instance(tvtk.Transform, (), transient=True)
    
    traits_view = View(VGroup(
                       Traceable.uigroup,
                       Item('diameter', editor=NumEditor),
                       Item('thickness', editor=NumEditor),
                       Item('offset', editor=NumEditor)
                        ),
                   )
    
    def _faces_default(self):
        return [CircularAbsorberFace(owner=self)]
    
    def make_step_shape(self):
        from raytrace.step_export import make_cylinder
        cyl = make_cylinder(self.centre, 
                             self.direction, 
                             self.diameter/2, 
                             self.thickness,
                             self.offset,
                             self.x_axis)
        return cyl, "green"
    
    @on_trait_change("diameter, thickness, offset")
    def config_pipeline(self):
        cyl = self.vtk_cylinder
        cyl.radius = self.diameter/2.
        thick = self.thickness
        cyl.height = thick
        
        t = self.cyl_trans
        t.identity()
        t.translate(self.offset,0,thick/2.)
        t.rotate_x(90.)
        
        self.update = True
    
    def _pipeline_default(self):
        cyl = self.vtk_cylinder
        norms = tvtk.PolyDataNormals(input=cyl.output)
        transF1 = tvtk.TransformFilter(input=norms.output, transform=self.cyl_trans)
        transF2 = tvtk.TransformFilter(input=transF1.output, transform=self.transform)
        self.config_pipeline()
        return transF2
    
    

class RectAbsorberFace(RectangularFace, AbsorberFace):
    pass


class RectAbsorber(BaseAbsorber):
    """
    A rectangular Absrober object
    """
    name = "rectangular absorber"
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
                             color = (0.,0.,0.))
    
    def _faces_default(self):
        return [RectAbsorberFace(owner=self)]

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
    from tracer import RayTraceModel, BuildRaySet, BuildConfocalRaySet
    
    mirror = AbsorberDisk(centre=(0,0,0),
                       orientation=0.0,
                       elevation=45.0)
    
#    input_rays = BuildRaySet(origin = (0,0,-20),
#                         direction = (0,0,1),
#                         radius=5.0,
#                         count=20)
    input_rays = BuildConfocalRaySet(focus=(0,0,5), 
                                     direction=(0,0,1),
                                     theta=20, 
                                     working_dist=20, 
                                     count=20)
    
    model = RayTraceModel(optics=[mirror], rays=input_rays)
    model.configure_traits()
