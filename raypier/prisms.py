#    Copyright 2009, Teraview Ltd., Bryan Cole
#
#    This file is part of Raypier.
#
#    Raypier is free software: you can redistribute it and/or modify
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


from traits.api import HasTraits, Array, Float, Complex,\
            Property, List, Instance, on_trait_change, Range, Any,\
            Tuple, Event, cached_property, Set, Int, Trait, Bool, observe
from traitsui.api import View, Item, ListEditor, VSplit,\
            RangeEditor, ScrubberEditor, HSplit, VGroup, Heading, InstanceEditor
from tvtk.api import tvtk
from tvtk.pyface.scene_model import SceneModel
from tvtk.pyface.scene_editor import SceneEditor
import numpy
from itertools import chain, islice, tee

#from raypier.tracer import Optic, VTKOptic, normaliseVector, RaySegment,\
#             Traceable, NumEditor, dotprod, transformPoints, transformNormals

from raypier.bases import Optic, Traceable
from raypier.core.cfaces import PolygonFace, ExtrudedPlanarFace, OrientedPolygonFace
from raypier.core.ctracer import FaceList
from raypier.materials import BaseOpticalMaterial, OpticalMaterial
from raypier.core.cmaterials import FullDielectricDispersiveMaterial


def pairwise(itr):
    a,b = tee(itr)
    fst = next(b)
    return zip(a, chain(b,[fst]))

class Extrusion(Optic):
    """a general flat-faced optic formed by extrusion 
    of a 2D polygon"""
    abstract=True
    profile = Array(shape=(None,2), dtype=numpy.double, transient=True)
    
    z_height_1 = Float(0.0)
    z_height_2 = Float(20.0)
    
    trace_ends = Bool(True, desc="include the end-faces in tracing")
    
    data_source = Instance(tvtk.ProgrammableSource, (), transient=True)
    
    extrude = Instance(tvtk.LinearExtrusionFilter, (), 
                       {'capping': True, 
                        'extrusion_type':'vector',
                        'vector': (0.,0.,1.)},
                        transient=True)
    
    def _faces_default(self):
        fl = FaceList(owner=self)
        fl.faces = self.make_faces()
        return fl
    
    def make_step_shape(self):
        from raypier.step_export import make_extruded_profile
        
        shape = make_extruded_profile(self.centre, 
                                      self.direction, 
                                      self.x_axis, 
                                      self.z_height_1, 
                                      self.z_height_2, 
                                      list(self.profile))
        return shape, "blue"
        
    
    def make_faces(self):
        z1 = self.z_height_1
        z2 = self.z_height_2
        profile = self.profile
        m = self.material
        #a plane is extruded differently from a 2d profile
        if profile.shape == (2,2):
            sides = [ExtrudedPlanarFace(owner=self, z1=z1, z2=z2, x1=profile[0,0], y1=profile[0,1], 
                        x2=profile[1,0], y2=profile[1,1], material=m)]
        else:
            sides = [ExtrudedPlanarFace(owner=self, z1=z1, z2=z2, x1=x1, y1=y1, 
                        x2=x2, y2=y2, material=m) for ((x2,y2),(x1,y1)) 
                                            in pairwise(profile)]
        if self.trace_ends:
            base = PolygonFace(owner=self, z_plane=z1,
                        xy_points=profile, material=m)
            top = PolygonFace(owner=self, z_plane=z2, material=m,
                        xy_points=profile, invert_normal=True)
            sides.extend([base, top])
        return sides
    
    @on_trait_change("z_height_1, z_height_2")
    def config_pipeline(self):
        self.extrude.scale_factor = self.z_height_2 - self.z_height_1
        self.faces.faces = self.make_faces()
        self.update=True
        
    def _trace_ends_changed(self):
        self.faces.faces = self.make_faces()
        self.update=True
        
    def _profile_changed(self):
        self.data_source.modified()
        self.faces.faces = self.make_faces()
        self.update=True
        
    def _vtk_profile(self):
        return self.profile
    
    def _pipeline_default(self):
        self.config_profile()
        self.config_pipeline()
        source = self.data_source
        def execute():
            xy = self._vtk_profile()
            z = numpy.ones(xy.shape[0]) * self.z_height_1
            points = numpy.column_stack((xy,z))
            
            cells = [list(range(len(z))),]
            
            output = source.poly_data_output
            output.points = points
            output.polys = cells
        source.set_execute_method(execute)
        print("Made PIPELINE")
        
        extrude = self.extrude
        extrude.input_connection = source.output_port
        
        t = self.transform
        transf = tvtk.TransformFilter(input_connection=extrude.output_port, 
                                      transform=t)
        return transf
    
    def config_profile(self):
        """abstract method to set the profile data"""
        pass


class Prism(Extrusion):
    name = "Prism"
    abstract = False
    height = Float(20.0) #distance from front (hypotenuse) face to apex
    width = Float(20.0) #width of front face
    
    traits_view = View(VGroup(
                       Traceable.uigroup,
                       Item('trace_ends'),
                       Item('n_inside'),
                       Item('z_height_1', label="Z top"),
                       Item('z_height_2', label="Z bottom"),
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
        

class RightAnglePrism(Optic):
    """
    A more convenient righ-angle prism optic which aligns the centre of the hypotenuse face as the origin and
    the normal of the hypotenuse is placed along the z-axis. The hypotenuse is aligned along the local x-axis.
    """
    abstract = False
    name = "Right-angle prism"
    ### The short side of the prism
    side_length = Float(10.0)
    ### height along the "extrusion axis"
    height = Float(10.0)
    
    dispersion = Instance(BaseOpticalMaterial)
    
    trace_ends = Bool(True)
    
    ###Visualisation objects
    data_source = Instance(tvtk.ProgrammableSource, (), transient=True)
    extrude = Instance(tvtk.LinearExtrusionFilter, (), 
                       {'capping': True, 
                        'extrusion_type':'vector',
                        'vector': (0.,1.,0.)},
                        transient=True)
    
    traits_view = View(VGroup(
                       Traceable.uigroup,
                       Item('trace_ends'),
                       Item('side_length'),
                       Item('height'),
                       Item('dispersion', style="custom", editor=InstanceEditor())
                       )
                       )
    
    def make_step_shape(self):
        from raypier.step_export import make_general_extrusion
        s_len = self.side_length
        h_len = numpy.sqrt(2*(s_len**2))/2
        y = self.height / 2.
        
        points = [(-h_len,-y,0.),(0.,-y,-h_len),(h_len,-y,0.)]
        length = self.height
        vector = (0.,1.,0.)
        
        shape = make_general_extrusion(self.centre, 
                                      self.direction, 
                                      self.x_axis, 
                                      points, length, vector)
        return shape, 'purple3'
    
    @on_trait_change("n_inside, n_outside")
    def n_changed(self):
        pass
    
    @observe("side_length, height, trace_ends")
    def on_param_changed(self, evt):
        self.faces.faces = self.make_faces()
    
    def _dispersion_default(self):
        return OpticalMaterial(glass_name="N-BK7")
    
    def _material_default(self):
        dispersion_curve = self.dispersion.dispersion_curve
        m = FullDielectricDispersiveMaterial(dispersion_inside=dispersion_curve)
        return m
    
    def _faces_default(self):
        fl = FaceList(owner=self)
        fl.faces = self.make_faces()
        return fl
    
    def make_faces(self):
        s_len = self.side_length
        h_len = numpy.sqrt(2*(s_len**2))/2
        y = self.height/2
        hypot = OrientedPolygonFace(normal=(0,0,1), 
                                    x_axis=(1,0,0), 
                                    origin=(0,0,0),
                                    xy_points=[(-h_len,-y),(-h_len,y),(h_len,y),(h_len,-y)])
        
        side1 = OrientedPolygonFace(normal=(1,0,-1), 
                                    x_axis=(-1,0,-1), 
                                    origin=(0,0,-h_len),
                                    xy_points=[(-s_len,-y),(-s_len,y),(0,y),(0,-y)])
        
        side2 = OrientedPolygonFace(normal=(-1,0,-1), 
                                    x_axis=(-1,0,1), 
                                    origin=(0,0,-h_len),
                                    xy_points=[(s_len,-y),(s_len,y),(0,y),(0,-y)])
        faces = [hypot, side1, side2]
        if self.trace_ends:
            pts = [(h_len,0),(0,h_len), (-h_len,0)]
            top = OrientedPolygonFace(normal=(0,1,0), 
                                    x_axis=(1,0,0), 
                                    origin=(0,y,0),
                                    xy_points=pts)
            bottom = OrientedPolygonFace(normal=(0,-1,0), 
                                    x_axis=(-1,0,0), 
                                    origin=(0,-y,0),
                                    xy_points=pts)
            faces.extend([top, bottom])
            
        for f in faces:
            f.material=self.material
            
        return faces
    
    @on_trait_change("side_length, height")
    def config_pipeline(self):
        self.extrude.scale_factor = self.height
        self.faces.faces = self.make_faces()
        self.update=True
        
    def _pipeline_default(self):
        self.config_pipeline()
        source = self.data_source
        def execute():
            s_len = self.side_length
            h_len = numpy.sqrt(2*(s_len**2))/2
            y = self.height / 2.
            
            points = numpy.array([(-h_len,-y,0.),(0.,-y,-h_len),(h_len,-y,0.)])            
            cells = [list(range(3)),]
            
            output = source.poly_data_output
            output.points = points
            output.polys = cells
        source.set_execute_method(execute)
        print("Made PIPELINE")
        
        extrude = self.extrude
        extrude.scale_factor = self.height
        extrude.input_connection = source.output_port
        
        t = self.transform
        transf = tvtk.TransformFilter(input_connection=extrude.output_port, 
                                      transform=t)
        return transf
        
        
class Rhomboid(Extrusion):
    name = "rhomboid"
    abstract = False
    height = Float(15.0) #distance between parallel faces
    width = Float(15.0) #width of parallel faces
    slant = Float(45.0) #angle of the oblique faces
    
    traits_view = View(VGroup(
                       Traceable.uigroup,
                       Item('trace_ends'),
                       Item('n_inside'),
                       Item('z_height_1'),
                       Item('z_height_2'),
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
        
        
class TruncatedRightanglePrism(Extrusion):
    name = "truncated right-angle prism"
    abstract = False
    depth = Float #distance from front face to apex
    width = Float #width of front face
    
    traits_view = View(VGroup(
                       Traceable.uigroup,
                       Item('trace_ends'),
                       Item('n_inside'),
                       Item('depth'),
                       Item('width'),
                       Item('z_height_1', label="Z top"),
                       Item('z_height_2', label="Z bottom"),
                       )
                       )

    @on_trait_change("depth, width")
    def config_profile(self):
        h = self.depth
        w = self.width/2.
        self.profile = [(w,0),
                        (w, h-w),
                        (0,h),
                        (-w, h-w),
                        (-w,0)]
        

class LDLF(Extrusion):
    name = "Linear Dialectric Light Funnle"
    abstract = False
    slat_width= Float #width of slats
    ap_width = Float #width of exit apperture
    slant = Float #angle of slanted sides
    
    traits_view = View(VGroup(
                       Traceable.uigroup,
                       Item('trace_ends'),
                       Item('n_inside'),
                       Item('length'),
                       Item('slat_width'),
                       Item('ap_width'),
                       Item('slant')
                       )
                       )

    @on_trait_change("slat_width, ap_width, slant")
    def config_profile(self):
        theta = self.slant*numpy.pi/180.
        l = self.slat_width
        h = l*numpy.sin(theta)
        dis = l*numpy.cos(theta)
        w = self.ap_width/2
        self.profile = [(-w,0),
                        (w,0),
                        (w+dis,h),
                        (-w-dis,h)]
                        
class Sheet(Extrusion):
    #just a wrapper for the extrudedplanarface
    name = "Extruded Sheet"
    abstract = False
    x1 = Float 
    y1 = Float
    x2 = Float
    y2 = Float
    trace_ends = Bool(False)
    
    traits_view = View(VGroup(
                       Traceable.uigroup,
                       Item('trace_ends'),
                       Item('n_inside'),
                       Item('x1'),
                       Item('x2'),
                       Item('y1'),
                       Item('y2')
                       )
                       )

    @on_trait_change("x1,x2,y1,y2")
    def config_profile(self):
        
        self.profile = [(self.x1,self.y1),
                        (self.x2,self.y2)]

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
    print("beamstop", beamstop.polydata)
    
    model = RayTraceModel(optics=[rhomboid, beamstop], rays=input_rays)
    model.configure_traits()
