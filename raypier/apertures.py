
from .bases import Traceable, NumEditor
from .core.cmaterials import CircularApertureMaterial, RectangularApertureMaterial
from .core.cfaces import CircularFace, RectangularFace
from .core.ctracer import FaceList

from traits.api import Float, on_trait_change, Instance, Property, Str, Bool, Constant
from traitsui.api import View, Item, VGroup
from tvtk.api import tvtk

import numpy
from scipy.special import erf


class BaseAperture(Traceable):
    edge_width = Float(2.0)
    invert = Bool(False)
    
    sc = Instance(tvtk.ProgrammableAttributeDataFilter, (),
                  transient=True)
    
    def _actors_default(self):
        pipeline = self.pipeline
        
        lut = tvtk.LookupTable()
        lut.range=(0.0,1.0)
        lut.ramp="linear"
        lut.alpha_range=(1.0,0.0)
        lut.hue_range=(0.65,0.65)
        lut.value_range=(0.4,0.4)
        lut.build()
        lut.below_range_color=(0.2,0.0,0.2,1.0)
        lut.above_range_color=(0.2,0.0,0.2,1.0)
        lut.use_below_range_color=True
        lut.use_above_range_color=True
        
        map = tvtk.PolyDataMapper(input_connection=pipeline.output_port,
                                  lookup_table=lut)
        map.scalar_range = (0.0,1.0)
        #map.scalar_visibility = True
        act = tvtk.Actor(mapper=map)
        #act.property = self.vtkproperty
        actors = tvtk.ActorCollection()
        actors.append(act)
        return actors


class CircularAperture(BaseAperture):
    abstract=False
    """
    This circular aperture is characterised by 3 diameters::
    
      outer_diameter : rays passing outside this diameter do not intersect
      inner_diameter : rays passing within the outer diameter but outside of the inner
                       diameter are blocked (creating no child rays).
                       Rays passing inside the inner_diameter intersect and generate
                       a child ray with amplitudes determined by the position 
                       relative to the aperture origin
      hole_diameter  : Rays passing inside the inner_diameter are attenuated if they 
                       pass outside the hole diameter and left unattenuated
                       if they are inside the hole. The boundary region is generates
                       a sigmoidal profile of width given by the edge_width parameter.
    """
    outer_diameter = Float(25.0)
    inner_diameter = Float(15.0)
    hole_diameter = Float(10.0)
    
    diameter = Property()
    offset = Float(0.0)
    
    name = Str("Circular aperture")
    
    vtkproperty = tvtk.Property(opacity = 0.7,
                             color = (0.0,0.2,0.8))
    
    vtk_disk = Instance(tvtk.DiskSource, (),
                            dict(circumferential_resolution=64,
                                 inner_radius=0), transient=True)
    vtk_disk2 = Instance(tvtk.DiskSource, (),
                            dict(circumferential_resolution=64,
                                 radial_resolution=50,
                                 inner_radius=0), transient=True)
    
    traits_view = View(VGroup(
                       Traceable.uigroup,
                       Item('outer_diameter', editor=NumEditor),
                       Item('inner_diameter', editor=NumEditor),
                       Item('hole_diameter', editor=NumEditor),
                       Item('edge_width', editor=NumEditor),
                       Item('invert'),
                        ),
                   )
    
    def _get_diameter(self):
        return self.outer_diameter
    
    def _material_default(self):
        m = self.make_material()
        return m
    
    def _faces_default(self):
        fl = FaceList(owner=self)
        fl.faces = self.make_faces()
        return fl
    
    def make_material(self):
        m = CircularApertureMaterial(origin=self.centre,
                                    radius = self.hole_diameter/2,
                                     outer_radius = self.inner_diameter/2,
                                    edge_width = self.edge_width/2,
                                    invert=self.invert)
        return m
    
    def make_faces(self):
        fl = [CircularFace(owner=self, diameter=self.outer_diameter,
                                material = self.material)]
        return fl
    
    @on_trait_change("inner_diameter, hole_diameter, edge_width, centre, invert")
    def on_material_params_changed(self):
        self.material = self.make_material()
        self.faces.faces = self.make_faces()
        self.vtk_disk.inner_radius = self.inner_diameter/2
        self.config_pipeline()
        self.update = True
        
    @on_trait_change("outer_diameter")
    def on_face_params_changed(self):
        self.faces.faces = self.make_faces()
        self.vtk_disk.outer_radius = self.outer_diameter/2
        self.config_pipeline()
        self.update = True
        
    def config_pipeline(self):
        disk = self.vtk_disk
        disk.inner_radius = self.inner_diameter/2
        disk.outer_radius = self.outer_diameter/2
        disk2 = self.vtk_disk2
        disk2.inner_radius = 0
        disk2.outer_radius = self.inner_diameter/2
        self.sc.modified()
        
    def _pipeline_default(self):
        disk1 = self.vtk_disk
        disk2 = self.vtk_disk2
        
        append = tvtk.AppendPolyData()
        append.add_input_connection(disk1.output_port)
        append.add_input_connection(disk2.output_port)
        
        sc = self.sc
        sc.input_connection = append.output_port

        def update( *args ):
            in_data = sc.get_input_data_object(0,0)
            points = in_data.points.to_array()
            centre = numpy.array(self.centre)
            r = numpy.sqrt(((points)**2).sum(axis=-1))
            width = self.edge_width
            r0 = self.hole_diameter
            alpha = 0.5 - 0.5*(erf(2*(r-r0)/width))
            alpha[r>=(0.99*self.inner_diameter/2)] = -0.01
            out = sc.get_output_data_object(0)
            if self.invert:
                alpha = 1 - alpha
            out.point_data.scalars=alpha
            
        sc.set_execute_method(update)
        
        trans = tvtk.TransformFilter(input_connection=sc.output_port, 
                                       transform=self.transform)
        self.config_pipeline()
        return trans
    
    
        
        
class RectangularAperture(BaseAperture):
    abstract=False
    name = Str("Rectangular aperture")
    
    hole_width = Float(2.0)
    hole_height = Float(5.0)
    outer_width = Float(20.0)
    outer_height = Float(15.0)
    inner_width = Float(10.0)
    inner_height = Float(7.0)
    
    offset = Constant(0.0)
    length = Property()
    width = Property()
    
    vtkplane = Instance(tvtk.PlaneSource, (), {'x_resolution':500, 'y_resolution':500})
    vtkbox = Instance(tvtk.Box, ())
    
    traits_view = View(VGroup(
                       Traceable.uigroup,
                       Item('hole_width', editor=NumEditor),
                       Item('hole_height', editor=NumEditor),
                       Item('inner_width', editor=NumEditor),
                       Item('inner_height', editor=NumEditor),
                       Item('outer_width', editor=NumEditor),
                       Item('outer_height', editor=NumEditor),
                       Item('edge_width', editor=NumEditor),
                       Item('invert'),
                        ),
                   )
    
    def _get_length(self):
        return self.outer_height
    
    def _get_width(self):
        return self.outer_width
    
    def _material_default(self):
        m = self.make_material()
        return m
    
    def _faces_default(self):
        fl = FaceList(owner=self)
        fl.faces = self.make_faces()
        return fl
    
    def make_material(self):
        m = RectangularApertureMaterial(origin=self.centre,
                                    width = self.hole_width,
                                    height = self.hole_height,
                                    outer_width = self.inner_width,
                                    outer_height = self.inner_height,
                                    edge_width = self.edge_width/2,
                                    invert=self.invert)
        return m
    
    def make_faces(self):
        fl = [RectangularFace(owner=self, length=self.outer_width,
                              width=self.outer_height, offset=0.0,
                              z_plane=0.0,
                                material = self.material)]
        return fl
    
    @on_trait_change("inner_width, inner_height, hole_width, hole_height, edge_width, centre, invert")
    def on_material_params_changed(self):
        self.material = self.make_material()
        self.faces.faces = self.make_faces()
        
        self.config_pipeline()
        self.update = True
        
    @on_trait_change("outer_diameter")
    def on_face_params_changed(self):
        self.faces.faces = self.make_faces()
        self.vtk_disk.outer_radius = self.outer_diameter/2
        self.config_pipeline()
        self.update = True
        
    def config_pipeline(self):
        plane = self.vtkplane
        w = self.outer_width/2.
        h = self.outer_height/2.
        plane.origin = (-w, -h, 0.0)
        plane.point1 = (-w, h, 0.0)
        plane.point2 = (w, -h, 0.0)
        
        box = self.vtkbox
        w = self.inner_width/2.
        h = self.inner_height/2.
        box.set_x_min(-w,-h,-1.0)
        box.set_x_max(w,h,1.0)
        
        self.sc.modified()
    
    def _pipeline_default(self):
        plane = self.vtkplane
        
        box = self.vtkbox
        
        clip = tvtk.ClipPolyData(input_connection=plane.output_port)
        clip.clip_function=box
        
        sc = self.sc
        sc.input_connection = plane.output_port

        def update( *args ):
            in_data = sc.get_input_data_object(0,0)
            points = in_data.points.to_array()
            centre = numpy.array(self.centre)
            
            width = self.edge_width
            x0 = self.hole_width/2
            y0 = self.hole_height/2
            x = points[:,0]
            y = points[:,1]
            
            alpha = 0.5 - 0.5*(erf((x-x0)/width))
            alpha *= 0.5 - 0.5*(erf(-(x+x0)/width))
            alpha *= 0.5 - 0.5*(erf((y-y0)/width))
            alpha *= 0.5 - 0.5*(erf(-(y+y0)/width))
            
            x1 = self.inner_width/2
            y1 = self.inner_height/2
            
            alpha[numpy.abs(x)>=x1] = -0.01
            alpha[numpy.abs(y)>=y1] = -0.01
            
            out = sc.get_output_data_object(0)
            if self.invert:
                alpha = 1 - alpha
            out.point_data.scalars=alpha
            
        sc.set_execute_method(update)
        
        trans = tvtk.TransformFilter(input_connection=sc.output_port, 
                                       transform=self.transform)
        self.config_pipeline()
        return trans