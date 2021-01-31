

from traits.api import Float, on_trait_change, Range, Instance

from traitsui.api import View, VGroup, Item

from raypier.bases import Optic, Traceable, NumEditor
from raypier.core.cfaces import CircularFace
from raypier.core.ctracer import FaceList
from raypier.core.cmaterials import WaveplateMaterial
            
import numpy
from tvtk.api import tvtk


class Waveplate(Optic):
    abstract=False
    name = "Waveplate"
    retardance = Range(0.0, 1.0, value=0.25)
    diameter = Float(10.0)
    offset = Float(0.0)
    
    traits_view = View(VGroup(
                       Traceable.uigroup,
                       Item('retardance', editor=NumEditor),
                       Item('diameter', editor=NumEditor),
                       Item('offset', editor=NumEditor),
                        ),
                   )
    
    vtk_disk = Instance(tvtk.DiskSource, (),
                            dict(circumferential_resolution=64,
                                 inner_radius=0), transient=True)
    cyl_trans = Instance(tvtk.Transform, (), transient=True)
    
    line = Instance(tvtk.LineSource, (), transient=True)
    
    def _retardance_changed(self, val):
        self.material.retardance = val
    
    def _faces_default(self):
        m = self.material
        fl = FaceList(owner=self)
        fl.faces = [CircularFace(owner=self, z_plane=0, material=m)]
        return fl
    
    def _material_default(self):
        m = WaveplateMaterial(retardance=self.retardance,
                                fast_axis = self.x_axis)
        return m
    
    @on_trait_change("_orientation, rotation")
    def update_material_fast_axis(self):
        self.material.fast_axis = self.x_axis
        
    @on_trait_change("diameter, offset")
    def config_pipeline(self):
        rad = self.diameter/2.
        cyl = self.vtk_disk
        cyl.outer_radius = rad
        
        t = self.cyl_trans
        t.identity()
        t.translate(self.offset,0,0)
        #t.rotate_x(90.)
        
        line = self.line
        line.point1 = [-rad, 0, 0]
        line.point2 = [rad, 0, 0]
        
        self.update = True
        
    def _pipeline_default(self):
        cyl = self.vtk_disk
        norms = tvtk.PolyDataNormals(input_connection=cyl.output_port)
        
        tube = tvtk.TubeFilter(input_connection=self.line.output_port)
        
        append = tvtk.AppendPolyData()
        append.add_input_connection(norms.output_port)
        append.add_input_connection(tube.output_port)
        
        transF1 = tvtk.TransformFilter(input_connection=append.output_port, 
                                       transform=self.cyl_trans)
        
        transF2 = tvtk.TransformFilter(input_connection=transF1.output_port, 
                                       transform=self.transform)
        self.config_pipeline()
        return transF2