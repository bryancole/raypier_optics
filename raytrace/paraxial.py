from enthought.traits.api import Float, Instance, on_trait_change

from enthought.traits.ui.api import View, Item, VGroup

from enthought.tvtk.api import tvtk

from raytrace.bases import Traceable, NumEditor
from raytrace.faces import ParaxialLensFace


class ParaxialLens(Traceable):
    name = "Paraxial Lens"
    diameter = Float(25.4)
    focal_length = Float(25.4)
    offset = Float(0.0)
    
    vtk_disk = Instance(tvtk.DiskSource, (),
                        dict(circumferential_resolution=32), 
                        transient=True)
    
    traits_view = View(VGroup(
                       Traceable.uigroup,
                       Item('diameter', editor=NumEditor),
                       Item('focal_length', editor=NumEditor),
                        ),
                   )
    
    def make_step_shape(self):
        from raytrace.step_export import make_cylinder
        cyl = make_cylinder(self.centre, 
                             self.direction, 
                             self.diameter/2, 
                             0.1,
                             self.offset,
                             self.x_axis)
        return cyl, "blue2"
    
    def _faces_default(self):
        return [ParaxialLensFace(owner=self)]
    
    def _vtkproperty_default(self):
        return tvtk.Property(opacity = 0.7,
                             color = (0.8,0.8,1.0))
        
    def _pipeline_default(self):
        transf = tvtk.TransformFilter(input=self.vtk_disk.output, 
                                      transform=self.transform)
        return transf
    
    @on_trait_change("diameter, focal_length")
    def config_pipeline(self):
        disk = self.vtk_disk
        disk.inner_radius = 0.0
        disk.outer_radius = self.diameter/2.        
        self.update = True