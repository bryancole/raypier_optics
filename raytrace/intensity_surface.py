

from traits.api import Float, Instance, Enum, Array

from traitsui.api import View, Item, VGroup, HGroup

from tvtk.api import tvtk
from tvtk.pyface.scene_model import SceneModel
from tvtk.pyface.scene_editor import SceneEditor

from raytrace.results import Result
from raytrace.fields import EFieldPlane


class ImageSource(object):
    def Initialise(self, vtkobj):
        vtkobj.SetNumberOfInputPorts(0)
        vtkobj.SetNumberOfOutputPorts(1)
        
    def FillInputPortInformation(self, vtkself, port, info):
        info.Set(tvtk.Algorithm().INPUT_REQUIRED_DATA_TYPE()._vtk_obj, "vtkDataSet")
        return 1
 
    def FillOutputPortInformation(self, vtkself, port, info):
        info.Set(tvtk.DataObject().DATA_TYPE_NAME()._vtk_obj, "vtkImageData")
        return 1
 
    def ProcessRequest(self, vtkself, request, inInfo, outInfo):
        if request.Has(tvtk.DemandDrivenPipeline().REQUEST_DATA()._vtk_obj):
            print( 'I am supposed to execute' )
        return 1


class IntensitySurface(Result):
    field_probe = Instance(EFieldPlane)
    
    display = Enum("Intensity", "E_x", "E_y", "E_z")
    
    intensity_data = Array()
    
    scale = Float(1.0)
    
    #scene = Instance(SceneModel, (), {'background':(0,0,0)}, transient=True)
    scene = Instance(SceneModel, transient=True)
    
    _data_source = Instance(tvtk.ProgrammableSource, ())
    
    traits_view = View(VGroup(
                        HGroup(
                            Item("display", show_label=False),
                            Item("scale", show_label=False)
                            ),
                        Item("scene", editor=SceneEditor(), show_label=False)
                        ))
        
    def _scene_default(self):
        scene = SceneModel(background=(0,0,0))
        
        src = self._data_source
        src.set_execute_method(self.create_grid)
        src.modified()
        
        geom = tvtk.ImageDataGeometryFilter(input_connection=src.output_port)
        
        mapper = tvtk.PolyDataMapper(input_connection=geom.output_port)
        
        act = tvtk.Actor(mapper=mapper)
        
        scene.add_actor(act)
        return scene
            
    def create_grid(self):
        source = self._data_source
        sp = source.structured_points_output
        w,h = self.intensity_data.shape
        sp.dimensions = (w,h,1)
        sp.whole_extent=(0,w,0,h,0,1)
        sp.origin = (-0.5, -0.5, 0)
        sp.spacing = (1./(w-1), 1./(h-1), 1)
        sp.set_update_extent_to_whole_extent()
        sp.point_data.scalars = self.intensity_data
        source.output = sp
    
    def calc_intensity(self, e_field):
        E = e_field
        mode = self.display
        if mode=="Intensity":
            U = (E.real**2).sum(axis=-1) + (E.imag**2).sum(axis=-1)
        else:
            idx = {"E_x":0, "E_y":1, "E_z":2}[mode]
            U = E[:,:,idx]
        self.intensity_data = U
        return U
    
    def _display_changed(self):
        self.on_field_changed(self.field_probe.E_field)
    
    def calc_result(self, model):
        pass
        #self.calc_intensity(self.field_probe.E_field)
        
    def _pan_tool_default(self):
        #return PanTool()
        return ProbePlanePanTool()
    
    def _field_probe_changed(self, old, new):
        if old is not None:
            old.on_trait_change(self.on_field_changed, "E_field", remove=True)
        new.on_trait_change(self.on_field_changed, "E_field")
        self.pan_tool.probe = new
        
    def on_field_changed(self, e_field):
        intensity_image = self.calc_intensity(e_field)
        
        probe = self.field_probe
        self.ds.set_data("img", intensity_image)
        side = probe.width
        yside = probe.height
        
        if self.plot is not None:
            a = probe.centre
            self.plot.title = f"Intensity @ ({a[0]:0.3f}, {a[1]:0.3f}, {a[2]:0.3f})"
            
            plot = self.image_plot
            #plot.
            plot.x_mapper.range.set_bounds(0,side)
            plot.y_mapper.range.set_bounds(0,yside)
            xdata = numpy.linspace(0,side,probe.size)
            ydata = numpy.linspace(0,yside, probe.size)
            plot.index.set_data(xdata,ydata)
            plot.request_redraw()
            
            
if __name__=="__main__":
    import numpy
    vm = IntensitySurface()
    
    _x = numpy.linspace(-1,1,50)
    _y = numpy.linspace(-1,1,50)
    x,y = numpy.meshgrid(_x,_y)
    
    z = numpy.exp(-(x**2 + y**2)/(0.3**2))
    
    vm.intensity_data = z
    
    vm.configure_traits()