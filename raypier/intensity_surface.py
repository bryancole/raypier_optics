

from traits.api import Float, Instance, Enum, Array, on_trait_change

from traitsui.api import View, Item, VGroup, HGroup

from tvtk.api import tvtk
from tvtk.pyface.scene_model import SceneModel
from tvtk.pyface.scene_editor import SceneEditor

from raypier.results import Result
from raypier.fields import EFieldPlane
from raypier.vtk_algorithms import NumpyImageSource
from .core.unwrap2d import unwrap2d

import numpy


class IntensitySurface(Result):
    name = "Field Plane Surface View"
    field_probe = Instance(EFieldPlane)
    
    display = Enum("Intensity", "E_x", "E_y", "E_z", "Phase_x", "Phase_y", "Phase_z")
    
    intensity_data = Array()
    
    scale = Float(1.0)
    
    scene = Instance(SceneModel, transient=True)
    
    #_ds = Instance(tvtk.DepthSortPolyData, transient=True)
    
    _data_source = Instance(NumpyImageSource, (), transient=True)
    _warp = Instance(tvtk.WarpScalar, transient=True)
    _map = Instance(tvtk.PolyDataMapper, transient=True)
    _axes = Instance(tvtk.CubeAxesActor2D, transient=True)
    
    traits_view = View(VGroup(
                        HGroup(
                            Item("display", show_label=False),
                            Item("scale", show_label=False)
                            ),
                        Item("scene", editor=SceneEditor(), show_label=False)
                        ))
    
    def _intensity_data_changed(self, data):
        sdata = data
        self._data_source.image_data = sdata
        map = self._map
        if map is not None:
            map.scalar_range = (sdata.min(), sdata.max())
        self.scene.render()
        
    @on_trait_change("scale, intensity_data")
    def _scale_changed(self, value):
        w = self._warp
        value = self.scale
        low = self.intensity_data.min()
        high = self.intensity_data.max()
        m=self._map
        if m is not None:
            m.scalar_range = (low, high)
        if w is not None:
            w.scale_factor = value/(high-low)
            w.modified()
        
    def _scene_default(self):
        scene = SceneModel(background=(0,0,0))
        
        src = self._data_source
        src.update()
        
        geom = tvtk.ImageDataGeometryFilter(input_connection=src.output_port)
        
        warp = tvtk.WarpScalar(input_connection=geom.output_port,
                               scale_factor=self.scale)
        self._warp = warp
        
        norm = tvtk.PolyDataNormals(input_connection=warp.output_port)
        
#         ###Doesn't work
#         ds = tvtk.DepthSortPolyData(input_connection=norm.output_port,
#                                     camera=tvtk.Camera())
#         self._ds = ds
        
        mapper = tvtk.PolyDataMapper(input_connection=norm.output_port)
        #mapper.scalar_range = (-2.,2.)
        self._map = mapper
        
        self._scale_changed(self.scale)
        
        act = tvtk.Actor(mapper=mapper)
        scene.add_actor(act)
        
        bar = tvtk.ScalarBarActor(lookup_table=mapper.lookup_table,
                                  number_of_labels=5)
        
        scene.add_actor(bar)
        
        axes = tvtk.CubeAxesActor2D(fly_mode="outer_edges",
                                    font_factor=1.5,
                                    scaling=True,
                                    z_axis_visibility=True,
                                    x_label="W",
                                    y_label="H",
                                    z_label="I")
        axes.set_input_connection(norm.output_port)
        self._axes = axes
        #axes.camera = scene.renderer.camera
        scene.add_actor(axes)
        scene.on_trait_change(self.on_activated, "activated")
        
#         def on_sync_camera(camera):
#             print( "CAM")
#             if camera is not None:
#                 ds.camera = camera
#                 
#         scene.on_trait_change(on_sync_camera, "_camera")
        
        return scene
    
    def on_activated(self, obj, val, old, new):
        a = self._axes
        a.camera = obj.camera
        a.use_ranges = True
        
    
    def calc_intensity(self, e_field):
        E = e_field
        mode = self.display
        select = {"x":0, "y":1, "z":2}
        if mode=="Intensity":
            U = (E.real**2).sum(axis=-1) + (E.imag**2).sum(axis=-1)
        elif mode.startswith("E"):
            idx = select[mode[-1]]
            U = E[:,:,idx].real
        elif mode.startswith("Phase"):
            idx = select[mode[-1]]
            e = E[:,:,idx]
            U = numpy.arctan2(e.imag, e.real)
            U, res = unwrap2d(U, anchor=(U.shape[0]//2,U.shape[1]//2))
        self.intensity_data = U
        return U
    
    def _display_changed(self):
        self.on_field_changed(self.field_probe.E_field)
    
    def calc_result(self, model):
        pass
        #self.calc_intensity(self.field_probe.E_field)
    
    def _field_probe_changed(self, old, new):
        if old is not None:
            old.on_trait_change(self.on_field_changed, "E_field", remove=True)
        new.on_trait_change(self.on_field_changed, "E_field")
        
    def on_field_changed(self, e_field):
        self.calc_intensity(e_field)
        p = self.field_probe
        data = self.intensity_data
        self._axes.ranges=[0, p.width*1000, 0, p.height*1000, data.min(), data.max()]
        
            
            
if __name__=="__main__":
    import numpy
    vm = IntensitySurface()
    
    _x = numpy.linspace(-1,1,50)
    _y = numpy.linspace(-1,1,50)
    x,y = numpy.meshgrid(_x,_y)
    
    z = numpy.exp(-(x**2 + y**2)/(0.5**2))
    
    vm.intensity_data = z
    
    vm.configure_traits()
    