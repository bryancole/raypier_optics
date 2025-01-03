"""
A module with classes to provide 2D image display using the vtkChart API , with 
traitsui integration.
"""


import numpy as np

from tvtk.api import tvtk

from traits.api import HasTraits, Instance, Array, Tuple, Float, observe
from traitsui.api import CustomEditor, View, Item, VGroup
from tvtk.pyface.ui.qt4.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor


def ctx_view_widget_factory(parent, editor, *args, **kwds):
    cview = editor.object._context_view
    widget = QVTKRenderWindowInteractor(parent=None)
    cview.render_window = widget.GetRenderWindow()
    #print("Interactor:", widget._Iren)
    return widget


CtxEditor = CustomEditor(factory=ctx_view_widget_factory)


class ImgDisplayModel(HasTraits):
    img_data = Array() ### n_dims=2, shape=(N,M)
    x_bounds = Tuple(Float, Float)
    y_bounds = Tuple(Float, Float)
    
    _chart = Instance(tvtk.ChartHistogram2D)
    
    _context_view = Instance(tvtk.ContextView)
    
    _img_data = Instance(tvtk.ImageData)
    
    traits_view = View(
            VGroup(
                Item("_context_view", style="custom", editor=CtxEditor)
                ),
            resizable=True
        )
    
    @observe("img_data")
    def on_img_changed(self, evt):
        data = evt.new
        nx,ny = data.shape
        d = tvtk.ImageData(dimensions=(nx,ny,1))
        d.point_data.scalars = data.ravel()
        d.dimensions = (nx,ny,1)
        d.origin = np.array([self.x_bounds[0],self.y_bounds[0],0.])
        d.extent = (0,nx-1,0,ny-1,0,0)
        x_extent = self.x_bounds[1]-self.x_bounds[0]
        y_extent = self.y_bounds[1]-self.y_bounds[0]
        sp = np.array([x_extent/(nx-1.), y_extent/(ny-1.), 1.0], dtype=np.double)
        print("spacing:", sp, x_extent, y_extent)
        d.spacing = sp
        self._img_data = d
        d.compute_bounds()
        print(d)
        print(data.shape, data)
        self._chart.set_input_data(d)
    
    def __chart_default(self):
        colors = tvtk.NamedColors()
        background_color = colors.get_color3d('SlateGray')
        title_color = colors.get_color3d('Orange')
        axis_title_color = colors.get_color3d('Orange')
        axis_label_color = colors.get_color3d('Black')
        legend_background_color = colors.get_color4ub('Tomato')
        # Define a chart
        chart = tvtk.ChartHistogram2D(title='2D Histogram')
        chart.title_properties.font_size = 36
        chart.title_properties.color = tuple(title_color)
    
        # Chart Axes.
        chart.get_axis(0).title_properties.font_size = 24
        chart.get_axis(0).title_properties.color = tuple(axis_title_color)
        chart.get_axis(0).label_properties.color = tuple(axis_label_color)
        chart.get_axis(0).label_properties.font_size = 18
    
        chart.get_axis(1).title_properties.font_size = 24
        chart.get_axis(1).title_properties.color = tuple(colors.get_color3d('orange'))
        chart.get_axis(1).label_properties.color = tuple(axis_label_color)
        chart.get_axis(1).label_properties.font_size = 18
    
        # Chart Legend.
        chart.legend.draw_border = True
        chart.legend.brush.color = tuple(legend_background_color)
        
        transfer_function = tvtk.ColorTransferFunction()
        transfer_function.add_hsv_segment(0.0, 0.0, 1.0, 1.0, 0.3333, 0.3333, 1.0, 1.0)
        transfer_function.add_hsv_segment(0.3333, 0.3333, 1.0, 1.0, 0.6666, 0.6666, 1.0, 1.0)
        transfer_function.add_hsv_segment(0.6666, 0.6666, 1.0, 1.0, 1.0, 0.2, 1.0, 0.3)
        transfer_function.build()
        chart.set_transfer_function(transfer_function)
        
        return chart
    
    def __context_view_default(self):
        view = tvtk.ContextView()
        view.scene.add_item(self._chart)
        return view
    
    
if __name__=="__main__":
    x,y = np.mgrid[0.0:25.0:500j,0.0:25.0:500j]
    z = np.sin(x)*np.cos(y)
    print(z)
    model = ImgDisplayModel(x_bounds=(-1.,1.), y_bounds=(0.0,5.0))
    model.img_data = z
    model.configure_traits()
    