from traits.api import Int, Str, HasTraits, Tuple, Array, Callable

from vtk.vtkCommonDataModel import vtkDataObject
from vtk.vtkCommonExecutionModel import vtkAlgorithm
from vtk.vtkCommonExecutionModel import vtkDemandDrivenPipeline
from vtk.vtkCommonExecutionModel import vtkStreamingDemandDrivenPipeline
from vtk.util import numpy_support

import numpy
import vtk
from tvtk.api import tvtk

to_tvtk = tvtk.to_tvtk


class VTKAlgorithm(HasTraits):
    """This is a superclass which can be derived to implement
    Python classes that work with vtkPythonAlgorithm. It implements
    Initialize(), ProcessRequest(), FillInputPortInformation() and
    FillOutputPortInformation().

    Initialize() sets the input and output ports based on data
    members.

    ProcessRequest() calls RequestXXX() methods to implement
    various pipeline passes.

    FillInputPortInformation() and FillOutputPortInformation() set
    the input and output types based on data members.
    """
    
    number_of_input_ports = Int(1)
    input_type = Str("vtkDataSet")
    number_of_output_ports = Int(1)
    output_type = Str("vtkPolyData")
    
#     def __init__(self, **traits):
#         klass = 
#         super(VTKAlgorithm,self).__init__(klass, **traits)

    def Initialize(self, vtkself):
        """Sets up number of input and output ports based on
        NumberOfInputPorts and NumberOfOutputPorts."""

        vtkself.SetNumberOfInputPorts(self.number_of_input_ports)
        vtkself.SetNumberOfOutputPorts(self.number_of_output_ports)

    def GetInputData(self, inInfo, i, j):
        """Convenience method that returns an input data object
        given a vector of information objects and two indices."""

        return inInfo[i].GetInformationObject(j).Get(vtkDataObject.DATA_OBJECT())

    def GetOutputData(self, outInfo, i):
        """Convenience method that returns an output data object
        given an information object and an index."""
        return outInfo.GetInformationObject(i).Get(vtkDataObject.DATA_OBJECT())

    def RequestDataObject(self, vtkself, request, inInfo, outInfo):
        """Overwritten by subclass to manage data object creation.
        There is not need to overwrite this class if the output can
        be created based on the OutputType data member."""
        return 1

    def RequestInformation(self, vtkself, request, inInfo, outInfo):
        """Overwritten by subclass to provide meta-data to downstream
        pipeline."""
        return 1

    def RequestUpdateExtent(self, vtkself, request, inInfo, outInfo):
        """Overwritten by subclass to modify data request going
        to upstream pipeline."""
        return 1

    def RequestData(self, vtkself, request, inInfo, outInfo):
        """Overwritten by subclass to execute the algorithm."""
        raise NotImplementedError('RequestData must be implemented')

    def ProcessRequest(self, vtkself, request, inInfo, outInfo):
        """Splits a request to RequestXXX() methods."""
        print(vtkself)
        print(request)
        print(inInfo)
        print(outInfo)
        if request.Has(vtkDemandDrivenPipeline.REQUEST_DATA_OBJECT()):
            return self.RequestDataObject(vtkself, request, inInfo, outInfo)
        elif request.Has(vtkDemandDrivenPipeline.REQUEST_INFORMATION()):
            return self.RequestInformation(vtkself, request, inInfo, outInfo)
        elif request.Has(vtkStreamingDemandDrivenPipeline.REQUEST_UPDATE_EXTENT()):
            return self.RequestUpdateExtent(vtkself, request, inInfo, outInfo)
        elif request.Has(vtkDemandDrivenPipeline.REQUEST_DATA()):
            return self.RequestData(vtkself, request, inInfo, outInfo)
        return 1

    def FillInputPortInformation(self, vtkself, port, info):
        """Sets the required input type to InputType."""
        info.Set(vtkAlgorithm.INPUT_REQUIRED_DATA_TYPE(), self.input_type)
        return 1

    def FillOutputPortInformation(self, vtkself, port, info):
        """Sets the default output type to OutputType."""
        info.Set(vtkDataObject.DATA_TYPE_NAME(), self.output_type)
        return 1


class PythonAlgorithmBase(tvtk.PythonAlgorithm):
    
    number_of_input_ports = Int(1)
    input_type = Str("vtkDataSet")
    number_of_output_ports = Int(1)
    output_type = Str("vtkPolyData")
    
    class InternalAlgorithm(object):
        "Internal class. Do not use."
        def Initialize(self, vtkself):
            pass

        def FillInputPortInformation(self, vtkself, port, info):
            return to_tvtk(vtkself).FillInputPortInformation(port, info)

        def FillOutputPortInformation(self, vtkself, port, info):
            return to_tvtk(vtkself).FillOutputPortInformation(port, info)

        def ProcessRequest(self, vtkself, request, inInfo, outInfo):
            return to_tvtk(vtkself).ProcessRequest(request, inInfo, outInfo)

    def __init__(self, *args, **kwds):
        """Sets up default NumberOfInputPorts, NumberOfOutputPorts,
        InputType and OutputType that are used by various methods.
        Make sure to call this method from any subclass' __init__"""

        number_of_input_ports = self.number_of_input_ports
        number_of_output_ports = self.number_of_output_ports
        super(PythonAlgorithmBase, self).__init__(*args, **kwds)
        vtk_obj = self._vtk_obj
        vtk_obj.SetPythonObject(PythonAlgorithmBase.InternalAlgorithm())
        self.number_of_input_ports = number_of_input_ports
        self.number_of_output_ports = number_of_output_ports
        

    def get_input_data(self, inInfo, i, j):
        """Convenience method that returns an input data object
        given a vector of information objects and two indices."""

        return inInfo[i].GetInformationObject(j).Get(vtkDataObject.DATA_OBJECT())

    def get_output_data(self, outInfo, i):
        """Convenience method that returns an output data object
        given an information object and an index."""
        return outInfo.GetInformationObject(i).Get(vtkDataObject.DATA_OBJECT())

    def FillInputPortInformation(self, port, info):
        """Sets the required input type to InputType."""
        info.Set(vtkAlgorithm.INPUT_REQUIRED_DATA_TYPE(), self.input_type)
        return 1

    def FillOutputPortInformation(self, port, info):
        """Sets the default output type to OutputType."""
        info.Set(vtkDataObject.DATA_TYPE_NAME(), self.output_type)
        return 1

    def ProcessRequest(self, request, inInfo, outInfo):
        """Splits a request to RequestXXX() methods."""
        if request.Has(vtkDemandDrivenPipeline.REQUEST_DATA_OBJECT()):
            return self.RequestDataObject(request, inInfo, outInfo)
        elif request.Has(vtkDemandDrivenPipeline.REQUEST_INFORMATION()):
            return self.RequestInformation(request, inInfo, outInfo)
        elif request.Has(vtkStreamingDemandDrivenPipeline.REQUEST_UPDATE_EXTENT()):
            return self.RequestUpdateExtent(request, inInfo, outInfo)
        elif request.Has(vtkDemandDrivenPipeline.REQUEST_DATA()):
            return self.RequestData(request, inInfo, outInfo)
        return 1

    def RequestDataObject(self, request, inInfo, outInfo):
        """Overwritten by subclass to manage data object creation.
        There is not need to overwrite this class if the output can
        be created based on the OutputType data member."""
        return 1

    def RequestInformation(self, request, inInfo, outInfo):
        """Overwritten by subclass to provide meta-data to downstream
        pipeline."""
        return 1

    def RequestUpdateExtent(self, request, inInfo, outInfo):
        """Overwritten by subclass to modify data request going
        to upstream pipeline."""
        return 1

    def RequestData(self, request, inInfo, outInfo):
        """Overwritten by subclass to execute the algorithm."""
        raise NotImplementedError('RequestData must be implemented')
        
        
class EmptyGridSource(PythonAlgorithmBase):
    number_of_input_ports = Int(0)
    output_type = Str("vtkStructuredPoints")
    
    dimensions = Tuple(1, 1, 1,)
    origin = Tuple(0.0, 0.0, 0.0)
    spacing = Tuple(0.1, 0.1, 0.1)
    
    def RequestInformation(self, request, inInfo, outInfo):
        dims = self.dimensions
        info = outInfo.GetInformationObject(0)
        info.Set(vtk.vtkStreamingDemandDrivenPipeline.WHOLE_EXTENT(),
            (0, dims[0]-1, 0, dims[1]-1, 0, dims[2]-1), 6)
        return 1
    
    def RequestData(self, request, inInfo, outInfo):
        output = vtk.vtkStructuredPoints.GetData(outInfo)
        dims = self.dimensions
        output.SetDimensions( dims )
        output.SetSpacing( *self.spacing)
        output.SetOrigin( *self.origin)
        return 1
    
    
class NumpyImageSource(EmptyGridSource):
    number_of_input_ports = Int(0)
    output_type = Str("vtkImageData")
    
    dimensions = Tuple(1, 1, 1,)
    origin = Tuple(0.0, 0.0, 0.0)
    spacing = Tuple(0.1, 0.1, 0.1)
    
    image_data = Array #Should be a 2D or 3D array
    
    def _image_data_changed(self, data):
        dims = list(data.shape)
        dims.extend([1,]*(3-len(dims)))
        self.dimensions = tuple(dims)
        self.modified()
        
    def RequestData(self, request, inInfo, outInfo):
        output = vtk.vtkImageData.GetData(outInfo)
        dims = self.dimensions
        output.SetDimensions( dims )
        output.SetSpacing( *self.spacing)
        output.SetOrigin( *self.origin)
        pd = output.GetPointData()
        data = numpy_support.numpy_to_vtk(self.image_data.ravel(), 1)
        pd.SetScalars(data)
        return 1
        
    
if __name__=="__main__":
    from tvtk.pyface.scene_model import SceneModel
    from tvtk.pyface.scene_editor import SceneEditor
    from traits.api import Instance
    from traitsui.api import View, Item
    import vtk
    
    class TestScene(HasTraits):
        scene = Instance(SceneModel, ())
        
        view = View(Item("scene", show_label=False, editor=SceneEditor()))
        
        def show(self):            
            grid = EmptyGridSource(dimensions=(20,30,40))
            src_output_port = grid.output_port
            src_output_port = grid._vtk_obj.GetOutputPort()
            assert src_output_port is not None
            extract = tvtk.ImageDataGeometryFilter(input_connection=src_output_port)
            mapper = tvtk.PolyDataMapper(input_connection=extract.output_port)
            act = tvtk.Actor(mapper=mapper)
            scene = self.scene
            scene.add_actor(act)
            self.configure_traits()
    
    scene = TestScene()
    scene.show()
