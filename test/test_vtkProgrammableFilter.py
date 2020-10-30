import vtk
from vtk.util.numpy_support import numpy_to_vtk
import numpy

x,y,z = numpy.ogrid[-10:10:50j,-10:10:50j,-10:10:50j]
data = (x**2 + y**2 + z**2)/(5**2)
data = data.ravel().astype("f")


src = vtk.vtkProgrammableSource()
###Initialise the OutputInformation on the ProgrammableSource
executive = src.GetExecutive()
outInfo = executive.GetOutputInformation(1)
outInfo.Set(executive.WHOLE_EXTENT(), 0, 49, 0, 49, 0, 49)

def make_grid():    
    sp = src.GetStructuredPointsOutput()
    sp.SetExtent(0,49,0,49,0,49)
    sp.SetOrigin(-10,-10,-10)
    sp.SetSpacing(0.4,0.4,0.4)
    sp.GetPointData().SetScalars(numpy_to_vtk(data))

src.SetExecuteMethod(make_grid)

func = vtk.vtkSphere()
func.SetCenter(0,0,0)
func.SetRadius(5.0)

clip = vtk.vtkImageMarchingCubes()
clip.SetInputConnection(0,src.GetOutputPort(1))
clip.SetValue(0,1.0)

map = vtk.vtkPolyDataMapper()
map.SetInputConnection(clip.GetOutputPort(0))
map.ScalarVisibilityOff()
 
surfaceActor = vtk.vtkActor()
surfaceActor.SetMapper(map)

ren = vtk.vtkRenderer()
renWin = vtk.vtkRenderWindow()
renWin.AddRenderer(ren)
iren = vtk.vtkRenderWindowInteractor()
iren.SetRenderWindow(renWin)
ren.AddActor(surfaceActor)
iren.Initialize()
renWin.Render()
iren.Start()