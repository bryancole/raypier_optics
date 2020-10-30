import vtk

source = vtk.vtkProgrammableSource()
def Execute():
    grid = source.GetStructuredGridOutput()
    points = vtk.vtkPoints()

    nx,ny = 25,25
    
    for i in range(nx):
        for j in range(ny):
            points.InsertNextPoint(i/5.0,j/5.0,0.0)
            
    grid.SetDimensions(nx,ny,1)
    grid.SetPoints(points)

source.SetExecuteMethod(Execute)
    
map = vtk.vtkDataSetMapper()
map.SetInputConnection(source.GetOutputPort(2))

act = vtk.vtkActor()
act.SetMapper(map)

ren = vtk.vtkRenderer()
ren.AddActor(act)

renwin = vtk.vtkRenderWindow()
renwin.AddRenderer(ren)

iren = vtk.vtkRenderWindowInteractor()
iren.SetRenderWindow(renwin)
iren.Start()