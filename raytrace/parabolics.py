from enthought.traits.api import Float, Instance, on_trait_change
from enthought.traits.ui.api import View, Item, ListEditor, VSplit,\
            RangeEditor, ScrubberEditor, HSplit, VGroup
from enthought.tvtk.api import tvtk

from raytrace.bases import NumEditor, Traceable
from raytrace.mirrors import BaseMirror
from raytrace.cfaces import OffAxisParabolicFace
from raytrace.ctracer import FaceList


class OffAxisParabloid(BaseMirror):
    name = "Off Axis Parabolic"
    abstract = False
    diameter = Float(25.0)
    EFL = Float(50.0, desc="Effective Focal Length")
    height = Float(50.0)
    
    max_length=Float(200.0) #what does this do?
    
    vtk_grid = Instance(tvtk.ProgrammableSource, ())
    vtk_cylinder = Instance(tvtk.Cylinder, ())
    vtk_quadric = Instance(tvtk.Quadric, ())
    
    traits_view = View(VGroup(
                        Traceable.uigroup,
                       Item('height', editor=NumEditor),
                       Item('diameter', editor=NumEditor),
                       Item('EFL', editor=NumEditor),
                       Item('max_length', editor=NumEditor),
                       ),
                       )
    
    def _faces_default(self):
        facelist = FaceList(owner=self)
        facelist.faces=[OffAxisParabolicFace(owner=self,
                                             EFL=self.EFL,
                                             diameter=self.diameter)]
        return facelist
    
    @on_trait_change("EFL, diameter")
    def config_pipeline(self):
        face = self.faces.faces[0]
        face.EFL = EFL = self.EFL
        face.diameter = self.diameter
        
        rad = self.diameter/2
        h = self.height
        
        self.vtk_grid.modified()
                      
        cyl = self.vtk_cylinder
        cyl.center = (EFL,0,0)
        cyl.radius = rad
        
        q = self.vtk_quadric
        A0 = A1 = -1/(2*EFL)
        A8 = 1 
        A9 = EFL/2
        q.coefficients = (A0, A1, 0, 
                          0, 0, 0, 
                          0, 0, A8, 
                          A9)
        self.update=True
        
    def create_grid(self):
        EFL = self.EFL
        r = self.diameter/2
        h = self.height
        l = self.max_length
        source = self.vtk_grid
        sp = source.structured_points_output
        size = 20
        spacing = 2*r / (size-1)
        lsize = int(l/spacing)
        sp.dimensions = (size,size,lsize)
        sp.whole_extent=(0,size,0,size,0,lsize)
        sp.origin = (EFL - r, -r, -h)
        sp.spacing = (spacing, spacing, spacing)
        sp.set_update_extent_to_whole_extent()
        
    def _pipeline_default(self):
        grid = self.vtk_grid
        grid.set_execute_method(self.create_grid)
        grid.modified()
        
        trans = tvtk.Transform()
        trans.rotate_x(90.)
        cyl = self.vtk_cylinder
        cyl.transform = trans
        
        clip1 = tvtk.ClipVolume(input=grid.structured_points_output,
                                 clip_function=self.vtk_cylinder,
                                 inside_out=1)
        
        clip2 = tvtk.ClipDataSet(input = clip1.output,
                                 clip_function=self.vtk_quadric,
                                 inside_out=1)
        
        topoly = tvtk.GeometryFilter(input=clip2.output)
        norms = tvtk.PolyDataNormals(input=topoly.output)
        
        transF = tvtk.TransformFilter(input=norms.output, transform=self.transform)
        self.config_pipeline()
        #clip1.update()
        grid.modified()
        return transF