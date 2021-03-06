
from traits.api import Float, Instance, \
        List, observe, Tuple, Str, ComparisonMode

from traitsui.api import View, Item, VGroup, Group, ListEditor

from tvtk.api import tvtk
import numpy
from operator import attrgetter

from raypier.bases import NumEditor, Traceable

from raypier.core.cmaterials import CoatedDispersiveMaterial, PECMaterial
from raypier.shapes import BaseShape
from raypier.lenses import BaseLens
from raypier.core.ctracer import FaceList
from . import faces
from .materials import OpticalMaterial, air


alist_editor = ListEditor(use_notebook=True,
                           deletable=False,
                           selected='selected',
                           export='DockWindowShell',
                           page_name='.name')


class GeneralLens(BaseLens):
    #: The name of the object in the model-tree
    name = Str("General Lens")
    
    #: An instance of a :py:class:`raypier.shapes.BaseShape` subclass, defining
    #: the 2D outline of the optic.
    shape = Instance(BaseShape)
    
    #: A list of :py:class:`raypier.faces.BaseFace` objects which defines the
    #: geometry of the surfaces comprising this optic.
    surfaces = List(faces.BaseFace, comparison_mode=ComparisonMode.identity)
    
    #: A list of :py:class:`raypier.materials.OpticalMaterial` objects describing 
    #: the materials properties of the dielectrics between the previously given
    #: surfaces. For a list of N surfaces, you should provide (N-1) materials.
    materials = List(OpticalMaterial, comparison_mode=ComparisonMode.identity)
    
    _grid = Instance(tvtk.PlaneSource, ())
    _clip = Instance(tvtk.ClipPolyData, ())
    
    #: A :py:class:`raypier.materials.OpticalMaterial` object to define the 
    #: characteristics of a single-layer dielectric coating applied to the outer
    #: surfaces of the object.
    coating_material = Instance(OpticalMaterial, ())
    
    #: The thickness, in microns, of the coating on the outer surfaces.
    coating_thickness = Float(0.25, desc="Thickness of the AR coating, in microns")
    
    
    grid_resolution = Tuple((200,200))
    glass_colour = Tuple((204,204,255,102))
    mirror_colour = Tuple((104,64,104,255))
    
    
    
    traits_view = View(Group(
                Traceable.uigroup,
                Group(
                    Item("shape", show_label=False, style="custom"),
                    label="Outline", dock="tab"
                    ),
                Group(
                    Item("surfaces", style="custom",
                         editor=alist_editor, 
                         show_label=False),
                    label="Faces", dock="tab"
                    ),
                Group(
                    VGroup(
                    Item("materials", style="custom", editor=alist_editor, show_label=False),
                    ),
                    Item("coating_material", style="custom"),
                    Item("coating_thickness", editor=NumEditor),
                    label="Materials", dock="tab"
                    ),
                layout="tabbed"
            )
        )
    
    def _actors_default(self):
        pipeline = self.pipeline
        
        map = tvtk.PolyDataMapper(input_connection=pipeline.output_port,
                                  scalar_mode="use_cell_data",
                                  color_mode="direct_scalars")
        map.scalar_visibility = True
        act = tvtk.Actor(mapper=map)
        actors = tvtk.ActorCollection()
        actors.append(act)
        return actors
    
    def _pipeline_default(self):
        nx,ny = self.grid_resolution
        xmin,xmax,ymin,ymax,zmin,zmax = self.grid_extent
        
        shape = self.shape
        func = shape.impl_func
        
        _grids = []
        
        append = tvtk.AppendPolyData()
        
        grid = self._grid
        grid.origin = (xmin,ymin,0)
        grid.point1=(xmin,ymax,0)
        grid.point2=(xmax,ymin,0)
        grid.x_resolution=nx
        grid.y_resolution=ny
    
        clip = self._clip
        clip.inside_out=True
        clip.input_connection = grid.output_port
        clip.clip_function = func
        
        bounds = tvtk.FeatureEdges(input_connection=clip.output_port)
        bounds.extract_all_edge_types_off()
        bounds.boundary_edges = True
        bounds.coloring = False
        
        for face in self.surfaces:
            attrb = tvtk.ProgrammableAttributeDataFilter(input_connection=clip.output_port)
            
            ###Beware, references to outer scope change in the loop. Capture state using kwd-args.
            def execute( *args , _face=face, _attrb=attrb):
                in_data = _attrb.get_input_data_object(0,0)
                points = in_data.points.to_array()
                z = _face.cface.eval_z_points(points.astype('d'))
                out = _attrb.get_output_data_object(0)
                out.point_data.scalars=z
                ncells = in_data.polys.number_of_cells
                colors = numpy.zeros((ncells,4), numpy.uint8)
                if _face.mirror:
                    colors[:] = self.mirror_colour
                else:
                    colors[:] = self.glass_colour
                out.cell_data.scalars = colors

            attrb.set_execute_method(execute)
            
            warp = tvtk.WarpScalar(input_connection=attrb.output_port,scale_factor=1.0)
            warp.normal=(0,0,1)
            warp.use_normal=True
            
            append.add_input_connection(warp.output_port)
                
        skirt = tvtk.ProgrammableFilter(input_connection=bounds.output_port)
        def calc_skirt(_skirt=skirt):
            in_data = _skirt.get_input_data_object(0,0)
            out = _skirt.get_output_data_object(0)
            points = in_data.points.to_array().astype('d')
            z_top = self.surfaces[0].cface.eval_z_points(points)
            z_bottom = self.surfaces[-1].cface.eval_z_points(points)
            
            size = points.shape[0]
            points_out = numpy.vstack([points,points])
            points_out[:size,2] = z_top
            points_out[size:,2] = z_bottom
            
            cell_data = in_data.lines.to_array().reshape(-1,3)
            ncells = cell_data.shape[0]
            cells_out = numpy.empty((ncells, 5))
            cells_out[:,0] = 4
            cells_out[:,1] = cell_data[:,1]
            cells_out[:,2] = cell_data[:,2]
            cells_out[:,3] = cell_data[:,2] + size
            cells_out[:,4] = cell_data[:,1] + size    
            out.points = points_out.reshape(-1,3)
            quads = tvtk.CellArray()
            quads.set_cells(5, cells_out)
            out.polys = quads
            out.lines = None

            colors = numpy.zeros((ncells*3,4), numpy.uint8)
            colors[:] = self.glass_colour
            out.cell_data.scalars = colors
            
        skirt.set_execute_method(calc_skirt)
        
        clean = tvtk.CleanPolyData(input_connection=skirt.output_port)
        
        append.add_input_connection(skirt.output_port)
        
        norms = tvtk.PolyDataNormals(input_connection=append.output_port)
        
        transF = tvtk.TransformFilter(input_connection=norms.output_port, 
                                      transform=self.transform)
        return transF
            
            
    def __init__(self, *args, **kwds):
        super().__init__(*args, **kwds)
        self.config_cfaces()
    
    def _faces_default(self):
        fl = FaceList(owner=self)
        return fl
    
    def config_cfaces(self):
        surfaces = sorted(self.surfaces, key=attrgetter("z_height"))
        cfaces = []
        mats = [air,] + self.materials + [air,]
        if len(mats) < len(surfaces) + 1:
            mats += [air,] * (len(surfaces)-len(mats) + 1)
        ###FIXME: need to iterate in order of z_height to get the
        ###the materials on the correct sides.
        for i, face in enumerate(surfaces):
            if face.trace:
                cface = face.cface
                cface.owner = self
                cface.shape = self.shape.cshape
                if face.mirror:
                    cface.material = PECMaterial()
                else:
                    mat = CoatedDispersiveMaterial()
                    mat.dispersion_inside = mats[i].dispersion_curve
                    mat.dispersion_outside = mats[i+1].dispersion_curve
                    if i in {0, len(surfaces)-1}:
                        mat.dispersion_coating = self.coating_material.dispersion_curve
                        mat.coating_thickness = self.coating_thickness
                    else:
                        mat.coating_thickness = 0.0
                    cface.material = mat
                cfaces.append(face.cface)
        self.faces.faces = cfaces
    
    @observe("materials.items.dispersion_curve, coating_material.dispersion_curve, coating_thickness")
    def on_material_change(self, evt):
        self.config_cfaces()
        self.update = True
        
    @observe("surfaces.items.updated, shape.updated")
    def on_surfaces_changed(self, evt):
        self.config_cfaces()
        print("Face changed", evt)
        self.on_bounds_change(None)
        self.update = True

    @observe("shape.impl_func, surfaces")
    def on_impl_func_change(self, evt):
        self._clip.clip_function = self.shape.impl_func
        self.render=True
    
    @observe("shape.bounds, surfaces")
    def on_bounds_change(self, evt):
        nx, ny = self.grid_resolution
        xmin, xmax, ymin, ymax = self.shape.bounds
        zmin, zmax = self.eval_z_extent(xmin, xmax, ymin, ymax)
        self.grid_extent = (xmin, xmax, ymin, ymax, zmin, zmax)
        g = self._grid
        g.origin=(xmin,ymin,0)
        g.point1=(xmin,ymax,0)
        g.point2=(xmax,ymin,0)
        g.x_resolution=nx
        g.y_resolution=ny
        g.modified()
    
    def eval_z_extent(self, xmin, xmax, ymin, ymax):
        nx, ny = self.grid_resolution
        try:
            top_face = self.surfaces[0]
            bottom_face = self.surfaces[-1]
        except IndexError:
            return 0,0
        x = numpy.linspace(xmin, xmax, nx)
        y = numpy.linspace(ymin, ymax, ny)
        tzmin, tzmax = top_face.cface.eval_z_extent(x,y)
        bzmin, bzmax = bottom_face.cface.eval_z_extent(x,y)
        return ( min(tzmin,bzmin), max(tzmax, bzmax) )
    
    