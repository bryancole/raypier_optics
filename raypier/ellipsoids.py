#    Copyright 2009, Teraview Ltd., Bryan Cole
#
#    This file is part of Raypier.
#
#    Raypier is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""
a module for parabolic optics e.g. OAPs
"""
from traits.api import HasTraits, Array, Float, Complex,\
            Property, List, Instance, on_trait_change, Range, \
            Tuple, Event, cached_property, Set, Int, Trait, Bool, \
            PrototypedFrom, BaseTuple
from traitsui.api import View, Item, ListEditor, VSplit,\
            RangeEditor, ScrubberEditor, HSplit, VGroup
from tvtk.api import tvtk
import numpy

from raypier.bases import Traceable, normaliseVector, NumEditor,\
            VectorEditor, transformPoints, transformNormals
#from raypier.sources import RayCollection
from raypier.mirrors import BaseMirror
from raypier.core.cfaces import EllipsoidalFace
from raypier.core.ctracer import FaceList

from raypier.vtk_algorithms import EmptyGridSource

class MinMax(BaseTuple):
    def validate(self, obj, name, value):
        valid = super(MinMax,self).validate(obj, name, value)
        a,b = valid
        if a>b:
            valid = (b,a)
        return valid
        


class Ellipsoid(BaseMirror):
    """
    An off-axis ellipsoidal mirror formed into a axis-aligned block.
    """    

    name = "Ellipsoid"
    abstract = False
    
    #: The position of the first focus, in local coordinates.
    focus1 = Tuple(-50.,0.,0.)
    
    #: The position of the second focus, in local coordinates.
    focus2 = Tuple(0., 50., 0.)
    
    #: twice the major axis length, or the distance from one
    #: focus to the ellipsoid edge to the other focus
    size = Float(100.0, desc="twice the major axis length, or the distance from one\
 focus to the ellipsoid edge to the other focus")
    
    #: A 2-tuple giving the lower and upper x-axis bounds of the mirror block. 
    X_bounds = MinMax(-25., 25.)
    
    #: A 2-tuple giving the lower and upper y-axis bounds of the mirror block. 
    Y_bounds = MinMax(-25., 25.)
    
    #: A 2-tuple giving the lower and upper z-axis bounds of the mirror block. 
    Z_bounds = MinMax(0., 50.)
    
    #: Render a glyph at the foci of the ellipsoid.
    show_foci = Bool(True)
    
    foci_Actors = Instance(tvtk.ActorCollection, ())
    
    axes = Property(Tuple, depends_on="focus1, focus2, size",
                 desc="(major, minor) axis lengths")
    
    ellipse_trans = Instance(tvtk.Transform, (), transient=True)
    
    combined_trans = Instance(tvtk.Transform, transient=True)
    
    f1_glyph = Instance(tvtk.SphereSource, (), transient=True)
    f2_glyph = Instance(tvtk.SphereSource, (), transient=True)
    
    f1_act = Instance(tvtk.Follower, (), transient=True)
    f2_act = Instance(tvtk.Follower, (), transient=True)
    
    vtk_grid = Instance(EmptyGridSource, (), transient=True)
    vtk_quadric = Instance(tvtk.Quadric, (), transient=True)
    
    vtkproperty = tvtk.Property(opacity = 0.7,
                             color = (0.8,0.8,0))
    
    traits_view = View(VGroup(
                        Traceable.uigroup,
                        Item("focus1", editor=VectorEditor),
                        Item("focus2", editor=VectorEditor),
                        Item("show_foci"),
                        Item("size", editor=NumEditor),
                        Item("X_bounds"),
                        Item("Y_bounds"),
                        Item("Z_bounds")
                       ),
                       )
    
    def _faces_default(self):
        facelist = FaceList(owner=self)
        facelist.faces=[EllipsoidalFace(owner=self,
                                x1=self.X_bounds[0],
                                x2=self.X_bounds[1],
                                y1=self.Y_bounds[0],
                                y2=self.Y_bounds[1],
                                z1=self.Z_bounds[0],
                                z2=self.Z_bounds[1],
                                major=self.axes[0],
                                minor=self.axes[1])]
        return facelist
    
    @on_trait_change("focus1, focus2, size")
    def config_trans(self, *args):
        f1 = numpy.asarray(self.focus1)
        f2 = numpy.asarray(self.focus2)
        size = self.size
        centre = (f2+f1)/2.
        ellipse_t = self.ellipse_trans
        ellipse_t.identity()
        #ellipse major axis along the X axis
        delta = f2 - f1
        ax = normaliseVector(delta)
        axy = numpy.sqrt(ax[0]**2 + ax[1]**2)
        ellipse_t.rotate_y(numpy.arctan2(ax[2], axy)*180/numpy.pi)
        
        ellipse_t.rotate_z(-numpy.arctan2(ax[1],ax[0])*180/numpy.pi)
        
        ellipse_t.translate(*-centre)
    
    @cached_property
    def _get_axes(self):
        f1 = numpy.asarray(self.focus1)
        f2 = numpy.asarray(self.focus2)
        delta = f2 - f1
        size = self.size
        h = numpy.sqrt((delta**2).sum())
        a = size/2.
        b = 0.5 * numpy.sqrt(size**2 - h**2)
        return (a,b)
    
    def _combined_trans_default(self):
        t = tvtk.Transform()
        t.concatenate(self.transform)
        t.concatenate(self.ellipse_trans.linear_inverse)
        return t
    
    def _show_foci_changed(self, val):
        for act in self.foci_Actors:
            act.visibility = val
        self.render = True
    
    def update_grid(self):
        xmin, xmax = self.X_bounds
        ymin, ymax = self.Y_bounds
        zmin, zmax = self.Z_bounds
        
        source = self.vtk_grid
        size = 20
        dx = (xmax-xmin) / (size-1)
        dy = (ymax-ymin) / (size-1)
        dz = (zmax-zmin) / (size-1)
        source.dimensions = (size,size,size)
        source.origin = (xmin, ymin, zmin)
        source.spacing = (dx, dy, dz)
        
    @on_trait_change("X_bounds, Y_bounds, Z_bounds")
    def change_bounds(self):
        self.update_grid()
        self.vtk_grid.modified()
        self.update=True
    
    @on_trait_change("focus1, focus2, size")
    def config_pipeline(self, *args):
        tp = self.transform.transform_point
        self.f1_act.position = tp(self.focus1)
        self.f2_act.position = tp(self.focus2)
        
        f1 = numpy.asarray(self.focus1)
        f2 = numpy.asarray(self.focus2)
        
        self.f1_glyph.center = f1
        self.f2_glyph.center = f2
        
        centre = (f2+f1)/2.
        ellipse_t = self.ellipse_trans
        ellipse_t.identity()
        #ellipse major axis along the X axis
        delta = f2 - f1
        ax = normaliseVector(delta)
        axy = numpy.sqrt(ax[0]**2 + ax[1]**2)
        ellipse_t.rotate_y(numpy.arctan2(ax[2], axy)*180/numpy.pi)
        
        ellipse_t.rotate_z(-numpy.arctan2(ax[1],ax[0])*180/numpy.pi)
        
        ellipse_t.translate(*-centre)
        
        a,b = self.axes
        #b = 0.5 * numpy.sqrt(size**2 - h**2)  ##not required
        
        q = self.vtk_quadric
        A1 = A2 = 1/(b**2)
        A0 = 1/(a**2)  
        #A8 = 1 
        A9 = -1
        q.coefficients = (A0, A1, A2, 
                          0, 0, 0, 
                          0, 0, 0, 
                          A9)
        self.update=True
        
                                         
    def _pipeline_default(self):
        grid = self.vtk_grid
        #grid.set_execute_method(self.create_grid)
        grid.modified()
        
        quad = self.vtk_quadric
        quad.transform = self.ellipse_trans
        
        clip = tvtk.ClipVolume(input_connection=grid.output_port,
                                 clip_function=quad,
                                 inside_out=0)
        
        topoly = tvtk.GeometryFilter(input_connection=clip.output_port)
        norm = tvtk.PolyDataNormals(input_connection=topoly.output_port)
        transF = tvtk.TransformFilter(input_connection=norm.output_port, transform=self.transform)
        self.config_pipeline()
        grid.modified()
        return transF
    
    def get_actors(self, scene):
        actors = []
        
        sList = [self.f1_glyph, self.f2_glyph]
        cList = [(0,1,0),(1,0,0)]
        
        for s,c in zip(sList, cList):
            s.radius = 1.0
            map = tvtk.PolyDataMapper(input_connection=s.output_port)
            act = tvtk.Actor(mapper=map, user_transform=self.transform)
            act.property.color=c
            actors.append(act)
        
        line = tvtk.LineSource(point1=(-100,0,0),
                               point2=(100,0,0))
        t_line = tvtk.TransformFilter(input_connection=line.output_port,
                                      transform=self.ellipse_trans.linear_inverse)
        map = tvtk.PolyDataMapper(input_connection=t_line.output_port)
        act = tvtk.Actor(mapper=map, user_transform=self.transform)
        act.property.color=(0,0,0)
        actors.append(act)
        
        l1 = tvtk.VectorText(text="F1")
        l2 = tvtk.VectorText(text="F2")
        m1 = tvtk.PolyDataMapper(input_connection=l1.output_port)
        m2 = tvtk.PolyDataMapper(input_connection=l2.output_port)
        act1 = self.f1_act
        act2 = self.f2_act
        
        act1.mapper = m1
        act2.mapper = m2
        
        scale = (5,5,5)
        act1.scale = scale
        act2.scale = scale

        act1.property.color=(0,0,0)
        act2.property.color=(0,0,0)
        
        act1.position = self.focus1
        act2.position = self.focus2
        
        def on_editor(new_ed):
            if new_ed is not None:
                act1.camera = new_ed._camera
                act2.camera = new_ed._camera
        
        scene.on_trait_change(on_editor, "scene_editor")
        
        actors.append(act1)
        actors.append(act2)
        
        for actor in actors:
            self.actors.append(actor)
            self.foci_Actors.append(actor)
            actor.visibility = self.show_foci
        
        return self.actors
    
    def make_step_shape(self):
        from raypier.step_export import make_ellipsoid_mirror
        ell = make_ellipsoid_mirror(self.focus1, 
                                     self.focus2, 
                                     self.size/2.,
                                     self.X_bounds, 
                                     self.Y_bounds, 
                                     self.Z_bounds, 
                                     self.centre, 
                                     self.direction,
                                     self.x_axis)
        return ell, "yellow"
