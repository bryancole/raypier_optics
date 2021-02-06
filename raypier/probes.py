#    Copyright 2009, Bryan Cole, Teraview Ltd.
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

from traits.api import on_trait_change, Float, Instance,Event, Int,\
        Property, Button, Str, Array, List, observe

from traitsui.api import View, Item, VGroup, DropEditor

from tvtk.api import tvtk
import numpy
import time

from raypier.bases import Probe, Traceable, NumEditor, Vector, BaseRayCollection
from raypier.sources import RayCollection, GaussletCollection
from raypier.core.ctracer import FaceList, select_ray_intersections, select_gausslet_intersections#detect_segment, detect_gausslet
from raypier.core.cfaces import RectangularFace

from raypier.utils import normaliseVector

try:
    from raypier.point_spread_func import PSF
except ImportError:
    PSF=None


class BaseCapturePlane(Probe):
    name = Str("Ray Capture Plane")
    width = Float(30.0)
    height = Float(20.0)
    
    face_list = Instance(FaceList)
    
    plane_src = Instance(tvtk.PlaneSource, ())
    
    captured = List(Instance(RayCollection))
    
    vtkproperty = Instance(tvtk.Property, (), {'opacity': 1.0, 'color': (0.1,0.1,0.1),
                                               'representation': 'wireframe'})
    _mtime = Float(0.0)
    
    ### A list of lists of rays. One list per source.
    _all_rays = List()
    
    
    traits_view = View(VGroup(
                       Traceable.uigroup,
                       Item('width', editor=NumEditor),
                       Item('height', editor=NumEditor),
                        ),
                   )
    
    @observe("centre, orientation, width, height")
    def config_pipeline(self, evt):
        hwidth = self.width/2.
        hheight = self.height/2.
        plane = self.plane_src
        plane.origin = (-hwidth, -hheight,0.0)
        plane.point1 = (hwidth, -hheight,0.0)
        plane.point2 = (-hwidth, +hheight,0.0)
        plane.modified()
        
        face = self.face_list.faces[0]
        face.length=self.height
        face.width = self.width
        
        self.update = True
    
    def _actors_default(self):
        plane = self.plane_src
        self.config_pipeline(None)
        
        trans = tvtk.TransformFilter(input_connection = plane.output_port,
                                     transform=self.transform)
        
        map = tvtk.PolyDataMapper(input_connection=trans.output_port)
        map.scalar_visibility = False
        act = tvtk.Actor(mapper=map)
        act.property = self.vtkproperty
        actors = tvtk.ActorCollection()
        actors.append(act)
        return actors
    
    def _face_list_default(self):
        face = RectangularFace(length=self.height, width=self.width, offset=0.0, z_plane=0.0)
        fl = FaceList(owner=self)
        fl.faces = [face,]
        return fl
    
    def evaluate(self, src_list):
        mtime = self._mtime
        if any(src._mtime > mtime for src in src_list):
            all_rays = []
            for src in src_list:
                all_rays.append(src.traced_rays)
            self._all_rays = all_rays
    
    @observe("centre, orientation, width, height, _all_rays")
    def on_input_changed(self, evt):
        self._evaluate()
    
    def _evaluate(self):
        raise NotImplementedError()
        
        
class RayCapturePlane(BaseCapturePlane):
    def _evaluate(self):
        all_rays = self._all_rays
        if not all_rays:
            return
        fl = self.face_list
        fl.sync_transforms()
        captured = [select_ray_intersections(fl, list(rc_list)) for rc_list in all_rays]
        self._mtime = time.monotonic()
        self.captured = captured
        
        
class GaussletCapturePlane(BaseCapturePlane):
    name = Str("Gausslet Capture Plane")
    captured = List(Instance(GaussletCollection))
    
    def _evaluate(self):
        all_rays = self._all_rays
        if not all_rays:
            return
        fl = self.face_list
        fl.sync_transforms()
        captured = [select_gausslet_intersections(fl, list(gc_list)) for gc_list in all_rays]
        self._mtime = time.monotonic()
        self.captured = captured
            

class PolarisationProbe(Probe):
    name = "Polarisation Probe"
    size = Float(25.0)
    
    update = Event() #request re-tracing
    
    plane_src = Instance(tvtk.PlaneSource, (), 
                    {"x_resolution":1, "y_resolution":1},
                    transient=True)
                    
    traits_view = View(VGroup(
                       Traceable.uigroup,
                       Item('size', editor=NumEditor),
                        ),
                   )
    
    @on_trait_change("size")
    def config_pipeline(self):
        src = self.plane_src
        size = self.size/2.
        
        src.origin = (-size,-size,0)
        src.point1 = (-size,size,0)
        src.point2 = (size,-size,0)
        
        self.update = True
        
    @on_trait_change("size, centre, direction")
    def on_change(self):
        self.update = True
    
    def _actors_default(self):
        source = self.plane_src
        trans_f = tvtk.TransformFilter(input_connection=source.output_port,
                        transform = self.transform)
        map = tvtk.PolyDataMapper(input_connection=trans_f.output_port)
        act = tvtk.Actor(mapper=map)
        actors = tvtk.ActorCollection()
        actors.append(act)
        return actors
        
    def get_actors(self, scene):
        return self.actors
        
    def find_intersections(self, ray_src):
        for rays in ray_src.TracedDetailRays:
            self.intersect_rays(rays)
            
    def intersect_rays(self, rays):
        """
        @param rays: a RayCollection instance
        """
        #raise NotImplementedError
        pass
      
if PSF is not None:
    class PointSpreadFunction(Probe):
        position = Vector(desc="centre point of the point grid on which the PSF is evaluated")
        direction = Vector(desc="normal vector of the grid plane")
        orientation = Vector(desc="x_axis direction of the grid plane. This is projected onto\
        the grid plane to get the actual axis direction")
        
        point_spacing = Float(0.05) #in mm
        size = Int(30)
        
        psf = Instance(PSF)
        
        wavelength = Float(100.0, desc="wavelength, in microns")
        
        ray_source = Instance(klass="raypier.sources.BaseRaySource")
        
        exit_pupil_offset = Float(100.0)
        aperture = Float(1.0)
        
        eval_btn = Button("calculate")
        
        def _psf_default(self):
            return PSF(owner=self) 
            
        def _eval_btn_fired(self):
            source = self.ray_source
            
            pass
            
    class FaceCenteredPSF(PointSpreadFunction):
        name = Property(Str, depends_on="target_face")
        position = Property(Vector, depends_on="target_face",
                            desc="centre point of the point grid on which the PSF is evaluated")
        direction = Property(Vector, depends_on="target_face",
                          desc="normal vector of the grid plane")
        orientation = Property(Vector, depends_on="target_face",
                               desc="x_axis direction of the grid plane. This is projected onto\
        the grid plane to get the actual axis direction")
        
        target_face = Instance(klass="raypier.faces.Face")
        
        target_points = Property(Array(shape=(None,3)),
                                  depends_on="position, direction, orientation, size, point_spacing")
        
        rays = Instance(klass="raypier.rays.RayCollection")
        
        source = Instance(tvtk.ProgrammableSource, ())
        
        traits_view = View(Item('eval_btn', show_label=False),
                           VGroup(
                               Item('point_spacing'),
                               Item('size'),
                               Item('exit_pupil_offset'),
                               Item('ray_source', editor=DropEditor()),
                               Item('target_face', editor=DropEditor())
                               ),
                            resizable=True)
        
        def _actors_default(self):
            source = self.source
            def execute():
                output = source.poly_data_output
                if self.rays is None:
                    return
                points = self.rays.origin
                cells = self.rays.cells
                output.points = points
                output.polys = cells.tolist()
                #print "PSF", cells, points
            source.set_execute_method(execute)
            
            map = tvtk.PolyDataMapper(input_connection=source.output_port)
            act = tvtk.Actor(mapper=map)
            act.property.representation="wireframe"
            actors = tvtk.ActorCollection()
            actors.append(act)
            return actors
        
        def _rays_changed(self):
            self.source.modified()
            self.render = True
        
        def _get_name(self):
            optic = self.target_face.owner
            idx = optic.faces.index(self.target_face)
            return "PSF on %s, face %d"%(optic.name, idx)
                           
        def _get_position(self):
            return numpy.asarray(self.target_face.owner.centre)
        
        def _get_direction(self):
            return -numpy.asarray(self.target_face.owner.direction)
        
        def _get_orientation(self):
            return numpy.asarray(self.target_face.owner.x_axis)
        
        def _get_target_points(self):
            length = (self.point_spacing * self.size)/2.
            size = self.size * 1j
            U,V = numpy.ogrid[-length:length:size, -length:length:size]
            
            axis1 = normaliseVector(numpy.cross(self.direction, self.orientation))
            axis2 = normaliseVector(numpy.cross(axis1, self.direction))
            
            newaxis = numpy.newaxis
            points = axis2[newaxis, newaxis,:] * U[:,:,newaxis] \
                    + axis1[newaxis, newaxis,:] * V[:,:,newaxis]
                    
            points += self.position
            return points
            
        
        def _eval_btn_fired(self):
            source = self.ray_source
            
            sequence = source.get_sequence_to_face(self.target_face)
            input_rays = source.InputDetailRays
            
            result = self.psf.do_trace(input_rays, sequence)
            
            self.rays = result
            
            target_points = self.target_points
            
            spot = self.psf.evaluate_scalar_amp(result, target_points)
            
            self.spot = spot
