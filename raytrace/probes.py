#    Copyright 2009, Bryan Cole, Teraview Ltd.
#
#    This file is part of Raytrace.
#
#    Raytrace is free software: you can redistribute it and/or modify
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

from enthought.traits.api import on_trait_change, Float, Instance,Event, Int,\
        Property, Button, Str, Array

from enthought.traits.ui.api import View, Item, VGroup, DropEditor

from enthought.tvtk.api import tvtk
import numpy

from raytrace.bases import Probe, Traceable, NumEditor, Vector
from raytrace.sources import RayCollection
from raytrace.rays import collectRays
from raytrace.utils import normaliseVector

try:
    from raytrace.point_spread_func import PSF
except ImportError:
    PSF=None


class PolarisationProbe(Probe):
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
        trans_f = tvtk.TransformFilter(input=source.output,
                        transform = self.transform)
        map = tvtk.PolyDataMapper(input=trans_f.output)
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
        raise NotImplementedError
      
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
        
        ray_source = Instance(klass="raytrace.sources.BaseRaySource")
        
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
        
        target_face = Instance(klass="raytrace.faces.Face")
        
        target_points = Property(Array(shape=(None,3)),
                                  depends_on="position, direction, orientation, size, point_spacing")
        
        rays = Instance(klass="raytrace.rays.RayCollection")
        
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
            
            map = tvtk.PolyDataMapper(input=source.output)
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
