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

from enthought.traits.api import on_trait_change, Float, Instance,Event

from enthought.traits.ui.api import View, Item, VGroup

from enthought.tvtk.api import tvtk
import numpy

from raytrace.tracer import Probe, Traceable, NumEditor
from raytrace.sources import RayCollection

class PolarisationProbe(Probe):
    size = Float(25.0)
    
    update = Event() #request re-tracing
    render = Event() #request rerendering
    
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
        