#    Copyright 2009, Teraview Ltd., Bryan Cole
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


import numpy
import itertools

from enthought.traits.api import HasTraits, Int, Float, \
     Bool, Property, Array, Event, List, cached_property, Str,\
     Instance, on_trait_change, Trait, Enum, Title

from enthought.traits.ui.api import View, Item, Tabbed, VGroup, Include, \
    Group

from enthought.tvtk.api import tvtk

from raytrace.rays import RayCollection
from raytrace.utils import normaliseVector, Range, TupleVector, Tuple
from raytrace.bases import Renderable
from raytrace.sources import BaseRaySource
from random import uniform

Vector = Array(shape=(3,))



class PlaneRaySource(BaseRaySource):
    """ A rectangular group of rays w/ a single direction."""
    origin = Tuple((0.,0.,0.))
    direction = Tuple((0.,0.,1.))
    X_width = Float(100.)
    X_number = Float(10.)
    Y_width = Float(30.)
    Y_number = Float(11.)
    randomness = Bool(True)
    
    view_ray_ids = numpy.arange(40)
    
    InputRays = Property(Instance(RayCollection), 
                         depends_on="origin, direction, resolution, radius, max_ray_len")
    
    geom_grp = VGroup(Group(Item('origin', show_label=False,resizable=True), 
                            show_border=True,
                            label="Origin position",
                            padding=0),
                       Group(Item('direction', show_label=False, resizable=True),
                            show_border=True,
                            label="Direction"),
                       Item('X_width'),
                       Item('X_number'),
                       Item('Y_width'),
                       Item('Y_number'),
                       label="Geometry")
    
    
    @on_trait_change("origin, direction, X_width, X_resolution, Y_width, Y_resolution, max_ray_len")
    def on_update(self):
        self.data_source.modified()
        self.update=True
    
    @cached_property
    def _get_InputRays(self):
        origin = numpy.array(self.origin)
        direction = numpy.array(self.direction)
        X_width = self.X_width
        X_num = self.X_number
        X_res = X_width/X_num
        Y_width = self.Y_width
        Y_num = self.Y_number
        Y_res = Y_width/Y_num
        randomness = self.randomness
        
        
        max_axis = numpy.abs(direction).argmax()
        if max_axis==0:
            v = numpy.array([0.,1.,0.])
        else:
            v = numpy.array([1.,0.,0.])
        d1 = numpy.cross(direction, v)
        d1 = normaliseVector(d1)
        d2 = numpy.cross(direction, d1)
        d2 = normaliseVector(d2)
        a = X_width/2
        b = Y_width/2
        X_range = numpy.linspace(-a, a, X_num)
        Y_range = numpy.linspace(-b, b, Y_num)
        offsets = [[]]*(X_range.size*Y_range.size)
        
        
        for i,x in enumerate(X_range):
            for j,y in enumerate(Y_range):
                if randomness:
                    x = x + uniform(-X_res/2, X_res/2)
                    y = y + uniform(-Y_res/2, Y_res/2)
                point = d1 * x + d2 * y
                block = i * Y_range.size
                offsets[block + j] = point
        

        
        origins = numpy.array([origin + offset for offset in offsets])
        directions = numpy.ones_like(origins) * direction
        rays = RayCollection(origin=origins, direction=directions,
                             max_length=self.max_ray_len)
        rays.set_polarisation(1, 0, 0)
        size = origins.shape[0],1
        rays.E1_amp = numpy.ones(size, dtype=numpy.complex128)
        rays.E2_amp = numpy.zeros(size, dtype=numpy.complex128)
        rays.refractive_index = numpy.ones(size, dtype=numpy.complex128)
        rays.offset_length = numpy.zeros_like(rays.length)
        rays.normals = numpy.zeros_like(origins)
        rays.normals[:,1]=1.0
        return rays
