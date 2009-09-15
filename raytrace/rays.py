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

from enthought.traits.api import HasTraits, Int, Float, Range,\
     Bool, Property, Array, Event, List, cached_property, Str,\
     Instance, Tuple, on_trait_change, Trait
     
from raytrace.utils import normaliseVector, dotprod

VectorArray = Array(shape=(None,3), dtype=numpy.double)

def area_of_triangles(p0,p1,p2):
    s1 = p1 - p0
    s2 = p1 - p2
    vp = numpy.cross(s1,s2)
    return 0.5 * numpy.sqrt((vp**2).sum(axis=-1))

def shift_cells(rayList):
    total = 0
    for rays in rayList:
        cells = rays.cells
        if cells is not None:
            yield rays.cells + total
        total += rays.number

def collectRays(*rayList):
    """Combines a sequence of RayCollections into a single larger one"""
    origin = numpy.vstack([r.origin for r in rayList])
    normals = numpy.vstack([r.normals for r in rayList])
    direction = numpy.vstack([r.direction for r in rayList])
    length = numpy.vstack([r.length for r in rayList])
    cum_length = numpy.vstack([r.cum_length for r in rayList])
    face = numpy.concatenate([r.face for r in rayList])
    end_face = numpy.concatenate([r.end_face for r in rayList])
    refractive_index = numpy.vstack([r.refractive_index for r in rayList])
    E_vector = numpy.vstack([r.E_vector for r in rayList])
    E1_amp = numpy.vstack([r.E1_amp for r in rayList])
    E2_amp = numpy.vstack([r.E2_amp for r in rayList])
    parent_ids = numpy.concatenate([r.parent_ids for r in rayList])
    
    #parent = rayList[0].parent
    max_length = max(r.max_length for r in rayList)
    
    cells = numpy.vstack(list(shift_cells(rayList)))
    
    newRays = RayCollection(origin=origin, direction=direction,
                            length=length, face=face,
                            refractive_index=refractive_index,
                            E_vector=E_vector, E1_amp=E1_amp,
                            E2_amp = E2_amp, parent_ids=parent_ids,
                            max_length=max_length,
                            cells = cells,
                            normals=normals)
    return newRays


class RayCollection(HasTraits):
    origin = VectorArray
    direction = Property(VectorArray)
    length = Array(shape=(None,1), dtype=numpy.double)
    offset_length = Array(shape=(None,1), dtype=numpy.double)
    cum_length = Property(Array(shape=(None,1), dtype=numpy.double),
                          depends_on="parent, parent_ids")
    termination = Property(depends_on="origin, direction, length")
    
    number = Property(Int, depends_on="origin")
    
    face = Array(shape=(None,), dtype=numpy.object,
                 desc="The face which _generates_ the ray")
    
    end_face = Array(shape=(None,), dtype=numpy.object,
                             desc="The face which terminates the ray (by intersection)")
    
    normals = VectorArray(desc="the normal vectors to the intersecting surface at the ray origin")
    
    refractive_index = Array(shape=(None,1), dtype=numpy.complex128)
    
    E_vector = VectorArray
    E1_amp = Array(shape=(None,1), dtype=numpy.complex128)
    E2_amp = Array(shape=(None,1), dtype=numpy.complex128)
    
    parent = Instance("RayCollection")
    parent_ids = Array(shape=(None,), dtype=numpy.uint32)
    
    max_length = Float
    
    cells = Array(shape=(None,None), dtype=numpy.uint32,
                  desc="ray topology (optional)")
    
    ray_areas = Property(Array(shape=(None)), depends_on="cells, origin")
    
    wavelength = Array(shape=(None,1), dtype=numpy.double)
    
    def __init__(self, *args, **kwds):
        super(RayCollection, self).__init__(*args, **kwds)
        if 'cells' not in kwds:
            self._evaluate_cells()
            
    def _get_cum_length(self):
        parent = self.parent
        if parent is None:
            cum_length = - self.offset_length
        else:
            pids = self.parent_ids
            cum_length = parent.cum_length[pids] + parent.length[pids] \
                         - self.offset_length
        return cum_length
            
    def _get_ray_areas(self):
        cells = self.cells
        if (cells==[[0]]).all():
            return None
        points = self.origin
        cl = [cells[:,i] for i in xrange(4)]
        p0, p1, p2, p3 = [points[c] for c in cl]
        
        cell_normal = normaliseVector(numpy.cross(p2-p0,p3-p1))
        cell_direction = normaliseVector(self.direction[cells].sum(axis=1))
        aspect = dotprod(cell_normal, cell_direction)
        
        area = area_of_triangles(p0, p1, p2)
        area += area_of_triangles(p2, p3, p0)
        area *= aspect[:,0]
        
        ray_areas = numpy.zeros(self.number, numpy.double)
        cell_count = numpy.zeros(self.number, numpy.int)
        
        for c in cl:
            ray_areas[c] += area
            cell_count[c] += 1
            
        mask = cell_count!=0
        ray_areas[mask] /= cell_count[mask]
        
        return ray_areas
    
    def _evaluate_cells(self):
        parent = self.parent
        if parent is None:
            return
        parent_cells = parent.cells
        if (parent_cells == [[0]]).all():
            return
        parent_ids = self.parent_ids
        forward_map = numpy.ones(parent.number, numpy.int) * -1
        forward_map[parent_ids] = numpy.arange(self.number)

        newcells = forward_map[parent_cells]
        #now remove invalid ones i.e. with id<0
        mask = (newcells>=0).all(axis=1)
        
        newcells = newcells[mask,:]
        self.cells = newcells
    
    def set_polarisation(self, Ex, Ey, Ez):
        E = numpy.array([Ex, Ey, Ez])
        
        orth = numpy.cross(E, self.direction)
        para = numpy.cross(orth, self.direction)
        para /= numpy.sqrt((para**2).sum(axis=1)).reshape(-1,1)
        
        self.E_vector = para
    
    def _length_default(self):
        return numpy.ones((self.origin.shape[0],1)) * numpy.Infinity
    
    def _get_termination(self):
        length = self.length.clip(0, self.max_length).reshape(-1,1)
        try:
            return self.origin + length*self.direction
        except:
            pass
    
    def _get_direction(self):
        return self._direction
    
    def _set_direction(self, d):
        self._direction = normaliseVector(d)
        
    def _get_number(self):
        return self.origin.shape[0]
