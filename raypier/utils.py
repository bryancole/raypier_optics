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

import numpy
from tvtk.api import tvtk
from traits.api import Range as _Range, Tuple as _Tuple,\
            BaseTuple
from traitsui.api import TupleEditor

### For backwards compatbility
from .core.utils import normaliseVector, dotprod, Convert_to_SP, rotation, z_rotation


class EditorTraits(object):
    def get_editor(self, *args, **kwds):
        e = super(EditorTraits, self).get_editor(*args, **kwds)
        editor_t = {'auto_set':False,
                    'enter_set':True}
        metadata = self._metadata
        if 'editor_traits' in metadata:
            editor_t.update(metadata['editor_traits'])
        e.set(**editor_t)
        return e

class Range(EditorTraits, _Range):
    pass
    
    
class Tuple(EditorTraits, _Tuple):
    pass


class UnitVectorTrait(EditorTraits, BaseTuple):
    def validate(self, object, name, value):
        value = super(UnitVectorTrait, self).validate(object,name,value)
        mag = numpy.sqrt(sum(a**2 for a in value))
        return tuple(float(a/mag) for a in value)


TupleVector = Tuple((0.,0.,0.), editor_traits={'cols':3,
                                             'labels':['x','y','z']})
                                             
UnitTupleVector = UnitVectorTrait((0.,0.,1.), editor_traits={'cols':3,
                                             'labels':['x','y','z']})



def transformPoints(t, pts):
    """Apply a vtkTransform to a numpy array of points"""
    out = tvtk.Points()
    t.transform_points(pts, out)
    return numpy.asarray(out)

def transformVectors(t, pts):
    """Apply a vtkTransform to a numpy array of vectors.
    I.e. ignore the translation component"""
    out = tvtk.DoubleArray(number_of_components=3)
    t.transform_vectors(pts, out)
    return numpy.asarray(out)

def transformNormals(t, pts):
    """Apply a vtkTransform to a numpy array of normals.
    I.e. ignore the translation component and normalise
    the resulting magnitudes"""
    out = tvtk.DoubleArray(number_of_components=3)
    t.transform_normals(pts, out)
    return numpy.asarray(out)

