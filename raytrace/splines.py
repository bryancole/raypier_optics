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


from enthought.traits.api import HasTraits, Array, Float, Complex,\
            Property, List, Instance, on_trait_change, Range, Any,\
            Tuple, Event, cached_property, Set, Int, Trait, Bool
from enthought.traits.ui.api import View, Item, ListEditor, VSplit,\
            RangeEditor, ScrubberEditor, HSplit, VGroup, Heading
from enthought.tvtk.api import tvtk
from enthought.tvtk.pyface.scene_model import SceneModel
from enthought.tvtk.pyface.scene_editor import SceneEditor
import numpy
from itertools import chain, izip, islice, tee

#from raytrace.tracer import Optic, VTKOptic, normaliseVector, RaySegment,\
#             Traceable, NumEditor, dotprod, transformPoints, transformNormals

from raytrace.bases import Optic, Traceable
from raytrace.cfaces import PolygonFace, ExtrudedBezier2Face, ExtrudedPlanarFace
from raytrace.ctracer import FaceList

import numpy as np
from scipy.interpolate import splprep, splev

def b_spline_to_bezier_series(tck, per = False):
  """Convert a parametric b-spline into a sequence of Bezier curves of the same degree.
 
  Inputs:
    tck : (t,c,k) tuple of b-spline knots, coefficients, and degree returned by splprep.
    per : if tck was created as a periodic spline, per *must* be true, else per *must* be false.
 
  Output:
    A list of Bezier curves of degree k that is equivalent to the input spline.
    Each Bezier curve is an array of shape (k+1,d) where d is the dimension of the
    space; thus the curve includes the starting point, the k-1 internal control
    points, and the endpoint, where each point is of d dimensions.


    from http://mail.scipy.org/pipermail/scipy-dev/2007-February/006651.html
    (actually, I got it from http://old.nabble.com/bezier-curve-through-set-of-2D-points-td27158642.html
    
  """
  from scipy.interpolate.fitpack import insert
  from numpy import asarray, unique, split, sum, transpose
  t,c,k = tck
  t = asarray(t)
  try:
    c[0][0]
  except:
    # I can't figure out a simple way to convert nonparametric splines to
    # parametric splines. Oh well.
    raise TypeError("Only parametric b-splines are supported.")
  new_tck = tck
  if per:
    # ignore the leading and trailing k knots that exist to enforce periodicity
    knots_to_consider = unique(t[k:-k])
  else:
    # the first and last k+1 knots are identical in the non-periodic case, so
    # no need to consider them when increasing the knot multiplicities below
    knots_to_consider = unique(t[k+1:-k-1])
  # For each unique knot, bring it's multiplicity up to the next multiple of k+1
  # This removes all continuity constraints between each of the original knots,
  # creating a set of independent Bezier curves.
  desired_multiplicity = k+1
  for x in knots_to_consider:
    current_multiplicity = sum(t == x)
    remainder = current_multiplicity%desired_multiplicity
    if remainder != 0:
      # add enough knots to bring the current multiplicity up to the desired multiplicity
      number_to_insert = desired_multiplicity - remainder
      new_tck = insert(x, new_tck, number_to_insert, per)
  tt,cc,kk = new_tck
  # strip off the last k+1 knots, as they are redundant after knot insertion
  bezier_points = transpose(cc)[:-desired_multiplicity]
  if per:
    # again, ignore the leading and trailing k knots
    bezier_points = bezier_points[k:-k]
  # group the points into the desired bezier curves
  return split(bezier_points, len(bezier_points) / desired_multiplicity, axis = 0)


class Extruded_interpolant(Optic):
    """a generalized extruded suface of a bezier spline approximating a surface"""
    
    #The profile is best with way too many points in it, because the spline
    #will determine the fewest number of ray tracing faces needed, but
    #but the end caps and visualization are polygon approximations
    #which are not speed critical use the profile directly.
    #
    #I haven't experiemnted to see what is excessive
    profile = Array(shape=(2,None), dtype=numpy.double)
    tck, uout = [None,None]      #the output of splprep 

    z_height_1 = Float(-30.0)   #must be smaller than z_height_2
    z_height_2 = Float(30.0)
    smoothness = Float(.0005)  #for splprep. found this number by trial and error. there are algorithms to guess better
    
    trace_ends = Bool(False, desc="include the end-faces in tracing")
    trace_top = Bool(True, desc="include the end-faces in tracing")

    data_source = Instance(tvtk.ProgrammableSource, ())
    
    extrude = Instance(tvtk.LinearExtrusionFilter, (), 
                       {'capping': False, 
                        'extrusion_type':'vector',
                        'vector': (0.,0.,1.)})
    
    def _faces_default(self):
        fl = FaceList(owner=self)
        fl.faces = self.make_faces()
        return fl
    
    def make_faces(self):
        z1 = self.z_height_1
        z2 = self.z_height_2
        profile = self.profile
        m = self.material
        
        #convert profile to a list of bezier curves defined by three 2D points (knots)
        tck, uout = splprep(profile, s=self.smoothness, k=2, per=False)
        self.tck, self.uout = [tck,uout]
        p2 = b_spline_to_bezier_series(tck)
        print "splprep used ",len(p2), " faces to make this spline"
        if len(p2) > 30:
            print "!!! thats alot of faces.  try adjusting smoothness. /nSpline degree 2 may be too small  !!!""" 
            tck, uout = splprep(profile, s=self.smoothness, k=3, per=False)
            print "3rd degree spline would use ",len(tck[0])-7," faces."
        curves = []
        for knots in p2:
            curves.append(ExtrudedBezier2Face(owner=self, X0=knots[0][0], Y0=knots[0][1], \
                                X1=knots[1][0], Y1=knots[1][1], X2=knots[2][0], Y2=knots[2][1], \
                                z_height_1 = self.z_height_1, z_height_2 = self.z_height_2, material=m))
                                
        
        if self.trace_ends:
            print "traced ends"
            """not perfect, just a polygon of the original profile"""
            prof =  np.column_stack((profile[:][0],profile[:][1]))
            base = PolygonFace(owner=self, z_plane=z1,
                        xy_points=prof, material=m)
            top = PolygonFace(owner=self, z_plane=z2, material=m,
                        xy_points=prof, invert_normal=True)
            curves.extend([base, top])

        if self.trace_top:
            curves.append( ExtrudedPlanarFace(owner=self, z1=z1, z2=z2, x1=profile[0][0], y1=profile[1][0], 
                    x2=profile[0][-1], y2=[1][-1], material=m) )

        return curves


            
    
    @on_trait_change("z_height_1, z_height_2")
    def config_pipeline(self):
        self.extrude.scale_factor = self.z_height_2 - self.z_height_1
        print "extrude factor: ",  self.extrude.scale_factor
        self.faces.faces = self.make_faces()
        self.update=True
        
    def _trace_ends_changed(self):
        self.faces.faces = self.make_faces()
        self.update=True
        
    def _smoothness_changed(self):
        self.faces.faces = self.make_faces()
        self.update=True
        
    def _profile_changed(self):
        self.data_source.modified()
        self.faces.faces = self.make_faces()
        self.update=True

    def get_real_profile(self):
        #return the evaluated profile of the spline
        x = np.linspace(0,1,100)    #splprep is defined of 0,1 of the parameter.  100 points
        profile = splev(x,self.tck)

        return profile

    def _pipeline_default(self):
        source = self.data_source
        def execute():
            profile = self.get_real_profile()
            xy = np.column_stack((profile[:][0],profile[:][1]))
            z = numpy.ones(xy.shape[0]) * self.z_height_1
            points = numpy.column_stack((xy,z))
            
            cells = [range(len(z)),]
            
            output = source.poly_data_output
            output.points = points
            output.polys = cells
        source.set_execute_method(execute)
        
        self.extrude.scale_factor = self.z_height_2 - self.z_height_1  #mm, put here because it wasn;t being initialized
        if self.trace_ends:
            print "drew ends"
            self.extrude.capping = True
        extrude = self.extrude
        extrude.input = source.output
        
        t = self.transform
        transf = tvtk.TransformFilter(input=extrude.output, transform=t)
        return transf

class test(Extruded_interpolant):
    name = "test"
    height = Float #distance from front face to apex
    width = Float #width of front face
    
    traits_view = View(VGroup(
                       Traceable.uigroup,
                       Item('trace_ends'),
                       Item('n_inside'),
                       Item('length'),
                       Item('height'),
                       Item('width'),
                       )
                       )

    @on_trait_change("height, width")
    def config_profile(self):
        h = self.height
        w = self.width/2
        self.profile = [(-w,0),
                        (w,0),
                        (0,h)]
            
'''
if __name__=="__main__":
    from ray_tracer import BuildRaySet, BeamStop, RayTraceModel
    
    input_rays = BuildRaySet(origin = (-7.42,15,0),
                             direction = (0,-1,0),
                             radius=1.0,
                             count=20)
    
    rhomboid = Rhomboid(n_inside = 1.764+0.0j,
                        orientation=0.0,
                        elevation=0.0,
                        centre=(0,0,0),
                        rotation=0,
                        length=10.0,
                        height=7.07,
                        width=14.0,
                        slant=45.0)
    
    beamstop = BeamStop(width=10,height=10,
                        centre=(7,-10,0),
                        direction=(0,1,0))
    
    #print "rhomboid", rhomboid.polydata
    print "beamstop", beamstop.polydata
    
    model = RayTraceModel(optics=[rhomboid, beamstop], rays=input_rays)
    model.configure_traits()
'''
