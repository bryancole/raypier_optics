#!/usr/bin/python
"""
An example of setting up a compound parabolic trough from parabolic troughs.
"""
from raytrace.test_sources import PlaneRaySource
from raytrace.tracer import RayTraceModel
from raytrace.troughs import TroughParabloid
from raytrace.constraints import BaseConstraint
from enthought.traits.api import Instance
from raytrace.has_queue import on_trait_change
from raytrace.absorbers import RectAbsorber

from raytrace.more_utils import transform_pts

import numpy


theta = 10      #acceptance angle in degrees
extreem_slope = numpy.tan(numpy.radians(theta))

len = 100
efl = 5
m1_centre = tuple([0,0,0])
right_tilt = tuple([extreem_slope, 1, 0])
left_tilt = tuple([-extreem_slope, 1, 0])

big_bounds = tuple([-100,100])

#m1 and m2 will be parabolas tilted towards each other.
# the focus of each parabola will lay on the surface of the other.
# the minimum of abs(x_bound) will be defined by these focii.
# the maximum of abs(x_bound) will be defined by the extreem ray that hits 
#    the focus (ie: must hit focus and not hit the back of the parabola first)

#to achieve this, first define m1 (right side) with extra big x_bounds so test rays 
# are sure to hit it.

m1 = TroughParabloid(X_bounds = big_bounds,
                      length = len,
                      EFL = efl,
                      centre = m1_centre,
                      direction = left_tilt)
                      
#set up a ray to find where the surface of m1 intersects the aperature (x=0)
#note: points for intersect method must be 2d, so use double square brakets
P1 = numpy.array([[0,0,-10]])
P2 = numpy.array([[60,0,-10]])
   
#where does m1's surface intersect apperature?
P = m1.intersect(P1, P2,100)
P['point'][0][2] = 0 #points had neg z inorder to hit, but center has z=0

#convert to tuple as m2's center
m2_centre = tuple([P['point'][0][0],P['point'][0][1],P['point'][0][2]])

#now set x bounds:
#lower bound is where m2 touches m1's focus
m2_xmin = m1_centre[0]

#upper bound is where the extreem ray that hits m2's focus would hit the back of m2 
#give ray a negative z value so that it is within the trough
extreem_ray_begin = numpy.array([[-100*extreem_slope, 100,-10]])

Q = m2_centre   #tuple
temp = numpy.zeros([1,3])   #numpy.array
temp[0] = Q[0], Q[1], Q[2]-10
temp2 = numpy.zeros([1,3])   #numpy.array
temp2[0] = Q[0]+extreem_ray_begin[0][0], Q[1]+extreem_ray_begin[0][1], Q[2]+extreem_ray_begin[0][2]

m2_temp = TroughParabloid(X_bounds = big_bounds,
                      length = len,
                      EFL = efl,
                      centre = m2_centre,
                      direction = right_tilt)
                      
pee = m2_temp.intersect(temp, temp2,400)

print "m2: ",pee['point'][0]
blah = transform_pts(m2_centre,right_tilt,numpy.array([[m2_xmin,0,0],pee['point'][0]]))

m2_bounds = tuple([blah[0][0],blah[1][0]])
print "m2_bounds: ", m2_bounds
#print "temp, tempt2, pee: ", temp, temp2, pee
#print "m2 actual center: ",m2_centre
       
#now that m2 is set up, we can use it to define the bounds on m1
m1_xmin = m2_centre[0]

extreem_ray_begin = numpy.array([[100*extreem_slope,100,0]])
#conversion from tuple to numpy.array needed
Q = m1_centre   #tuple
temp = numpy.zeros([1,3])   #numpy.array
temp[0] = Q[0], Q[1], Q[2]
temp2 = numpy.zeros([1,3])   #numpy.array
temp2[0] = Q[0]+extreem_ray_begin[0][0], Q[1]+extreem_ray_begin[0][1], Q[2]+extreem_ray_begin[0][2]
Pm1 = m1.intersect(temp, temp2,100)

print "???: ",temp[0],temp2[0],Pm1['point']
blah = transform_pts(m1_centre,left_tilt,numpy.array([[m1_xmin,0,0],Pm1['point'][0]]))

m1_bounds = tuple([blah[0][0],blah[1][0]])

"""
s1 = ParallelRaySource(origin=(0,130,-40),
                            direction=(-.17,-1,-.5),
                            working_dist = 100.,
                            number=20,
                            radius=3,
                            display = "wires")
                            
s2 = ParallelRaySource(origin=(0,130,-60),
                            direction=(.17,-1,-.2),
                            working_dist = 100.,
                            number=20,
                            radius=4,
                            display = "wires")
"""                            
s3 = PlaneRaySource(origin=(0,100,0),
                            direction=(0,-1,0),
                            X_width = 100.,
                            X_resolution = 10.,
                            Y_width = 10.,
                            Y_resolution = 2.,
                            display = "wires") 

s4 = PlaneRaySource(origin=(0,100,0),
                            direction=(.17,-1,0),
                            X_width = 100.,
                            X_resolution = 10.,
                            Y_width = 10.,
                            Y_resolution = 2.,
                            display = "wires") 

right_side = TroughParabloid(X_bounds = m1_bounds,
                      length = len,
                      EFL = efl,
                      centre = m1_centre,
                      direction = left_tilt)

left_side = TroughParabloid(X_bounds = m2_bounds,
                      length = len,
                      EFL = efl,
                      centre = m2_centre,
                      direction = right_tilt)

sorber = RectAbsorber(width = m1_xmin, length=len, direction = (0,1,0),
                        centre = (m1_xmin/2,0,0))
                        
print "M1 bounds: ",m1_bounds
print "M1 center: ",m1_centre
print "M2 bounds: ",m2_bounds
print "M2 center: ",m2_centre

model = RayTraceModel(optics=[left_side, right_side],
                        sources= [s3,s4])



#s1 = ParallelRaySource(origin=(0,30,30),
#                            direction=(0,-1,0),
#                            working_dist = 50.,
#                            #number=20,
#                            radius=14.)
                            
                            
#print source.InputDetailRays.origin.shape
                

                
#model = RayTraceModel(optics=[m1,],
#                    sources=[s1,],
#                    probes=[])
 
#model.trace_detail_async()
#import time
#start = time.clock()
#model.trace_all()
#end = time.clock()
#print "traced in", end - start                   

model.configure_traits()
