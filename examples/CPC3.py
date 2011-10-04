#!/usr/bin/python
"""
An example of setting up a compound parabolic trough from parabolic troughs.
"""
from raytrace.splines import Extruded_bezier
from raytrace.beamstop import RectTarget
from raytrace.sources import RectRaySource
from raytrace.tracer import RayTraceModel
from raytrace.cmaterials import OpaqueMaterial, PECMaterial,DielectricMaterial
from raytrace.results import RayPaths
import numpy as np

#from raytrace.more_utils import transform_pts

import numpy as np

system_length = 100     #trough length

theta_p = np.radians(10)      #acceptance half angle (degrees becomes radians)
theta_n = np.radians(10)
len_exit_ap = 4.
exit_ap =[[-len_exit_ap/2.,0],[len_exit_ap/2,0]]
#start with an untilted parabola (directrix parallel to x axis), phi = 0 is -y
#in polar: r = 2a/(1+cos(phi))
#determine a, knowing that at phi = 90-theta+, r = width of entry aperature

lea = len_exit_ap
theta_D = np.pi/2. - theta_p
a = (lea*(1+np.cos(theta_D)))/2.

#next, determine the point where the other extreem ray from the orgin hits the parabola
#this is where phi = 180 - (theta+ + theta-)
#determine r there. call this point E

theta_E = np.pi - theta_p - theta_n
rE = (2.*a)/(1+np.cos(theta_E))

#determine the (x,y) of E

E = np.array([rE*np.sin(theta_E),-rE*np.cos(theta_E)])

#call the edge of the entry aperature D, and determine its (x,y)

D = np.array([lea*np.sin(theta_D),-lea*np.cos(theta_D)])

#determine quadratic control points:
#http://alecmce.com/as3/parabolas-and-quadratic-bezier-curves

#directrix is parallel to x axis, so line perp to directrix is paralell to y
A = np.array([D[0],-2*a])
B = np.array([E[0],-2*a])

#C is at origin, so halfway points are easy to calculate
H = A/2.
J = B/2.

#find the intersection of DH and EJ, which is control point K

try:
 m1 = (D[1]-H[1])/(D[0]-H[0])
 m2 = (E[1]-J[1])/(E[0]-J[0])
except ZeroDivisionError:
 raise ZeroDivisionError("Zero divison: this is not a real CPC")

b1 = D[1] - m1*D[0]
b2 = E[1] - m2*E[0]

Kx = (b2-b1)/(m1-m2)
K = np.array([Kx,m1*Kx+b1]) 
#D,K and E are control points of our quadratic, but we have a cubic bezier to work with.
#the formula for control points c1 and c2 of cubic bezier curve that is equal to quadratic bezier curve with control point K is:
#
c11 = 2./3.*K + 1./3.*D
c12 = 2./3.*K + 1./3.*E

#now rotate these points around the origin by theta+, and translate to the left by 1/2 the entry aperatures length (for symetry)


rot_matrix = np.array([[np.cos(theta_p),-np.sin(theta_p)],[np.sin(theta_p),np.cos(theta_p)]])
trans_vector = np.array([lea/2.,0])

D1 = np.dot(rot_matrix,D)-trans_vector
c11 = np.dot(rot_matrix,c11)-trans_vector
c12 = np.dot(rot_matrix,c12)-trans_vector
E1 = np.dot(rot_matrix,E)-trans_vector

#D should now be at (.5ea,0)

#print "D:",D
#print "c1:",c1
#print "c2:",c2
#print "E:",E
right_arm = Extruded_bezier(name = "right arm",control_points = np.array([[D1,c11,c12,E1]]),z_height_1 = 0, z_height_2 = system_length, material = PECMaterial())

#do it for the other arm
theta_D = np.pi/2. - theta_n
a = (lea*(1+np.cos(theta_D)))/2.

theta_E = np.pi - theta_p - theta_n
rE = (2.*a)/(1+np.cos(theta_E))

E2 = np.array([-rE*np.sin(theta_E),-rE*np.cos(theta_E)])

D2 = np.array([-lea*np.sin(theta_D),-lea*np.cos(theta_D)])

#determine quadratic control points:
#http://alecmce.com/as3/parabolas-and-quadratic-bezier-curves

#directrix is parallel to x axis, so line perp to directrix is paralell to y
A = np.array([D2[0],-2*a])
B = np.array([E2[0],-2*a])

#C is at origin, so halfway points are easy to calculate
H = A/2.
J = B/2.

#find the intersection of DH and EJ, which is control point K

try:
 m1 = (D2[1]-H[1])/(D2[0]-H[0])
 m2 = (E2[1]-J[1])/(E2[0]-J[0])
except ZeroDivisionError:
 raise ZeroDivisionError("Zero divison: this is not a real CPC")

b1 = D2[1] - m1*D2[0]
b2 = E2[1] - m2*E2[0]

Kx = (b2-b1)/(m1-m2)
K = np.array([Kx,m1*Kx+b1]) 
#D,K and E are control points of our quadratic, but we have a cubic bezier to work with.
#the formula for control points c1 and c2 of cubic bezier curve that is equal to quadratic bezier curve with control point K is:
#
c21 = 2./3.*K + 1./3.*D2
c22 = 2./3.*K + 1./3.*E2

#now rotate these points around the origin by - theta_n and translate to the right by 1/2 the entry aperatures length (for symetry)


rot_matrix = np.array([[np.cos(-theta_n),-np.sin(-theta_n)],[np.sin(-theta_n),np.cos(-theta_n)]])
trans_vector = np.array([-lea/2.,0])

D2 = np.dot(rot_matrix,D2)-trans_vector
c21 = np.dot(rot_matrix,c21)-trans_vector
c22 = np.dot(rot_matrix,c22)-trans_vector
E2 = np.dot(rot_matrix,E2)-trans_vector

left_arm = Extruded_bezier(name = "left arm",control_points = np.array([[D2,c21,c22,E2]]),z_height_1 = 0, z_height_2 = system_length, material = PECMaterial())



target = RectTarget(width = len_exit_ap, length = system_length, elevation=90., rotation =90., centre = (0,0,system_length/2.))

source = RectRaySource(origin=(0,100,system_length/4.),
                            direction=(0,-1,0),
                            working_dist = 200.,
                            number=10,
                            length = system_length/4.,
            			    width = 10,   #lea+2*E[0],
            			    randomness = True,
                            theta = 10.00000001)

raypaths = RayPaths()

model = RayTraceModel(sources=[source],
                      optics=[target, right_arm,left_arm], results = [raypaths], recursion_limit=10)


model.configure_traits()


