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

def compute_length(start, end):
    #takes [n,3] vectors and returns [1,n] array of distances
    a = start - end
    d = a**2
    e = sum(d.T)
    distance = numpy.sqrt(e)
    
    mask = distance < .0001
    distance[mask]=numpy.Infinity
    return distance 
    
def interpolate_z(first, threeD, twoD):
        # takes an 3d origin (x,y,z), a 3d point on the line, and intrpolates the third
        # dimesnsion of another point on that line, for which x and y are given 
        # and z is 0.
        # pass numpy.fliplr(P1[:,1:]).copy() to interpolate y instead,
        # and then flip the result back
        
        # len2d1/len2d2 = LenInZ1/LenInZ2
    
    len2d1 = compute_length(first[:,:2], threeD[:,:2])
    len2d2 = compute_length(first[:,:2], twoD[:,:2])
       
    k = len2d2 / len2d1
        
    #Zf = Z1 - k*ChangeinZ         
    z = first[:,2] - (first[:,2]-threeD[:,2])*k
        
    return  z