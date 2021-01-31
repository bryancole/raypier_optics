
import numpy


def normaliseVector(a):
    """normalise a (3,) vector or a (n,3) array of vectors"""
    a= numpy.asarray(a)
    mag = numpy.sqrt((a**2).sum(axis=-1))[...,numpy.newaxis]
    return a/mag

def dotprod(a,b):
    """dot-product along last axis"""
    return (a*b).sum(axis=-1)[...,numpy.newaxis]

def Convert_to_SP(input_v, normal_v, E1_vector, E1_amp, E2_amp):
    """
    All inputs are 2D arrays
    """
    cross = numpy.cross
    
    E2_vector = cross(input_v, E1_vector)
    
    v = cross(input_v, normal_v)
    S_vector = numpy.where(numpy.all(v==0, axis=1).reshape(-1,1)*numpy.ones(3), 
                           normaliseVector(E1_vector),
                           normaliseVector(v) )
    
    v = cross(input_v, S_vector)
    P_vector = normaliseVector(v)
    
    S_amp = E1_amp*dotprod(E1_vector,S_vector) + E2_amp*dotprod(E2_vector, S_vector)
                
    P_amp = E1_amp*dotprod(E1_vector,P_vector) + E2_amp*dotprod(E2_vector, P_vector)
    
    return S_amp, P_amp, S_vector, P_vector

def rotation(theta):
    #this is maybe deprecated. added z_rotation() because this was rotating in all three dimensions
    #by the same angle.  why would one want to do that?
    tx=ty=tz = theta
    cos = numpy.cos
    sin = numpy.sin
    Rx = numpy.array([[1,0,0], [0, cos(tx), -sin(tx)], [0, sin(tx), cos(tx)]])
    Ry = numpy.array([[cos(ty), 0, -sin(ty)], [0, 1, 0], [sin(ty), 0, cos(ty)]])
    Rz = numpy.array([[cos(tz), -sin(tz), 0], [sin(tz), cos(tz), 0], [0,0,1]])

    return numpy.dot(Rx, numpy.dot(Ry, Rz))

def z_rotation(theta):
    tz = theta
    cos = numpy.cos
    sin = numpy.sin
    Rz = numpy.array([[cos(tz), -sin(tz), 0], [sin(tz), cos(tz), 0], [0,0,1]])
    return Rz

