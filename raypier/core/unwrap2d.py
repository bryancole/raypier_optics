"""
Functions for performing 2D phase unwrapping.
"""

import numpy


def unwrap2d(phase_array, anchor=(0,0)):
    """
    A basic phase unwrapper using numpy.unwrap and the Itoh method.
    
    phase_array - a (N,M) shaped array of values in the range -pi .. +pi
    """
    ###dPhi/dx
    ddx = numpy.diff(phase_array, axis=0) #shape=(N-1,M)
    ddy = numpy.diff(phase_array, axis=1) #shape=(N,M-1)
    
    ### Calc path-integral of phase delta over each 2x2 pixel cell
    residuals = ddx[:,:-1] - ddx[:,1:] - ddy[:-1,:] + ddy[1:,:]
    
    int1 = numpy.unwrap(phase_array, axis=0)
    int2 = numpy.unwrap(int1, axis=1)
    
    a,b = anchor
    offset = int((int2[a,b] - phase_array[a,b])/(2*numpy.pi))
    int2 -= offset*numpy.pi*2
    return int2, residuals
    
    
    
if __name__=="__main__":
    x_ = numpy.linspace(-0.1,0.1,100)
    
    x,y = numpy.meshgrid(x_,x_)
    z = 1 - numpy.sqrt(1 - x**2 - y**2)
    z *= 5000.
    
    noise = numpy.random.standard_normal(z.shape)*(numpy.pi/10)
    z += noise
    
    z_w = (z + numpy.pi)%(2*numpy.pi) - numpy.pi
    
    
    from matplotlib import pyplot as pp
    
    pp.imshow(z_w)
    pp.show()
    
    residuals, unwrapped = unwrap2d(z_w, anchor=(50,50))
    
    pp.imshow(residuals)
    pp.show()
    
    pp.imshow(unwrapped)
    pp.show()
    