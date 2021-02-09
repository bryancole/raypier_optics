"""
Utility module
"""

from numpy.linalg import solve
import numpy as np
from raypier.core.ctracer import RayCollection, GaussletCollection


def find_ray_focus(ray_collection):
    if isinstance(ray_collection, GaussletCollection):
        data = ray_collection.copy_as_array()
        ray = data['base_ray']
        E1 = ray['E1_amp']
        E2 = ray['E2_amp']
        origin = ray['origin']
        direction = ray['direction']
    else:
        if not isinstance(ray_collection, RayCollection):
            raise TypeError("Expecting a RayCollection object")
        E1 = ray_collection.E1_amp
        E2 = ray_collection.E2_amp
        origin = ray_collection.origin
        direction = ray_collection.direction
    weights = E1.real**2 + E1.imag**2 + E2.real**2 + E2.imag**2
    weights /= weights.sum()
    return find_focus(origin, direction, weights=weights)


def find_focus(ray_origins, ray_directions, weights=None):
    """
    givens Nx3 arrays for ray origin and directions, return the 
    3D point of closest intersection, found in the least-squares 
    sense.
    
    See https://math.stackexchange.com/q/1762491
    """
    
    o = ray_origins
    n = ray_directions / (ray_directions**2).sum(axis=1, keepdims=True)
    
    A = np.matmul(n[:,:,np.newaxis], n[:,np.newaxis,:]) - np.eye(3)[np.newaxis,:,:]
    # A.shape = (N,3,3), for N rays
    if weights is not None:
        A *= weights[:,None,None]
        
    AA = A.sum(axis=0) #AA is a (3,3) shape matrix
    
    b = np.matmul(A,o[:,:,np.newaxis]).sum(axis=0)
    
    out = solve(AA,b[:,0])
    return out


if __name__=="__main__":
    directions = np.array([[ 0.00000000e+00,  1.00000000e+00,  3.58332104e-17],
       [ 2.00168533e-01,  9.79761480e-01,  3.04549358e-17],
       [ 1.22567877e-17,  9.79761480e-01,  2.00168533e-01],
       [-2.00168533e-01,  9.79761480e-01,  5.49685112e-17],
       [-3.67703631e-17,  9.79761480e-01, -2.00168533e-01]])
    origins = np.array([[ 0.00000000e+00, -1.40000000e+01,  1.02258724e-15],
       [-9.83472949e+00, -1.51142336e+01,  8.82453980e-16],
       [-6.02203500e-16, -1.51142336e+01, -9.83472949e+00],
       [ 9.83472949e+00, -1.51142336e+01, -3.21953019e-16],
       [ 1.80661050e-15, -1.51142336e+01,  9.83472949e+00]])
    
    print(directions.shape, origins.shape)
    out = find_focus(origins, directions)
    print(out)
    
    ###Now check numerically
    def get_errs(y_coord):
        y = np.array([[0.0,1.0,0.0]])
        y_ = y_coord - origins[:,1] 
        det = np.dot(directions, y[0])
        a = y_/det
        intersection = origins + (a.reshape(-1,1)*directions)
        centroid = intersection.mean(axis=0)
        delta = ((intersection - centroid)**2).sum()
        return delta
    
    #for i in np.linspace(20,40,100):
    #    print i, get_errs(i)
        
    from matplotlib import pyplot as pp
    
    a=np.linspace(0.0,100.0,200)
    
    for i in range(5):
        x = origins[i,0] + a * directions[i,0]
        y = origins[i,1] + a * directions[i,1]
    
        pp.plot(y,x)
    pp.show()
    
    
