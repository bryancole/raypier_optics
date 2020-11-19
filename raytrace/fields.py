"""
Functions relating to the evaluation of the optical E-field by
summation of General Astigmatic Gaussian Beams
"""


def project_to_sphere(rays, centre=(0,0,0), radius=10.0):
    """project the given set of rays back to their intercept 
    with a sphere at the given centre and radius.
    Returns a new RayCollection"""
    
    
def evaluate_modes(rays, neighbours):
    """For each ray in rays, use its nearest neighbours
    to compute the best fit for the Astigmatic Gaussian Beam
    parameters.
    For N rays, return a Nx3 complex array of coeffs"""
    ### Do linear least squares on each ray and neighbours
    
    
def evaluate__E(rays, mode_coeffs, points):
    """Evaluate the E-field at each position in the
    _points_ array.
    """