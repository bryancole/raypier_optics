"""routines used in the evaluation of the point-spread-function probes"""
from enthought.traits.api import HasTraits, Instance, Property, Float

from raytrace.bases import normaliseVector, dotprod

from scipy.interpolate import interp2d
from scikits.delaunay import NNInterpolator, Triangulation
import numpy

class PSF(HasTraits):
    owner = Instance(klass="raytrace.bases.Probe")
    
    sample_spacing = Property(Float, depends_on="owner")
    
    def do_trace(self, input_rays, face_sequence):
        rays = input_rays
        for faceList in face_sequence:
            if len(faceList) > 1:
                print "rays more than 1"
                rayList = [f.trace_rays(rays) for f in faceList]
                rays = collectRays(*rayList)
            else:
                rays = faceList[0].trace_rays(rays)
                
        self.project_onto_exit_pupil(rays)
        return rays
    
    def project_onto_exit_pupil(self, rays):
        owner = self.owner
        #mask = rays.length != numpy.Infinity
        directions = normaliseVector(rays.direction)#[mask]
        E1 = rays.E1_amp
        E2 = rays.E2_amp
        amp = numpy.sqrt(E1*E1.conjugate() + E2*E2.conjugate()).real
        mean_dir = normaliseVector((directions*amp).mean(axis=0))
               
        exit_pupil_offset = owner.exit_pupil_offset
        aperture = exit_pupil_offset /(owner.aperture*2.)
        pupil_centre = owner.position - exit_pupil_offset * mean_dir
        
        #project onto pupil plane
        delta = rays.origin - pupil_centre
        height = dotprod(delta, mean_dir)
        aspect = dotprod(directions, mean_dir)
        new_origin = rays.origin - (height/aspect) * directions
        
        rays.origin = new_origin
        return
    
    def resample_exit_pupil(self):
        sample_spacing = self.sample_spacing #FIXME
        axes = numpy.array([[1,0,0],[0,1,0],[0,0,1]])
        proj = dotprod(axes,mean_dir).argmax()
        best_ax = axes[proj]
        axis1 = normaliseVector(numpy.cross(best_ax, mean_dir))
        axis2 = normaliseVector(numpy.cross(axis1, mean_dir))
        
        grid_size = int(aperture/sample_spacing)
        
        A,B = numpy.ogrid[-aperture:aperture:grid_size*1j,
                          -aperture:aperture:grid_size*1j]
        
        #convert points to UV coords
        delta = new_origin - pupil_centre
        U = dotprod(delta, axis1)
        V = dotprod(delta, axis2)
        
        cells = rays.cells
        for cell in cells:
            pass
    
    def sum_wavelets(rays):
        pass
            
    