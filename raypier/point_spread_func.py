"""routines used in the evaluation of the point-spread-function probes"""
from traits.api import HasTraits, Instance, Property, Float

from raypier.bases import normaliseVector, dotprod
from raypier.rays import collectRays

from scipy.interpolate import interp2d
from scikits.delaunay import NNInterpolator, Triangulation
import numpy

class PSF(HasTraits):
    owner = Instance(klass="raypier.bases.Probe")
    
    sample_spacing = Property(Float, depends_on="owner")
    
    def do_trace(self, input_rays, face_sequence):
        rays = input_rays
        for faceList in face_sequence:
            if len(faceList) > 1:
                print("rays more than 1")
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
        rays.offset_length = height/aspect
        return
    
    def evaluate_scalar_amp(self, rays, target):
        newaxis = numpy.newaxis
        direction = rays.direction[newaxis,newaxis,...]
        origin = rays.origin[newaxis,newaxis,...]
        optical_path = rays.cum_length[newaxis,newaxis,:,0]
        target = target[:,:,newaxis,:]
        
        delta = target - origin
        z = dotprod(direction, delta)
        r_bar = delta - (z*direction)
        r = numpy.sqrt(dotprod(r_bar,r_bar))[...,0]
        z = z[...,0]
        
        wavelen = self.owner.wavelength / 1000.0 #convert to mm
        
        ray_areas = rays.ray_areas[newaxis,newaxis,:]
        w0 = numpy.sqrt(ray_areas/numpy.pi) #beam waists
        z_r = ray_areas / wavelen #Rayleigh length
        w = w0 * numpy.sqrt(1 + (z/z_r)**2)
        R = z * ( 1 + (z_r/z)**2 )
        
        q = 1./( 1/R - (1j*wavelen)/(numpy.pi * w *w) )
        
        ik = -2j*numpy.pi/wavelen
        
        E = numpy.exp(ik*optical_path + (2*ik*r**2)/q)
        
        E_total = E.sum(axis=2)
        return E_total
    
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
            
    