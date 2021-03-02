"""
Functions relating to the creation and manipulations of gausslets
"""

import numpy
import time

from .core.utils import normaliseVector
from .core.ctracer import FaceList
from .bases import Traceable
from .core.cmaterials import ResampleGaussletMaterial
from .core.cfaces import CircularFace
from .core.fields import eval_Efield_from_gausslets
from .core.gausslets import make_hexagonal_grid, decompose_angle

from traits.api import Range, Float, Array, Property, Instance, observe, Str
from traitsui.api import View, Group, Item

from tvtk.api import tvtk
from raypier.editors import NumEditor, IntEditor
from raypier.core.gausslets import decompose_position




class BaseDecompositionPlane(Traceable):
    ### Sets the geometry of the intersection face
    diameter = Float(25.0)
    offset = Float(0.0)
    
    material = Instance(ResampleGaussletMaterial)
    
    def _material_default(self):
        m = ResampleGaussletMaterial(eval_func=self.evaluate_decomposed_rays)
        return m
    
    def _faces_default(self):
        fl = FaceList(owner=self)
        fl.faces =  [CircularFace(owner=self, z_plane=0.0, material=self.material)]
        return fl
        
    def _pipeline_default(self):
        radius = self.radius
        circle = tvtk.ArcSource(use_normal_and_angle=True,
                                center=[0.0,0.0,0.0],
                                polar_vector=[0.0, radius, 0.0],
                                normal=[0.0,0.0,1.0],
                                angle=360.0,
                                resolution=100
                                )
        
        trans = tvtk.TransformFilter(input_connection=circle.output_port,
                                     transform=self.transform)
        return trans
    
    def evaluate_decomposed_rays(self, input_rays):
        raise NotImplementedError()


class PositionDecompositionPlane(BaseDecompositionPlane):
    abstract=False
    name = Str("Position Decomposition")
    
    radius = Property()
    
    ### Wavefront curvature estimate, as a radius from the focus to decomposition plane centre
    curvature = Float(0.0)
    
    resolution = Float(10.0)
    
    blending = Float(1.2)
    
    ### Keep the sample points and sampled fields around 
    ### for debugging purposes
    _sample_points = Array()
    _unwrapped_phase = Array()
    _E_field = Array()
    
    traits_view = View(Group(
                    Traceable.uigroup,
                    Item("radius", editor=NumEditor),
                    Item("resolution", editor=NumEditor),
                    Item("curvature", editor=NumEditor),
                    Item("blending", editor=NumEditor),
                    ))
    
    def _get_radius(self):
        return self.diameter/2.
    
    def _set_radius(self, val):
        self.diamter = val*2.0
        
    @observe("radius, curvature, resolution, blending")
    def _on_changed(self, evt):
        self.update=True
    
    def evaluate_decomposed_rays(self, input_rays):
        start_time = time.monotonic()
        radius = self.diameter/2.
        curvature = self.curvature
        resolution = self.resolution
        blending = self.blending
        
        origin = numpy.asarray(self.centre)
        direction = normaliseVector(numpy.asarray(self.direction))
        axis1 = normaliseVector(numpy.asarray(self.x_axis))
        
        if curvature==0.0:
            curvature=None
        
        gausslets, E_field, uphase = decompose_position(input_rays, origin, direction, axis1, radius, \
                                                        resolution, curvature=curvature, blending=blending)
        
        self._unwrapped_phase = uphase
        self._E_field = E_field
        end_time = time.monotonic()
        print(f"Decomposition in {end_time-start_time} seconds")
        return gausslets


class AngleDecompositionPlane(BaseDecompositionPlane):
    abstract=False
    name = Str("Angle Decomposition")
    
    sample_spacing = Float(1.0) #in microns
    
    width = Range(0,512,64)
    height = Range(0,512,64)
    
    _mask = Instance(numpy.ndarray)
    mask = Property()
    
    ### Maximum ray ange in degrees
    max_angle = Float(30.0)
    
    _E_field = Array() #store the E_field at the aperture for debugging purposes
    _k_field = Array()
    
    traits_view = View(Group(
                    Traceable.uigroup,
                    Item("sample_spacing", editor=NumEditor),
                    Item("width", editor=IntEditor),
                    Item("height", editor=IntEditor),
                    Item("max_angle", editor=NumEditor)
                    ))
    
    def evaluate_decomposed_rays(self, input_rays):
        origin = numpy.asarray(self.centre)
        direction = numpy.asarray(self.direction)
        axis1 = numpy.asarray(self.x_axis)
        axis2 = numpy.cross(axis1, direction)
        spacing = self.sample_spacing
        
        wavelengths = input_rays.wavelengths
        
        x = numpy.arange(self.width)*spacing/1000.0
        x -= x.mean()
        y = numpy.arange(self.height)*spacing/1000.0
        y -= y.mean()
        
        points = x[:,None,None]*axis1[None,None,:] + y[None,:,None]*axis2[None,None,:]
        points += origin[None,None,:]
        
        flat_points = points.reshape(-1,3)
        E_field = eval_Efield_from_gausslets(input_rays, flat_points).reshape(*points.shape)
        self._E_field = E_field
        mask = self._mask
        if mask is not None:
            E_field *= mask[:,:,None]
        max_angle = self.max_angle
        
        pwr = (E_field.real**2 + E_field.imag**2).sum(axis=-1)
        imax = pwr.argmax()
        E_max = E_field.reshape(-1,3)[imax,:]
        p_max = flat_points[imax,:]
        
        new_rays, debug_k = decompose_angle(origin, direction, axis1, E_field, spacing, 
                                   max_angle, wavelengths[0], oversample=4, E_max=E_max, pos_max=p_max)
        print("Decomposed rays count:", len(new_rays))
        self._k_field = debug_k
        return new_rays
    
    def _vtkproperty_default(self):
        return tvtk.Property(color=(0.1,0.1,0.1) )
    
    def _get_mask(self):
        return self._mask
    
    def _set_mask(self, obj):
        w,h = obj.shape
        self._mask = numpy.clip(obj.astype('d'), 0.0, 1.0)
        self.width = w
        self.height = h
        
    def set_circular_mask(self, diameter):
        w = self.width
        h = self.height
        x = numpy.arange(w)*self.sample_spacing
        x -= x.mean()
        y = numpy.arange(h)*self.sample_spacing
        y -= y.mean()
        mask = numpy.zeros((w,h))
        X,Y = numpy.meshgrid(x,y)
        inside = (X**2) + (Y**2) <= ((diameter/2)**2)
        mask[inside] = 1.0
        self.mask = mask
        self.update = True
        
    def set_rectangular_mask(self, width, height):
        w = self.width
        h = self.height
        x = numpy.arange(w)*self.sample_spacing
        x -= x.mean()
        y = numpy.arange(h)*self.sample_spacing
        y -= y.mean()
        X,Y = numpy.meshgrid(x,y)
        inside = numpy.maximum(numpy.abs(X)-(width/2), numpy.abs(Y)-(height/2)) <= 0.0
        mask = numpy.zeros((w,h))
        mask[inside] = 1.0
        self.mask = mask
        self.update = True
        
    @observe("width, height")
    def _chk_dims(self, evt):
        mask = self._mask
        if mask is None:
            return
        if (self.width, self.height) != mask.shape:
            raise ValueError("Width and height must match shape of mask.")
    
    
