

from traits.api import on_trait_change, Float, Instance,Event, Int,\
        Property, Str, Array, cached_property, List, Bool, observe, Button

from traitsui.api import View, Item, VGroup, Tabbed

from tvtk.api import tvtk

from .bases import Probe, Traceable, NumEditor
from .probes import BaseCapturePlane
from .sources import RayCollection, BaseRaySource

from .core.find_focus import find_ray_focus
from .core.utils import normaliseVector, dotprod
from .core.ctracer import GaussletCollection, RayCollection
from .core.fields import eval_Efield_from_gausslets, eval_Efield_from_rays
from .editors import IntEditor

import numpy
    
import time
import traceback
    
    
class EFieldPlane(Probe):
    #: The name of the object in the model-tree
    name = Str("E-Field Probe")
    
    #: An instance of a :py:class:`raypier.probes.BaseCapturePlane` which provides the
    #: incident rays for which the E-field should be evaluated.
    detector = Instance(BaseCapturePlane)
    
    #: If True, the detector object should be positioned at the centre position for this EFieldPlane
    #: whenever the object moves. 
    align_detector = Bool(True)
    
    #: If no detector-object is given, select the input rays by indexing the source traced_rays list. 
    gen_idx = Int(-1)
    
    #: Width of the EFieldPlane
    width = Float(0.5) #in mm
    
    #: Height of the EField Plane
    height = Float(0.5)
    
    #: The resolution of the EFieldPlane. The same value is used for both x- and y-sides.
    size = Int(30)
            
    #: For RayFields (not Gausslets), be can back-project the indicent rays on to a plane 
    #: and evaluate the Gaussia modes from this "exit pupil" plane. This is ignored for 
    #: the E-field evaluation of Gausslets.
    exit_pupil_offset = Float(10.0) #in mm
    
    #: Controls the overlap of the Gaussian modes at the point of evaluation. It's best to
    #: leave this at unity, otherwise one can get confused between this and the blending factor
    #: offered by the gausslet_source classes.
    blending = Float(1.0)
    
    #: When assigned to (triggered) the EField plane will be repositioned on the geometric focus
    #: of the input rays (calculated as the point of closest approach of the input rays).
    centre_on_focus_btn = Button()
    
    #: The output of the probe. A (size,size,3)-shaped complex array.
    E_field = Array()
    
    #: property - The "|E-field|**2" calculated from the E-field.
    intensity = Property(Array, depends_on="E_field")
    
    #: property - The sum of the intensity array.
    total_power = Property(depends_on="intensity, width, height, size")
    
    #: The mean real part of the refractive index for the material where the probe plane lies
    #: Get get this from the selected/captured rays.
    refractive_index = Float(1.0)
        
    _mtime = Float(0.0)
    _src_list = List()
    _plane_src = Instance(tvtk.PlaneSource, (), 
                    {"x_resolution":1, "y_resolution":1},
                    transient=True)
    
    _attrib = Instance(tvtk.ProgrammableAttributeDataFilter,(), transient=True)
                    
    traits_view = View(Tabbed(
                        VGroup(
                       Traceable.uigroup,
                       Item("total_power", style="readonly"),
                       Item('size', editor=IntEditor),
                       Item('width', editor=NumEditor),
                       Item('height', editor=NumEditor),
                       Item('exit_pupil_offset', editor=NumEditor),
                       Item('blending', editor=NumEditor),
                       Item('gen_idx', editor=IntEditor),
                       Item('centre_on_focus_btn', show_label=False, label="Centre on focus")
                   )))
    
    @observe("centre")
    def on_move(self, evt):
        detector = self.detector
        if detector is not None and self.align_detector:
            detector.centre = evt.new
        self._mtime = 0.0
        self.on_change()
    
    @on_trait_change("orientation, size, width, height, exit_pupil_offset, blending, gen_idx")
    def config_pipeline(self):
        src = self._plane_src
        size = self.size
        side = self.width/2.
        yside = self.height/2
        
        src.origin = (-side,-yside,0)
        src.point1 = (-side,yside,0)
        src.point2 = (side,-yside,0)
        
        src.x_resolution = size
        src.y_resolution = size
        
        self._mtime = 0.0
        self.on_change()
        
    def on_change(self):
        try:
            self.evaluate(None)
        except:
            traceback.print_exc()
        self.update=True
    
    def _actors_default(self):
        source = self._plane_src
        attr = self._attrib
        attr.input_connection = source.output_port
        
        trans_f = tvtk.TransformFilter(input_connection=attr.output_port,
                        transform = self.transform)
        map = tvtk.PolyDataMapper(input_connection=trans_f.output_port)
        
        def execute():
            dataobj = attr.poly_data_output
            data = self.intensity.T
            dataobj.cell_data.scalars = data.ravel()
            map.scalar_range = (data.min(), data.max()) 
        attr.set_execute_method(execute)
        
        act = tvtk.Actor(mapper=map)
        actors = tvtk.ActorCollection()
        actors.append(act)
        return actors
    
    @cached_property
    def _get_intensity(self):        
        E = self.E_field
        return (E.real**2).sum(axis=-1) + (E.imag**2).sum(axis=-1)
    
    def _get_total_power(self):
        power = self.intensity.sum()*(self.width*self.height/(self.size**2))
        power *= self.refractive_index
        return power
        
    def evaluate(self, src_list):
        mtime = self._mtime
        detector = self.detector
        if src_list is None:
            src_list = self._src_list
        else:
            self._src_list = src_list
        
        start = time.monotonic()
            
        if detector is not None:
            detector.evaluate(src_list)
            if detector.captured is None:
                return
            if mtime > detector._mtime:
                return
            ray_list = detector.captured
        else:
            idx = self.gen_idx
            ray_list = [src.traced_rays[idx] for src in src_list if src.traced_rays]
        
        n_list = []
        for ct, rays in enumerate(ray_list):
            wavelengths = rays.wavelengths
            
            size = self.size
            side = self.width/2.
            yside = self.height/2
            px = numpy.linspace(-side,side,size)
            py = numpy.linspace(-yside, yside, size)
            
            centre = numpy.asarray(self.centre)
            axis2 = numpy.cross(self.direction, self.x_axis)
            axis1 = numpy.cross(axis2, self.direction)
            
            points = centre[None,None,:] + px[None,:,None]*axis1 + py[:,None,None]*axis2
            points2 = points.reshape(-1,3)
            
            if isinstance(rays, GaussletCollection):
                n_list.append(rays.base_rays.refractive_index.real)
                E = eval_Efield_from_gausslets(rays, points2,
                                               blending=self.blending)
            else:
                n_list.append(rays.refractive_index.real)
                E = eval_Efield_from_rays(rays, points2, wavelengths, 
                                          blending=self.blending,
                                          exit_pupil_offset=self.exit_pupil_offset,
                                          exit_pupil_centre=self.centre)
            
            if ct==0:
                E_field = E.reshape(self.size, self.size, 3)
            else:
                E_field += E.reshape(self.size, self.size, 3)
                
        if not n_list:
            return
                
        self.refractive_index = numpy.concatenate(n_list).mean()
        self.E_field = E_field
        
        self._attrib.modified()
        end = time.monotonic()
        self._mtime = end
        print(f"Field calculation took: {end-start} s")
        
            
    def intersect_plane(self, rays):
        """
        @param rays: a numpy array of ray_t dtype
        
        intersect rays with plane, returning the indices of the intersecting rays
        """
        #raise NotImplementedError
        pass
    
    @observe("centre_on_focus_btn")
    def do_centre_on_focus(self, evt):
        self.centre_on_focus()
    
    def centre_on_focus(self, idx=-1, offset=(0,0,0)):
        """
        Locate the centre of the field-probe on the point of closest approach 
        for the indicated RayCollection (by default, the last one).
        """
        detector = self.detector
        if detector is None:
            src = self.source
            rays = src.traced_rays[idx]
        else:
            rays = detector.captured[0]
        self.centre = tuple(a+b for a,b in zip(find_ray_focus(rays), offset))
    
    
    