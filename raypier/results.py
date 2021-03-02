"""
A module for Results subclasses
"""


from traits.api import Instance, Float, on_trait_change,\
            Button, DictStrFloat, Str

from traitsui.api import View, Item, DropEditor

from raypier.bases import Result, Traceable
from raypier.sources import BaseRaySource
#from raypier.tracer import RayTraceModel
from raypier.core.ctracer import Face
from raypier.dispersion import FusedSilica

from traitsui.editors.drop_editor import DropEditor
from traitsui.editors.title_editor import TitleEditor

import numpy
import itertools
import traceback


class TargetResult(Result):
    """A base-class for result calculated from rays terminating on a chosen face"""
    source = Instance(BaseRaySource)
    target = Instance(Face)
    name = Str("TargetResult")
    result = Float(label="Result", transient=True)
    
    traits_view = View(Item("name", style="readonly", editor=TitleEditor()), 
                       Item('result', style="readonly"),
                       Item('source', editor=DropEditor()),
                       Item('target', editor=DropEditor()),
                       title="Optical path length",
                       resizable=True,
                       )
    
    def calc_result(self, tracer):
        try:
            self._calc_result()
        except:
            traceback.print_exc()
        
    def _calc_result(self):
        raise NotImplementedError()
        
    @on_trait_change("source, target")
    def update(self):
        if self.target is None:
            return
        if self.source is None:
            return
        if self.source.traced_rays:
            try:
                self._calc_result()
            except:
                traceback.print_exc()
        

class MeanOpticalPathLength(TargetResult):
    name = "Mean optical path length"
    abstract = False
    
    result = Float(label="Path length", transient=True)
    
    traits_view = View(Item('result', style="readonly"),
                       Item('source', editor=DropEditor()),
                       Item('target', editor=DropEditor()),
                       title="Optical path length",
                       resizable=True,
                       )
    
    @on_trait_change("source, target")
    def update(self):
        if self.target is None:
            return
        if self.source is None:
            return
        if self.source.traced_rays:
            self._calc_result()
    
    def _calc_result(self):
        all_rays = [r.copy_as_array() for r in reversed(self.source.traced_rays)]
        idx = self.target.idx
        last = all_rays[0]
        selected_idx = numpy.argwhere(last['end_face_idx']==idx).ravel()
        total = numpy.zeros(len(selected_idx), 'd')
        for ray in all_rays:
            selected = ray[selected_idx]
            total += selected['length'] * selected['refractive_index'].real
            selected_idx = selected['parent_idx']
        self.result = total.mean()
        
        
def evaluate_phase(all_wavelengths, traced_rays, target_face,
                   glass_length=0.0, glass_dispersion=FusedSilica()):
    """
    
    all_wavelengths - a numpy array giving the wavelengths, in microns
    traced_rays - the list of RayCollections from the source object after tracing
    target_face - the Face object at which we terminate the tracing
    glass_length - an additional length of glass added to the computation
    glass_dispersion - the DispersionCurve for the extra glass (default Fused Silica)
    
    returns - (freq, phase) #freq in THz
    """
    c = 2.99792458e8 * 1e-9 #convert to mm/ps
    idx = target_face.idx
    last = traced_rays[-1].copy_as_array() #all_rays[0]
    selected_idx = numpy.argwhere(last['end_face_idx']==idx).ravel()
    wavelengths = all_wavelengths[last['wavelength_idx'][selected_idx]]
    sort_idx = numpy.argsort(wavelengths)[::-1]
    wavelengths = wavelengths[sort_idx]
    selected_idx = selected_idx[sort_idx]
    
    idx = len(selected_idx)//2
    phase = last['phase'][selected_idx].copy()
    phase -= phase.mean()
    
    total = last['accumulated_path'][selected_idx]

    if len(total) < 6:
        raise ValueError("Not enough rays to evaluate 2nd and 3rd derivatives")
    ave_path = total.mean()
    total -= ave_path #because we don't care about the overall offset
    f = 1000.0*c/wavelengths #in THz
    
    n_fs = glass_dispersion.evaluate_n(wavelengths).real
    fs_total = n_fs*glass_length*1000.0 #convert path length to mm
    fs_total -= fs_total.mean()
    total += fs_total
    
    phase += total*f*((2*numpy.pi)/c)
    return (f, phase)
    

class GroupVelocityDispersion(MeanOpticalPathLength):
    name = "Group Velocity Dispersion"
    wavelength = Float(label="Wavelength /um", transient=True)
    result = Float(label="GDD /fs^2", transient=True)
    tod = Float(label="TOD", transient=True)
    glass_path = Float(0.0)
    
    _fs = Instance(FusedSilica, ())
    
    traits_view = View(Item("wavelength", style="readonly"),
                       Item('result', style="readonly"),
                       Item('tod', style="readonly"),
                       Item("glass_path"),
                       Item('source', editor=DropEditor()),
                       Item('target', editor=DropEditor()),
                       title="Optical path length",
                       resizable=True,
                       )
    
    @on_trait_change("glass_path")
    def _glass_path_update(self):
        self.update()
    
    def _calc_result(self):
        all_wavelengths = numpy.asarray(self.source.wavelength_list)
        traced_rays = self.source.traced_rays
        target_face = self.target
        glass_length = self.glass_path
        glass_dispersion = self._fs
        f, phase = evaluate_phase(all_wavelengths, traced_rays, target_face, glass_length, glass_dispersion)
        
        omega = f*(2*numpy.pi)
        dw = numpy.diff(omega)
        dw_mean = dw.mean()
        #dw_sd = numpy.std(dw)
        #print("uniformity in f:", dw_sd, dw_mean)
        #print "FREQ:", f
        second_deriv = numpy.diff(phase,2)/(dw_mean**2)
        third_deriv = numpy.diff(phase,3)/(dw_mean**3)
        second_deriv *= 1e6 #convert to fs^2
        third_deriv *= 1e9 #convert to fs^3
        
        c = 2.99792458e8 * 1e-9 #convert to mm/ps
        
        idx = target_face.idx
        self.result = second_deriv[idx]
        self.wavelength = c*1000.0/numpy.median(f) #convert to nm
        self.tod = third_deriv[idx]
        
        
class FocalPoint(TargetResult):
    pass #TODO


class Divergence(TargetResult):
    abstract = False
    name = "Beam Divergence"
    
    traits_view = View(Item('result', style="readonly"),
                       Item('source', editor=DropEditor()),
                       Item('target', editor=DropEditor()),
                       title="Optical path length",
                       resizable=True,
                       )
    
    def _calc_result(self):
        all_rays = [r.copy_as_array() for r in reversed(self.source.traced_rays)]
        idx = self.target.idx
        last = all_rays[0]
        selected_idx = numpy.argwhere(last['end_face_idx']==idx).ravel()
        selected_rays = last[selected_idx]
        directions = selected_rays['direction']
        ave_direction = directions.mean(axis=0,keepdims=True)
        cp_vectors = numpy.cross(ave_direction,directions)
        rms_cp = (cp_vectors**2).sum(axis=-1).mean()
        self.result = rms_cp 
        
        
def get_total_intersections(raysList, face):
    all_rays = itertools.chain(*raysList)
    '''#debug
    print face.idx, " all_rays: "
    for ray in all_rays:
    print ray.end_face_idx  #'''
    idx = face.idx
    return sum(1 for ray in all_rays if ray.end_face_idx==idx)


def get_total_power(raysList, face):
    all_rays = itertools.chain(*raysList)
    idx = face.idx
    return sum(ray.power for ray in all_rays if ray.end_face_idx==idx)


class RayPaths(Result):
    ''' keep a tally of how many rays follow a particular path from face to face.  result 
    is a dictionary where the keys are lists of face indicies like [0,3,2,3,463772762] and 
    the values are the percentage of rays that encountered faces in that order.  The Huge
    face_idx is infinity and is always the last "face" 

    check result.keys() to see which paths were taken.'''

    name = "ray paths"
    result = DictStrFloat(label="Dictionary of ray paths")

    _tracer = Instance("raypier.tracer.RayTraceModel") #to cache the tracer instance
    
    traits_view = View(Item('result', style="readonly"),
                       title="Paths taken by rays",
                       resizable=True,
                       )

    @on_trait_change("_tracer")
    def update(self):
        if self._tracer is not None:
            self._calc_result()
        
    def calc_result(self, tracer):
        self._tracer = tracer
        self._calc_result()
    
    def _calc_result(self):
                
        result = {}    
        for source in self._tracer.sources:
            #a list of RayCollection instances
            rays = source.get_ray_list_by_id()
            #print "how many rays?",len(rays)
            for ray in rays:
                path = []
                for segment in ray:
                    if segment['end_face_idx'] < 4294967200:
                        c = segment['end_face_idx']
                    else: c = numpy.Infinity
                    path.append(c)
                #get_ray_list_by_id goes from last ray to first, which isn't as pretty as forwards in time.
                path.reverse()      
                if str(path) in list(result.keys()):
                    result[str(path)] += 1
                else:
                    result[str(path)] = 1
        total = float(sum(result.values()))
        for key in list(result.keys()):
            result[key] = result[key]/total
        self.result = result   


class Ratio(Result):
    name = "a ratio"
    abstract = False
    nominator = Instance(Face)
    denominator = Instance(Face)
    
    #because I always get them the wrong way round!
    switch_faces = Button("Switch faces")
    
    result = Float(label="Ratio of intersections")
    
    _tracer = Instance("raypier.tracer.RayTraceModel") #to cache the tracer instance
    
    traits_view = View(Item('result', style="readonly"),
                       Item('nominator', editor=DropEditor()),
                       Item('denominator', editor=DropEditor()),
                       Item('switch_faces', show_label=False),
                       title="Face intersection ratio",
                       resizable=True,
                       )
    
    def _switch_faces_changed(self):
        self.nominator, self.denominator = \
            self.denominator, self.nominator
        self._calc_result()
    
    @on_trait_change("nominator, denominator, _tracer")
    def update(self):
        if self._tracer is not None:
            self._calc_result()
        
    def calc_result(self, tracer):
        self._tracer = tracer
        self._calc_result()
    
    def _calc_result(self):
        nom = self.nominator
        denom = self.denominator
        if not all((nom, denom)):
            return
        
        nom_count = 0
        denom_count = 0
        #just sum result from multiple sources?
        #maybe a dictionary or something would be better?
        for source in self._tracer.sources:
            #a list of RayCollection instances
            raysList = source.traced_rays
        
            nom_count = nom_count + get_total_intersections(raysList, nom)
            denom_count = denom_count + get_total_intersections(raysList, denom)
	#print "nom and denom counts", nom_count, denom_count
        try:
            self.result = float(nom_count)/float(denom_count)
        except ZeroDivisionError:
            self.result = numpy.Infinity
        
        
class IncidentPower(Ratio):
    name = "Incident power"
    abstract = False
    object = Instance(Traceable)
    
    traits_view = View(Item('result', style="readonly"),
                       Item('object', editor=DropEditor()),
                       title="Incident power on object",
                       resizable=True,
                       )
    
    def _calc_result(self):
        nom = self.object
        if not nom:
            return
        
        nom_count = 0.0
        power_in = 0.0
        #just sum result from multiple sources?
        #maybe a dictionary or something would be better?
        for source in self._tracer.sources:
            #a list of RayCollection instances
            raysList = source.traced_rays
            input = source.traced_rays[0]
            E1_amp = input.E1_amp
            E2_amp = input.E2_amp
            n = input.refractive_index.real
            power_in += (n*E1_amp*E1_amp.conjugate()).real.sum() + (n*E2_amp*E2_amp.conjugate()).real.sum()
            
            for f in nom.faces.faces:
                for rc in raysList:
                    rca = rc.copy_as_array()
                    selected = rca[rca['end_face_idx']==f.idx]
                    E1_amp = selected['E1_amp']
                    E2_amp = selected['E2_amp']
                    n = selected['refractive_index'].real
                    nom_count += (n*E1_amp*E1_amp.conjugate()).real.sum() +\
                                (n*E2_amp*E2_amp.conjugate()).real.sum()
            
    #print "nom and denom counts", nom_count, denom_count
        try:
            self.result = float(nom_count)/float(power_in)
        except ZeroDivisionError:
            self.result = numpy.Infinity
