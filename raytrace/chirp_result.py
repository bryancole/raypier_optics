

from raytrace.results import TargetResult
from raytrace.dispersion import FusedSilica

from traits.api import Float, Instance, Bool
from traitsui.api import View, Item, VGroup
from chaco.api import Plot, ArrayPlotData
from enable.api import ComponentEditor
from chaco.tools.api import PanTool, ZoomTool

import numpy as np
from numpy import fft
from numpy.core._umath_tests import euclidean_pdist
from traits.has_traits import on_trait_change
from posix import pwrite



class ChirpResult(TargetResult):
    plot = Instance(Plot)
    
    pulse_width = Float(0.2) #in picoseconds
    glass_path = Float(0.0) #glass path in metres
    ACF = Bool(False)
    fwhm = Float(0.0)
    
    wavelength = Float(0.8) #in microns
    
    _fs = Instance(FusedSilica, ())
    
    _pd = Instance(ArrayPlotData, ())
    
    traits_view = View(VGroup(Item("pulse_width"),
                              Item("glass_path"),
                              Item("ACF"),
                              Item("fwhm", label="FWHM", style="readonly", 
                                   tooltip="femtoseconds")
                              ),
                       Item("plot", editor=ComponentEditor(), show_label=False))
    
    def _plot_default(self):
        x = np.linspace(-14, 14, 100)
        y = np.sin(x) * x**3
        plotdata = self._pd
        plotdata.set_data("x", x)
        plotdata.set_data("y", y)

        plot = Plot(plotdata)
        plot.plot(("x", "y"), type="line", color="blue")
        #plot.title = "sin(x) * x^3"
        plot.tools.append(PanTool(plot, drag_button="right"))
        plot.overlays.append(ZoomTool(plot))
        
        return plot
    
    @on_trait_change("glass_path, ACF, pulse_width")
    def do_update(self):
        self.update()
    
    def _calc_result(self):
        c = 2.99792458e8 * 1e-9 #convert to mm/ps
        all_wavelengths = np.asarray(self.source.wavelength_list) #in microns
        all_rays = [r.copy_as_array() for r in reversed(self.source.TracedRays)]
        idx = self.target.idx
        last = all_rays[0]
        selected_idx = np.argwhere(last['end_face_idx']==idx).ravel()
        wavelength_idx = last['wavelength_idx'][selected_idx]
        wavelengths = all_wavelengths[wavelength_idx]
        sort_idx = np.argsort(wavelengths)[::-1]
        wavelengths = wavelengths[sort_idx]
        selected_idx = selected_idx[sort_idx]
        wavelength_idx = last['wavelength_idx'][selected_idx]
        
        phase = last['phase'][selected_idx].copy() #The Grating Phase
        phase -= phase.mean()
        total = np.zeros(len(selected_idx), 'd')
        for ray in all_rays:
            selected = ray[selected_idx]
            wl_idx = selected['wavelength_idx']
            assert np.array_equal(wl_idx, wavelength_idx)
            total += selected['length'] * selected['refractive_index'].real
            selected_idx = selected['parent_idx']
            
        #print "Phase:", phase
        if len(total) < 6:
            return
        total -= total.mean() #because we don't care about the overall offset
        f = 1000.0*c/wavelengths #in THz
        omega = 2 * np.pi * f
        
        n_fs = self._fs.evaluate_n(wavelengths).real
        fs_total = n_fs*self.glass_path*1000.0 #convert path length to mm
        fs_total -= fs_total.mean()
        total += fs_total
        
        hsize = len(phase)//2
        
        phase += total*omega/c
        
        dphi = np.diff(phase)[hsize-hsize//5:hsize + hsize//5].mean()
        phase -= np.arange(len(phase))*dphi
        phase -= phase[hsize]
        
        f0 = c/(self.wavelength*0.001) #in THz
        dt = 0.2*(1/f0)/3.0 #in ps
        size = int(20.0* self.pulse_width/dt)
        #print("SIZE:", size, "F0", f0)
        t = np.arange(size)*dt - 10*self.pulse_width
        omega0 = 2*np.pi*f0
        pulse_in = np.exp(-(t**2)/((self.pulse_width/2)**2))*np.exp(1.0j*omega0*t)
        
        spec = fft.fft(pulse_in)
        f_max = 1./(2*dt)
        f2 = np.roll(np.linspace(-f_max,f_max,len(spec)), len(spec)//2)
        
        phase2 = np.interp(f2, f, phase, 0.0, 0.0) + np.interp(-f2, f, phase, 0.0, 0.0)
        
        bad = np.logical_or(np.abs(f2)<f.min(), np.abs(f2)>f.max())
        spec *= np.exp(1.0j*phase2/1.0)
        spec[bad] = 0.0
        
        pulse_out = fft.ifft(spec)
        
        pwr = pulse_out.real**2 + pulse_out.imag**2
        
        if self.ACF:
            psp = fft.rfft(pwr)
            pwr = fft.irfft(psp * psp.conjugate())
            pwr = np.roll(pwr, len(pwr)//2)
        
        mask = np.arange(len(pwr))[(pwr > (pwr.max()/2))]
        fwhm = t[mask.max()+1] - t[mask.min()]
        self.fwhm = fwhm*1000 #convert to fs
        
        self._pd.set_data('x', t)
        self._pd.set_data('y', pwr)
        self.plot.request_redraw()
        