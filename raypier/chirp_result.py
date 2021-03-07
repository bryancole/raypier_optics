

from raypier.results import TargetResult, evaluate_phase
from raypier.dispersion import FusedSilica

from traits.api import Float, Instance, Bool
from traitsui.api import View, Item, VGroup, HGroup
from chaco.api import Plot, ArrayPlotData
from enable.api import ComponentEditor
from chaco.tools.api import PanTool, ZoomTool

import numpy as np
from numpy import fft
from numpy.core._umath_tests import euclidean_pdist
from traits.has_traits import on_trait_change
from posix import pwrite



class ChirpResult(TargetResult):
    name = "Pulse Intensity Envelope"
    plot = Instance(Plot)
    
    #: The input pulse width, in picoseconds
    pulse_width = Float(0.2) #in picoseconds
    
    #: The length of additional fibre (fused silica) in the optical path, in metres
    glass_path = Float(0.0) #glass path in metres
    
    #: Display the autocorrelation function
    ACF = Bool(False)
    
    #: Force the output pulse peak to be at the centre of the array
    centre = Bool(True)
    
    #: Calculated FWHM for the output pulse width
    fwhm = Float(0.0)
    
    #: Centre wavelength, in microns.
    wavelength = Float(0.780) #in microns
    
    _fs = Instance(FusedSilica, ())
    
    _pd = Instance(ArrayPlotData, ())
    
    traits_view = View(VGroup(HGroup(Item("pulse_width", tooltip="input pulse width (picoseconds)."), 
                                     Item("glass_path", tooltip="Fused silica glass path (metres)")),
                              HGroup(Item("ACF", tooltip="Autocorrelation Function"), Item("centre"),
                              Item("fwhm", label="FWHM", style="readonly", 
                                   tooltip="femtoseconds"))
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
        self._calc_result(redraw=False)
        return plot
    
    @on_trait_change("glass_path, ACF, pulse_width")
    def do_update(self):
        self.update()
    
    def _calc_result(self, redraw=True):
        all_wavelengths = np.asarray(self.source.wavelength_list)
        traced_rays = self.source.traced_rays
        target_face = self.target
        glass_length = self.glass_path
        glass_dispersion = self._fs
        f, phase = evaluate_phase(all_wavelengths, traced_rays, target_face, glass_length, glass_dispersion)
        
        c = 2.99792458e8 * 1e-9 #convert to mm/ps
        hsize = len(phase)//2
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
            
        if self.centre:
            imax = np.argmax(pwr)
            middle = len(pwr)//2
            pwr = np.roll(pwr, middle-imax)
        
        mask = np.arange(len(pwr))[(pwr > (pwr.max()/2))]
        fwhm = t[mask.max()+1] - t[mask.min()]
        self.fwhm = fwhm*1000 #convert to fs
        
        self._pd.set_data('x', t)
        self._pd.set_data('y', pwr)
        if redraw: #prevent a recursive loop.
            self.plot.request_redraw()
        