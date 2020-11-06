
from raytrace.aspherics import Thorlabs355440_B, Thorlabs355440, Thorlabs352440
from raytrace.tracer import RayTraceModel
from raytrace.sources import ConfocalRaySource
from raytrace.mirrors import PlanarDispersiveWindow, PlanarWindow

from raytrace.dispersion import NamedDispersionCurve

from raytrace.find_focus import find_ray_focus

import numpy
import time

lens = Thorlabs355440_B(centre=(0,0,0),
                    direction=(1,0,0))
#lens.centre = (lens.CT,0,0)

window = PlanarWindow(centre=(-1.0,0,0),
                        direction=(1,0,0),
                        thickness=0.25,
                        diameter=5.0,
                        n_inside=1.5112)

source = ConfocalRaySource(focus = (-7.09,0,0),
                           direction = (1,0,0),
                           working_dist= 10.0,
                           theta=12.0,
                           wavelen=0.78)
 
model = RayTraceModel(optics=[lens],
                      sources=[source])

model.configure_traits()

x1 = numpy.linspace(2.0,10.0,500)
x2 = []
spr = []

def find_rms_spread(ray_set, x_offset):
    s = ray_set.origin
    d = ray_set.direction
    y_vals = s[:,1] + (x_offset-s[:,0])*d[:,1]/d[:,0]
    z_vals = s[:,2] + (x_offset-s[:,0])*d[:,2]/d[:,0]
    y_vals -= y_vals.mean()
    z_vals -= z_vals.mean()
    return numpy.sqrt((y_vals**2 + z_vals**2).mean())

model.trace_all()
start = time.clock()
for x in x1:
    source.focus = (-x,0,0)
    #window.centre = (-x + 0.5, 0.0, 0.0)
    source.theta = 6.0 * numpy.arctan2(2.,x)/0.2
    model.trace_all()
    last_rays = source.TracedRays[-1]
    focus = find_ray_focus(last_rays)
    x_offset = focus[0]
    x2.append( x_offset - lens.CT )
    
    rms = find_rms_spread(last_rays, x_offset)
    spr.append( rms )
    
    #print(x, source.theta,  x2[-1] , rms)
end = time.clock()
print("Took:", end-start)

from matplotlib import pyplot as pp
pp.plot(x1,x2,'b-', label="focus")
pp.grid(True)
ax = pp.gca()
ax.set_title(f"{lens.name} - reversed - with 0.25mm laser window")
#ax.set_title(f"{lens.name} - no laser window")
ax.set_xlabel("S1 Working dist /mm")
ax.set_ylabel("S2 Working dist /mm")


ax2 = ax.twinx()
ax2.set_ylabel("Point spread RMS /microns")
ax2.plot(x1, 1000*numpy.array(spr), 'r-', label="point spread")

lines, labels = ax.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
ax2.legend(lines + lines2, labels + labels2, loc=0)
pp.show()

model.configure_traits()

