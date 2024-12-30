

from raypier.aspherics import Thorlabs_A230_A
from raypier.sources import ParallelRaySource
from raypier.tracer import RayTraceModel
from raypier.results import FocalPoint
from raypier.windows import CircularWindow

lens1 = Thorlabs_A230_A()

win = CircularWindow(centre=(0.,0.,4.0), direction=(0.,0.,1.),
                     thickness=0.25, diameter=5.0)

src = ParallelRaySource(origin=(0.,0.,-10.), 
                        direction=(0.,0.,1.),
                        radius=2.5,
                        rings=7,
                        wavelength=0.78,
                        max_ray_len=20.0,
                        show_normals=True,
                        display="wires",
                        scale_factor=0.02)

focus = FocalPoint(source=src)

model = RayTraceModel(optics=[lens1, win], sources=[src], results=[focus])
model.configure_traits()

"""
In our model the focal point comes out as 5.8mm from origin. 
The PDF datasheet for the A230 optic gives a focal point of 5.85mm.
I'm not sure where the discrepancy comes from.
"""

