import numpy as np

from raypier.tracer import RayTraceModel
from raypier.achromats import EdmundOptic45805
from raypier.diffraction_gratings import RectangularGrating
from raypier.mirrors import PECMirror, RectMirror
from raypier.sources import BroadbandRaySource
from raypier.beamstop import BeamStop
from raypier.results import GroupVelocityDispersion
from raypier.constraints import BaseConstraint
from raypier.chirp_result import ChirpResult

from traits.api import Range, on_trait_change
from traitsui.api import View, Item


direction = np.array([0.0,0.0,2.5]) - np.array([-29.685957,124.73850,1.5])

source = BroadbandRaySource(origin=(-29.685957,124.73850,1.5),
                            direction=tuple(direction),
                            wavelength_start=0.76,
                            wavelength_end = 0.80,
                            number=260,
                            uniform_deltaf=True,
                            max_ray_len=300.0)


grating = RectangularGrating(centre=(0,0,0),
                             direction=(0,1,0),
                             lines_per_mm=1400,
                             order=-1)

grating_init_rotation = 41.0
grating.orientation = grating_init_rotation

lens = EdmundOptic45805(centre=(0,75,0),
                        direction=(0,-1,0))
init_lens_rotation = lens.orientation

mir = PECMirror(centre=(0,150,0),
                direction=(0,1,0),
                thickness=5)

mir2 = RectMirror(centre=(-10.413843, 52.149718, -3.0),
                  direction=(0.23150872699334185, -0.9727849212359605, -0.009654343160915735),
                  width=6.0,
                  length=10.0)

bs = BeamStop(centre=(-35,150,15),
              direction=(-0.2314983 ,  0.97131408,  0.05438228))

gvd = GroupVelocityDispersion(source=source, target=bs.faces.faces[0])


class GVDConstraint(BaseConstraint):
    name = "Dispersion Compensator Adjustment"
    focus_adjust = Range(70.0,80.0, value=73.5)
    gdd_adjust = Range(-20.0,40.0, value=0.0)
    lens_rotation = Range(-20.0, 20.0, value=0.0)
    grating_angle = Range(-5.0,5.0, value=0.0)
    
    traits_view = View(Item("focus_adjust"),
                       Item("gdd_adjust"),
                       Item("lens_rotation"),
                       Item("grating_angle"),
                       resizable=True)
    
    @on_trait_change("focus_adjust, gdd_adjust")
    def on_new_values(self):
        lens.centre = (0.0, 75.0 + self.gdd_adjust, 0.0)
        mir.centre = (0.0, 75.0 + self.gdd_adjust + self.focus_adjust, 0.0)
        self.update = True
        
    def _lens_rotation_changed(self):
        lens.orientation = (init_lens_rotation + self.lens_rotation+180.0)%360.0 - 180.0
        self.update = True
        
    def _grating_angle_changed(self):
        grating.orientation = grating_init_rotation + self.grating_angle
        self.update = True
    
        
gvd_cstr = GVDConstraint()
    

chirp = ChirpResult(source=source, target=bs.faces.faces[0])

model = RayTraceModel(optics=[grating, lens, mir, mir2, bs], 
                      sources=[source,],
                      results=[gvd, chirp],
                      constraints=[gvd_cstr,])

model.configure_traits(kind="live")
