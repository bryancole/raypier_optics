
from raytrace.tracer import RayTraceModel
from raytrace.sources import HexagonalRayFieldSource
from raytrace.lenses import PlanoConvexLens
from raytrace.fields import EFieldPlane
from raytrace.constraints import BaseConstraint

from traits.api import Range, on_trait_change
from traitsui.api import View, Item


lens = PlanoConvexLens(centre=(0,0,20),
                       direction=(0,0,1),
                       diameter=25.0,
                       thickness=6.0,
                       curvature=40.0,
                       n_inside=1.5)

src = HexagonalRayFieldSource(spacing=7.0, direction=(0,0,1))

probe = EFieldPlane(source=src,
                    centre=(0,0,70),
                    direction=(0,0,1),
                    exit_pupil_offset=100.,
                    width=1,
                    size=100)


class FocalPlane(BaseConstraint):
    z_pos = Range(20.0,100.0)
    
    traits_view = View(Item("z_pos"))
    
    @on_trait_change("z_pos")
    def on_change_z_pos(self):
        probe.centre = (0,0,self.z_pos)
        
        

model = RayTraceModel(sources=[src], optics=[lens],
                      probes=[probe], constraints=[FocalPlane()])



model.configure_traits()

for item in src.iter_neighbours():
    print(item)