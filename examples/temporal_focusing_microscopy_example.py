

from raypier.beamsplitters import UnpolarisingBeamsplitterCube
from raypier.lenses import PlanoConvexLens
from raypier.diffraction_gratings import RectangularGrating
from raypier.gausslet_sources import BroadbandGaussletSource, SingleGaussletSource
from raypier.tracer import RayTraceModel


src = BroadbandGaussletSource(
    origin = (0,0,0),
    direction=(1,0,0.001),
    working_dist=0.0,
    wavelength = 1.0,
    wavelength_extent = 0.1,
    beam_waist = 1000.0,
    display='wires',
    show_paras=True
    )

grating = RectangularGrating(centre=(80.0,0.,0.),
                             direction=(-1,1,0.),
                             length=15.,
                             width=20.0,
                             thickness=3.0,
                             lines_per_mm=1400.0)

lens1 = PlanoConvexLens(centre=(40.0,0.0,0.0),
                        direction=(-1,0,0),
                        diameter=25.0,
                        CT=6.0,
                        n_inside=1.6,
                        curvature=25.0)

lens2 = PlanoConvexLens(centre=(10.0,-30.0,0.0),
                        direction=(0,1,0),
                        diameter=25.0,
                        CT=6.0,
                        n_inside=1.6,
                        curvature=25.0)

# src = SingleGaussletSource(origin = (0,0,0),
#     direction=(1,0,0)
#     )


bs = UnpolarisingBeamsplitterCube(centre = (10.0, 0., 0.),
                                  size=15.0,
                                  )



model = RayTraceModel(optics = [bs, grating, lens1, lens2],
                      sources = [src])

model.configure_traits()

