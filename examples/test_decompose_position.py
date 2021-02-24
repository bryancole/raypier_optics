
from raypier.api import CollimatedGaussletSource, PositionDecompositionPlane, RayTraceModel,\
        PlanoConvexLens, GaussletCapturePlane, EFieldPlane, IntensitySurface
        
from raypier.intensity_image import IntensityImageView
        
        
src = CollimatedGaussletSource(origin=(0,0,0),
                               direction=(0,0,1),
                               E_vector=(1,0,0),
                               E1_amp=1.0,
                               E2_amp=0.0,
                               radius=10.0,
                               resolution=6.0,
                               beam_waist=8.0,
                               wavelength=1.0,
                               display="wires",
                               opacity=0.1)

lens = PlanoConvexLens(centre=(0,0,100),
                       direction=(0,0,1),
                       n_inside=1.5,
                       curvature=25.0,
                       diameter=25.0,
                       CT=5.0)

decomp = PositionDecompositionPlane(centre=(0,0,95), direction=(0,0,1),
                                    radius=10.0, resolution=16., curvature=0.0, blending=1.2)

cap = GaussletCapturePlane(width=30.0,height=30.0)

probe = EFieldPlane(detector=cap,
                    centre=(0,0,110),
                    direction=(0,0,1),
                    width=20.0,
                    height=20.0,
                    size=100)

img = IntensityImageView(field_probe=probe)
surf = IntensitySurface(field_probe=probe)


model = RayTraceModel(optics=[lens, decomp], sources=[src], probes=[cap,probe], results=[img, surf])
model.configure_traits()
