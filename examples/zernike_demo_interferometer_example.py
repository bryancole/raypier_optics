"""
A modification of the original Michelson interferometer example where the spherical face has been
replaced by a plane-flat with a Zernike polynomial distortion applied.

Any number of Zernike poly terms can be applied. The ANSI standard J-indexing scheme is used.
Use the "Faces" tab of the Propterties window of the Distorted Mirror. Zernike coefficient
amplitudes can be added, updated or removed "live" in the GUI.

"""

from raypier.api import RayTraceModel, UnpolarisingBeamsplitterCube, CollimatedGaussletSource, CircleShape, GeneralLens,\
        PlanarFace, GaussletCapturePlane, EFieldPlane, IntensitySurface
from raypier.intensity_image import IntensityImageView
from raypier.distortions import ZernikeSeries
from raypier.faces import DistortionFace
        

src=CollimatedGaussletSource(origin=(-30,0,0),
                             direction=(1,0,0),
                             radius=5.0,
                             beam_waist=10.0,
                             resolution=10.0,
                             E_vector=(0,1,0),
                             wavelength=1.0,
                             display="wires",
                             opacity=0.05,
                             max_ray_len=50.0)
        
        
bs = UnpolarisingBeamsplitterCube(centre=(0,0,0),
                                  size=10.0)

shape = CircleShape(radius=10.0)

f1 = PlanarFace(mirror=True)

dist = ZernikeSeries(unit_radius=5.0, coefficients=[(3,0.0005),(7,0.0001)])
f2 = PlanarFace(mirror=True) #SphericalFace(curvature=2000.0, mirror=True)
f2 = DistortionFace(base_face=f2, mirror=True, distortion=dist)

m1 = GeneralLens(name="Flat Mirror",
                 shape=shape,
                 surfaces=[f1],
                 centre=(0,20,0),
                 direction=(0,-1,0))

m2 = GeneralLens(name="Distorted Mirror",
                 shape=shape,
                 surfaces=[f2],
                 centre=(20,0,0),
                 direction=(-1,0,0))

cap = GaussletCapturePlane(centre=(0,-20,0),
                           direction=(0,1,0))

field = EFieldPlane(centre=(0,-20,0),
                    direction=(0,1,0),
                    detector=cap,
                    align_detector=True,
                    width=10.0,
                    height=10.0,
                    size=100)

image = IntensityImageView(field_probe=field)
surf = IntensitySurface(field_probe=field)


model = RayTraceModel(optics=[bs, m1, m2],
                      sources=[src],
                      probes=[field, cap],
                      results=[image, surf])

model.configure_traits()
