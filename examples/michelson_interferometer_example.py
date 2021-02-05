

from raypier.api import RayTraceModel, UnpolarisingBeamsplitterCube, CollimatedGaussletSource, CircleShape, GeneralLens,\
        SphericalFace, PlanarFace, OpticalMaterial, PolarisingBeamsplitterCube, ParallelRaySource
        
        

src=CollimatedGaussletSource(origin=(-30,0,0),
                             direction=(1,0,0),
                             radius=1.0,
                             resolution=1.5,
                             E_vector=(0,1,0),
                             wavelength=1.0)

# src = ParallelRaySource(origin=(-30,0,0),
#                              direction=(1,0,0),
#                              radius=1.0,
#                              rings=0,
#                              E_vector=(0,1,0),
#                              wavelength=1.0)
        
        
bs = UnpolarisingBeamsplitterCube(centre=(0,0,0),
                                  size=10.0)
# bs = PolarisingBeamsplitterCube(centre=(0,0,0),
#                                   size=10.0)


model = RayTraceModel(optics=[bs],
                      sources=[src])

model.configure_traits()
