from raypier.sources import ParallelRaySource
from raypier.tracer import RayTraceModel
from raypier.corner_cubes import HollowRetroreflector, SolidRetroreflector


source = ParallelRaySource(origin=(-10,0.1,50),
                           direction=(0,0,-1),
                           radius=4,
                           scale=0.1,
                           number=20)

cc = SolidRetroreflector() #
#cc = HollowRetroreflector()

data = """
- !ParallelRaySource
  __traits_version__: 3.3.0
  detail_resolution: 32
  direction: !!python/tuple [0.0, 0.0, -1.0]
  display: pipes
  export_pipes: false
  max_ray_len: 200.0
  name: Ray Source
  number: 20
  origin: !!python/tuple [-10.0, 0.10000000000000001, 50.0]
  radius: 4.0
  rings: 3
  scale: 0.10000000000000001
  scale_factor: 0.20000000000000001
  show_direction: false
  show_normals: false
  show_polarisation: false
  show_start: true
- !SolidRetroreflector
  __traits_version__: 3.3.0
  _orientation: !!python/tuple [0.0, 0.0]
  all_rays: false
  centre: !!python/tuple [0.0, 0.0, 0.0]
  diameter: 25.399999999999999
  direction: !!python/tuple [0.0, 0.0, 1.0]
  display: shaded
  elevation: 0.0
  material: !!python/object:raypier.core.cmaterials.DielectricMaterial {}
  n_inside: !!python/complex '1.5'
  n_outside: !!python/complex '1.0'
  name: Solid Retroreflector
  orientation: 0.0
  rotation: 0.0
  thickness: 20.0
"""

import yaml

print(yaml.dump([source,cc]))

source, cc = yaml.load(data)

model = RayTraceModel(sources=[source],
                      optics=[cc])

model.configure_traits()
#model.update = True
#print source.traced_rays