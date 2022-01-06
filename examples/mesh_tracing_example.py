

from raypier.tracer import RayTraceModel
from raypier.meshes import STLFileMesh
from raypier.sources import ParallelRaySource


src = ParallelRaySource(origin=(0.,0.,-50.0)
    )

m = STLFileMesh(file_name = "../experiments/monkey.stl", scale_factor=20.0)

model = RayTraceModel(sources=[src], optics=[m])
model.configure_traits()

