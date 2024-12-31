

from raypier.tracer import RayTraceModel
from raypier.meshes import STLFileMesh
from raypier.sources import ParallelRaySource


src = ParallelRaySource(origin=(0.,0.,-50.0),
                        show_normals=True
    )

m = STLFileMesh(file_name = "../experiments/monkey.stl", scale_factor=20.0,
                n_inside=1.5)

model = RayTraceModel(sources=[src], optics=[m])
model.configure_traits()

