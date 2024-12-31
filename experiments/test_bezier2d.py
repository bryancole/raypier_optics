
from raypier.tracer import RayTraceModel
from raypier.bezier2d import BezierSinglePatch, BSplineSinglePatch

from raypier.sources import ParallelRaySource
from raypier.gausslet_sources import CollimatedGaussletSource


ctrl_pts = [[[-9,-9,0], [-3,-9,0], [3,-9,0], [9, -9, 0]],
                [[-9,-3,0], [-3,-3,9], [3,-3,9], [9, -3, 0]],
                [[-9,3,0], [-3,3,9], [3,3,9], [9, 3, 0]],
                [[-9,9,0], [-3,9,0], [3,9,0], [9, 9, 0]],
                ]

patch = BezierSinglePatch(control_points=ctrl_pts)

patch2 = BSplineSinglePatch()


if True:
    src = ParallelRaySource(origin=(0,0,-50.0),
                            direction=(0,0,1),
                            radius = 3.0,
                            rings=2,
                            display="pipes",
                            show_normals=True,
                            opacity=0.4)
else:
    src = CollimatedGaussletSource(origin=(0,0,-50.0),
                            direction=(0,0,1),
                            radius = 3.0,
                            resolution=0.5,
                            display="pipes",
                            show_normals=True,
                            opacity=0.4)

model = RayTraceModel(optics=[patch], sources=[src], recursion_limit=2)
#model = RayTraceModel(optics=[patch2], sources=[], recursion_limit=2)
model.configure_traits()
 