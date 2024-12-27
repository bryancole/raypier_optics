
from raypier.tracer import RayTraceModel
from raypier.bezier2d import BSplinePatch

from raypier.sources import ParallelRaySource


ctrl_pts = [[[-9,-9,0], [-3,-9,0], [3,-9,0], [9, -9, 0]],
                [[-9,-3,0], [-3,-3,9], [3,-3,9], [9, -3, 0]],
                [[-9,3,0], [-3,3,9], [3,3,9], [9, 3, 0]],
                [[-9,9,0], [-3,9,0], [3,9,0], [9, 9, 0]],
                ]

patch = BSplinePatch(control_points=ctrl_pts)

src = ParallelRaySource(origin=(0,0,-50.0),
                        direction=(0,0,1),
                        radius = 3.0,
                        rings=2,
                        display="wires",
                        show_normals=True,
                        opacity=0.1)

model = RayTraceModel(optics=[patch], sources=[src], recursion_limit=2)
model.configure_traits()
 