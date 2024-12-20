
from raypier.tracer import RayTraceModel
from raypier.bezier2d import BSplinePatch


ctrl_pts = [[[-9,-9,0], [-3,-9,0], [3,-9,0], [9, -9, 0]],
                [[-9,-3,0], [-3,-3,3], [3,-3,3], [9, -3, 0]],
                [[-9,3,0], [-3,3,3], [3,3,3], [9, 3, 0]],
                [[-9,9,0], [-3,9,0], [3,9,0], [9, 9, 0]],
                ]

patch = BSplinePatch(control_points=ctrl_pts)

model = RayTraceModel(optics=[patch])
model.configure_traits()
 