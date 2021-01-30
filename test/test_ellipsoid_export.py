from raypier.step_export import make_ellipsoid_mirror, export_shapes

f1 = (22.5, 45.1, 50.0)
f2 = (92.5, 45.1, 20.0)

major_axis=52.0

x_bounds = (0.0, 45.0)
y_bounds = (0.0, 45.0)
z_bounds = (0.0, 40.0)

x_axis = (1.0,0.0,0.0)
centre = (0.0,0.0,0.0)
direction = (0.0,0.0,1.0)

e = make_ellipsoid_mirror(f1, f2, major_axis, x_bounds, y_bounds, z_bounds, centre, direction, x_axis)

export_shapes([e], "/home/bryan/EllipsoidMirror.step")
