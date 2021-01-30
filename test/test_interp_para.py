import os
os.environ['CSF_GraphicShr'] = r"/usr/local/lib/libTKOpenGl.so"

from OCC.Display.wxSamplesGui import display, start_display

from raypier.step_export import make_interp_parabola

para = make_interp_parabola(12.5, 5.0, 15.0, segments=20)

display.DisplayShape(para)

start_display()