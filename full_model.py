import pickle as pickle
from ray_tracer import  BuildRaySet, RayTraceModel
from prisms import Prism, Rhomboid

input_rays = BuildRaySet(origin = (-7.42,15,0),
                         direction = (0,-1,0),
                         radius=1.0,
                         count=20)

prism = Prism(n_inside = 1.764+0.0j,
                    orientation=0.0,
                    elevation=0.0,
                    centre=(0,0,0),
                    rotation=0,
                    length=10.0,
                    height=7.07,
                    width=14.0,
                    )

rhomboid = Rhomboid(n_inside = 1.764+0.0j,
                    orientation=0.0,
                    elevation=0.0,
                    centre=(0,0,0),
                    rotation=0,
                    length=10.0,
                    height=7.07,
                    width=14.0,
                    slant=45.0)

#beamstop = BeamStop(width=10,height=10,
#                    centre=(7,0,-20))

#print "rhomboid", rhomboid.polydata
#print "beamstop", beamstop.polydata

#fname = "RayTracer.pickle"
#try:
#    raise IOError("pass")
#    fobj = file(fname, 'rb')
#    model = pickle.load(fobj)
#    fobj.close()
#except (IOError, EOFError, AttributeError):
#    model = RayTraceModel(optics=[prism, rhomboid], rays=input_rays)
    
model = RayTraceModel(optics=[prism, rhomboid], rays=input_rays)
model.configure_traits()

#fobj = file("RayTracer.pickle",'wb')
#pickle.dump(model, fobj)
#fobj.close()