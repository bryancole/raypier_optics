import pydot
import types
import itertools

from raytrace import (bases, ellipsoids, lenses, 
                      mirrors, beamstop,
                      sources, tracer, splines,
                      beamstop, prisms, dielectrictroughs,
                      results, corner_cubes, cfaces, cmaterials, ctracer)

modules = (bases, ellipsoids, lenses, 
                      mirrors, beamstop,
                      sources, tracer, splines,
                      beamstop, prisms, dielectrictroughs,
                      results, corner_cubes, cfaces, cmaterials, ctracer)

roots = (bases.Renderable, sources.BaseRaySource,
         tracer.RayTraceModel, ctracer.Face, ctracer.InterfaceMaterial)

def iterclasses(m):
    for name in dir(m):
        obj = getattr(m, name)
        try:
            if issubclass(obj, roots):
                yield obj
        except TypeError:
            pass
            
all_classes = set(itertools.chain(*[iterclasses(m) for m in modules]))

def gen_pairs(cls):
    for base in cls.__bases__:
        for pair in gen_pairs(base):
            yield pair
        yield (cls, base)
        
all_pairs = itertools.chain(*[gen_pairs(c) for c in all_classes])

def make_name(cls):
    return "%s.%s"%(cls.__module__, cls.__name__)

all_names = ((make_name(b), make_name(a)) for a,b in all_pairs)

edges = set(all_names)
for a in edges:
    print(a)

g = pydot.graph_from_edges(edges)
g.aspect = (8, 5)
g.set_rankdir("LR")
g.set_ranksep("0.1")
#g.set_ratio("1.0")
g.write_png("class_heirarchy.png", prog='dot')

