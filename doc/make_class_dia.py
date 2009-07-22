import pydot
import types
import itertools

from raytrace import (bases, faces, ellipsoids, lenses, 
                      mirrors, parabolics, paraxial, rays,
                      sources, tracer)

modules = (bases, faces, ellipsoids, lenses, 
                      mirrors, parabolics, paraxial, rays,
                      sources, tracer)

roots = (bases.Renderable, faces.Face, sources.BaseRaySource,
         tracer.RayTraceModel)

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
    print a

g = pydot.graph_from_edges(edges)
g.aspect = (1.0, 5)
g.write_png("class_heirarchy.png", prog='dot')

