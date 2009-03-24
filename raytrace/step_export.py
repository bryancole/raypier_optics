#    Copyright 2009, Teraview Ltd., Bryan Cole
#
#    This file is part of Raytrace.
#
#    Raytrace is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.

from OCC import BRepPrimAPI, BRep, STEPControl, TopoDS, gp, \
        BRepBuilderAPI, BRepAlgoAPI, BRepOffsetAPI, Geom
from itertools import izip, tee
import itertools
import numpy

def pairs(itr):
    a,b = tee(itr)
    b.next()
    return izip(a,b)

def MakeVertex(pt): 
    vt = BRepBuilderAPI.BRepBuilderAPI_MakeVertex(gp.gp_Pnt(*pt))
    #print "make vertex", pt, vt
    return vt

def MakeEdge(v1, v2):
    e = BRepBuilderAPI.BRepBuilderAPI_MakeEdge(v1.Vertex(), v2.Vertex())
    #print "make edge", v1, v2, e
    return e

def toshape(builder):
    shape = builder.Shape()
    shape._builder = builder
    return shape

def make_cylinder(position, direction, radius, length):
    cyl = BRepPrimAPI.BRepPrimAPI_MakeCylinder(radius, length)
    ax = gp.gp_Ax3(gp.gp_Pnt(*position), gp.gp_Dir(*direction))
    trans = gp.gp_Trsf()
    trans.SetTransformation(ax, gp.gp_Ax3())
    t_cyl = BRepBuilderAPI.BRepBuilderAPI_Transform(cyl.Shape(), trans)
    print position, direction, radius, length
    return toshape(t_cyl)

def make_cylinder_2(start, end, radius):
    start = numpy.asarray(start)
    end = numpy.asarray(end)
    direction = end-start
    length = numpy.sqrt((direction**2).sum())
    cyl = BRepPrimAPI.BRepPrimAPI_MakeCylinder(radius, length)
    ax = gp.gp_Ax3(gp.gp_Pnt(*start), gp.gp_Dir(*direction))
    trans = gp.gp_Trsf()
    trans.SetTransformation(ax, gp.gp_Ax3())
    t_cyl = BRepBuilderAPI.BRepBuilderAPI_Transform(cyl.Shape(), trans)
    return toshape(t_cyl)

def scale(sx, sy, sz):
    t = gp.gp_GTrsf()
    t.SetValue(1,1,sx)
    t.SetValue(2,2,sy)
    t.SetValue(3,3,sz)
    return t

def make_ellipsoid(focus1, focus2, major_axis):
    """
    @param focus1: length 3 sequence giving first focus location
    @param focus2: length 3 sequence giving second focus location
    @param path_length: major axis length
    """
    f1 = numpy.asarray(focus1)
    f2 = numpy.asarray(focus2)
    direction = -(f1 - f2)
    centre = (f1 + f2)/2.
    sep = numpy.sqrt((direction**2).sum())
    minor_axis = numpy.sqrt( major_axis**2 - (sep/2.)**2 )
    
    sphere = BRepPrimAPI.BRepPrimAPI_MakeSphere(minor_axis)
    
    scale = gp.gp_GTrsf()
    scale.SetValue(3,3, major_axis/minor_axis)
    
    ellipse = BRepBuilderAPI.BRepBuilderAPI_GTransform(sphere.Shape(), scale)
    
    loc = gp.gp_Ax3()
    loc.SetLocation(gp.gp_Pnt(*centre))
    loc.SetDirection(gp.gp_Dir(*direction))
    
    tr = gp.gp_Trsf()
    tr.SetTransformation(loc, gp.gp_Ax3())
    
    trans = BRepBuilderAPI.BRepBuilderAPI_Transform(ellipse.Shape(), tr)
    shape = toshape(trans)
    return shape

def make_ellipsoid_2(focus1, focus2, major_axis):
    f1 = numpy.asarray(focus1)
    f2 = numpy.asarray(focus2)
    direction = -(f1 - f2)
    centre = (f1 + f2)/2.
    sep = numpy.sqrt((direction**2).sum())
    minor_axis = numpy.sqrt( major_axis**2 - (sep/2.)**2 )
    
    el = gp.gp_Elips(gp.gp_Ax2(gp.gp_Pnt(0,0,0), gp.gp_Dir(1,0,0)), 
                          major_axis, minor_axis)
    edge1 = BRepBuilderAPI.BRepBuilderAPI_MakeEdge(el, 
                                                  gp.gp_Pnt(0,0,major_axis),
                                                  gp.gp_Pnt(0,0,-major_axis))
    edge2 = BRepBuilderAPI.BRepBuilderAPI_MakeEdge(gp.gp_Pnt(0,0,-major_axis),
                                                  gp.gp_Pnt(0,0,major_axis))
    
    wire = BRepBuilderAPI.BRepBuilderAPI_MakeWire()
    wire.Add(edge1.Edge())
    wire.Add(edge2.Edge())
    
    face = BRepBuilderAPI.BRepBuilderAPI_MakeFace(wire.Wire())
    
    el = BRepPrimAPI.BRepPrimAPI_MakeRevol(face.Shape(), 
                                           gp.gp_Ax1(gp.gp_Pnt(0,0,0),
                                                     gp.gp_Dir(0,0,1)),
                                           numpy.pi*2)
    
    loc = gp.gp_Ax3()
    loc.SetLocation(gp.gp_Pnt(*centre))
    loc.SetDirection(gp.gp_Dir(*direction))
    
    tr = gp.gp_Trsf()
    tr.SetTransformation(loc, gp.gp_Ax3())
    
    trans = BRepBuilderAPI.BRepBuilderAPI_Transform(el.Shape(), tr)
    shape = toshape(trans)
    return shape

def make_ellipsoid_mirror(f1, f2, major_axis, 
                          x_bounds, y_bounds, z_bounds,
                          centre, direction, x_axis):
    ellipsoid = make_ellipsoid(f1, f2, major_axis)
    P1 = gp.gp_Pnt(x_bounds[0], y_bounds[0], z_bounds[0])
    P2 = gp.gp_Pnt(x_bounds[1], y_bounds[1], z_bounds[1])
    block = BRepPrimAPI.BRepPrimAPI_MakeBox(P1, P2)
    
    #cut = BRepAlgoAPI.BRepAlgoAPI_Cut(block.Shape(), ellipsoid)
    cut = make_compound([block.Shape(), ellipsoid])
    
    print "ell", centre, direction, x_axis
    ax = gp.gp_Ax3()
    ax.SetLocation(gp.gp_Pnt(*centre))
    ax.SetDirection(gp.gp_Dir(*direction))
    ax.SetXDirection(gp.gp_Dir(*x_axis))
    
    tr = gp.gp_Trsf()
    tr.SetTransformation(ax, gp.gp_Ax3())
    
    trans = BRepBuilderAPI.BRepBuilderAPI_Transform(cut, tr)
    shape = toshape(trans)
    return shape

def make_wire(listOfPoints, close=False):
    vertices = [BRepBuilderAPI.BRepBuilderAPI_MakeVertex(gp.gp_Pnt(*p))
                for p in listOfPoints]
    edges = [BRepBuilderAPI.BRepBuilderAPI_MakeEdge(v1.Vertex(),v2.Vertex())
         for v1,v2 in pairs(vertices)]
    if close:
        edges.append(BRepBuilderAPI.BRepBuilderAPI_MakeEdge(vertices[-1].Vertex(),
                                                            vertices[0].Vertex()))
    wire = BRepBuilderAPI.BRepBuilderAPI_MakeWire()
    for e in edges:
        wire.Add(e.Edge())
    return toshape(wire)

def make_compound(listOfShapes):
    aRes = TopoDS.TopoDS_Compound()
    aBuilder = BRep.BRep_Builder()
    aBuilder.MakeCompound(aRes)
    for item in listOfShapes:
        aBuilder.Add(aRes, item)
    return aRes

def make_rays_both(listOfRays, scale):
    wires = make_rays(listOfRays, scale)
    pipes = make_rays_4(listOfRays, scale)
    return make_compound([wires, pipes])

def make_rays(listOfRays, scale=None):
    raysItr = iter(listOfRays)
    first = raysItr.next()
    def MakeVertex(pt): 
        vt = BRepBuilderAPI.BRepBuilderAPI_MakeVertex(gp.gp_Pnt(*pt))
        #print "make vertex", pt, vt
        return vt
    def MakeEdge(v1, v2):
        e = BRepBuilderAPI.BRepBuilderAPI_MakeEdge(v1.Vertex(), v2.Vertex())
        #print "make edge", v1, v2, e
        return e
    wires = [BRepBuilderAPI.BRepBuilderAPI_MakeWire() 
             for i in xrange(first.origin.shape[0])]
    v_start = [MakeVertex(pt) for pt in first.origin]
    v_end = [MakeVertex(pt) for pt in first.termination]
    first_edges = [MakeEdge(v1,v2) for v1, v2 in izip(v_start, v_end)]
    for edge, wire in izip(first_edges, wires):
        wire.Add(edge.Edge())
    id_map = range(len(wires))
    for rays in raysItr:
        id_map = [id_map[pid] for pid in rays.parent_ids]
        v_start = [v_end[pid] for pid in rays.parent_ids]
        v_end = [MakeVertex(pt) for pt in rays.termination]
        edges = [MakeEdge(v1,v2) for v1, v2 in izip(v_start, v_end)]
        for edge, w_id in izip(edges, id_map):
            wires[w_id].Add(edge.Edge())
    return make_compound([w.Shape() for w in wires])

def make_rays_2(listOfRays, radius):
    raysItr = iter(listOfRays)
    first = raysItr.next()
    sections = []
    v_start = list(first.origin)
    v_end = list(first.termination)
    sections = zip(v_start, v_end)
    id_map = range(len(sections))
    for rays in raysItr:
        id_map = [id_map[pid] for pid in rays.parent_ids]
        v_start = [v_end[pid] for pid in rays.parent_ids]
        v_end = list(rays.termination)
        sections.extend(zip(v_start, v_end))
    cylinders = [make_cylinder_2(s, e, radius) for s,e in sections]
    return make_compound(cylinders)

def make_rays_3(listOfRays):
    sections = []
    for rays in listOfRays:
        section = BRepOffsetAPI.BRepOffsetAPI_ThruSections(True)
        wire1 = make_wire(list(rays.origin), close=True)._builder
        wire2 = make_wire(list(rays.termination), close=True)._builder
        section.AddWire(wire1.Wire())
        section.AddWire(wire2.Wire())
        sections.append(toshape(section))
    return make_compound(sections)

def make_rays_4(listOfRays, radius=0.1):
    itrRays = iter(listOfRays)
    first = itrRays.next()
    v_lists = [[s,e] for s,e in izip(first.origin, first.termination)]
    id_map = range(len(v_lists))
    for rays in itrRays:
        id_map = [id_map[pid] for pid in rays.parent_ids]
        v_end = list(rays.termination)
        for w_id, v in izip(id_map, v_end):
            v_lists[w_id].append(v)
        
    shapes=[]    
    for vlist in v_lists:
        cylinders = [make_cylinder_2(s, e, radius) for s,e in pairs(vlist)]
        shapes.append(fuse_shapes(cylinders))
        #shapes.extend(cylinders)
    return make_compound(shapes)

def fuse_shapes(shapeList):
    s1 = shapeList[0]
    _builders = []
    for s2 in shapeList[1:]:
        fuse = BRepAlgoAPI.BRepAlgoAPI_Fuse(s1, s2)
        _builders.append(fuse)
        s1 = fuse.Shape()
    s1._builders = _builders
    return s1

def export_shapes(shapeList, filename):
    print shapeList
    step_export = STEPControl.STEPControl_Writer()
    for shape in shapeList:
        step_export.Transfer(shape,STEPControl.STEPControl_AsIs)
    step_export.Write(str(filename))
    
if __name__=="__main__":
    cyl = make_cylinder((1,2,3),(7,6,5),5,2)
    export_shapes([cyl]*5, "test_step_export.stp")