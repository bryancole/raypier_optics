#    Copyright 2009, Teraview Ltd., Bryan Cole
#
#    This file is part of Raypier.
#
#    Raypier is free software: you can redistribute it and/or modify
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
        BRepBuilderAPI, BRepAlgoAPI, BRepOffsetAPI, Geom, TopAbs, TopExp,\
        XCAFApp, STEPCAFControl, TDocStd, TCollection,\
        XCAFDoc, Quantity, TopLoc
from itertools import tee
import itertools
import numpy

def pairs(itr):
    a,b = tee(itr)
    next(b)
    return zip(a,b)

def sign(val):
    return val/abs(val)

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

def position_shape(shape, centre, direction, x_axis):
    ax = gp.gp_Ax3()
    ax.SetLocation(gp.gp_Pnt(*centre))
    ax.SetDirection(gp.gp_Dir(*direction))
    ax.SetXDirection(gp.gp_Dir(*x_axis))
    
    tr = gp.gp_Trsf()
    tr.SetTransformation(ax, gp.gp_Ax3())
    
    loc = TopLoc.TopLoc_Location(tr)
    
    shape.Location(loc)
    #trans = BRepBuilderAPI.BRepBuilderAPI_Transform(shape, tr)
    #return toshape(trans)
    return shape

def make_box(position, direction, x_axis, dx, dy, dz):
    box = BRepPrimAPI.BRepPrimAPI_MakeBox(gp.gp_Pnt(-dx/2., -dy/2., 0.0),
                                          dx, dy, -dz)
    ax = gp.gp_Ax2(gp.gp_Pnt(*position), 
                   gp.gp_Dir(*direction),
                   gp.gp_Dir(*x_axis))
    ax3 = gp.gp_Ax3()
    trans = gp.gp_Trsf()
    trans.SetTransformation(gp.gp_Ax3(ax), ax3)
    t_box = BRepBuilderAPI.BRepBuilderAPI_Transform(box.Shape(), trans)
    return toshape(t_box)

def make_cylinder(position, direction, radius, length, offset, x_axis):
    cyl_ax = gp.gp_Ax2(gp.gp_Pnt(offset,0,0), 
                       gp.gp_Dir(0,0,1), 
                       gp.gp_Dir(1,0,0))
    cyl = BRepPrimAPI.BRepPrimAPI_MakeCylinder(cyl_ax, radius, length)
    ax = gp.gp_Ax2(gp.gp_Pnt(*position), 
                   gp.gp_Dir(*direction),
                   gp.gp_Dir(*x_axis))
    ax3 = gp.gp_Ax3()
    trans = gp.gp_Trsf()
    trans.SetTransformation(gp.gp_Ax3(ax), ax3)
    t_cyl = BRepBuilderAPI.BRepBuilderAPI_Transform(cyl.Shape(), trans)
    print(position, direction, radius, length)
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

def make_sphere(position, radius):
    sph = BRepPrimAPI.BRepPrimAPI_MakeSphere(radius)
    trans = gp.gp_Trsf()
    trans.SetTranslation(gp.gp_Vec(*position))
    t_sph = BRepBuilderAPI.BRepBuilderAPI_Transform(sph.Shape(), trans)
    return toshape(t_sph)

# def make_spherical_lens(position, direction, x_axis, 
#                         diameter, curvature1, curvature2, 
#                         centre_height1, centre_height2, offset):
#     radius = diameter/2.
#     lower_z = min(centre_height1+curvature1, centre_height2+curvature2,
#                   centre_height1, centre_height2)
#     upper_z = max(centre_height1+curvature1, centre_height2+curvature2,
#                   centre_height1, centre_height2)
#     length = upper_z-lower_z
#     cyl_ax = gp.gp_Ax2(gp.gp_Pnt(offset,0,lower_z), 
#                        gp.gp_Dir(0,0,1), 
#                        gp.gp_Dir(1,0,0))
#     cyl = BRepPrimAPI.BRepPrimAPI_MakeCylinder(cyl_ax, radius, length)
#     
#     
    

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
    
    ##comment these out to reinstate the cut
    #t_block = position_shape(toshape(block), centre, direction, x_axis)
    #t_ellipse = position_shape(ellipsoid, centre, direction, x_axis)
    #return make_compound([toshape(block), ellipsoid])
    
    #cut = toshape(BRepAlgoAPI.BRepAlgoAPI_Cut(toshape(block), ellipsoid))
    cut = make_compound([toshape(block), ellipsoid])
    return position_shape(cut, centre, direction, x_axis)

def make_interp_parabola(FL, rmin, rmax, segments=50):
    A = 1./(4*FL)
    x = numpy.linspace(rmin, rmax, segments)
    y = (A * x**2) - FL
    
    points = [(X,0,Z) for X,Z in zip(x,y)]
    points.append((x[0],0,y[-1]))
    
    def pairs(itr):
        a,b = itertools.tee(itr)
        next(b)
        return zip(a,b)
    
    edges = (BRepBuilderAPI.BRepBuilderAPI_MakeEdge(
                    gp.gp_Pnt(*p1), gp.gp_Pnt(*p2))
                    for p1, p2 in pairs(points))
    last_edge = BRepBuilderAPI.BRepBuilderAPI_MakeEdge(
                    gp.gp_Pnt(*points[-1]), 
                    gp.gp_Pnt(*points[0]))
    
    wire = BRepBuilderAPI.BRepBuilderAPI_MakeWire()
    for e in edges:
        wire.Add(e.Edge())
    wire.Add(last_edge.Edge())
    
    face = BRepBuilderAPI.BRepBuilderAPI_MakeFace(wire.Wire())
    
    ax = gp.gp_Ax1(gp.gp_Pnt(0,0,0),
                   gp.gp_Dir(0,0,1))
    
    revol = BRepPrimAPI.BRepPrimAPI_MakeRevol(face.Shape(), ax)
    
    return revol.Shape()

def make_true_para(FL, rmin, rmax):
    """
    makes a parabloid, with rotation axis along Z and focus
    at the origin
    """
    ax = gp.gp_Ax2(gp.gp_Pnt(0.,0.,-FL), #origin
                   gp.gp_Dir(1.,0.,0.), #main direction is Z
                   gp.gp_Dir(0.,0.,1.)) #X Direction is X
    para = Geom.Geom_Parabola(ax, FL)
    h_para = Geom.Handle_Geom_Parabola(para)
    
    ax2 = gp.gp_Ax2(gp.gp_Pnt(0,0,0), #origin
                   gp.gp_Dir(0.,0.,1.), #main direction is Z
                   gp.gp_Dir(1.,0.,0.)) #X Direction is X
    
    pbl_shape = BRepPrimAPI.BRepPrimAPI_MakeRevolution(ax2, h_para,
                                                       rmin, rmax)
    return pbl_shape.Shape()

def make_OAP(EFL, diameter, height, centre, direction, x_axis):
    FL = EFL/2. #focal length
    radius = diameter/2.
    outside = EFL + radius
    inside = EFL - radius
    length = (outside**2)/(4.*FL) - FL + height
    
    #pbl_shape = make_true_para(FL, 5.0, outside+1)
    pbl_shape = make_interp_parabola(FL, inside-1, outside+1)
    
    ax3 = gp.gp_Ax2(gp.gp_Pnt(EFL,0,-height), #origin
                   gp.gp_Dir(0.,0.,1.), #main direction is X
                   gp.gp_Dir(0.,1.,0.)) #X Direction is Y
    cyl_solid = BRepPrimAPI.BRepPrimAPI_MakeCylinder(ax3, radius, length)
    
    
    cut = BRepAlgoAPI.BRepAlgoAPI_Cut(cyl_solid.Shape(), 
                                      pbl_shape)
        
    loc= position_shape(toshape(cut), centre, 
                          direction, x_axis)
    #nurbs = BRepBuilderAPI.BRepBuilderAPI_NurbsConvert(loc)
    return loc #toshape(nurbs)

def make_spherical_lens(CT, diameter, curvature, 
                        centre, direction, x_axis):
    cax = gp.gp_Ax2(gp.gp_Pnt(0,0,CT-curvature),
                    gp.gp_Dir(0,1,0),
                    gp.gp_Dir(1,0,0))
    circ = Geom.Geom_Circle(cax, curvature)
    h_circ = Geom.Handle_Geom_Circle(circ)
    
    r = diameter/2.
    h2 = CT - curvature + numpy.sqrt(curvature**2 - r**2)
    p1 = gp.gp_Pnt(0,0,CT)
    p2 = gp.gp_Pnt(r,0,h2)
    p3 = gp.gp_Pnt(r,0,0)
    p4 = gp.gp_Pnt(0,0,0)
    
    #ps = p1,p2,p3,p4
    #vs = [BRepBuilderAPI.BRepBuilderAPI_MakeVertex(p) for p in ps]
    
    e1 = BRepBuilderAPI.BRepBuilderAPI_MakeEdge(h_circ, p1,
                                                p2)    
    e2 = BRepBuilderAPI.BRepBuilderAPI_MakeEdge(p2,p3)
    e3 = BRepBuilderAPI.BRepBuilderAPI_MakeEdge(p3,p4)
    e4 = BRepBuilderAPI.BRepBuilderAPI_MakeEdge(p4,p1)
    
    wire = BRepBuilderAPI.BRepBuilderAPI_MakeWire()
    for e in (e1,e2,e3,e4):
        print(e)
        wire.Add(e.Edge())
    
    face = BRepBuilderAPI.BRepBuilderAPI_MakeFace(wire.Wire())
    
    ax = gp.gp_Ax1(gp.gp_Pnt(0,0,0),
                   gp.gp_Dir(0,0,1))
    solid = BRepPrimAPI.BRepPrimAPI_MakeRevol(face.Shape(), ax)
    
    return position_shape(toshape(solid), centre, direction, x_axis)


def make_spherical_lens2(CT1, CT2, diameter, 
                         curvature1, curvature2, 
                        centre, direction, x_axis):
    cax = gp.gp_Ax2(gp.gp_Pnt(0,0,CT1-curvature1),
                    gp.gp_Dir(0,sign(curvature1),0),
                    gp.gp_Dir(1,0,0))
    circ = Geom.Geom_Circle(cax, abs(curvature1))
    h_circ = Geom.Handle_Geom_Circle(circ)
    
    cax2 = gp.gp_Ax2(gp.gp_Pnt(0,0,CT2-curvature2),
                    gp.gp_Dir(0,-sign(curvature2),0),
                    gp.gp_Dir(1,0,0))
    circ2 = Geom.Geom_Circle(cax2, abs(curvature2))
    h_circ2 = Geom.Handle_Geom_Circle(circ2)
    
    r = diameter/2.
    
    h2 = CT1 - curvature1 + numpy.sqrt(curvature1**2 - r**2)*sign(curvature1)
    h3 = CT2 - curvature2 + numpy.sqrt(curvature2**2 - r**2)*sign(curvature2)
    p1 = gp.gp_Pnt(0,0,CT1)
    p2 = gp.gp_Pnt(r,0,h2)
    p3 = gp.gp_Pnt(r,0,h3)
    p4 = gp.gp_Pnt(0,0,CT2)
    
    e1 = BRepBuilderAPI.BRepBuilderAPI_MakeEdge(h_circ, p1, p2)    
    e2 = BRepBuilderAPI.BRepBuilderAPI_MakeEdge(p2,p3)
    e3 = BRepBuilderAPI.BRepBuilderAPI_MakeEdge(h_circ2, p3,p4)
    e4 = BRepBuilderAPI.BRepBuilderAPI_MakeEdge(p4,p1)
    
    wire = BRepBuilderAPI.BRepBuilderAPI_MakeWire()
    for e in (e1,e2,e3,e4):
        print(e)
        wire.Add(e.Edge())
    
    face = BRepBuilderAPI.BRepBuilderAPI_MakeFace(wire.Wire())
    
    ax = gp.gp_Ax1(gp.gp_Pnt(0,0,0),
                   gp.gp_Dir(0,0,1))
    solid = BRepPrimAPI.BRepPrimAPI_MakeRevol(face.Shape(), ax)
    
    return position_shape(toshape(solid), centre, direction, x_axis)
    

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
    aBuilder._InputShapes = listOfShapes
    aRes._builder = aBuilder
    return aRes    

def make_rays_both(listOfRays, scale):
    wires = make_rays_wires(listOfRays, scale)
    pipes = make_rays_pipes(listOfRays, scale)
    return make_compound([wires, pipes])

def make_rays_wires(listOfRays, scale=None):
    raysItr = iter(listOfRays)
    first = next(raysItr)
    def MakeVertex(pt): 
        vt = BRepBuilderAPI.BRepBuilderAPI_MakeVertex(gp.gp_Pnt(*pt))
        #print "make vertex", pt, vt
        return vt
    def MakeEdge(v1, v2):
        e = BRepBuilderAPI.BRepBuilderAPI_MakeEdge(v1.Vertex(), v2.Vertex())
        #print "make edge", v1, v2, e
        return e
    wires = [BRepBuilderAPI.BRepBuilderAPI_MakeWire() 
             for i in range(first.origin.shape[0])]
    v_start = [MakeVertex(pt) for pt in first.origin]
    v_end = [MakeVertex(pt) for pt in first.termination]
    first_edges = [MakeEdge(v1,v2) for v1, v2 in zip(v_start, v_end)]
    for edge, wire in zip(first_edges, wires):
        wire.Add(edge.Edge())
    id_map = list(range(len(wires)))
    for rays in raysItr:
        id_map = [id_map[pid] for pid in rays.parent_idx]
        v_start = [v_end[pid] for pid in rays.parent_idx]
        v_end = [MakeVertex(pt) for pt in rays.termination]
        edges = [MakeEdge(v1,v2) for v1, v2 in zip(v_start, v_end)]
        for edge, w_id in zip(edges, id_map):
            wires[w_id].Add(edge.Edge())
    return make_compound([w.Shape() for w in wires])

def make_rays_pipes(listOfRays, radius=0.1):
    itrRays = iter(listOfRays)
    first = next(itrRays)
    v_lists = [[s,e] for s,e in zip(first.origin, first.termination)]
    id_map = list(range(len(v_lists)))
    for rays in itrRays:
        id_map = [id_map[pid] for pid in rays.parent_idx]
        v_end = list(rays.termination)
        for w_id, v in zip(id_map, v_end):
            v_lists[w_id].append(v)
        
    ### Fusing shape is finicky
    #_fuse = fuse_shapes
    _fuse = make_compound
    shapes=[]    
    for vlist in v_lists:
        cylinders = (make_cylinder_2(s, e, radius) for s,e in pairs(vlist))
        spheres = (make_sphere(p, radius) for p in vlist[1:-1])
        interleave = list(next(it) for it in itertools.cycle([cylinders,
                                                          spheres]))
        shapes.append(_fuse(interleave))
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
    print(shapeList)
    step_export = STEPControl.STEPControl_Writer()
    for shape in shapeList:
        step_export.Transfer(shape,STEPControl.STEPControl_AsIs)
    step_export.Write(str(filename))
    
def get_color_map():
    c = [a for a in dir(Quantity) if a.startswith("Quantity_NOC_")]
    names = [a.split('_')[-1].lower() for a in c]
    vals = [getattr(Quantity, n) for n in c]
    return dict(list(zip(names,vals)))
    
def export_shapes2(shapeList, filename, colorList=[]):
    h_doc = TDocStd.Handle_TDocStd_Document()
    app = XCAFApp.XCAFApp_Application_GetApplication().GetObject()
    app.NewDocument(TCollection.TCollection_ExtendedString("MDTV-CAF"),h_doc)
    doc = h_doc.GetObject()
    h_shape_tool = XCAFDoc.XCAFDoc_DocumentTool().ShapeTool(doc.Main())
    h_color_tool = XCAFDoc.XCAFDoc_DocumentTool().ColorTool(doc.Main())
    shape_tool = h_shape_tool.GetObject()
    color_tool = h_color_tool.GetObject()
    top_label = shape_tool.NewShape()
    
    colorList.extend(["brown"]*(len(shapeList) - len(colorList)))
    cmap = get_color_map()
    
    tr = gp.gp_Trsf()
    loc = TopLoc.TopLoc_Location(tr)
    print(colorList)
    print(shapeList)
    colorMap = dict((c, Quantity.Quantity_Color(cmap.get(c, 133))) for c in colorList)
    
    for shape, color in zip(shapeList, colorList):
        if not shape:
            continue
        print("color:", color, "shape", shape)
        label = shape_tool.AddShape(shape, False, True)
        ref_label = shape_tool.AddComponent(top_label, label, loc)
        c = colorMap[color]
        color_tool.SetColor(ref_label, c, XCAFDoc.XCAFDoc_ColorGen)
        
    mode = STEPControl.STEPControl_AsIs
    writer = STEPCAFControl.STEPCAFControl_Writer()
    writer.Transfer(h_doc, mode)
    writer.Write(str(filename))
    
if __name__=="__main__":
    cyl = make_cylinder((1,2,3),(7,6,5),5,2, 0, (1,1,1) )
    export_shapes2([cyl]*5, "test_step_export.stp")