from OCC import BRepPrimAPI, BRep, STEPControl, TopoDS, gp, BRepBuilderAPI
from itertools import izip, tee

def pairs(itr):
    a,b = tee(itr)
    b.next()
    return izip(a,b)

def make_cylinder(position, direction, radius, length):
    cyl = BRepPrimAPI.BRepPrimAPI_MakeCylinder(radius, length)
    ax = gp.gp_Ax3(gp.gp_Pnt(*position), gp.gp_Dir(*direction))
    trans = gp.gp_Trsf()
    trans.SetTransformation(ax, gp.gp_Ax3())
    t_cyl = BRepBuilderAPI.BRepBuilderAPI_Transform(cyl.Shape(), trans)
    print position, direction, radius, length
    return t_cyl

def make_wire(listOfPoints):
    vertices = [BRepBuilderAPI.BRepBuilderAPI_MakeVertex(gp.gp_Pnt(*p))
                for p in listOfPoints]
    edges = [BRepBuilderAPI.BRepBuilderAPI_MakeEdge(v1.Vertex(),v2.Vertex())
         for v1,v2 in pairs(vertices)]
    wire = BRepBuilderAPI.BRepBuilderAPI_MakeWire()
    for e in edges:
        wire.Add(e.Edge())
    return wire

def make_compound(listOfShapes):
    aRes = TopoDS.TopoDS_Compound()
    aBuilder = BRep.BRep_Builder()
    aBuilder.MakeCompound(aRes)
    for item in listOfShapes:
        aBuilder.Add(aRes, item.Shape())
    return aRes

def export_shapes(shapeList, filename):
    print shapeList
    step_export = STEPControl.STEPControl_Writer()
    for shape in shapeList:
        step_export.Transfer(shape.Shape(),STEPControl.STEPControl_AsIs)
    step_export.Write(str(filename))
    
if __name__=="__main__":
    cyl = make_cylinder((1,2,3),(7,6,5),5,2)
    export_shapes([cyl]*5, "test_step_export.stp")