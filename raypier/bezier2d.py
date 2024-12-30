
from raypier.core.cbezier import BezierPatch, BezierPatchFace, BSplinePatch
from raypier.core.ctracer import FaceList
from raypier.bases import Optic, Traceable

from traits.api import Instance, Array, Str
from traitsui.api import View, Item, VGroup
from tvtk.api import tvtk

import numpy as np


class BaseUVPatch(Optic):
    control_points = Array(shape=(None,None,3), dtype=np.float64)
    mesh_src = Instance(tvtk.ProgrammableSource, (), transient=True)
    ctrl_pt_src = Instance(tvtk.ProgrammableSource, (), transient=True)
    
    traits_view = View(VGroup(
                       Traceable.uigroup,
                       )
                       )
    
    def _pipeline_default(self):
        mesh_src = self.mesh_src
        def make_mesh():
            pts, cells, uvs = self.bezier.get_mesh(50,50)
            out = mesh_src.poly_data_output
            out.points = pts
            out.polys = cells
        mesh_src.set_execute_method(make_mesh)
        
        t = self.transform
        transf = tvtk.TransformFilter(input_connection=mesh_src.output_port, 
                                      transform=t)
        return transf
    
    def _actors_default(self):
        actors = super()._actors_default()
        
        src = self.ctrl_pt_src
        def make_dataset():
            ctrl_pts = self.control_points
            pts = ctrl_pts.copy().reshape(-1,3)
            n,m = ctrl_pts.shape[:2]
            ids = np.arange(n*m).reshape(n,m)
            lines = []
            for i in range(n-1):
                for j in range(m-1):
                    lines.append( (ids[i,j], ids[i+1,j]) )
                    lines.append( (ids[i,j], ids[i,j+1]) )
            for i in range(n-1):
                lines.append( (ids[i,m-1], ids[i+1,m-1]) )
            for j in range(m-1):
                lines.append( (ids[n-1,j], ids[n-1,j+1]) )
            out = src.poly_data_output
            out.points = pts
            out.lines = np.array(lines)
        src.set_execute_method(make_dataset)
        
        sph_src = tvtk.SphereSource(radius=0.5)
        gly = tvtk.Glyph3D(input_connection=src.output_port)
        gly.set_source_connection(sph_src.output_port)
        
        apd = tvtk.AppendPolyData()
        apd.add_input_connection(0, src.output_port)
        apd.add_input_connection(0, gly.output_port)
        
        t = self.transform
        transf = tvtk.TransformFilter(input_connection=apd.output_port, 
                                      transform=t)
        
        map = tvtk.PolyDataMapper(input_connection=transf.output_port)
        map.scalar_visibility = False
        act = tvtk.Actor(mapper=map)
        actors.append(act)
        
        return actors


class BezierSinglePatch(BaseUVPatch):
    """
    The default material is a PECMaterial, i.e. mirror
    """
    name = Str("Bezier Surface Patch")
    bezier = Instance(BezierPatch)
    
    def _faces_default(self):
        #We need quite a lax tolerance because the initial intersection test using
        #the bezier mesh has a more significant deviation from the true surface
        cfaces = [BezierPatchFace(bezier=self.bezier, invert_normals=True, tolerance=0.01,
                                  u_res=50, v_res=50, atol=1e-10)]
        fl = FaceList(owner=self)
        fl.faces = cfaces
        return fl
    
    def _bezier_default(self):
        ctrl_pts = self.control_points
        N,M, *args = ctrl_pts.shape
        patch = BezierPatch(N-1,M-1)
        patch.control_pts = ctrl_pts
        return patch
        
    def _control_points_default(self):
        data = [[[-9,-9,0], [-3,-9,0], [3,-9,0], [9, -9, 0]],
                [[-9,-3,0], [-3,-3,3], [3,-3,3], [9, -3, 0]],
                [[-9,3,0], [-3,3,3], [3,3,3], [9, 3, 0]],
                [[-9,9,0], [-3,9,0], [3,9,0], [9, 9, 0]],
                ]
        return np.array(data, dtype=np.float64)
    
    
class BSplineSinglePatch(BaseUVPatch):
    name = Str("BSpline Surface Patch")
    bspline = Instance(BSplinePatch)
    
    u_knots = Array()
    v_knots = Array()
    
    # def _faces_default(self):
    #     #We need quite a lax tolerance because the initial intersection test using
    #     #the bezier mesh has a more significant deviation from the true surface
    #     cfaces = [BezierPatchFace(bezier=self.bezier, invert_normals=True, tolerance=0.01,
    #                               u_res=50, v_res=50, atol=1e-10)]
    #     fl = FaceList(owner=self)
    #     fl.faces = cfaces
    #     return fl
    
    def _bspline_default(self):
        ctrl_pts = self.control_points
        N,M, *args = ctrl_pts.shape
        patch = BSplinePatch(N-1,M-1)
        patch.control_pts = ctrl_pts
        patch.u_knots = self.u_knots
        patch.v_knots = self.v_knots
        return patch
        
    def _control_points_default(self):
        data = [[(-44.3783965, -13.956608000000003, 0.0), (-44.3783965, -13.956608000000003, 50.0)], 
                [(-51.80648887000841, 9.784359970257904, 0.0), (-26.902931958404587, -12.498303172462593, 50.0)], 
                [(-31.628220000000002, 31.385614, 0.0), (-28.27639032948056, 37.47853386793429, 50.0)], 
                [(28.903247680789157, 28.39060884033357, 0.0), (-6.514713, 3.249176, 50.0)], 
                [(6.906227945225628, -23.726293952654803, 0.0), (11.917509383107664, -24.67958413884696, 50.0)], 
                [(23.633918710105423, -40.34000411166008, 0.0), (30.655313796701574, 10.287457120427124, 50.0)], 
                [(53.319759, -34.146309, 0.0), (53.319759, -34.146309, 50.0)]]
        return np.array(data, dtype=np.float64)
    
    def _u_knots_default(self):
        data = [0.0, 0.0, 0.0, 0.0, 46.86875896639466, 86.8150208656012, 124.91998295572799, 166.76909466075705, 
                166.76909466075705, 166.76909466075705, 166.76909466075705]
        return np.array(data, dtype=np.float64)
    
    def _v_knots_default(self):
        data = [0.0, 0.0, 1.0, 1.0]
        return np.array(data, dtype=np.float64)
    
    