
from raypier.core.cbezier import BezierPatch, BezierPatchFace
from raypier.core.ctracer import FaceList
from raypier.bases import Optic, Traceable

from traits.api import Instance, Array, Str
from traitsui.api import View, Item, VGroup
from tvtk.api import tvtk

import numpy as np


class BSplinePatch(Optic):
    name = Str("Bezier Surface Patch")
    control_points = Array(shape=(None,None,3), dtype=np.float64)
    
    bezier = Instance(BezierPatch)
    
    mesh_src = Instance(tvtk.ProgrammableSource, (), transient=True)
    ctrl_pt_src = Instance(tvtk.ProgrammableSource, (), transient=True)
    
    traits_view = View(VGroup(
                       Traceable.uigroup,
                       )
                       )
    
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
            