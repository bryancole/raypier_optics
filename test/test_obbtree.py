
import unittest

from raypier.core.obbtree import jacobi, OBBTree

import numpy as np

from tvtk.api import tvtk


def view_obb(obb, in_points=None):
    corner = obb.corner
    axes = obb.axes
    points = [corner,
              corner + axes[0],
              corner + axes[1],
              corner + axes[2],
              corner + axes[0] + axes[1],
              corner + axes[0] + axes[2],
              corner + axes[1] + axes[2],
              corner + axes[0] + axes[1] + axes[2]]
    points = np.array(points)
    lines = np.array([(0,1),(0,2),(0,3),(2,6),(3,6),(6,7),(1,4),(1,5),(4,7),(5,7),(2,4),(3,5)])
    pd = tvtk.PolyData(points=points, lines=lines)
    
    mapper = tvtk.PolyDataMapper(color_mode=2)
    mapper.set_input_data(pd)
    act = tvtk.Actor(mapper=mapper)
    ren = tvtk.Renderer()
    ren.add_actor(act)
    
    if in_points is not None:
        pd2 = tvtk.PolyData(points=np.asarray(in_points))
        dot = tvtk.SphereSource(radius=0.5)
        gly = tvtk.Glyph3D()
        gly.set_source_connection(dot.output_port)
        gly.set_input_data(pd2)
        
        map2 = tvtk.PolyDataMapper(input_connection=gly.output_port)
        act2 = tvtk.Actor(mapper=map2)
        ren.add_actor(act2)
    
    
    renwin = tvtk.RenderWindow()
    renwin.add_renderer(ren)
    iren = tvtk.RenderWindowInteractor()
    iren._set_render_window(renwin)
    iren.start()
    


class TestComputeOBBCells(unittest.TestCase):
    def test_compute_obb_cells(self):
        src = tvtk.ArrowSource()
        tris = tvtk.TriangleFilter(input_connection=src.output_port)
        tris.update()
        pd = tris.output
        points = np.asarray(pd.points)
        cells = pd.polys.to_array().reshape(-1,4)[:,1:]
        cells = np.ascontiguousarray(cells, dtype=np.int32)
        print(points.shape)
        print(points)
        print(cells)
        obbtree = OBBTree(points, cells)
        
        celllist = np.arange(len(cells)).astype(np.int32)
        obb = obbtree.compute_obb_cells(celllist)
        
        print(obb.corner)
        print(obb.axes)
        
        view_obb(obb, in_points=points)


class TestComputeOBB(unittest.TestCase):
    def test_compute_point_obb(self):
        pts = np.array([(0.0,0.0,0.0),
                        (1.0,0.0,0.0),
                        (0.0,2.0,0.0),
                        (0.0,0.0,3.0),
                        (1.0,2.0,3.0),
                        (0.0,2.0,3.0),
                        (1.0,0.0,3.0),
                        (1.0,2.0,0.0)])
        cells = np.array([(0,1,2),
                          (3,4,5),
                          (6,7,0)], np.int32)
        obbtree = OBBTree(pts, cells)
        
        obb = obbtree.compute_obb(cells.ravel())
        
        print(obb.corner)
        print(obb.axes)
        
    def test_compute_point_obb2(self):
        ax = np.array([[1.,2.,3.]])**2
        pts = np.random.normal(0.0,1.0,300).reshape(-1,3)
        pts *= ax
        cells = np.array([(0,1,2),
                          (3,4,5),
                          (6,7,0)], np.int32)
        obbtree = OBBTree(pts, cells)
        
        pt_ids = np.arange(pts.shape[0], dtype=np.int32)
        obb = obbtree.compute_obb(pt_ids)
        
        print(obb.corner)
        print(obb.axes)
        
        view_obb(obb, in_points=pts)
        


class TestJacobiIteration(unittest.TestCase):
    def test_jacobi2(self):
        A = np.array([[1,2],[2,4]], dtype='d')
        vals, vects = jacobi(A)
        self.assertEqual(vals[0],0.0)
        self.assertEqual(vals[1],5.0)
        
        residual = (vects[:,0]*vects[:,1]).sum()
        self.assertEqual(residual, 0.0)
        print(vals,vects, residual)
        
    def test_jacobi3(self):
        rt2 = np.sqrt(2)
        A = np.array([[1., rt2, 2],
                      [rt2, 3. ,rt2],
                      [2, rt2, 1]])
        vals, vects = jacobi(A)
        residual1 = (vects[:,0]*vects[:,1]).sum()
        residual2 = (vects[:,0]*vects[:,2]).sum()
        residual3 = (vects[:,1]*vects[:,2]).sum()
        for r in (residual1, residual2, residual3):
            self.assertAlmostEqual(r, 0.0)
        print(vals, vects)
        
    def test_jacobi4(self):
        rt2 = np.sqrt(2)
        A = np.array([[1., rt2, 2],
                      [rt2, 3. ,rt2],
                      [2, rt2, 1]])
        A *= 0.3
        vals, vects = jacobi(A)
        for i in range(3):
            print(np.sqrt((vects[:,i]**2).sum()))
            