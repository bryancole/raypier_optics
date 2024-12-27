import unittest

from raypier.bezier2d import BSplinePatch
from raypier.core.obbtree import OBBTree

import numpy as np


class TestBezierMeshIntersect(unittest.TestCase):
    def test_intersect(self):
        patch = BSplinePatch()
        
        bezier = patch.bezier
        
        print(bezier.control_pts)
        pts, cells, uvs = bezier.get_mesh(20,20)
        
        tree = OBBTree(pts, np.ascontiguousarray(cells, np.int32))
        tree.max_level = 100
        tree.number_of_cells_per_node = 2
        
        if tree.level <= 0:
            tree.build_tree()
        workspace = np.zeros(tree.level+1, np.int64)
        
        
        p1 = (0.,0.,-50.)
        p2 = (0.,0., 150.)
        it = tree.intersect_with_line(p1,p2)
        
        print("intersection:", it)