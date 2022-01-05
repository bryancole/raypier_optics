
from tvtk.api import tvtk

import numpy as np

from raypier.core.obbtree import OBBTree


def get_monkey_actor():
    reader = tvtk.STLReader(file_name="monkey.stl")
    m = tvtk.PolyDataMapper(input_connection=reader.output_port)
    a = tvtk.Actor(mapper=m)
    return a


def show_monkey():
    a = get_monkey_actor()
        
    ren = tvtk.Renderer()
    ren.add_actor(a)
    
    renwin = tvtk.RenderWindow()
    renwin.add_renderer(ren)
    
    iren = tvtk.RenderWindowInteractor(render_window=renwin)
    iren.start()
    
def get_monkey_mesh():
    reader = tvtk.STLReader(file_name="monkey.stl")
    tris = tvtk.TriangleFilter(input_connection=reader.output_port)
    tris.update()
    pd = tris.output
    points = np.asarray(pd.points)
    cells = pd.polys.to_array().reshape(-1,4)[:,1:]
    cells = np.ascontiguousarray(cells, dtype=np.int32)
    return points.copy(), cells.copy()


def get_monkey_obbtree():
    points, cells = get_monkey_mesh()
    print(cells.dtype, cells.shape, points.shape)
    obbtree = OBBTree(points, cells)
    obbtree.max_level = 20
    obbtree.number_of_cells_per_node = 8
    obbtree.build_tree()
    return obbtree
    
    
def append_obb_to_dataset(obb, point_list, line_list=None, poly_list=None):
    corner = obb.corner
    axes = obb.axes
    points = [corner, #0
              corner + axes[0], #1
              corner + axes[1], #2
              corner + axes[2], #3
              corner + axes[0] + axes[1], #4
              corner + axes[0] + axes[2], #5
              corner + axes[1] + axes[2], #6
              corner + axes[0] + axes[1] + axes[2]] #7
    if any(any(not np.isfinite(a) for a in pt) for pt in points):
        return
    start = len(point_list)
    point_list.extend(points)
    if line_list is not None:
        lines = [(0,1),(0,2),(0,3),(2,6),(3,6),(6,7),(1,4),(1,5),(4,7),(5,7),(2,4),(3,5)]
        for i, line in enumerate(lines):
            this = (start + line[0], start + line[1])
            line_list.append(this)

    if poly_list is not None:
        tris = [
            (0,1,2), (1,2,4),
            (0,1,3), (1,3,5),
            (0,2,3), (2,3,6),
            (7,3,5), (7,6,3),
            (7,5,1), (7,1,4),
            (7,6,2), (7,2,4) 
            ]
        for i, tri in enumerate(tris):
            this = tuple(start + vtx for vtx in tri)
            poly_list.append(this)
            
            
def add_subtree(obb, points, lines, polys, level=-1):
    append_obb_to_dataset(obb, points, lines, polys)
    if level == 0:
        return
    level -= 1
    if obb.child1 is not None:
        add_subtree(obb.child1, points, lines, polys, level=level)
    if obb.child2 is not None:
        add_subtree(obb.child2, points, lines, polys, level=level)
            
            
def show_monkey_tree():
    obbtree = get_monkey_obbtree()
    
    print(obbtree.intersect_with_line([0.,0.,0.], [10.,10.,10.]))
    
    print("Level", obbtree.level)
    
    points = []
    lines = []
    polys=[]
    
    add_subtree(obbtree.root, points, lines, polys, level=4)
    
    points = np.array(points)
    lines = np.array(lines)
    polys = np.array(polys)
    
    #print(points)
    #print(polys)
    
    pd = tvtk.PolyData(points=points, lines=lines, polys=polys)
    
    mapper = tvtk.PolyDataMapper(color_mode=2)
    mapper.set_input_data(pd)
    act = tvtk.Actor(mapper=mapper)
    act.property.opacity=0.1
    ren = tvtk.Renderer()
    ren.add_actor(act)
    
    mk = get_monkey_actor()
    ren.add_actor(mk)    
    
    renwin = tvtk.RenderWindow()
    renwin.add_renderer(ren)
    iren = tvtk.RenderWindowInteractor()
    iren._set_render_window(renwin)
    iren.start()
    
    
if __name__=="__main__":
    show_monkey_tree()
    