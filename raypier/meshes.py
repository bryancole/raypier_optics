

from raypier.bases import Traceable, Optic
from raypier.core.obbtree import OBBTree, OBBTreeFace
from raypier.core.ctracer import FaceList, InterfaceMaterial
from raypier.core.cmaterials import OpaqueMaterial

from tvtk.api import tvtk

from traits.api import File, Instance, observe, Int, Float
from traitsui.api import View, VGroup, Item

import numpy as np


class STLFileMesh(Optic):
    file_name = File()
    
    obbtree = Instance(OBBTree)
    
    scale_factor = Float(10.0)
    max_level = Int(100)
    cells_per_node = Int(6)
    
    _reader = Instance(tvtk.TriangleFilter)
    
    traits_view = View(VGroup(
                       Traceable.uigroup,
                       Item('n_inside'),
                       ))
    
    #material = Instance(InterfaceMaterial)
    
    # def _vtkproperty_default(self):
    #     return tvtk.Property(opacity=0.4)
    #
    # def _material_default(self):
    #     print("FOO")
    #     return OpaqueMaterial()
    
    def _faces_default(self):
        #self.apply_materials()
        fl = FaceList(owner=self)
        fl.faces = [ 
                  OBBTreeFace(tree = self.obbtree, material=self.material)
                ]
        return fl
    
    @observe(["file_name", "scale_factor"])
    def _on_fname_changed(self, evt):
        src = tvtk.STLReader(file_name=self.file_name)
        scale = tvtk.Transform()
        scale.scale(self.scale_factor, self.scale_factor, self.scale_factor)
        f1 = tvtk.TransformFilter(input_connection=src.output_port,
                             transform=scale)
        tris = tvtk.TriangleFilter(input_connection=f1.output_port)
        tris.update()
        self._reader = tris
        pd = tris.output
        points = np.asarray(pd.points)
        cells = pd.polys.to_array().reshape(-1,4)[:,1:]
        cells = np.ascontiguousarray(cells, dtype=np.int32)
        obbtree = OBBTree(points.copy(), 
                          cells.copy())
        obbtree.max_level = self.max_level
        obbtree.number_of_cells_per_node = self.cells_per_node
        self.obbtree = obbtree

    def _pipeline_default(self):
        src = self._reader
        
        transF = tvtk.TransformFilter(input_connection=src.output_port, 
                                      transform=self.transform)
        return transF