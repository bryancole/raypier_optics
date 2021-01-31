from raypier.mirrors import PECMirror, RectMirror
from raypier.core.cmaterials import OpaqueMaterial
from tvtk.api import tvtk


class BeamStop(PECMirror):
    name = "Beam Stop"
    thickness = 0.5
    
    def _material_default(self):
        return OpaqueMaterial()

class RectTarget(RectMirror):
    name = "Beam Stop"
    thickness = 0.5

    vtkproperty = tvtk.Property(opacity = 0.4,
                             color = (0,0,0))

    def _material_default(self):
        return OpaqueMaterial()
    
