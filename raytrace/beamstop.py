from raytrace.mirrors import PECMirror
from raytrace.cmaterials import OpaqueMaterial



class BeamStop(PECMirror):
    name = "Beam Stop"
    thickness = 0.5
    
    def _material_default(self):
        return OpaqueMaterial()
    