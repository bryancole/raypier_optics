from raytrace.faces import CircularFace
from raytrace.mirrors import PECMirror

class BeamStopFace(CircularFace):
    pass


class BeamStop(PECMirror):
    name = "Beam Stop"
    thickness = 0.5
    
    def _faces_default(self):
        return [BeamStopFace(owner=self)]
    