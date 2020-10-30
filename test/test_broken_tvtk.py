from traits.api import HasTraits, Instance
from traitsui.api import View

from tvtk.api import tvtk
from tvtk.pyface.scene_model import SceneModel
from tvtk.pyface.scene_editor import SceneEditor


class Tester(HasTraits):
    scene = Instance(SceneModel, (), {'background':(1,1,0.8)}, transient=True)
    
    traits_view = View