
from traitsui.api import TextEditor, TupleEditor


NumEditor = TextEditor(auto_set=False, enter_set=True, evaluate=float)
IntEditor = TextEditor(auto_set=False, enter_set=True, evaluate=int)
ComplexEditor = TextEditor(auto_set=False, enter_set=True, 
                           evaluate=float)

ROField = TextEditor(auto_set=False, enter_set=True, evaluate=float)
VectorEditor = TupleEditor(labels=['x','y','z'], auto_set=False, enter_set=True)