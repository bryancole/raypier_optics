Jupyter Integration
===================

Sometimes the live GUI isn't what you want. Sometimes you want to manipulate or
analyse a model in a jupyter notebooks environment. If you want to display the model 
view in a notebook, you simply call the 
:py:meth:`~raypier.tracer.RayTraceModel.ipython_view` method of the
model object to show the view. The arguments `width` and `height` specify the size of the 
displayed image. The optional argument `view` is a dictionary describing the view camera parameters. 

.. py:class:: raypier.tracer.RayTraceModel

.. py:method:: raypier.tracer.RayTraceModel.ipython_view(width : int, height : int, view={}) -> view dict

    Renders the 3D view of the model as a bitmap to display in a jupyter notebook. 
    The image is displayed along with a set of ipywidget buttons to allow the 
    user to rotate, pan and zoom the view.
    
    The method returns a `dict` containing the view pan/zoom/rotation state. This dict
    can then be passed to subsequent calls to `ipython_view()` as the `view` keyword argument,
    to create another view with similar state.
          
An example of jupyter notebook integration can be found here: `Jupyter Notebook Example <_static/Jupyter_Notebook_Example.html>`_
The original .ipynb document is included in the examples folder.

Note, the view-control buttons are not present in the html-output of the link above. You need to run the notebook
in a jupyter notebook server to get the widgets.
