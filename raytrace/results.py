"""
A module for Results subclasses
"""

from enthought.traits.api import Instance, Float, on_trait_change,\
            Button

from enthought.traits.ui.api import View, Item, DropEditor

from raytrace.bases import Result
from raytrace.faces import Face
from raytrace.tracer import RayTraceModel
from enthought.traits.ui.editors.drop_editor import DropEditor

import numpy

def get_total_intersections(raysList, face):
    return sum((rays.end_face==face).sum() for rays in raysList)


class Ratio(Result):
    name = "a ratio"
    nominator = Instance(Face)
    denominator = Instance(Face)
    
    #because I always get them the wrong way round!
    switch_faces = Button("Switch faces")
    
    result = Float(label="Ratio of intersections")
    
    _tracer = Instance(RayTraceModel) #to cache the tracer instance
    
    traits_view = View(Item('result', style="readonly"),
                       Item('nominator', editor=DropEditor()),
                       Item('denominator', editor=DropEditor()),
                       Item('switch_faces', show_label=False),
                       title="Face intersection ratio",
                       resizable=True,
                       )
    
    def _switch_faces_changed(self):
        self.nominator, self.denominator = \
            self.denominator, self.nominator
        self._calc_result()
    
    @on_trait_change("nominator, denominator, _tracer")
    def update(self):
        if self._tracer is not None:
            self._calc_result()
        
    def calc_result(self, tracer):
        self._tracer = tracer
        self._calc_result()
    
    def _calc_result(self):
        nom = self.nominator
        denom = self.denominator
        if not all((nom, denom)):
            return
        
        #maybe the object needs a source trait
        #instead, just take the first source found
        source = self._tracer.sources[0]
        
        #a list of RayCollections
        raysList = source.TracedRays
        
        nom_count = get_total_intersections(raysList, nom)
        denom_count = get_total_intersections(raysList, denom)
        
        try:
            self.result = float(nom_count)/float(denom_count)
        except ZeroDivisionError:
            self.result = numpy.Infinity
        