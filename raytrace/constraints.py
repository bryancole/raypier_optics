from enthought.traits.api import HasTraits, Event, Str

from raytrace.has_queue import HasQueue, on_trait_change

class BaseConstraint(HasQueue):
    update = Event()
    name = Str("a constraint")