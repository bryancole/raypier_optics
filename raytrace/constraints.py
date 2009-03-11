from enthought.traits.api import HasTraits, Event, Str


class BaseConstraint(HasTraits):
    update = Event()
    name = Str("a constraint")