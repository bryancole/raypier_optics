from enthought.traits.has_traits import HasTraits, \
        on_trait_change as _on_trait_change
        
from enthought.traits.trait_notifiers import TraitChangeNotifyWrapper

from collections import deque, defaultdict
import threading
from functools import wraps


@wraps(_on_trait_change)
def on_trait_change ( name, post_init = False, dispatch="same", *names ):
    """ overrides the has_traits.on_trait_change decorator to
    add some extra metadata to the trait handler function to
    indicate the dispatch-method
    """
    def decorator( function ):
        function = _on_trait_change(name, post_init, *names)( function )
        function.dispatch = dispatch
        return function
    return decorator


class QueuedTraitChangeNotifyWrapper(TraitChangeNotifyWrapper):
    def dispatch(self, handler, *args):
        this_thd = threading.currentThread()
        notification_queue = self.object()._notification_queue
        this_queue = notification_queue[this_thd]
        this_queue.append((handler, args))
            

class HasQueue(HasTraits):    
    """
    This subclass of HasTraits add a new method of trait notification
    dispatch: "queued". The notification queue is per-instance and per-thread. When 
    any attribute on the instance is set, the queue is processed (although
    setting non-trait attributes will not place anything on the queue, so
    no actions will take place).
    
    Use the on_trait_change decorator defined in this module to define the dispatch
    method (the has_trait.on_trait_change decorator doesn't allow a dispatch
    keyword).
    
    "queued" trait-change-handlers are effectively processed after all 
    other ("same") handlers. 
    
    multiple firing of handlers will still occur (i.e. setting a handler to "queued"
    doesn't change the number of times it fires, only the ordering). If re-triggering is 
    not desired, this can easily be caught in the model logic.
    
    N.B. If *all* handlers are set to "queued" this is equivalent to a 
    breadth-first traversal of the dependancy tree (whereas all "same" gives
    a depth-first traversal).
    """
    def __init__(self, *args, **kwds):
        this = super(HasQueue, self)
        this.__init__(*args, **kwds)
        this.__setattr__("_notification_queue", 
                            defaultdict(deque) )
        
    def _on_trait_change( self, handler, name = None, remove = False,
                                 dispatch = 'same', priority = False ):
        dispatch = getattr(handler, "dispatch", dispatch)
        super(HasQueue, self)._on_trait_change( handler, name=name, 
                                                remove=remove,
                                                dispatch=dispatch, 
                                                priority=priority)
    
    def __setattr__(self, name, val):
        """is there a better way to intercept trait-assignment?"""
        q = self._notification_queue
        thd = threading.currentThread()
        if thd not in q:
            q[thd] = this_q = deque()
            try:
                super(HasQueue, self).__setattr__(name, val)
                while True:
                    try:
                        handler, args = this_q.popleft()
                        handler(*args)
                    except IndexError:
                        break
            finally:
                del q[thd]
        else:
            super(HasQueue, self).__setattr__(name, val)

HasQueue.set_trait_dispatch_handler("queued", QueuedTraitChangeNotifyWrapper)




if __name__=="__main__":
    from enthought.traits.api import Float, Bool
    class test(HasQueue):
        """
        A simple test object with a diamond-shape trait dependency
        tree. You can control the ordering of trait handler execution
        by changing the various dispatch methods to "same" or "queued".
        
        trait 'e' includes some logic to prevent re-triggering, as an example.
        """
        
        a = Float(0.0)
        b = Float(0.0)
        c = Float(0.0)
        d = Float(0.0)
        e = Float(0.0)
        
        e_fired = Bool(False)
        
        @on_trait_change("a", dispatch="same")
        def change_a(self, vnew):
            print "a changed to", vnew
            self.b = vnew + 1
            self.c = vnew + 2
            
        @on_trait_change("b", dispatch="same")
        def change_b(self, vnew):
            print "b changed to", vnew
            self.d = vnew + 3
            
        @on_trait_change("c", dispatch="queued")
        def change_c(self, vnew):
            print "c changed to", vnew
            self.d = vnew + 4
            
        @on_trait_change("d", dispatch="same")
        def change_d(self, vnew):
            print "d changed to", vnew
            if not self.e_fired:
                self.e_fired = True
                self.e = vnew + 5
            
        @on_trait_change("e", dispatch="queued")
        def change_e(self, vnew):
            print "e changed to", vnew, self.d + 5
            self.e_fired = False
            
    t = test()
    t.a = 1.0