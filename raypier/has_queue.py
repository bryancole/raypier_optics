from traits.has_traits import HasTraits, \
        on_trait_change as _on_trait_change
        
from traits.trait_notifiers import TraitChangeNotifyWrapper

from collections import deque
from threading import currentThread
from functools import wraps


@wraps(_on_trait_change)
def on_trait_change ( name, post_init = False, 
                      dispatch="same", 
                      retrigger="all",
                      *names ):
    """ overrides the has_traits.on_trait_change decorator to
    add some extra metadata to the trait handler function to
    indicate the dispatch-method
    """
    def decorator( function ):
        function = _on_trait_change(name, post_init, *names)( function )
        function.dispatch = dispatch
        function.retrigger = retrigger
        return function
    return decorator


class QueuedTraitChangeNotifyWrapper(TraitChangeNotifyWrapper):
    def __init__(self, handler, owner, target=None):
        retrigger = getattr(handler, 'retrigger', 'all')
        self.dispatch = getattr(self, self._policy_map[retrigger])
        TraitChangeNotifyWrapper.__init__(self, handler, owner)
    
    def _dispatch_all(self, handler, *args):
        this_thd = currentThread()
        notification_queue = self.object().__notification_queue__
        this_queue = notification_queue[this_thd]
        this_queue.append((handler, args))
        
    def _dispatch_replace(self, handler, *args):
        this_thd = currentThread()
        notification_queue = self.object().__notification_queue__
        this_queue = notification_queue[this_thd]
        handlers = [h for h,a in this_queue]
        count = handlers.count(handler)
        for i in range(count):
            idx = [h for h,a in this_queue].index(handler)
            del this_queue[idx]
        this_queue.append((handler, args))
        
    _policy_map = {'all': '_dispatch_all',
                   'replace': '_dispatch_replace'}
            
            
def wrap_queue( func ):
    def wrapped( self, *args, **kwds ):
        s = super(HasQueue, self)
        q = self.__notification_queue__
        thd = currentThread()
        if thd not in q:
            q[thd] = this_q = deque()
            try:
                func(self, *args, **kwds)
                while True:
                    try:
                        handler, args = this_q.popleft()
                    except IndexError:
                        break
                    handler(*args)
            finally:
                del q[thd]
        else:
            func(self, *args, **kwds)
    return wrapped


class HasQueue(HasTraits):    
    """
    This subclass of HasTraits add a new method of trait notification
    dispatch: "queued". By default the scope of the queue is global across all
    subclasses of HasQueue but you can limit this by adding a __notification_queue__={}
    class attribute to specific subclasses, or even instances. 
    
    When any trait on the instance is set, the queue is processed after all
    synchronous trait handlers (dispatch="same") have completed.
    
    Use the on_trait_change decorator defined in this module to define the dispatch
    method (the has_trait.on_trait_change decorator doesn't allow a dispatch
    keyword).
    
    multiple firing of handlers will still occur (i.e. setting a handler to "queued"
    doesn't change the number of times it fires, only the ordering). If re-triggering is 
    not desired, this can easily be caught in the model logic.
    
    N.B. If *all* handlers are set to "queued" this is equivalent to a 
    breadth-first traversal of the dependancy tree (whereas all "same" gives
    a depth-first traversal).
    """
    __notification_queue__ = {}
    
    __init__ = wrap_queue( HasTraits.__init__ )
    
    trait_set = wrap_queue( HasTraits.trait_set )
        
    def _on_trait_change( self, handler, name = None, remove = False,
                                 dispatch = 'same', priority = False, target = None):
        dispatch = getattr(handler, "dispatch", dispatch)
        super(HasQueue, self)._on_trait_change( handler, name=name, 
                                                remove=remove,
                                                dispatch=dispatch, 
                                                priority=priority)
    
    def __setattr__(self, name, val):
        """is there a better way to intercept trait-assignment?"""
        s = super(HasQueue, self)
        if name in self.trait_names():
            q = self.__notification_queue__
            thd = currentThread()
            if thd not in q:
                q[thd] = this_q = deque()
                try:
                    s.__setattr__(name, val)
                    while True:
                        try:
                            handler, args = this_q.popleft()
                        except IndexError:
                            break
                        handler(*args)
                finally:
                    del q[thd]
            else:
                s.__setattr__(name, val)
        else:
            s.__setattr__(name, val)

HasQueue.set_trait_dispatch_handler("queued", QueuedTraitChangeNotifyWrapper)




if __name__=="__main__":
    from traits.api import Float, Bool
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
        
        @on_trait_change("a", dispatch="same")
        def change_a(self, vnew):
            print("a changed to", vnew)
            self.c = vnew + 2
            self.b = vnew + 1
            
        @on_trait_change("b", dispatch="same")
        def change_b(self, vnew):
            print("b changed to", vnew)
            self.d = vnew + 3
            
        @on_trait_change("c", dispatch="same")
        def change_c(self, vnew):
            print("c changed to", vnew)
            self.d = vnew + 4
            
        @on_trait_change("d", dispatch="queued", retrigger="replace")
        def change_d(self, vnew):
            print("d changed to", vnew)
            self.e = vnew + 5
            
    t = test()
    t.a = 1.0