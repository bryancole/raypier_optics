===========
Ray Sources
===========

Ray sources generate the input rays for the model. The Ray-source objects also hold the results of the trace, in the 
:py:attr:`traced_rays` attribute, as a list of :py:class:`RayCollection` objects. Each item in this list represents one "generation"
of rays. 

Ray source classes are subclasses of :py:class:`raypier.sources.BaseRaySource`. Besides generating the :py:attr:`input_rays` (a :py:class:`RayCollection`
instance) attribute as the input to the model, and holding the trace results in the :py:attr:`traced_rays` member, the source
objects also control the visualisation of the rays. The following visualisation-related attributes are available.


