.. SpeX Prism Library Analysis Toolkit documentation master file, created by
   sphinx-quickstart on Sat Jul 11 20:07:28 2015.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Main SPLAT module
===============================================================

.. toctree
   :maxdepth: 3


Examples
---------------------------


Classes
---------------------------
.. autoclass:: splat.Spectrum
  :members:


Routines
---------------------------

Data Access
^^^^^^^^^^^^^^
.. autofunction:: splat.getSpectrum
.. autofunction:: splat.loadSpectrum
.. autofunction:: splat.readSpectrum
.. autofunction:: splat.searchLibrary
.. autofunction:: splat.fetchDatabase
.. autofunction:: splat.keySource
.. autofunction:: splat.keySpectrum


Spectral Classification
^^^^^^^^^^^^^^
.. autofunction:: splat.classifyByIndex
.. autofunction:: splat.classifyByStandard
.. autofunction:: splat.classifyByTemplate
.. autofunction:: splat.classifyGravity
.. autofunction:: splat.metallicity


Other Spectral Analyses
^^^^^^^^^^^^^^
.. autofunction:: splat.compareSpectra
.. autofunction:: splat.filterMag
.. autofunction:: splat.measureIndex
.. autofunction:: splat.measureIndexSet
.. autofunction:: splat.measureEW
.. autofunction:: splat.measureEWSet


Source Analysis
^^^^^^^^^^^^^^
.. autofunction:: splat.typeToMag
.. autofunction:: splat.typeToTeff
.. autofunction:: splat.estimateDistance


Conversion Routines
^^^^^^^^^^^^^^
   
.. autofunction:: splat.caldateToDate
.. autofunction:: splat.dateToCaldate
.. autofunction:: splat.coordinateToDesignation
.. autofunction:: splat.designationToCoordinate
.. autofunction:: splat.designationToShortName
.. autofunction:: splat.typeToNum
.. autofunction:: splat.properCoordinates
.. autofunction:: splat.isNumber


I/O Routines
^^^^^^^^^^^^^^

.. autofunction:: splat.checkOnline
.. autofunction:: splat.checkAccess
.. autofunction:: splat.checkFile
.. autofunction:: splat.checkLocalFile

Other Routines
^^^^^^^^^^^^^^
.. autofunction:: splat.test
.. autofunction:: splat.weightMeanVar




* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

