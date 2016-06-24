API
===============================================

SPLAT Classes
-------------

Spectrum
^^^^^^^^
.. autoclass:: splat.Spectrum
	:members:


SPLAT Routines
--------------

Data Access
^^^^^^^^^^^
.. autofunction:: splat.getSpectrum
.. autofunction:: splat.getStandard
.. autofunction:: splat_db.getPhotometry
.. autofunction:: splat_db.searchLibrary
.. autofunction:: splat_db.fetchDatabase
.. autofunction:: splat_db.keySource
.. autofunction:: splat_db.keySpectrum


Spectral Classification
^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: splat.classifyByIndex
.. autofunction:: splat.classifyByStandard
.. autofunction:: splat.classifyByTemplate
.. autofunction:: splat.classifyGravity
.. autofunction:: splat.initiateStandards
.. autofunction:: splat.metallicity


Spectrophotometry
^^^^^^^^^^^^^^^^^^
.. autofunction:: splat.filterInfo
.. autofunction:: splat.filterProperties
.. autofunction:: splat.filterMag


Other Spectral Analyses
^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: splat.compareSpectra
.. autofunction:: splat.measureIndex
.. autofunction:: splat.measureIndexSet
.. autofunction:: splat.measureEW
.. autofunction:: splat.measureEWSet


Source Analysis
^^^^^^^^^^^^^^^
.. autofunction:: splat.typeToMag
.. autofunction:: splat.typeToTeff
.. autofunction:: splat.estimateDistance


Conversion Routines
^^^^^^^^^^^^^^^^^^^
   
.. autofunction:: splat.caldateToDate
.. autofunction:: splat.dateToCaldate
.. autofunction:: splat.coordinateToDesignation
.. autofunction:: splat.designationToCoordinate
.. autofunction:: splat.designationToShortName
.. autofunction:: splat.typeToNum
.. autofunction:: splat.properCoordinates
.. autofunction:: splat.isNumber


Plotting Routines
^^^^^^^^^^^^^^^^^
.. autofunction:: splat_plot.plotSpectrum
.. autofunction:: splat_plot.plotBatch
.. autofunction:: splat_plot.plotSequence
.. autofunction:: splat_plot.plotSED
.. autofunction:: splat_plot.plotIndices


Modeling Routines
^^^^^^^^^^^^^^^^^
.. autofunction:: splat_model.getModel
.. autofunction:: splat_model.loadModel
.. autofunction:: splat_model.fitGrid
.. autofunction:: splat_model.fitMCMC



Evolutionary Model Routines
^^^^^^^^^^^^^^^^^^^^^^^^^^^

EUCLID Analysis Routines
^^^^^^^^^^^^^^^^^
.. autofunction:: splat_euclid.spexToEuclid
.. autofunction:: splat_euclid.addEuclidNoise


Other Calculation Routines
^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: splat.test
.. autofunction:: splat.weightMeanVar
.. autofunction:: splat_model.distributionStats


BibTeX Routines
^^^^^^^^^^^^^^^

.. autofunction:: splat_db.getBibTex


I/O Routines
^^^^^^^^^^^^

.. autofunction:: splat_db.checkOnline
.. autofunction:: splat_db.checkAccess
.. autofunction:: splat_db.checkFile
.. autofunction:: splat_db.checkLocalFile





*Search*


* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

