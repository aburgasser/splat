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
^^^^^^^^^^^^^^^^^
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


Utility Routines
^^^^^^^^^^^^^^^^
.. autofunction:: splat.test
.. autofunction:: splat.weightMeanVar
.. autofunction:: splat.distributionStats


I/O Routines
^^^^^^^^^^^^
.. autofunction:: splat_db.checkOnline
.. autofunction:: splat_db.checkAccess
.. autofunction:: splat_db.checkFile
.. autofunction:: splat_db.checkLocalFile


Modeling Routines
-----------------

Spectral Modeling Routines
^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: splat_model.loadModel
.. autofunction:: splat_model.getModel
.. autofunction:: splat_model.loadInterpolatedModel
.. autofunction:: splat_model.loadModelParameters
.. autofunction:: splat_model.modelFitGrid
.. autofunction:: splat_model.modelFitMCMC
.. autofunction:: splat_model.modelFitEMCEE
.. autofunction:: splat_model.modelFitEMCEE_lnlikelihood
.. autofunction:: splat_model.modelFitEMCEE_lnprior_limits
.. autofunction:: splat_model.modelFitEMCEE_lnprior_normal
.. autofunction:: splat_model.modelFitEMCEE_lnprob
.. autofunction:: splat_model.modelFitEMCEE_bestparameters
.. autofunction:: splat_model.modelFitEMCEE_plotchains
.. autofunction:: splat_model.modelFitEMCEE_plotcomparison
.. autofunction:: splat_model.modelFitEMCEE_plotbestcomparison
.. autofunction:: splat_model.modelFitEMCEE_plotcorner
.. autofunction:: splat_model.modelFitEMCEE_summary
.. autofunction:: splat_model.reportModelFitResults


Evolutionary Model Routines
^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: splat_evolve.loadEvolModel
.. autofunction:: splat_evolve.modelParameters
.. autofunction:: splat_evolve.modelParametersSingle
.. autofunction:: splat_evolve.plotModelParameters


Specialty Packages
------------------


EUCLID Analysis Routines
^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: splat_euclid.spexToEuclid
.. autofunction:: splat_euclid.addEuclidNoise



BibTeX Routines
^^^^^^^^^^^^^^^

.. autofunction:: splat_db.getBibTex







*Search*


* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

