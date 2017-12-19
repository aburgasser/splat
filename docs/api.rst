API
===============================================

SPLAT Classes
-------------

Spectrum
^^^^^^^^
.. autoclass:: splat.core.Spectrum
	:members:


SPLAT Routines
--------------

Database Access
^^^^^^^^^^^^^^^
.. autofunction:: splat.utilities.checkFile
.. autofunction:: splat.utilities.checkAccess
.. autofunction:: splat.core.getSpectrum
.. autofunction:: splat.core.getStandard
.. autofunction:: splat.core.initiateStandards
.. autofunction:: splat.core.searchLibrary
.. autofunction:: splat.core.addUserData
.. autofunction:: splat.core.keySource
.. autofunction:: splat.core.keySpectrum


Spectral Comparison
^^^^^^^^^^^^^^^^^^^
.. autofunction:: splat.core.compareSpectra
.. autofunction:: splat.core.stitch
.. autofunction:: splat.core.generateMask


Spectral Classification and Characterization
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: splat.core.measureIndex
.. autofunction:: splat.core.measureIndexSet
.. autofunction:: splat.core.measureEW
.. autofunction:: splat.core.measureEWElement
.. autofunction:: splat.core.measureEWSet
.. autofunction:: splat.core.classifyByIndex
.. autofunction:: splat.core.classifyByStandard
.. autofunction:: splat.core.classifyByTemplate
.. autofunction:: splat.core.classifyGravity
.. autofunction:: splat.core.metallicity
.. autofunction:: splat.empirical.redden


Spectrophotometry
^^^^^^^^^^^^^^^^^
.. autofunction:: splat.photometry.checkFilter
.. autofunction:: splat.photometry.filterInfo
.. autofunction:: splat.photometry.filterMag
.. autofunction:: splat.photometry.filterProfile
.. autofunction:: splat.photometry.filterProperties
.. autofunction:: splat.photometry.magToFlux
.. autofunction:: splat.photometry.visualizeFilter


High-resolution Spectroscopic Functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: splat.utilities.baryVel
.. autofunction:: splat.utilities.lsfRotation


Empirical Relationships
^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: splat.empirical.estimateDistance
.. autofunction:: splat.empirical.typeToColor
.. autofunction:: splat.empirical.typeToMag
.. autofunction:: splat.empirical.typeToTeff
.. autofunction:: splat.empirical.typeToLuminosity
.. autofunction:: splat.empirical.typeToBC


Conversion Routines
^^^^^^^^^^^^^^^^^^^
.. autofunction:: splat.utilities.caldateToDate
.. autofunction:: splat.utilities.dateToCaldate
.. autofunction:: splat.utilities.coordinateToDesignation
.. autofunction:: splat.utilities.designationToCoordinate
.. autofunction:: splat.utilities.designationToShortName
.. autofunction:: splat.utilities.typeToNum
.. autofunction:: splat.utilities.properCoordinates
.. autofunction:: splat.utilities.properDate


Plotting Routines
^^^^^^^^^^^^^^^^^
.. autofunction:: splat.plot.plotMap
.. autofunction:: splat.plot.plotSpectrum
.. autofunction:: splat.plot.plotBatch
.. autofunction:: splat.plot.plotSequence
.. autofunction:: splat.plot.plotSED
.. autofunction:: splat.plot.plotIndices


Utility Routines
^^^^^^^^^^^^^^^^
.. autofunction:: splat.utilities.checkInstrument
.. autofunction:: splat.utilities.checkKeys
.. autofunction:: splat.utilities.distributionStats
.. autofunction:: splat.utilities.isNumber
.. autofunction:: splat.core.readSpectrum
.. autofunction:: splat.utilities.weightedMeanVar


I/O Routines
^^^^^^^^^^^^
.. autofunction:: splat.utilities.checkOnline
.. autofunction:: splat.utilities.checkOnlineFile
.. autofunction:: splat.utilities.checkAccess
.. autofunction:: splat.utilities.checkFile
.. autofunction:: splat.utilities.checkLocal


Modeling Routines
-----------------

Spectral Modeling Routines
^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: splat.utilities.checkSpectralModelName
.. autofunction:: splat.model.checkModelName
.. autofunction:: splat.model.loadModel
.. autofunction:: splat.model.getModel
.. autofunction:: splat.model.loadInterpolatedModel
.. autofunction:: splat.model.blackbody
.. autofunction:: splat.model.modelFitGrid
.. autofunction:: splat.model.modelFitMCMC
.. autofunction:: splat.model.modelFitEMCEE



Evolutionary Model Routines
^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: splat.utilities.checkEvolutionaryModelName
.. autofunction:: splat.evolve.loadEvolModel
.. autofunction:: splat.evolve.modelParameters
.. autofunction:: splat.evolve.plotModelParameters
.. autofunction:: splat.evolve.simulatePopulation
.. autofunction:: splat.evolve.simulateAges
.. autofunction:: splat.evolve.simulateMasses
.. autofunction:: splat.evolve.simulateMassRatios
.. autofunction:: splat.evolve.simulateSpatialDistribution
.. autofunction:: splat.evolve.simulateBinaryOrbits
.. autofunction:: splat.evolve.simulateGalacticOrbits
.. autofunction:: splat.evolve.simulateKinematics
.. autofunction:: splat.evolve.simulatePhotometry


Specialty Packages
------------------

EUCLID Analysis Routines
^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: splat_euclid.spexToEuclid
.. autofunction:: splat_euclid.addEuclidNoise


BibTeX Routines
^^^^^^^^^^^^^^^
.. autofunction:: splat.citations.bibTexParser
.. autofunction:: splat.citations.citeURL
.. autofunction:: splat.citations.getBibTex
.. autofunction:: splat.citations.getBibTexOnline
.. autofunction:: splat.citations.veryShortRef
.. autofunction:: splat.citations.shortRef
.. autofunction:: splat.citations.longRef
.. autofunction:: splat.citations.veryLongRef


Astroquery-based Data Access
^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: splat.database.getPhotometry
.. autofunction:: splat.database.querySimbad
.. autofunction:: splat.database.queryVizier



*Search*


* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

