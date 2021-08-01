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

Data and Database Access
^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: splat.core.getSpectrum
.. autofunction:: splat.core.getStandard
.. autofunction:: splat.core.searchLibrary
.. autofunction:: splat.core.keySource
.. autofunction:: splat.core.keySpectrum
.. autofunction:: splat.core.initiateStandards
.. autofunction:: splat.core.addUserData
.. autofunction:: splat.core.Spectrum.export
.. autofunction:: splat.core.Spectrum.save
.. autofunction:: splat.core.readSpectrum
.. autofunction:: splat.utilities.checkInstrument
.. autofunction:: splat.utilities.checkLocation
.. autofunction:: splat.utilities.checkTelescope
.. autofunction:: splat.utilities.checkKeys
.. autofunction:: splat.core.Spectrum.updateSourceInfor



Spectral Manipulation
^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: splat.core.Spectrum.toFluxUnit
.. autofunction:: splat.core.Spectrum.toWaveUnit
.. autofunction:: splat.core.Spectrum.toFnu
.. autofunction:: splat.core.Spectrum.toFlam
.. autofunction:: splat.core.Spectrum.toSED
.. autofunction:: splat.core.Spectrum.toSurface
.. autofunction:: splat.core.Spectrum.toBrightnessTemperature
.. autofunction:: splat.core.Spectrum.toTemperature
.. autofunction:: splat.core.Spectrum.mapto
.. autofunction:: splat.core.Spectrum.toWavelengths
.. autofunction:: splat.core.Spectrum.toInstrument
.. autofunction:: splat.core.Spectrum.redden
.. autofunction:: splat.core.Spectrum.normalize
.. autofunction:: splat.core.Spectrum.scale
.. autofunction:: splat.core.Spectrum.fluxCalibrate
.. autofunction:: splat.core.Spectrum.smooth
.. autofunction:: splat.core.Spectrum._smoothToResolution
.. autofunction:: splat.core.Spectrum._smoothToSlitWidth
.. autofunction:: splat.core.Spectrum._smoothToSlitPixelWidth
.. autofunction:: splat.core.Spectrum.maskFlux
.. autofunction:: splat.core.Spectrum.remove
.. autofunction:: splat.core.Spectrum.trim
.. autofunction:: splat.core.stitch
.. autofunction:: splat.core.Spectrum.addNoise



Spectral Classification and Characterization
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: splat.core.compareSpectra
.. autofunction:: splat.core.measureIndex
.. autofunction:: splat.core.measureIndexSet
.. autofunction:: splat.plot.visualizeIndices
.. autofunction:: splat.core.measureEW
.. autofunction:: splat.core.measureEWElement
.. autofunction:: splat.core.measureEWSet
.. autofunction:: splat.core.classifyByIndex
.. autofunction:: splat.core.classifyByStandard
.. autofunction:: splat.core.classifyByTemplate
.. autofunction:: splat.core.classifyGravity
.. autofunction:: splat.core.metallicity
.. autofunction:: splat.empirical.redden
.. autofunction:: splat.core.generateMask


Spectrophotometry
^^^^^^^^^^^^^^^^^
.. autofunction:: splat.photometry.checkFilter
.. autofunction:: splat.photometry.filterInfo
.. autofunction:: splat.photometry.filterMag
.. autofunction:: splat.core.Spectrum.filterMag
.. autofunction:: splat.photometry.filterProfile
.. autofunction:: splat.photometry.filterProperties
.. autofunction:: splat.photometry.magToFlux
.. autofunction:: splat.photometry.visualizeFilter


High-resolution Spectroscopic Functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: splat.utilities.baryVel
.. autofunction:: splat.utilities.lsfRotation
.. autofunction:: splat.core.Spectrum.rvShift
.. autofunction:: splat.core.Spectrum.broaden
.. autofunction:: splat.core.Spectrum.rotate


Empirical Relationships
^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: splat.empirical.estimateDistance
.. autofunction:: splat.empirical.redden
.. autofunction:: splat.empirical.typeToColor
.. autofunction:: splat.empirical.typeToMag
.. autofunction:: splat.empirical.typeToTeff
.. autofunction:: splat.empirical.typeToLuminosity
.. autofunction:: splat.empirical.typeToBC
.. autofunction:: splat.utilities.checkAbsMag
.. autofunction:: splat.utilities.checkBC
.. autofunction:: splat.utilities.checkLbol



Conversion Routines
^^^^^^^^^^^^^^^^^^^
.. autofunction:: splat.utilities.caldateToDate
.. autofunction:: splat.utilities.dateToCaldate
.. autofunction:: splat.utilities.coordinateToDesignation
.. autofunction:: splat.utilities.designationToCoordinate
.. autofunction:: splat.utilities.designationToCoordinateString
.. autofunction:: splat.utilities.designationToShortName
.. autofunction:: splat.utilities.typeToNum
.. autofunction:: splat.utilities.properCoordinates
.. autofunction:: splat.utilities.properDate
.. autofunction:: splat.utilities.UVW
.. autofunction:: splat.utilities.XYZ



Plotting Routines
^^^^^^^^^^^^^^^^^
.. autofunction:: splat.plot.plotMap
.. autofunction:: splat.plot.plotSpectrum
.. autofunction:: splat.plot.plotBatch
.. autofunction:: splat.plot.plotSequence
.. autofunction:: splat.plot.plotIndices


Utility Routines
^^^^^^^^^^^^^^^^
.. autofunction:: splat.utilities.distributionStats
.. autofunction:: splat.utilities.integralResample
.. autofunction:: splat.utilities.reMap
.. autofunction:: splat.core.readSpectrum
.. autofunction:: splat.utilities.weightedMeanVar
.. autofunction:: splat.utilities.isNumber
.. autofunction:: splat.utilities.isUnit
.. autofunction:: splat.utilities.numberList
.. autofunction:: splat.utilities.randomSphereAngles


I/O Routines
^^^^^^^^^^^^
.. autofunction:: splat.utilities.checkOnline
.. autofunction:: splat.utilities.checkOnlineFile
.. autofunction:: splat.utilities.checkAccess
.. autofunction:: splat.utilities.checkFile
.. autofunction:: splat.utilities.checkLocal


Modeling and Simulation Routines
--------------------------------

Spectral Modeling Routines
^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: splat.model.getModel
.. autofunction:: splat.model.loadModel
.. autofunction:: splat.model.loadOriginalModel
.. autofunction:: splat.model.loadTelluric
.. autofunction:: splat.model.blackbody
.. autofunction:: splat.model.modelFitGrid
.. autofunction:: splat.model.modelFitMCMC
.. autofunction:: splat.model.modelFitEMCEE
.. autofunction:: splat.model.makeForwardModel
.. autofunction:: splat.model.mcmcForwardModelFit
.. autofunction:: splat.model.mcmcForwardModelReport
.. autofunction:: splat.model.addUserModels
.. autofunction:: splat.model.processModelsToInstrument
.. autofunction:: splat.utilities.checkSpectralModelName
.. autofunction:: splat.model.info


Evolutionary Model Routines
^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: splat.utilities.checkEvolutionaryModelName
.. autofunction:: splat.evolve.loadEvolModel
.. autofunction:: splat.evolve.modelParameters
.. autofunction:: splat.evolve.plotModelParameters


Population Simulation Routines
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: splat.evolve.simulateAges
.. autofunction:: splat.evolve.simulateMasses
.. autofunction:: splat.evolve.simulateMassRatios
.. autofunction:: splat.evolve.simulateDistances
.. autofunction:: splat.evolve.simulateUVW
.. autofunction:: splat.evolve.simulatePopulation
.. autofunction:: splat.simulate.galactic_density_juric
.. autofunction:: splat.simulate.UVWpopulation
.. autofunction:: splat.simulate.volumeCorrection
.. autofunction:: splat.simulate.galacticPotential

Specialty Packages
------------------

EUCLID Analysis Routines
^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: splat_euclid.spexToEuclid
.. autofunction:: splat_euclid.addEuclidNoise


BibTeX Routines
^^^^^^^^^^^^^^^
.. autofunction:: splat.citations.getBibTex
.. autofunction:: splat.citations.getBibTexOnline
.. autofunction:: splat.citations.veryShortRef
.. autofunction:: splat.citations.shortRef
.. autofunction:: splat.citations.longRef
.. autofunction:: splat.citations.veryLongRef
.. autofunction:: splat.citations.processBiblibrary
.. autofunction:: splat.citations.bibTexParser
.. autofunction:: splat.citations.citeURL


Astroquery-based Data Access
^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: splat.database.getPhotometry
.. autofunction:: splat.database.querySimbad
.. autofunction:: splat.database.queryVizier
.. autofunction:: splat.database.queryNist
.. autofunction:: splat.database.queryXMatch



*Search*


* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

