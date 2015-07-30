.. SpeX Prism Library Analysis Toolkit documentation master file, created by
   sphinx-quickstart on Sat Jul 11 20:07:28 2015.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.



Known Bugs
===============================================================

.. toctree
   :maxdepth: 3

As of 29 July 2015

High priority 
---------------------------
Currently broken or doesn't exist and should:

 * Documentation incomplete
 * MCMC fitting code in splat_model.py currently not functional (temporary code in place)
 * ``searchLibrary`` needs to be ditched in exchange for SQL
 * calling ``searchLibrary(spt = sp, logic = “and”)`` returns the full database
 * no ``plotIndices`` function
 * no ``spectralBinary`` function
 * no tabular output of data
 * need a reddening function

Medium priority 
---------------------------
Annoying but functional:

 * ``plotSpectrum``: add inset plots
 * ``plotSpectrum``: if you set don't set "multipage" then layout doesn't do anything
 * ``plotSpectrum``: font scaling is not consistent for all devices - add a fontScale parameter
 * dividing two spectra cannot be plotted - not sure why
 * need a way to store/read in indices
 * evolutionary models may not be returning accurate values - this needs to be validated
 * if an error occurs or you stop a spectrum download midway, that (blank) file persists and causes problems - deleting in advance does not seem to work
 * not entirely sure EW program is returning correct values


Low priority 
---------------------------
Would be nice if these worked better, or some convenient things to add:

 * ``compareSpectra`` and ``plotSpectrum``: would be nice if this plotted the difference spectrum
 * ``plotSpectrum``: add in alternate feature labeling
 * ``plotSpectrum``: add more atomic/molecular features
 * add in subdwarf and gravity standards
 * add in mean templates (from Kelle Cruz)
 * model fitting that is grid style (i.e., not MCMC)
 * visualization for evolutionary model parameters




* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

