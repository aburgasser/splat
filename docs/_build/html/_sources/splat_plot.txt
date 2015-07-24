.. SpeX Prism Library Analysis Toolkit documentation master file, created by
   sphinx-quickstart on Sat Jul 11 20:07:28 2015.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

SPLAT Plotting Routines
===============================================================

.. toctree
   :maxdepth: 3


Examples
-----------

    **Example 1: A simple view of a random spectrum**
      This example shows various ways of displaying a random spectrum in the library
    
       >>> import splat
       >>> spc = splat.getSpectrum(spt = 'T5', lucky=True)[0]	# select random spectrum
       >>> spc.plot()                   				# this automatically generates a "quicklook" plot
       >>> splat.plotSpectrum(spc)      				# does the same thing
       >>> splat.plotSpectrum(spc,uncertainty=True,tdwarf=True)     # show the spectrum uncertainty and T dwarf absorption features
    
      The last plot should look like the following:
    .. image:: _images/plot_example1.png
      :width: 400
      :align: center
    
    **Example 2: Compare two spectra**
      Optimally scale and compare two spectra. 
    
       >>> import splat
       >>> spc = splat.getSpectrum(spt = 'T5', lucky=True)[0]	# select random spectrum
       >>> spc2 = splat.getSpectrum(spt = 'T4', lucky=True)[0]	# read in another random spectrum
       >>> comp = splat.compareSpectra(spc,spc2)	# compare spectra to get optimal scaling
       >>> spc2.scale(comp[1])			# apply optimal scaling
       >>> splat.plotSpectrum(spc,spc2,colors=['black','red'],labels=[spc.name,spc2.name])     # show the spectrum uncertainty and T dwarf absorption features
    
    .. image:: _images/plot_example2.png
      :width: 400
      :align: center
     
    
    **Example 3: Compare several spectra for a given object**
      In this case we'll look at all of the spectra of TWA 30B in the library, sorted by year and compare each to the first epoch data. This is an example of using both multiplot and multipage.
    
       >>> splist = splat.getSpectrum(name = 'TWA 30B')         # get all spectra of TWA 30B
       >>> junk = [sp.normalize() for sp in splist]             # normalize the spectra
       >>> dates = [sp.date for sp in splist]                   # observation dates
       >>> spsort = [s for (d,s) in sorted(zip(dates,splist))]   # sort spectra by dates
       >>> dates.sort()                                         # don't forget to sort dates!
       >>> splat.plotSpectrum(spsort,multiplot=True,layout=[2,2],multipage=True,\   # here's our plot statement
           comparison=spsort[0],uncertainty=True,mdwarf=True,telluric=True,legends=dates,\
           legendLocation='lower left',output='TWA30B.pdf')
    
      Here is the first page of the resulting 5 page pdf file
    .. image:: _images/plot_example3.png
      :width: 500
      :align: center
           
    **Example 4: Display the spectra sequence of L dwarfs**
            This example uses the list of standard files contained in SPLAT, and illustrates the stack feature
    
       >>> spt = [splat.typeToNum(i+20) for i in range(10)]     # generate list of L spectral types
       >>> files = [splat.spex_stdfiles[s] for s in spt]        # get the standard files
       >>> splist = [splat.Spectrum(f) for f in files]          # read in list of Spectrum objects
       >>> junk = [sp.normalize() for sp in splist]             # normalize the spectra
       >>> labels = [sp.shortname for sp in splist]             # set labels to be names
       >>> splat.plotSpectrum(splist,figsize=[10,20],labels=labels,stack=0.5,\  # here's our plot statement
           colorScheme='copper',legendLocation='outside',telluric=True,output='lstandards.pdf')
    
    .. image:: _images/plot_example4.png
      :width: 400
      :align: center


Routines
---------- 

.. autofunction:: splat.plotSpectrum


* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

