.. SpeX Prism Library Analysis Toolkit documentation master file, created by
   sphinx-quickstart on Sat Jul 11 20:07:28 2015.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.


Empirical Relations
===============================================================

.. toctree
   :maxdepth: 3


There are a number of empirical relations contained in the SPLAT code to help convert between common parameters. 
These are all contained in the `splat.empirical <splat_empirical.html>`_ module, which should be imported separately:
    
    >>> import splat.empirical as spemp

Examples include:

.. _`typeToColor()` : api.html#splat.empirical.typeToColor
.. _`typeToTeff()` : api.html#splat.empirical.typeToTeff
.. _`typeToLbol()` : api.html#splat.empirical.typeToLuminosity
.. _`typeToLuminosity()` : api.html#splat.empirical.typeToLuminosity
.. _`typeToBC()` : api.html#splat.empirical.typeToBC
.. _`typeToMag()` : api.html#splat.empirical.typeToMag
.. _`estimateDistance()` : api.html#splat.estimateDistance

* `typeToColor()`_
	Takes a spectral type and optionally a ``color`` (string) and returns the typical color of the source, currently using only the `Skryzpek et al. (2015) <http://adsabs.harvard.edu/abs/2015A%26A...574A..78S>`_ : relations.

    Example:

    >>> import splat
    >>> print splat.typeToColor('L3', 'J-K')
        (1.46, nan)
    >>> print splat.typeToMag('M5', 'i-z', ref = 'skrzypek', unc=0.5)
        (0.91, 0.57797809947624645)
    >>> print splat.typeToMag('M0', 'i-z', ref = 'skrzypek')
        Spectral type M0.0 is outside the range for reference set Skrzypek et al. (2015)
        (nan, nan)


* `typeToMag()`_
	Takes a spectral type and a filter, and returns absolute magnitude, using one of several empirical relationships; you can see a list by typing ``print(splat.SPT_ABSMAG_SETS.keys())``.

    Example:

    >>> import splat
    >>> print splat.typeToMag('L3', '2MASS J')
        (12.730064813273996, 0.4)
    >>> print splat.typeToMag(21, 'MKO K', ref = 'burgasser')
        (10.705292820099999, 0.26)
    >>> print splat.typeToMag(24, '2MASS J', ref = 'faherty')
        Invalid filter given for Abs Mag/SpT relation from Faherty et al. (2012)
        (nan, nan)
    >>> print splat.typeToMag('M0', '2MASS H', ref = 'dupuy')
        Spectral Type is out of range for Abs Mag/SpT relation from Dupuy & Liu (2012) Abs Mag/SpT relation
        (nan, nan)


* `typeToTeff()`_
	Returns an effective temperature (Teff) and its uncertainty for a given spectral type, using one of the following empirical relationships (set by ``ref`` keyword):

        - `Golimowski et al. (2004) <http://adsabs.harvard.edu/abs/2004AJ....127.3516G>`_ : Allowed spectral type range is M6 to T8  (``ref`` = 'golimowski')
        - `Looper et al. (2008) <http://adsabs.harvard.edu/abs/2008ApJ...685.1183L>`_ : Allowed spectral type range is L0 to T8  (``ref`` = 'looper')
        - `Stephens et al. (2009) <http://adsabs.harvard.edu/abs/2009ApJ...702..154S>`_ : Allowed spectral type range is M6 to T8 and uses alternate coefficients for L3 to T8.  (``ref`` = 'stephens')
        - `Marocco et al. (2013) <http://adsabs.harvard.edu/abs/2013AJ....146..161M>`_ : Allowed spectral type range is M7 to T8  (``ref`` = 'marocco')
        - `Filippazzo et al. (2015). <http://adsabs.harvard.edu/abs/2015ApJ...810..158F>`_ : Allowed spectral type range is M6 to T9 (``ref`` = 'filippazzo')

    Example:

    >>> import splat
    >>> print splat.typeToTeff(20)
        (2233.4796740905499, 100.00007874571999)
    >>> print splat.typeToTeff(20, unc = 0.3, ref = 'golimowski')
        (2305.7500497902788, 127.62548366132124)


* `estimateDistance()`_
	Takes the apparent magnitude in a given ``filter`` and either takes or determines the absolute magnitude from empirical relations, then uses the absolute magnitude/distance relation to estimate the distance to the object in parsecs. Returns estimated distance and uncertainty in parsecs. If given only a spectrum object, this routine will measure the apparent magnitude, classify the spectrum, estimate the absolute magnitude, and estimate the distance; any additional inputs such as ``mag`` (for apparent magnitude), ``spt`` (for spectral type), ``absmag`` (for absolute magnitude) and their uncertainties will reduce dependence on the Spectrum object. With all three parameters, this routine operates without a Spectrum object.

    Example:
    >>> import splat
    >>> sp = splat.getSpectrum(shortname='1555+0954')[0]
    Retrieving 2 files
    >>> splat.estimateDistance(sp)
        Please specify the filter used to determine the apparent magnitude
        (nan, nan)
    >>> splat.estimateDistance(sp, filter='2MASS J')
    	(212.20546914411625, 50.593458481040173)
    >>> sp.fluxCalibrate('2MASS J',12.83)
    >>> splat.estimateDistance(sp, filter='2MASS J')
    	(6.8967647325911665, 1.7439740983732679)
    >>> splat.estimateDistance(sp, filter='2MASS J', mag=12.83)
    	(6.4528658994336521, 1.6853855848066823)
    >>> splat.estimateDistance(sp, filter='2MASS J', mag=12.83, mag_e=0.03)
    	(6.1292809243737336, 1.4986946706101478)
    >>> splat.estimateDistance(sp, filter='2MASS J', mag=12.83, mag_e=0.03, spt='L5')
    	(6.9954039276140554, 1.1679437846129084)
    >>> splat.estimateDistance(filter='2MASS J', mag=12.83, mag_e=0.03, spt='L5', absmag=13.56, absmag_e = 0.2)
    	(7.1788501442275461, 0.74878521889450711)


Useful Program Constants
------------------------

``splat.SPT_BC_SETS``, ``splat.SPT_LBOL_SETS``, ``splat.SPT_ABSMAG_SETS``
    Dictionaries containing information on the empirical relations contained in the SPLAT code




* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

