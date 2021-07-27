from __future__ import print_function, division

"""
.. note::
         These are the plotting functions for the SPLAT code 
"""


# imports: internal
import copy
import glob
import os

# imports: external
from astropy.coordinates import Angle,SkyCoord,Galactic,BarycentricTrueEcliptic      # coordinate conversion
import astropy.units as u
import matplotlib
#matplotlib.use('agg')
import matplotlib.cm as cm
import matplotlib.colors as colmap
import matplotlib.patches as patches
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy
from scipy.interpolate import interp1d 
from scipy import ndimage

# imports: splat
from splat.initialize import *
from splat.utilities import *
import splat.core as splat


def plotMap(*args,**kwargs):
    '''
    :Purpose: Plot coordinates onto an equatorial map grid

    :Input
    One or more coordinates, specified as astropy.coordinates.SkyCoord objects, two-element arrays 
    (RA and declination in decimal degrees), or string coordinate designation; the latter two are 
    converted to SkyCoord variables using `splat.properCoordinates()`

    :Parameters:
    projection = 'mollweide'
        projection type; see the `Basemap documentation <https://matplotlib.org/basemap/users/mapsetup.html>`_
    reference_declination or rdec
        a single or array of declinations to indicates as references, useful for denoting pointing limits

    galactic = False
        Plot in galactic coordinates (NOT YET IMPLEMENTED)
    ecliptic = False
        Plot in ecliptic coordinates (NOT YET IMPLEMENTED)
    extrgalactic or supergalactic = False
        Plot in supergalactic coordinates (NOT YET IMPLEMENTED)

    marker or markers or symbol or symbols = 'o'
        symbol type; see the `matplotlib documentation <https://matplotlib.org/api/markers_api.html>`_
    size or sizes or symsize or symsizes = 10
        symbol size
    color or colors = 'k'
        color of plot symbols
    alpha or alphas = 0.5
        transparency of plot symbols

    legend, legends, label or labels = None
        list of strings providing legend-style labels for each spectrum plotted
    grid = True:
        add a grid

    file or filename or output:
        filename or filename base for output
    figsize = (8,6)
        set the figure size; set to default size if not indicated
    tight = True:
        makes a ``tight'' box plot to elimate extra whitespace
        
    :Example: 
       >>> import splat
       >>> import splat.plot as splot
       >>> s = splat.searchLibrary(young=True)
       >>> c = [splat.properCoordinates(x) for x in s['DESIGNATION']]
       >>> splot.plotMap(c)
    '''
    
    projection = kwargs.get('projection','mollweide')
    figsize = kwargs.get('figsize',(8,6))
    colors = kwargs.get('colors',['k' for i in range(len(args))])
    colors = kwargs.get('color',colors)
    alphas = kwargs.get('alphas',[0.5 for i in range(len(args))])
    alphas = kwargs.get('alpha',alphas)
    markers = kwargs.get('symbols',['o' for i in range(len(args))])
    markers = kwargs.get('symbol',markers)
    markers = kwargs.get('markers',markers)
    markers = kwargs.get('marker',markers)
    symsizes = kwargs.get('size',[10 for i in range(len(args))])
    symsizes = kwargs.get('sizes',symsizes)
    symsizes = kwargs.get('symsize',symsizes)
    symsizes = kwargs.get('symsizes',symsizes)
    rdec = kwargs.get('reference_declination',None)
    rdec = kwargs.get('rdec',rdec)
    file = kwargs.get('file','')
    file = kwargs.get('output',file)
    file = kwargs.get('filename',file)
    fontsize = kwargs.get('fontsize',14)
    fontsize = kwargs.get('charsize',fontsize)

# process markers
#    if isinstance(symsizes,float) or isinstance(symsizes,int):
#        symsizes = [symsizes for i in range(len(args))]
    if not isinstance(colors,list): colors = [colors]
    if not isinstance(alphas,list): alphas = [alphas]
    if not isinstance(symsizes,list): symsizes = [symsizes]
    if not isinstance(markers,list): markers = [markers]
    while len(colors) < len(args): colors.append(colors[-1])
    while len(alphas) < len(args): alphas.append(alphas[-1])
    while len(symsizes) < len(args): symsizes.append(symsizes[-1])
    while len(markers) < len(args): markers.append(markers[-1])


    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(111, projection=projection)
    ax.set_xticklabels(['14h','16h','18h','20h','22h','0h','2h','4h','6h','8h','10h'],fontsize=fontsize)
    ax.set_yticklabels([r'-75$\degree$',r'-60$\degree$',r'-45$\degree$',r'-30$\degree$',r'-15$\degree$','0$\degree$','15$\degree$','30$\degree$','45$\degree$','60$\degree$','75$\degree$'],fontsize=fontsize)
    ax.grid(kwargs.get('grid',True))

    for i,pcoords in enumerate(args):
# sense what the input values are
        if not isinstance(pcoords,list):
            pcoords = [pcoords]
        if isinstance(pcoords[0],float) and len(pcoords) == 2:
            pcoords = [pcoords]
        if isinstance(pcoords[0],str):
            pcoords = [splat.properCoordinates(c) for c in pcoords]
        if isinstance(pcoords[0],list):
            if isinstance(pcoords[0][0],float):
                pcoords = [splat.properCoordinates(c) for c in pcoords]
        if not isinstance(pcoords[0],SkyCoord):
            raise ValueError('\nCould not parse coordinate input {}'.format(pcoords))

# convert coordinates and scatter plot
        ra = [c.ra for c in pcoords]
        ra = [c.wrap_at(180*u.degree) for c in ra]
        ra = [c.radian for c in ra]
        dec = [c.dec.radian for c in pcoords]

        p = ax.scatter(ra, dec,color=colors[i],alpha=alphas[i],s=symsizes[i],marker=markers[i])

# declination reference
    if rdec != None:
        raref = Angle(numpy.arange(-180,180.,1.)*u.degree)
        raref.wrap_at(180.*u.degree)
        if not isinstance(rdec,list): rdec = [rdec]
        for d in rdec:
            pref = ax.plot(raref.radian,Angle([d]*len(raref)*u.degree).radian,'k--')

# plot galactic plane
    if kwargs.get('galactic',False) != False:
        lng = Angle(numpy.arange(-180.,180.,1.)*u.degree)
        lat = Angle(numpy.zeros(len(lng))*u.degree)
        s = SkyCoord(Galactic,l=lng,b=lat)
        mra = s.icrs.ra
        mra = [c.wrap_at(180*u.degree) for c in mra]
        mra = [c.radian for c in mra]
        mdec = s.icrs.dec.radian
        mdecp = [x for y, x in sorted(zip(mra, mdec))]
        mrap = sorted(mra)
        p = ax.plot(mrap, mdecp,color='k',alpha=1,ls='--')

# plot ecliptic plane
    if kwargs.get('ecliptic',False) != False:
        lng = Angle(numpy.arange(-180.,180.,1.)*u.degree)
        lat = Angle(numpy.zeros(len(lng))*u.degree)
        s = SkyCoord(BarycentricTrueEcliptic,lon=lng,lat=lat)
        mra = s.icrs.ra
        mra = [c.wrap_at(180*u.degree) for c in mra]
        mra = [c.radian for c in mra]
        mdec = s.icrs.dec.radian
        mdecp = [x for y, x in sorted(zip(mra, mdec))]
        mrap = sorted(mra)
        p = ax.plot(mrap, mdecp,color='k',alpha=1,ls=':')

# plot legend
    if kwargs.get('legend',None) != None:
        plt.legend(kwargs['legend'],bbox_to_anchor=(1, 1),bbox_transform=plt.gcf().transFigure,fontsize=fontsize)

# plot legend
    if file != None:
        if kwargs.get('tight',True) == True:
            plt.savefig(file, bbox_inches='tight')
        else:
            plt.savefig(file)
    
    return fig 



def plotSpectrum(*args, **kwargs):
    '''
    Purpose
    -------
        Primary plotting program for splat Spectrum objects.

    Parameters
    ----------
    input(s) : Spectrum objects, either sequentially, in a list, or in list of lists
        These are the spectra to be plotted; the input is flexible:
            * `Spec1, Spec2, ...`: plot multiple spectra together, or separately if multiplot = True
            * `[Spec1, Spec2, ...]`: plot multiple spectra together, or separately if multiplot = True
            * `[[Spec1, Spec2], [Spec3, Spec4], ..]`: plot multiple sets of spectra (multiplot forced to be True)
    
    xrange = [0.85,2.42]: list of 2 floats or unitted astropy quantities (optional)
        plot range for wavelength axis
    
    yrange = [-0.02,1.2]*fluxMax: : list of 2 floats or unitted astropy quantities (optional)
        plot range for flux axis
    
    xlabel = wave.unit: string (optional)
        wavelength axis label; by default set by wlabel keywords and wave.unit in first spectrum object
    
    ylabel = flux.unit: string (optional)
        flux axis label; by default set by fscale, flabel and flux.unit in first spectrum object
    
    xlog, ylog = False : boolean (optional)
        set the x (wavelength) or y (flux) axis to plot as a log scale
    
    features = [] : list of strings (optional)
        a list of strings indicating chemical features to label on the spectra
        options include H2O, CH4, CO, TiO, VO, FeH, H2, HI, KI, NaI, SB (for spectral binary)
    
    mdwarf, ldwarf, tdwarf, young, binary = False : boolean (optional)
        add in pre-defined features characteristic of these classes

    telluric = False : boolean (optional)
        indicate telluric absorption features

    band(s) = [] : list of 2-element float arrays
        a single or array of 2-element arrays that indicate bands that you want to specifically shade in
    
    bandcolor(s) = 'k' : single or array of strings
        a single or array of colors to shade the bands (default = 'k')
    
    bandalpha(s) = 0.2 : single or array of floats
        a single or array of alphas to shade the bands (default = 0.2)
    
    bandlabel(s) = '' : single or array of strings
        a single or array of labels to annotate the bands (default = '')
    
    bandlabelposition(s) = 'bottom' : single or array of strings
        a single or array of strings indicating the position of the labels; 
        can be 'bottom', 'middle', or 'top' (default = 'bottom')
    
    legend = [] : array of strings
        list of strings providing legend-style labels for each spectrum plotted
        Alternate notation: legends, label, labels
    
    legendLocation = 'upper right' : string
        place of legend; options are 'upper left', 'center middle', 'lower right' (variations thereof) 
        and 'outside'
        Alternate notation: labelLocation 
    
    legendfontscale = 1 : float
        sets the scale factor for the legend fontsize (defaults to fontscale)
    
    grid = False : boolean
        set to True to add a grid

    stack = 0:
        set to a numerical offset to stack spectra on top of each other
    
    zeropoint = [0,...]:
        list of offsets for each spectrum, giving finer control than stack
    
    showZero = True:
        plot the zeropoint(s) of the spectra

    comparison:
        a comparison Spectrum to compare in each plot, useful for common reference standard
    
    noise, showNoise or uncertainty = False:
        plot the uncertainty for each spectrum
    
    residual = False:
        plots the residual between two spectra
    
    colorComparison:
        color of comparison source plot lines; by default all grey

    color or colors:
        color of plot lines; by default all black
    
    colorUnc or colorsUnc:
        color of uncertainty lines; by default same as line color but reduced opacity
    
    colorScheme or colorMap:
        color map to apply based on matplotlib colormaps; 
        see http://matplotlib.org/api/pyplot_summary.html?highlight=colormaps#matplotlib.pyplot.colormaps
    
    linestyle or linestyles:
        line style of plot lines; by default all solid
    
    fontscale = 1:
        sets a scale factor for the fontsize

    inset = False:
        place an inset panel showing a close up region of the spectral data
    
    inset_xrange = False:
        wavelength range for inset panel
    
    inset_yrange = False:
        flux range for inset panel (otherwise set by flux)
    
    inset_position = [0.65,0.60,0.20,0.20]
        position of inset planet in normalized units, in order left, bottom, width, height
    
    inset_features = False
        list of features to label in inset plot

    file or filename or output:
        filename or filename base for output
    
    filetype = 'pdf':
        output filetype, generally determined from filename
    
    multiplot = False: 
        creates multiple plots, depending on format of input (optional)
    
    multipage = False: 
        spreads plots across multiple pages; output file format must be PDF
        if not set and plots span multiple pages, these pages are output sequentially as separate files
    
    layout or multilayout = [1,1]:
        defines how multiple plots are laid out on a page
    
    figsize:
        set the figure size; set to default size if not indicated
    
    interactive = False:
        if plotting to window, set this to make window interactive
    
    tight = True:
        makes a ``tight'' box plot to elimate extra whitespace

        
    :Example 1: A simple view of a random spectrum
       >>> import splat
       >>> import splat.plot as splot
       >>> spc = splat.getSpectrum(spt = 'T5', lucky=True)[0]
       >>> spc.plot()                       # this automatically generates a "quicklook" plot
       >>> splot.plotSpectrum(spc)          # does the same thing
       >>> splot.plotSpectrum(spc,uncertainty=True,tdwarf=True)     # show the spectrum uncertainty and T dwarf absorption features

    :Example 2: Viewing a set of spectra for a given object
        In this case we'll look at all of the spectra of TWA 30B in the library, sorted by year and compared to the first epoch data
        This is an example of using multiplot and multipage

       >>> import splat
       >>> import splat.plot as splot
       >>> splist = splat.getSpectrum(name = 'TWA 30B')         # get all spectra of TWA 30B
       >>> junk = [sp.normalize() for sp in splist]             # normalize the spectra
       >>> dates = [sp.date for sp in splist]                   # observation dates
       >>> spsort = [s for (s,d) in sorted(zip(dates,splis))]   # sort spectra by dates
       >>> dates.sort()                                         # don't forget to sort dates!
       >>> splot.plotSpectrum(spsort,multiplot=True,layout=[2,2],multipage=True,\   # here's our plot statement
           comparison=spsort[0],uncertainty=True,mdwarf=True,telluric=True,legends=dates,\
           legendLocation='lower left',output='TWA30B.pdf')
       
    :Example 3: Display the spectra sequence of L dwarfs
        This example uses the list of standard files contained in SPLAT, and illustrates the stack feature

       >>> import splat
       >>> import splat.plot as splot
       >>> spt = [splat.typeToNum(i+20) for i in range(10)] # generate list of L spectral types
       >>> splat.initiateStandards()                        # initiate standards
       >>> splist = [splat.SPEX_STDS[s] for s in spt]       # extact just L dwarfs
       >>> junk = [sp.normalize() for sp in splist]         # normalize the spectra
       >>> labels = [sp.shortname for sp in splist]         # set labels to be names
       >>> splot.plotSpectrum(splist,figsize=[10,20],labels=labels,stack=0.5,\  # here's our plot statement
           colorScheme='copper',legendLocation='outside',telluric=True,output='lstandards.pdf')
       
    '''

# keyword parameters
#    from .splat import Spectrum
    nsamples = kwargs.get('nsamples',1000)
    multiplot = kwargs.get('multiplot',False)           # create multiple plots
    multipage = kwargs.get('multipage',False)           # create a multiple page sequence of plots
    multilayout = kwargs.get('multilayout',[1,1])       # layout of multiple plots, [# horizontal, # vertical]
    multilayout = kwargs.get('layout',multilayout)      
    stack = kwargs.get('stack',0)                   # stack spectra on top of each other
    grid = kwargs.get('grid',False)                 # plot internal grid lines
    filename = kwargs.get('filename','')            # output filename
    filename = kwargs.get('file',filename)
    filename = kwargs.get('output',filename)
    if filename == False:
        filename = ''
    title = kwargs.get('title','')
    fontscale = kwargs.get('fontscale',1)
    legendfontscale = kwargs.get('legendfontscale',0.8*fontscale)
    fb = filename.split('.')[:-1]               # filebase for multiple files
    if filename != '' and len(fb) > 1:
        filebase = fb[0]
        for x in fb[1:]:
            filebase+='.'+x
    elif filename != '' and len(fb) == 1:
        filebase = fb[0]
    else:
        filebase = filename
    filetype = kwargs.get('format',filename.split('.')[-1])
    filetype.lower()
    if filetype == '':
        filetype = 'pdf'
    comparison = kwargs.get('comparison',False)
    if not isinstance(comparison,splat.Spectrum):
        comparison = False
    residual = kwargs.get('residual',False)
    inset = kwargs.get('inset',False)
    inset_xrange = kwargs.get('inset_xrange',[])
    inset_yrange = kwargs.get('inset_yrange',[])
    inset_position = kwargs.get('inset_position',[.65, .6, .2, .2])
#    inset_color = kwargs.get('inset_color','k')
    inset_features = kwargs.get('inset_features',False)
    

#    mask = kwargs.get('mask',False)                # not yet implemented

# features to label on spectra
    feature_labels = { \
        'h2o': {'label': r'H$_2$O', 'type': 'band', 'wavelengths': [[0.925,0.95],[1.08,1.20],[1.325,1.550],[1.72,2.14]]}, \
        'ch4': {'label': r'CH$_4$', 'type': 'band', 'wavelengths': [[1.1,1.24],[1.28,1.44],[1.6,1.76],[2.2,2.35]]}, \
        'co': {'label': r'CO', 'type': 'band', 'wavelengths': [[2.29,2.39]]}, \
        'tio': {'label': r'TiO', 'type': 'band', 'wavelengths': [[0.6569,0.6852],[0.705,0.727],[0.76,0.80],[0.825,0.831],[0.845,0.86]]}, \
#        'tio': {'label': r'TiO', 'type': 'band', 'wavelengths': [[0.76,0.80],[0.825,0.831]]}, \
        'vo': {'label': r'VO', 'type': 'band', 'wavelengths': [[1.04,1.08]]}, \
        'vo': {'label': r'VO', 'type': 'band', 'wavelengths': [[1.04,1.08]]}, \
        'young vo': {'label': r'VO', 'type': 'band', 'wavelengths': [[1.17,1.20]]}, \
#        'feh': {'label': r'FeH', 'type': 'band', 'wavelengths': [[0.86,0.90],[0.98,1.03],[1.19,1.25],[1.57,1.64]]}, \
        'cah': {'label': r'CaH', 'type': 'band', 'wavelengths': [[0.6346,0.639],[0.675,0.705]]}, \
        'crh': {'label': r'CrH', 'type': 'band', 'wavelengths': [[0.8611,0.8681]]}, \
        'feh': {'label': r'FeH', 'type': 'band', 'wavelengths': [[0.8692,0.875],[0.98,1.03],[1.19,1.25],[1.57,1.64]]}, \
        'h2': {'label': r'H$_2$', 'type': 'band', 'wavelengths': [[1.5,2.4]]}, \
        'sb': {'label': r'*', 'type': 'band', 'wavelengths': [[1.6,1.64]]}, \
        'h': {'label': r'H I', 'type': 'line', 'wavelengths': [[1.004,1.005],[1.093,1.094],[1.281,1.282],[1.944,1.945],[2.166,2.166]]},\
        'hi': {'label': r'H I', 'type': 'line', 'wavelengths': [[1.004,1.005],[1.093,1.094],[1.281,1.282],[1.944,1.945],[2.166,2.166]]},\
        'h1': {'label': r'H I', 'type': 'line', 'wavelengths': [[1.004,1.005],[1.093,1.094],[1.281,1.282],[1.944,1.945],[2.166,2.166]]},\
        'na': {'label': r'Na I', 'type': 'line', 'wavelengths': [[0.8186,0.8195],[1.136,1.137],[2.206,2.209]]}, \
        'nai': {'label': r'Na I', 'type': 'line', 'wavelengths': [[0.8186,0.8195],[1.136,1.137],[2.206,2.209]]}, \
        'na1': {'label': r'Na I', 'type': 'line', 'wavelengths': [[0.8186,0.8195],[1.136,1.137],[2.206,2.209]]}, \
        'cs': {'label': r'Cs I', 'type': 'line', 'wavelengths': [[0.8521,0.8521],[0.8943,0.8943]]}, \
        'csi': {'label': r'Cs I', 'type': 'line', 'wavelengths': [[0.8521,0.8521],[0.8943,0.8943]]}, \
        'cs1': {'label': r'Cs I', 'type': 'line', 'wavelengths': [[0.8521,0.8521],[0.8943,0.8943]]}, \
        'rb': {'label': r'Rb I', 'type': 'line', 'wavelengths': [[0.78,0.78],[0.7948,0.7948]]}, \
        'rbi': {'label': r'Rb I', 'type': 'line', 'wavelengths': [[0.78,0.78],[0.7948,0.7948]]}, \
        'rb1': {'label': r'Rb I', 'type': 'line', 'wavelengths': [[0.78,0.78],[0.7948,0.7948]]}, \
        'mg': {'label': r'Mg I', 'type': 'line', 'wavelengths': [[1.7113336,1.7113336],[1.5745017,1.5770150],[1.4881595,1.4881847,1.5029098,1.5044356],[1.1831422,1.2086969],]}, \
        'mgi': {'label': r'Mg I', 'type': 'line', 'wavelengths': [[1.7113336,1.7113336],[1.5745017,1.5770150],[1.4881595,1.4881847,1.5029098,1.5044356],[1.1831422,1.2086969],]}, \
        'mg1': {'label': r'Mg I', 'type': 'line', 'wavelengths': [[1.7113336,1.7113336],[1.5745017,1.5770150],[1.4881595,1.4881847,1.5029098,1.5044356],[1.1831422,1.2086969],]}, \
        'ca': {'label': r'Ca I', 'type': 'line', 'wavelengths': [[0.6573,0.6573],[2.263110,2.265741],[1.978219,1.985852,1.986764],[1.931447,1.945830,1.951105]]}, \
        'cai': {'label': r'Ca I', 'type': 'line', 'wavelengths': [[0.6573,0.6573],[2.263110,2.265741],[1.978219,1.985852,1.986764],[1.931447,1.945830,1.951105]]}, \
        'ca1': {'label': r'Ca I', 'type': 'line', 'wavelengths': [[0.6573,0.6573],[2.263110,2.265741],[1.978219,1.985852,1.986764],[1.931447,1.945830,1.951105]]}, \
        'caii': {'label': r'Ca II', 'type': 'line', 'wavelengths': [[1.184224,1.195301],[0.985746,0.993409]]}, \
        'ca2': {'label': r'Ca II', 'type': 'line', 'wavelengths': [[1.184224,1.195301],[0.985746,0.993409]]}, \
        'al': {'label': r'Al I', 'type': 'line', 'wavelengths': [[1.672351,1.675511],[1.3127006,1.3154345]]}, \
        'ali': {'label': r'Al I', 'type': 'line', 'wavelengths': [[1.672351,1.675511],[1.3127006,1.3154345]]}, \
        'al1': {'label': r'Al I', 'type': 'line', 'wavelengths': [[1.672351,1.675511],[1.3127006,1.3154345]]}, \
        'fe': {'label': r'Fe I', 'type': 'line', 'wavelengths': [[1.5081407,1.5494570],[1.25604314,1.28832892],[1.14254467,1.15967616,1.16107501,1.16414462,1.16931726,1.18860965,1.18873357,1.19763233]]}, \
        'fei': {'label': r'Fe I', 'type': 'line', 'wavelengths': [[1.5081407,1.5494570],[1.25604314,1.28832892],[1.14254467,1.15967616,1.16107501,1.16414462,1.16931726,1.18860965,1.18873357,1.19763233]]}, \
        'fe1': {'label': r'Fe I', 'type': 'line', 'wavelengths': [[1.5081407,1.5494570],[1.25604314,1.28832892],[1.14254467,1.15967616,1.16107501,1.16414462,1.16931726,1.18860965,1.18873357,1.19763233]]}, \
        'k': {'label': r'K I', 'type': 'line', 'wavelengths': [[0.7699,0.7665],[1.169,1.177],[1.244,1.252]]}, \
        'ki': {'label': r'K I', 'type': 'line', 'wavelengths': [[0.7699,0.7665],[1.169,1.177],[1.244,1.252]]}, \
        'k1': {'label': r'K I', 'type': 'line', 'wavelengths': [[0.7699,0.7665],[1.169,1.177],[1.244,1.252]]}}

    features = kwargs.get('features',[])
    if not isinstance(features,list):
        features = [features]
    if (kwargs.get('ldwarf',False) or kwargs.get('mdwarf',False)):
        features.extend(['k','na','feh','tio','co','h2o','h2'])
    if (kwargs.get('tdwarf',False)):
        features.extend(['k','ch4','h2o','h2'])
    if (kwargs.get('young',False)):
        features.extend(['vo'])
    if (kwargs.get('binary',False)):
        features.extend(['sb'])
# clean repeats while maintaining order - set does not do this
    fea = []
    for i in features:
        if i not in fea:
            fea.append(i)
    features = fea

# error check - make sure you're plotting something
    if (len(args) < 1):
        print('plotSpectrum needs at least one Spectrum object to plot')
        return

# if a list is passed, use this list
    elif (len(args) == 1 and isinstance(args[0],list)):
        splist = args[0]
    
# if a set of objects is passed, turn into a list
    else:
        splist = []
        for a in args:
            if isinstance(a,splat.Spectrum):      # a spectrum object
                splist.append(a)
            elif isinstance(a,list):
                splist.append(a)
            else:
                print('\nplotSpectrum: Ignoring input object {} as it is neither a Spectrum object nor a list\n\n'.format(a))

# set up for multiplot
    if (len(splist) == 1):
        multiplot = False
    
# array of lists => force multiplot
    elif (len(splist) > 1 and isinstance(splist[0],list)):
        multiplot = True

# reformat array of spectra of multiplot is used (i.e., user forgot to set)
    if multiplot == True and isinstance(splist[0],splat.Spectrum):
        splist = [[s] for s in splist]

    elif multiplot == False and isinstance(splist[0],splat.Spectrum):
        splist = [splist]
        
# flatten array if multiplot is not set
    elif multiplot == False and isinstance(splist[0],list) and len(splist) > 1:
        splist = [[item for sublist in splist for item in sublist]]       # flatten

    tot_sp = len([item for sublist in splist for item in sublist])    # Total number of spectra
    
# prep legend
    legend = kwargs.get('legend',[str() for x in numpy.arange(tot_sp)])
    legend = kwargs.get('legends',legend)
    legend = kwargs.get('label',legend)
    legend = kwargs.get('labels',legend)
    if not isinstance(legend,list):
        legend = [legend]
    if(len(legend) < tot_sp):
        legend.extend([str() for x in numpy.arange(tot_sp-len(legend))])
    legendLocation = kwargs.get('legendLocation','upper right')       # sets legend location
    legendLocation = kwargs.get('labelLocation',legendLocation)       # sets legend location

# now run a loop through the input subarrays
    plt.close('all')

# set up here for multiple file output
    nplot = 1
    if multipage == True or multiplot == True:
        nplot = multilayout[0]*multilayout[1]
        numpages = int(len(splist) / nplot) + 1
        if (len(splist) % nplot == 0):
               numpages -= 1
        fig = []
        
    if multipage == True and filetype == 'pdf':
        pdf_pages = PdfPages(filename)
        
    if multipage == False:
        if len(splist) > 1:
            files = [filebase+'{}.'.format(i+1)+filetype for i in numpy.arange(len(splist))]
        else:
            files = [filebase+'.'+filetype]

    pg_n = 0        # page counter
    plt_n = 0       # plot per page counter
    lg_n = 0        # legend per plot counter
    for plts,sp in enumerate(splist):
# set specific plot parameters
        if not isinstance(sp[0],splat.Spectrum):
            raise ValueError('\nInput to plotSpectrum has wrong format:\n\n{}\n\n'.format(sp[0]))
        zeropoint = kwargs.get('zeropoint',numpy.zeros(len(sp)))

# settings that work if the spectrum was read in as legitmate Spectrum object
        try:
            xlabel = kwargs.get('xlabel','{} ({})'.format(sp[0].wlabel,sp[0].wave.unit))
            ylabel = kwargs.get('ylabel','{} {} ({})'.format(sp[0].fscale,sp[0].flabel,sp[0].flux.unit))
        except:
            xlabel = kwargs.get('xlabel','Wavelength')
            ylabel = kwargs.get('ylabel','Flux')
        xrng = kwargs.get('xrange',[numpy.nanmin(sp[0].wave.value),numpy.nanmax(sp[0].wave.value)])
        if isUnit(xrng[0]): xrng = [x.value for x in xrng]
        bound = []
        bound.extend(xrng)
#        ymax = [s.fluxMax().value for s in sp]
        ymax = [numpy.nanquantile(s.flux.value,0.98) for s in sp]
        yrng = kwargs.get('yrange',numpy.array([-0.02,1.2])*numpy.nanmax(ymax)+numpy.nanmax(zeropoint))
        if isUnit(yrng[0]): yrng = [x.value for x in yrng]
        bound.extend(yrng)
        linestyle = kwargs.get('linestyle',['-' for x in numpy.arange(len(sp))])
        linestyle = kwargs.get('linestyles',linestyle)
        if (len(linestyle) < len(sp)):
            linestyle.extend(['-' for x in numpy.arange(len(sp)-len(linestyle))])

# colors
# by default all black lines
        colors = kwargs.get('colors',['k' for x in numpy.arange(len(sp))])
        colors = kwargs.get('color',colors)
        if not isinstance(colors,list):
            colors = [colors]
        if (len(colors) < len(sp)):
            while len(colors) < len(sp):
                colors.append(colors[-1])
        colorScheme = kwargs.get('colorScheme',None)
        colorScheme = kwargs.get('colorMap',colorScheme)
        if (colorScheme != None):
            values = numpy.arange(len(sp))
            color_map = plt.get_cmap(colorScheme)
            norm  = colmap.Normalize(vmin=0, vmax=1.0*values[-1])
            scalarMap = cm.ScalarMappable(norm=norm, cmap=color_map)
            for i in numpy.arange(len(sp)):
                colors[i] = scalarMap.to_rgba(values[i])
        colorsUnc = kwargs.get('colorsUnc',colors)
        colorsUnc = kwargs.get('colorUnc',colorsUnc)
        if (len(colorsUnc) < len(sp)):
            while len(colorsUnc) < len(sp):
                colorsUnc.append(colors[-1])
        linewidths = kwargs.get('linewidths',1.5)
        linewidths = kwargs.get('linewidth',linewidths)
        linewidths = kwargs.get('lw',linewidths)
        if not isinstance(linewidths,list):
            linewidths = [linewidths]
        if (len(linewidths) < len(sp)):
            while len(linewidths) < len(sp):
                linewidths.append(linewidths[-1])


# show uncertainties
        showNoise = kwargs.get('showNoise',[False for x in numpy.arange(len(sp))])
        showNoise = kwargs.get('noise',showNoise)
        showNoise = kwargs.get('uncertainty',showNoise)
        if not isinstance(showNoise, list):
            showNoise = [showNoise]
        if (len(showNoise) < len(sp)):
            showNoise.extend([True for x in numpy.arange(len(sp)-len(showNoise))])

# zero points - by default true
        showZero = kwargs.get('showZero',[True for x in numpy.arange(len(sp))])
        if not isinstance(showZero, list):
            showZero = [showZero]
        if (len(showZero) < len(sp)):
            while len(showZero) < len(sp):
                showZero.extend(showZero[-1])

# GENERATE PLOTS
        if (multiplot == True or multipage == True):
            plt_n = plts % nplot
            if (plt_n == 0):# and plts != len(splist)):
#                ax = range(nplot)
#                t = tuple([tuple([i+b*multilayout[1] for i in range(multilayout[1])]) for b in range(multilayout[0])])
#                fig[pg_n], ax = plt.subplots(multilayout[0], multilayout[1], sharex = True, sharey = True)
#
# NOTE THE FOLLOWING LINE IS HAVING PROBLEMS IN PYTHON3
#
#
                fig.append(plt.figure())
                pg_n += 1
            ax = fig[pg_n-1].add_subplot(multilayout[0], multilayout[1], plt_n+1)
            
# plotting a single plot with all spectra
        else:
            plt.close('all')
#            ax = range(1)
            plt_n = 0
            fig = []
            if (kwargs.get('figsize') != None):
                fig.append(plt.figure(figsize = kwargs.get('figsize')))
            else:
                fig.append(plt.figure())
            ax = fig[0].add_subplot(111)
        
        for ii, a in enumerate(sp):
            flx = [i+zeropoint[ii] for i in a.flux.value]
#stack
            if stack > 0:
                flx = [f + (len(sp)-ii-1)*stack for f in flx]
#                if kwargs.get('yrange') == None:
#                    bound[3] = bound[3] + stack

            ax.plot(a.wave.value,flx,color=colors[ii],linestyle=linestyle[ii], lw=linewidths[ii], zorder = 10, label = legend[lg_n])  

# add comparison
            if comparison != False:
                colorComparison = kwargs.get('colorComparison',colors[0])
                linestyleComparison = kwargs.get('linestyleComparison',linestyle[0])
                linewidthComparison = kwargs.get('linewidthComparison',linewidths[0])
                cflx = [i+zeropoint[ii] for i in comparison.flux.value]

                if stack > 0:
                    cflx = [f + (len(sp)-ii-1)*stack for f in cflx]

                ax.plot(comparison.wave.value,cflx,color=colorComparison,linestyle=linestyleComparison, lw=linewidthComparison, alpha=0.5, zorder = 10)
    
# add residual
            if residual == True and len(sp) == 2:
                # Save flux values from first spectrum
                if ii == 0:
                    flx0 = [f - (len(sp)-ii-1)*stack for f in flx]
                    
                # Subtract fluxes and plot
                elif ii == 1:
                    res = [flx0[f_n] - f for f_n, f in enumerate(flx)]
                    ax.plot(a.wave.value, res, alpha = 0.3, color = 'g')
                    
                    # Fix bound[2] if res goes below 0
                    if min(res) < 0:
                        b0 = numpy.argmax(a.wave.value > bound[0])
                        b1 = numpy.argmin(a.wave.value < bound[1])
                        bound[2] = bound[2] + min(res[b0:b1])

# noise
            if (showNoise[ii]):
                ns = [i+zeropoint[ii] for i in a.noise.value]
                ax.plot(a.wave.value,ns,color=colorsUnc[ii],linestyle=linestyle[ii],alpha=0.3, zorder = 10)


# zeropoint
            if (showZero[ii]):
                ze = numpy.ones(len(a.flux))*zeropoint[ii]
                ax.plot(a.wave.value,ze,color=colors[ii],linestyle=':',alpha=0.3, zorder = 10)

# determine maximum flux for all spectra
            f = interp1d(a.wave,flx,bounds_error=False,fill_value=0.)
            if (ii == 0):
                wvmax = numpy.arange(bound[0],bound[1],0.001)
                flxmax = f(wvmax)
            else:
                flxmax = numpy.maximum(flxmax,f(wvmax))

# legend counter
            lg_n = lg_n + 1 # Increment lg_n


# label features
# THIS NEEDS TO BE FIXED WITH GRETEL'S STUFF
        yoff = 0.02*(bound[3]-bound[2])
        dbx = 0.01*(bound[1]-bound[0])
        fontsize = 10-numpy.min([(multilayout[0]*multilayout[1]-1),6])
        for ftr in features:
            ftr = ftr.lower()
            if ftr in feature_labels:
                for ii,waveRng in enumerate(feature_labels[ftr]['wavelengths']):
                    wRng = ((waveRng*u.micron).to(sp[0].wave.unit)).value
                    if (numpy.min(wRng) > bound[0] and numpy.max(wRng) < bound[1]):
                        dwrng = numpy.nanmax(wRng)-numpy.nanmin(wRng)
#                        if dwrng == 0.: wRng = [wRng[0]-0.01*(bound[1]-bound[0]),wRng[0]+0.01*(bound[1]-bound[0])]
                        if dwrng < dbx: dwrng = dbx
                        x = (numpy.arange(0,nsamples+1.0)/nsamples)*dwrng*1.1+numpy.nanmin(wRng)-dwrng*0.05
#                        print(ftr,wRng,numpy.nanmin(x),numpy.nanmax(x))
                        wfeature = numpy.where(numpy.logical_and(wvmax >= x[0],wvmax <= x[-1]))
                        f = interp1d(numpy.array(wvmax)[wfeature],numpy.array(flxmax)[wfeature],bounds_error=False,fill_value=0.)
                        y = numpy.nanmax(f(x))+1.*yoff

                        if feature_labels[ftr]['type'] == 'band':
                            ax.plot(wRng,[y+yoff]*2,color='k',linestyle='-')
                            ax.plot([wRng[0]]*2,[y,y+yoff],color='k',linestyle='-')
                            ax.text(numpy.mean(wRng),y+1.5*yoff,feature_labels[ftr]['label'],horizontalalignment='center',fontsize=fontsize)
                        else:
                            for w in wRng: ax.plot([w]*2,[y,y+yoff],color='k',linestyle='-')
                            ax.text(numpy.mean(wRng),y+1.5*yoff,feature_labels[ftr]['label'],horizontalalignment='center',fontsize=fontsize)
#                            wRng = [wRng[0]-0.02,wRng[1]+0.02]   # for overlap
#                        print(ftr,y,y+yoff,len(wfeature))

                        foff = numpy.zeros(len(flxmax))
                        foff[wfeature] = 3.*yoff
                        flxmax = list(numpy.array(flxmax)+foff)
        bound[3] = numpy.nanmax([numpy.nanmax(flxmax)+2.*yoff,bound[3]])


# grid
        if (grid):
            ax.grid()            

# axis labels 
        fontsize = (numpy.round(numpy.max([13./((multilayout[0]*multilayout[1])**0.33),5]))) * fontscale        # Added in fontscale
#        print(fontsize)
        legendfontsize = (13-numpy.min([(multilayout[0]*multilayout[1]-1),8])) * legendfontscale        # Added in fontscale
        ax.set_xlabel(xlabel, fontsize = fontsize)
        ax.set_ylabel(ylabel, fontsize = fontsize)
        ax.tick_params(axis='x', labelsize=fontsize)
        ax.tick_params(axis='y', labelsize=fontsize)

# log scale?
        if kwargs.get('xlog',False):
            ax.set_xscale('log',nonposx='clip')
        if kwargs.get('ylog',False):
            ax.set_yscale('log',nonposy='clip')

# place legend
        if len(legend) > 0:
            if legendLocation == 'outside':
                box = ax.get_position()
                ax.set_position([box.x0, box.y0 + box.height * 0.15, box.width * 0.7, box.height * 0.7])
                ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), prop={'size':legendfontsize})
            else:
                ax.legend(loc=legendLocation, prop={'size':legendfontsize})
                bound[3] = bound[3]+0.1*(bound[3]-bound[2])     # extend axis for in-plot legends

# overplot telluric absorption
        if (kwargs.get('telluric',False) == True):
            twv = [[1.1,1.2],[1.3,1.5],[1.75,2.0]]
            for waveRng in twv:
                wR = ((waveRng*u.micron).to(sp[0].wave.unit)).value
                rect = patches.Rectangle((wR[0],bound[2]),wR[1]-wR[0],bound[3]-bound[2],facecolor='0.95',alpha=0.2,color='0.95')
                ax.add_patch(rect)
                ax.text(numpy.mean(wR),bound[2]+3*yoff,r'$\oplus$',horizontalalignment='center',fontsize=fontsize)

# overplot color swaths for pre-specified bands
        bands = kwargs.get('bands',[])
        bands = kwargs.get('band',bands)
        if len(bands) > 0:
            if len(bands) == 2 and isinstance(bands[0],list) == False: bands = [bands]
            bandcolors = kwargs.get('bandcolors',[])
            bandcolor = kwargs.get('bandcolor',bandcolors)
            if not isinstance(bandcolors,list): bandcolors = [bandcolors]*len(bands)
            if len(bandcolors) < len(bands): 
                for i in range(len(bands)-len(bandcolors)): bandcolors.append('k')
            bandalphas = kwargs.get('bandalphas',[])
            bandalpha = kwargs.get('bandalpha',bandalphas)
            if not isinstance(bandalphas,list): bandalphas = [bandalphas]*len(bands)
            if len(bandalphas) < len(bands): 
                for i in range(len(bands)-len(bandalphas)): bandalphas.append(0.2)
            bandlabels = kwargs.get('bandlabels',[])
            bandlabel = kwargs.get('bandlabel',bandlabels)
            if not isinstance(bandlabels,list): bandlabels = [bandlabels]*len(bands)
            if len(bandlabels) < len(bands): 
                for i in range(len(bands)-len(bandlabels)): bandlabels.append('')
            bandlabelpositions = kwargs.get('bandlabelpositions',[])
            bandlabelposition = kwargs.get('bandlabelposition',bandlabelpositions)
            if not isinstance(bandlabelpositions,list): bandlabelpositions = [bandlabelpositions]*len(bands)
            if len(bandlabelpositions) < len(bands): 
                for i in range(len(bands)-len(bandlabelpositions)): bandlabelpositions.append('bottom')
            bandwidth = kwargs.get('bandwidth',0.1)
            for i,b in enumerate(bands):
                if not isinstance(b,list): 
                    try:
                        b = [float(b)-0.5*bandwidth,float(b)+0.5*bandwidth]
                    except:
                        print('\nWarning: plotSpectrum bands variables should be array of 2-element arrays; you passed {}'.format(bands))
                        b = [0.,0.]
                rect = patches.Rectangle((b[0],bound[2]),b[1]-b[0],bound[3]-bound[2],facecolor=bandcolors[i],color=bandcolors[i],alpha=bandalphas[i])
                ax.add_patch(rect)
                if bandlabelpositions[i].lower() == 'top':
                    ax.text(numpy.mean(b),bound[3]-3*yoff,bandlabels[i],horizontalalignment='center',fontsize=fontsize)
                elif bandlabelpositions[i].lower() == 'middle':
                    ax.text(numpy.mean(b),0.5*(bound[2]+bound[3]),bandlabels[i],horizontalalignment='center',fontsize=fontsize)
                else:
                    ax.text(numpy.mean(b),bound[2]+3*yoff,bandlabels[i],horizontalalignment='center',fontsize=fontsize)

# place inset - RIGHT NOW ONLY SETTING LIMITS WITH FIRST SPECTRUM IN LIST
        if not isinstance(inset,bool):
            inset_xrange = kwargs.get('inset_xrange',inset)
            inset = True
        if inset == True and len(inset_xrange) == 2:
            ax_inset = fig[pg_n-1].add_axes(inset_position) #, axisbg='white')
            bound2 = inset_xrange
            if len(inset_yrange) == 0:
                b0 = numpy.argmax(sp[0].wave.value > bound2[0])
                b1 = numpy.argmin(sp[0].wave.value < bound2[1])
                inset_yrange = [min(sp[0].flux.value[b0:b1]),max(sp[0].flux.value[b0:b1])]
            bound2.extend(inset_yrange)
            db = (bound2[3]-bound2[2])
            bound2[2] = bound2[2]-0.05*db
            bound2[3] = bound2[3]+0.05*db
            ax_inset.axis(bound2)
            inset_fontsize = fontsize*0.7

            for ii,a in enumerate(sp):
                flx = [i+zeropoint[ii] for i in a.flux.value]
                ax_inset.plot(a.wave.value,flx,color=colors[ii],linestyle=linestyle[ii])            
                ax_inset.set_xlabel('')
                ax_inset.set_ylabel('')
                ax_inset.tick_params(axis='x', labelsize=inset_fontsize)
                ax_inset.tick_params(axis='y', labelsize=inset_fontsize)
#                ax_inset.legend()

# inset feature labels
            if inset_features != False:
                f = interp1d(sp[0].wave,flx,bounds_error=False,fill_value=0.)
                wvmax = numpy.arange(bound2[0],bound2[1],0.001)
                flxmax = f(wvmax)
                yoff = 0.05*(bound2[3]-bound2[2])
                for ftr in inset_features:
                    ftr = ftr.lower()
                    if ftr in feature_labels:
                        for ii,waveRng in enumerate(feature_labels[ftr]['wavelengths']):
                            wR = ((waveRng*u.micron).to(sp[0].wave.unit)).value
                            if (numpy.min(wR) > bound2[0] and numpy.max(wR) < bound2[1]):
                                x = (numpy.arange(0,nsamples+1.0)/nsamples)* \
                                    (numpy.nanmax(wR)-numpy.nanmin(wR)+0.04)+numpy.nanmin(wR)-0.02
                                f = interp1d(wvmax,flxmax,bounds_error=False,fill_value=0.)
                                y = numpy.nanmax(f(x))+0.5*yoff
        
                                if feature_labels[ftr]['type'] == 'band':
                                    ax_inset.plot(wR,[y+yoff]*2,color='k',linestyle='-')
                                    ax_inset.plot([wR[0]]*2,[y,y+yoff],color='k',linestyle='-')
                                    ax_inset.text(numpy.mean(wR),y+2*yoff,feature_labels[ftr]['label'],horizontalalignment='center',fontsize=inset_fontsize)
                                else:
                                    for w in waveRng:
                                        ax_inset.plot([w]*2,[y,y+yoff],color='k',linestyle='-')
                                    ax_inset.text(numpy.mean(wR),y+2*yoff,feature_labels[ftr]['label'],horizontalalignment='center',fontsize=inset_fontsize)
                                    waveRng = [wR[0]-0.02,wR[1]+0.02]   # for overlap
        
# update offset
                                foff = [y+3*yoff if (w >= wR[0] and w <= wR[1]) else 0 for w in wvmax]
                                flxmax = [numpy.max([xx,yy]) for xx, yy in zip(flxmax, foff)]
                bound2[3] = numpy.max([bound2[3],numpy.max(flxmax)+5.*yoff])
                ax_inset.axis(bound2)

# finalize bounding
#        if kwargs.get('xrange',None) != None:
        bound[0:2] = xrng
#        if kwargs.get('yrange',None) != None:
        bound[2:4] = yrng
        for i,b in enumerate(bound):
            if isUnit(b): bound[i]=b.value
        ax.axis(bound)

    
# save to file or display
        if multipage == False:
            if filebase != '' and (plts % nplot == 3 or plts == len(splist)-1):
                if kwargs.get('tight',True) == True: 
                    plt.savefig(files[plts], format=filetype, bbox_inches='tight')
                else:
                    plt.savefig(files[plts], format=filetype)
    if filename == '' and not kwargs.get('web',False):
        plt.show()
        if (kwargs.get('interactive',False) != False):
            plt.ion()        # make window interactive 
        else:
            plt.ioff()


# save figures in multipage format and write off pdf file
    if (multipage == True):    
        for pg_n in numpy.arange(numpages):
#            fig[pg_n].text(0.5, 0.04, xlabel, ha = 'center', va = 'center')
#            fig[pg_n].text(0.06, 0.5, ylabel, ha = 'center', va = 'center', rotation = 'vertical')
            fig[pg_n].tight_layout
            fig[pg_n].suptitle(title, fontsize = 14, fontweight = 'bold')
            pdf_pages.savefig(fig[pg_n])
        if filetype == 'pdf':
            pdf_pages.close()

    return fig


def plotBatch(inp,output='spectra_plot.pdf',comparisons=None,classify=False,normalize=False,basecolors=['k','r'],legend=[],fontscale=0.7,layout=[2,2],normrange=[0.9,1.4],classify_kwargs={},plot_kwargs={},**kwargs):
    '''
    Purpose 
    -------
        Plots a batch of spectra into a 2x2 set of PDF files, with options of overplotting 
        comparison spectra, including best-match spectral standards. 

    Parameters 
    ----------
    input : Spectrum, string, or array array of these
        A single or list of Spectrum objects or filenames, or the glob search string 
        for a set of files (e.g., '/Data/myspectra/*.fits'). 

    output = 'spectra_plot.pdf' : string (optional), default = 'spectra_plot.pdf'
        Filename for PDF file output; full path should be included if not saving to local directory

    comparisons = None : Spectrum, string, or array array of these (optional) 
        list of Spectrum objects or filenames for comparison spectra. If comparisons list is shorter 
        than source list, then the last comparison source will be repeated. If the comparisons list is 
        longer, the list will be truncated. 
    
    classify = False : boolean (optional)
        Set to True to classify sources based on comparison to MLT spectral standards 
        following the method of `Kirkpatrick et al. (2010) <http://adsabs.harvard.edu/abs/2010ApJS..190..100K>`_. 
        This option normalizes the spectra by default

    normalize = False : boolean (optional)
        Set to True to normalize source and (if passed) comparison spectra.

    normrange = [0.9,1.4] : list (optional)
        Range in micron to compute the normalization

    layout = [2,2] : list (optional)
        Layout of plots on page (# per row x # per column)

    fontscale = 0.7 : float (optional)
        Set to list of legends for plots. The number of legends should equal the number of 

    legend = [] : list (optional)
        Set to list of legends for plots. The number of legends should equal the number of 
        sources and comparisons (if passed) in an alternating sequence. The default is to display
        the file name for each source and object name for comparison (if passed)

    classify_kwargs = {'method':'kirkpatrick'} : dictionary (optional)
        Dictionary of keywords for `splat.classifyByStandard()`_ , currently set to just use method='kirkpatrick'

    plot_kwargs = {} : dictionary (optional)
        Dictionary of additional keywords for `splat.plot.plotSpectrum()`_ 

    Outputs 
    -------
    matplotlib figure object showing individual panels of each source with optional comparison star

    Example
    -------
       >>> import glob, splat
       >>> import splat.plot as splot
       >>> files = glob.glob('/home/mydata/*.fits')
       >>> sp = splot.plotBatch(files,classify=True,output='comparison.pdf')
       >>> sp = splot.plotBatch('/home/mydata/*.fits',classify=True,output='comparison.pdf')
       >>> sp = splot.plotBatch([splat.Spectrum(file=f) for f in files],classify=True,output='comparison.pdf')
       
       All three of these commands produce the same result

    Dependencies
    ------------
    copy
    glob
    numpy
    os
    `splat.classifyByStandard()`_
    `splat.Spectrum()`_
    `splat.plot.plotSpectrum()`_ 

    .. _`splat.plot.plotSpectrum()` : api.html#splat.plot.plotSpectrum
    .. _`splat.classifyByStandard()` : api.html#splat.classifyByStandard
    .. _`splat.Spectrum()` : api.html#splat.Spectrum

    '''

# alt keyword check
    for k in ['file','filename']: 
        if kwargs.get(k,'') != '': output = kwargs[k]
    for k in ['legends','labels']: 
        if kwargs.get(k,'') != '': legend = kwargs[k]

# force input into a list
    if isinstance(inp,list): inputlist = copy.deepcopy(inp)
    else: inputlist = [inp]

# if input is a string of filenames, read in each file to a spectrum object
    if isinstance(inputlist[0],str):
# try a glob search string  
        files = glob.glob(os.path.normpath(inputlist[0]))
        if len(files) > 1 or (len(files) == 1 and inputlist[0].find('*') != -1):
            inputlist = files
# try reading in files into Spectrum object
        try:
            splist = [splat.Spectrum(file = f) for f in inputlist]
        except:
            raise ValueError('\nCould not read in list of files {} - make sure the full path is specified and the files are correctly formatted'.format(inputlist))

# if filenames, read in each file to a spectrum object
    elif isinstance(inputlist[0],splat.Spectrum):
        splist = copy.deepcopy(inputlist)
    else:
        raise ValueError('\nInput should be list of splat.Spectrum objects or filenames')

# normalize if desired
    if normalize==True:
        tmp = [sp.normalize(normrange) for sp in splist]

# comparison files are present
    complist = []
    if comparisons != None:
        comp = copy.deepcopy(comparisons)
        if not isinstance(comp,list): comp = [comp]
        if isinstance(comp[0],str):
            try:
                complist = [splat.Spectrum(file = f) for f in comp]
            except:
                print('\nCould not read in comparison files: ignoring comparisons')
        if isinstance(comp[0],splat.Spectrum):
            complist = comp
        if len(complist) < len(splist):
            while len(complist) < len(splist):
                complist.append(complist[-1])
# normalize
        if normalize==True:
            tmp = [sp.normalize(normrange) for sp in complist]

# set comparison files to be standards for spectral classification
# overrules input comparison sample
    if classify == True:
        complist = []
        base_kwargs={
        'return_standard': True,
        'method': 'kirkpatrick',
        }
        base_kwargs.update(classify_kwargs)
        for sp in splist:
            complist.append(splat.classifyByStandard(sp,**base_kwargs))

# prep for plotting
    plotlist = []
    clist = []
    for i,sp in enumerate(splist):
        if len(complist) == len(splist):
            plotlist.append([sp,complist[i]])
            clist.extend(basecolors)
        else:
            plotlist.append([sp])
            clist.extend(basecolors[0])

# manage legends
    if len(legend) != 0:
        if not isinstance(legend,list): legend = [legend]
        if len(legend) < (len(splist)+len(complist)):
# guess: just left out the comparison legends            
            if len(complist) > 0 and len(legend) == len(splist):
                legtmp = []
                for i,l in enumerate(legend):
                    legtmp.extend([l,'{}'.format(complist[i].name)])
                legend = legtmp
            else:
# otherwise: pad the remaining legends with the last legend (pairs)           
                while len(legend) < (len(splist)+len(complist)):
                    if len(complist)>0:
                        legend.extend([legend[-2],legend[-1]])
                    else:
                        legend.extend([legend[-1]])
        if len(legend) > (len(splist)+len(complist)):
            legend = legend[0:(len(splist)+len(complist))]
    else:
        legend = []
        for i,sp in enumerate(splist):
            l = []
            if 'name' in list(sp.__dict__.keys()): l.append(sp.name)
            else: l.append(os.path.basename(sp.filename))
            if len(complist)>0:
                if 'name' in list(complist[i].__dict__.keys()): l.append(complist[i].name)
                else: l.append(os.path.basename(complist[i].filename))
            legend.extend(l)

# generate plot
    base_kwargs={
    'multiplot': True,
    'multipage': True,
    'legends': legend,
    'colors': clist,
    'layout': layout,
    'fontscale': fontscale,
    'output': output,
    }
    base_kwargs.update(plot_kwargs)
    fig = plotSpectrum(plotlist,**base_kwargs)

    return fig




def visualizeIndices(sp,indices,**kwargs):
    '''
    :Purpose: ``Plot index-index plots.``

    indices should be a dictionary of {'index_name': {'ranges': [[],[]], 'value': #, 'unc': #}, ...}
    Not currently implemented
    '''

# check inputs
    if not isinstance(sp,splat.Spectrum): raise ValueError('\nFirst argument in visualizeIndices must be a Spectrum object')
    if not isinstance(indices, dict): raise ValueError("\nSecond argument in visualizeIndices must be a dictionary of the form {'index_name': {'ranges': [[],[]], 'value': #, 'unc': #}, ...}")
    file = kwargs.get('output',None)
    file = kwargs.get('file',file)
    file = kwargs.get('filename',file)

    mkwargs = copy.deepcopy(kwargs)
    spc = copy.deepcopy(sp)
    spc.normalize()

    if mkwargs.get('filename', False) != False: del mkwargs['filename']
    if mkwargs.get('file', False) != False: del mkwargs['file']
    if mkwargs.get('output', False) != False: del mkwargs['output']

    fig = plotSpectrum(sp,**mkwargs)
    for iname in list(indices.keys()):
        if isinstance(indices[iname],list):
            ranges = indices[iname]
        elif isinstance(indices[iname],dict):
            if 'ranges' in list(indices[iname].keys()):
                ranges = indices[iname]['ranges']
            else: raise ValueError('\nCannot find index wavelength ranges in indices variable: {}'.format(indices))
        else: raise ValueError('\nCannot find index wavelength ranges in indices variable: {}'.format(indices))
        if not isinstance(ranges[0],list):
            ranges = [ranges]
        vals = []
        for r in ranges:
            v,u = splat.measureIndex(sp,[r],method='single',sample='median')
            rect = patches.Rectangle((r[0],v-u),r[1]-r[0],2.*u,facecolor='0.95',alpha=0.2,color='0.95')
            ax.add_patch(rect)
            vals.append(v+u)
        if len(ranges) == 1: 
            ax.text(numpy.mean(ranges[0]),vals[0],r'{}\n'.format(iname),horizontalalignment='center')
        else:
            rflat = [x for y in ranges for x in y]
            rconn = []
            for i,r in enumerate(ranges):
                rconn.append(numpy.mean(r))
                plt.plot([rconn[-1]]*2,[vals[i],numpy.max(vals)+0.05*sp.fluxMax().value])
            plt.plot([numpy.min(rconn),numpy.max(rconn)],[numpy.max(vals)+0.05*sp.fluxMax().value]*2)
            plt.plot([numpy.mean(rconn)]*2,[numpy.max(vals)+0.05*sp.fluxMax().value,numpy.max(vals)+0.07*sp.fluxMax().value])
            ax.text(numpy.mean(rconn),numpy.max(vals)+0.07*sp.fluxMax().value,r'{}\n'.format(iname),horizontalalignment='center')

    if file != None:
        plt.savefig(file)

    return fig
    

def plotSED(*args, **kwargs):
    '''
    :Purpose: ``Plot SED photometry with SpeX spectrum.``

    Not currently implemented
    '''
    pass
    return



def plotSequence(spec,type_range=2, std_class='dwarf', spt='', output='', verbose=False, **kwargs):
    '''
    :Purpose: 

        Visualze compares a spectrum to a sequence of standards laid out vertically. The standards are chosen to be some number of
        type about the best-guess classification, either passed as a parameter or determined with `classifyByStandard()`_.
        More than one spectrum can be provided, in which case multiple plots are returned

    :Required Inputs: 

        :param: spec: A single or series of Spectrum objects or filenames, or the glob search string for a set of files to read in (e.g., '/Data/myspectra/*.fits'). At least one input must be provided

    :Optional Inputs: 

        :param: spt = '': Default spectral type for source; this input skips `classifyByStandard()`_
        :param: type_range = 2: Number of subtypes to consider above and below best-fit spectral type 
        :param: std_type = 'dwarf': Type of standards to compare to. Should be one of the following: 'dwarf' (default), 'sd', 'dsd', 'esd', 'vlg', 'intg'. These can also be defined by setting the equalivalent keyword parameter; e.g., plotSequence(spec,dwarf=True).
        :param: output = '': Filename for output; full path should be include if not saving to current directory. If blank, plot is shown on screen
        :param: verbose = False: Set to True to provide additional feedback

        In addition, relevant parameters for `plotSpectrum()`_ and `classifyByStandard()`_ may be provided

    :Outputs: 

        A matplotlib figure object containing the plot  of the spectrum(a) compared to sequence of standards on screen or saved to file

    :Example:
       >>> import splat
       >>> import splat.plot as splot
       >>> sp = splat.getSpectrum(lucky=True)[0]
       >>> fig = splat.plotSequence(sp,output='classify.pdf')

.. _`plotSpectrum()` : api.html#splat.plot.plotSpectrum
.. _`classifyByStandard()` : api.html#splat.classifyByStandard

    '''

# check inputs

#    from .splat import classifyByStandard, Spectrum
#    parameters = ['type_range']
#    checkKeys(kwargs,parameters,forcekey=False)
#    type_range =kwargs.get('type_range',2)

# some taste preferences
    kwargs['stack']=kwargs.get('stack',0.5)
#    kwargs['legendLocation']=kwargs.get('legendLocation','outside')
    kwargs['telluric']=kwargs.get('telluric',True)
    kwargs['method']=kwargs.get('method','kirkpatrick')
#    kwargs['color']=kwargs.get('color','r')
    kwargs['colorComparison']=kwargs.get('colorComparison','k')
    kwargs['colorScheme']=kwargs.get('colorScheme','winter')
    kwargs['figsize']=kwargs.get('figsize',[10,10*type_range*kwargs['stack']])
    kwargs['fontscale']=kwargs.get('fontscale',1.5)

# process input
    if isinstance(spec,str):
        if len(glob.glob(os.path.normpath(spec))) == 0:
            raise ValueError('\nCannot find input file {} - make sure full path is included'.format(spec))
        try:
            sp = splat.Spectrum(file = spec)
        except:
            raise ValueError('\nCould not read in file {} - make sure the file is correctly formatted'.format(spec))
    elif isinstance(spec,splat.Spectrum):
        sp = copy.deepcopy(spec)
    else:
        raise ValueError('\nInput should be a Spectrum object or filename')
    sp.normalize()

# choose the standard class set
    allowed_classes = ['dwarf','sd','dsd','esd','vlg','intg','subdwarf','lowg']
    for a in allowed_classes:
        if kwargs.get(a,False) == True: std_class = a
    std_class = std_class.lower()
    if std_class not in allowed_classes: 
        if verbose == True: print('\nDo not recognize class {}; defaulting to dwarf'.format(allowed_classes))
        std_class = 'dwarf'

    if verbose==True: print('Using {} class standards'.format(std_class))

    if std_class == 'dwarf': std_ref = splat.STDS_DWARF_SPEX
    elif std_class == 'sd' or std_class == 'subdwarf': std_ref = splat.STDS_SD_SPEX
    elif std_class == 'dsd': std_ref = splat.STDS_DSD_SPEX
    elif std_class == 'esd': std_ref = splat.STDS_ESD_SPEX
    elif std_class == 'vlg' or std_class == 'lowg': std_ref = splat.STDS_VLG_SPEX
    elif std_class == 'intg': std_ref = splat.STDS_INTG_SPEX
    else:
        raise ValueError('\nUnknown class type {}'.format(std_class))

# classify by comparison to standards
    spt = kwargs.get('spt',splat.classifyByStandard(sp,std_class=std_class,**kwargs)[0])
    if not isinstance(spt,str):
        spt = typeToNum(spt,subclass=std_class)

# produce range of standards for plot
    std_ref_sptn = [typeToNum(s) for s in std_ref]
    std_ref_sptn.sort()
    try:
        ref_ind = std_ref_sptn.index(typeToNum(spt))
    except:
        std_ref_diff = [numpy.absolute(s-splat.typeToNum('sdL0.0')) for s in std_ref_sptn]
        ref_ind = numpy.argmin(numpy.array(std_ref_diff))
    ref_range = [int(ref_ind-type_range),int(ref_ind+type_range)+1]
    if ref_range[0]<0: ref_range[0] = 0
    if ref_range[-1]>len(std_ref_sptn): ref_range[-1] = -1
    std_spts = [splat.typeToNum(s,subclass=std_class) for s in std_ref_sptn[ref_range[0]:ref_range[1]]]

    stds = []
    stdlabels = []
    for s in std_spts:
        stds.append(std_ref[s])
        if s == spt: stdlabels.append('{} Standard (Best)'.format(s))
        else: stdlabels.append('{} Standard'.format(s))
    if len(stds) == 0:
        raise ValueError('\nCould not find any standards between {} and {} in class {}; try a wide range or different standard class'.format(stdspt[0],stdspt[1],std_class))
    kwargs['yrange']=kwargs.get('yrange',[0,len(stds)*kwargs['stack']+1.])
    plotlist = []
#    labels = []
#    colors = []
    for i,stdsp in enumerate(stds):
        plotlist.append([stdsp,sp])
#        labels.extend([])
    fig = plotSpectrum(stds,comparison=sp,labels=stdlabels,**kwargs)

    return fig

# main testing of program
if __name__ == '__main__':
    out_folder = '/Users/adam/projects/splat/code/testing/'
    def test_plotBatch():
        data_folder = '/Users/adam/projects/splat/adddata/done/daniella/'
        files = glob.glob(os.path.normpath(data_folder+'*.fits'))
        plotBatch(files,classify=True,output=out_folder+'plotBatch_test1.pdf',telluric=True)
        plotBatch(data_folder+'*.fits',classify=True,output=out_folder+'plotBatch_test2.pdf',noise=True)
        splist = [splat.Spectrum(file=f) for f in files]
        plotBatch(splist,classify=True,output=out_folder+'plotBatch_test3.pdf',features=['h2o','feh','co'],legend=[s.name for s in splist])
        return

    def test_plotSequence():
        sp = splat.Spectrum(10001)
        plotSequence(sp,output=out_folder+'plotSequence_test1.pdf')
        plotSequence(sp,type_range=3,output=out_folder+'plotSequence_test2.png')
        data_folder = '/Users/adam/projects/splat/adddata/done/daniella/'
        files = glob.glob(os.path.normpath(data_folder+'*.fits'))
        plotSequence(files[0],telluric=True,stack=0.7,spt='M0',output=out_folder+'plotSequence_test3.eps')
        sp = splat.getSpectrum(shortname='0415-0935')[0]
        plotSequence(sp,telluric=True,stack=0.3,output=out_folder+'plotSequence_test4.eps')
        plotSequence(sp,telluric=True,stack=0.3)
        return

    test_plotBatch()
#    test_plotSequence()


    
