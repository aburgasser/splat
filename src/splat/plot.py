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
import splat


def plotMap(*args,**kwargs):
    """
    :Purpose: Plot coordinates onto an equatorial map grid

    :Input:
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
    """
    
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
        s = SkyCoord(l=lng,b=lat,frame='galactic')
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
        s = SkyCoord(lon=lng,lat=lat,frame='barycentrictrueecliptic')
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



def plotSpectrum(inp,xrng=[],yrng=[],xlabel='',ylabel='',xlog=False,ylog=False,grid=False,
    legend=[],legend_location='upper right',fontscale=1,legend_fontscale=1,title='',
    color='k',colormap=None,linestyle='-',linewidth=1.5,alpha=1.,
    show_noise=True,color_noise='k',linestyle_noise='-',linewidth_noise=1.5,alpha_noise=0.5,
    comparison=None,color_comparison='grey',linestyle_comparison='-',linewidth_comparison=1.5,alpha_comparison=1,
    residual=False,color_residual='m',linestyle_residual='-',linewidth_residual=1.5,alpha_residual=0.5,
    telluric=False,color_telluric='grey',linestyle_telluric='-',linewidth_telluric=1.5,alpha_telluric=0.2,
    features=[],mdwarf=False,ldwarf=False,tdwarf=False,young=False,binary=False,nsamples=100,
    band=[],band_color='k',band_alpha=0.2,band_label='',band_label_position='bottom',band_width=0.1,
    show_zero=True,stack=0.,zeropoint=0.,color_zero='k',linestyle_zero=':',linewidth_zero=1.5,alpha_zero=0.3,
    inset=False,inset_xrange=[],inset_yrange=[],inset_position=[0.65,0.60,0.20,0.20],inset_features=False,
    output='',multiplot=False,multipage=False,layout=[1,1],figsize=[],tight=True,
    interactive=False,**kwargs):
    """
    Primary plotting program for splat Spectrum objects.

    Parameters
    ----------
    inp : Spectrum objects or array of Spectrum objects or nested array of Spectrum objects
        These are the spectra to be plotted; the input is flexible:
            * `Spec`: plot the spectrum Spec
            * `[Spec1, Spec2, ...]`: plot multiple spectra together, or separately if multiplot = True
            * `[[Spec1, Spec2], [Spec3, Spec4], ..]`: plot multiple sets of spectra (multiplot forced to be True)
    
    xrange: array of two floats or two unitted astropy quantities, default = [0.85,2.42]
        plot range for wavelength axis
    
    yrange : array of two floats or two unitted astropy quantities, default = [-0.02,1.2] times `fluxMax()`_
        plot range for flux axis
    
    xlabel : string, default = wave.unit
        wavelength axis label; by default set by wave.unit in first Spectrum object
    
    ylabel : string, default = flux.unit
        flux axis label; by default set by flux.unit in first Spectrum object
    
    xlog, ylog : bool, default = False
        set the x (wavelength) or y (flux) axis to plot as a log scale
    
    grid : bool, default = False
        set to True to add a grid

    telluric : bool, default = False
        indicate telluric absorption features

    features : array of strings, default = []
        A list of strings indicating chemical features to label on the spectra
        options include H2O, CH4, CO, TiO, VO, FeH, H2, HI, KI, NaI, SB (for spectral binary)
    
    mdwarf, ldwarf, tdwarf, young, binary : boolean, default = False
        Set to True to add pre-defined features characteristic of these classes

    legend : array of strings, default = [] 
        list of strings providing legend-style labels for each spectrum plotted
    
    legend_location : string, default = 'upper right'
        place of legend; options are 'upper left', 'center middle', 'lower right' (variations thereof) 
        and 'outside'
    
    legend_fontscale: float, default = 1 
        sets the scale factor for the legend fontsize (defaults to fontscale)
    
    band : array of two floats or array of arrays of two floatS, default = []
        a single or array of 2-element arrays that indicate bands that you want to specifically shade in
    
    band_color : string or array of strings, default = 'k'
        a single or array of colors to shade the bands
    
    band_alpha : float or array of floats, default = 0.2 
        a single or array of alphas to shade the bands (default = 0.2)
    
    band_label : string or array of strings, default = ''
        a single or array of labels to annotate the bands (default = '')
    
    band_label_position : string or array of strings, default = 'bottom'
        a single or array of strings indicating the position of the labels; 
        can be 'bottom', 'middle', or 'top' (default = 'bottom')
    
    stack : float, default = 0.
        numerical offset to stack spectra on top of each other
    
    zeropoint : float or array of floats, default = 0.
        list of offsets for each spectrum, giving finer control than stack
    
    show_zero : book, default = True
        plot the zeropoint(s) of the spectra

    comparison:
        a comparison Spectrum to compare in each plot, useful for common reference standard
    
    show_noise : boolean, default = False
        set to True to plot the uncertainty for each spectrum
    
    residual = False:
        plots the residual between two spectra
    
    color_comparison:
        color of comparison source plot lines; by default all grey

    color:
        color of plot lines; by default all black
    
    color_noise:
        color of uncertainty lines; by default same as line color but reduced opacity
    
    colormap:
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

    output:
        filename or filename base for output
        
    multiplot = False: 
        creates multiple plots, depending on format of input (optional)
    
    multipage = False: 
        spreads plots across multiple pages; output file format must be PDF
        if not set and plots span multiple pages, these pages are output sequentially as separate files
    
    layout = [1,1]:
        defines how multiple plots are laid out on a page
    
    figsize:
        set the figure size; set to default size if not indicated
    
    interactive = False:
        if plotting to window, set this to make window interactive
    
    tight = True:
        makes a ``tight'' box plot to elimate extra whitespace

    Returns
    -------
    array of matplotlib figure objects
        Each array element is one of the plots

    Examples
    --------
    Case 1: A simple view of a random spectrum.

       >>> import splat
       >>> import splat.plot as splot
       >>> spc = splat.getSpectrum(spt='T5',lucky=True)[0]  # grab one T5 spectrum
       >>> spc.plot()                                       # this automatically generates a "quicklook" plot
       >>> splot.plotSpectrum(spc)                          # does the same thing
       >>> splot.plotSpectrum(spc,show_noise=True,tdwarf=True)  # shows the spectrum uncertainty and T dwarf absorption features

    Case 2: Viewing a set of spectra for a given object
        In this case we'll look at all of the spectra of TWA 30B in the library, 
        sorted by year and compared to the first epoch data
        This is an example of using multiplot and multipage.

       >>> import splat
       >>> import splat.plot as splot
       >>> splist = splat.getSpectrum(name = 'TWA 30B')         # get all spectra of TWA 30B
       >>> junk = [sp.normalize([0.9,1.2]) for sp in splist]    # normalize the spectra
       >>> dates = [sp.observation_date for sp in splist]       # observation dates
       >>> spsort = [s for (d,s) in sorted(zip(dates,splist))]  # sort spectra by dates
       >>> dates.sort()                                         # don't forget to sort dates!
       >>> splot.plotSpectrum(spsort,multiplot=True,layout=[2,2],multipage=True,
       ... comparison=spsort[0],show_noise=True,mdwarf=True,telluric=True,legend=dates,
       ... yrange=[0,1.2],legend_location='lower left',color_comparison='m',alpha_comparison=0.5,
       ... output='TWA30B.pdf')

       
    Case 3: Display a spectra sequence of L dwarfs
        This example uses the list of standard L dwarf spectra contained in SPLAT, 
        and illustrates the stack feature.

       >>> import splat
       >>> import splat.plot as splot
       >>> spt = [splat.typeToNum(i+20) for i in range(10)]     # generate list of L spectral types
       >>> splat.initiateStandards()                            # initiate standards
       >>> splist = [splat.STDS_DWARF_SPEX[s] for s in spt]     # extact just L dwarfs
       >>> junk = [sp.normalize([0.9,1.3]) for sp in splist]    # normalize the spectra
       >>> legend = [sp.name for sp in splist]                  # set labels to be names
       >>> splot.plotSpectrum(splist,figsize=[10,20],legend=legend,stack=0.5,  # here's our plot statement
       ... colormap='copper',legend_location='outside',telluric=True,
       ... xrange=[0.8,2.45],output='lstandards.pdf')
    
    See Also
    --------
    plotBatch : plots an array of spectra into a grid
    plotSequence : plots a sequence of spectra vertically

    """

# keyword parameters (for backward compatability)
    for k in ['showZero','showzero']: show_zero=kwargs.get(k,show_zero)
    for k in ['showNoise','noise','uncertainty','shownoise','showuncertainty','show_uncertainty']: show_noise=kwargs.get(k,show_noise)

    for k in ['line_style','lineStyle','ls','linestyles','line_styles']: linestyle=kwargs.get(k,linestyle)
    for k in ['line_width','lineWidth','width','lw','linewidths','line_widths']: linewidth=kwargs.get(k,linewidth)
    for k in ['colors','colour','colours']: color=kwargs.get(k,color)
    for k in ['colorScheme','color_scheme','colorscheme','colorMap','color_map']: colormap=kwargs.get(k,colormap)

    for k in ['colornoise','colorNoise','colorUnc','coloruncertainty','color_uncertainty','colorUncertainty']: color_noise=kwargs.get(k,color_noise)
    for k in ['linestylenoise','line_style_noise','linestyleNoise']: linestyle_noise=kwargs.get(k,linestyle_noise)
    for k in ['linewidthnoise','linewidthNoise','line_width_noise']: linewidth_noise=kwargs.get(k,linewidth_noise)
    for k in ['alphanoise','alphaNoise']: alpha_noise=kwargs.get(k,alpha_noise)

    for k in ['colorzero','colorZero']: color_zero=kwargs.get(k,color_zero)
    for k in ['linestylezero','line_style_zero','linestyleZero']: linestyle_zero=kwargs.get(k,linestyle_zero)
    for k in ['linewidthzero','linewidthZero','line_width_zero']: linewidth_zero=kwargs.get(k,linewidth_zero)
    for k in ['alphazero','alphaZero']: alpha_zero=kwargs.get(k,alpha_zero)

    for k in ['colorcomparison','colorComparison']: color_comparison=kwargs.get(k,color_comparison)
    for k in ['linestyleComparison','line_style_comparison','linestylecomparison']: linestyle_comparison=kwargs.get(k,linestyle_comparison)
    for k in ['linewidthcomparison','linewidthComparison','line_width_comparison']: linewidth_comparison=kwargs.get(k,linewidth_comparison)
    for k in ['alphacomparison','alphaComparison']: alpha_comparison=kwargs.get(k,alpha_comparison)

    for k in ['colorresidual','colorResidual']: color_residual=kwargs.get(k,color_residual)
    for k in ['linestyleresidual','line_style_residual','linestyleResidual']: linestyle_residual=kwargs.get(k,linestyle_residual)
    for k in ['linewidthresidual','linewidthResidual','line_width_residual']: linewidth_residual=kwargs.get(k,linewidth_residual)
    for k in ['alpharesidual','alphaResidual']: alpha_residual=kwargs.get(k,alpha_residual)

    for k in ['bands']: band=kwargs.get(k,band)
    if len(band) == 2 and isinstance(band[0],list) == False: band = [band]
    for k in ['bandcolors','bandcolor','band_colors']: band_color=kwargs.get(k,band_color)
    for k in ['bandalphas','band_alphas','bandalpha']: band_alpha=kwargs.get(k,band_alpha)
    for k in ['band_labels','bandlabel','bandlabels']: band_label=kwargs.get(k,band_label)
    for k in ['band_label_positions','bandlabelposition','bandlabelpositions']: band_label_position=kwargs.get(k,band_label_position)
    for k in ['bandwidth','bandwidths','band_widths']: band_width=kwargs.get(k,band_width)
    for par in [band_color,band_alpha,band_label,band_label_position,band_width]:
        if not isinstance(par,list): par = [par]*len(band)
        if len(par) < len(band): par.extend([par[-1] for x in range(len(band)-len(par))])

    for k in ['legends','label','labels']: legend=kwargs.get(k,legend)
    if not isinstance(legend,list): legend = [legend]
    for k in ['legendfontscale','legendFontscale']: legend_fontscale=kwargs.get(k,legend_fontscale)
    legend_fontscale=legend_fontscale*fontscale
    for k in ['legendLocation','legendlocation','labelLocation','labellocation','label_location']: legend_location=kwargs.get(k,legend_location)

    for k in ['xrange','x_range','wave_range','wrange','wrng']: xrng=kwargs.get(k,xrng)
    if not isinstance(xrng,list): xrng = [xrng]
    for k in ['yrange','y_range','flux_range','frange','frng']: yrng=kwargs.get(k,yrng)
    if not isinstance(yrng,list): yrng = [yrng]

    for k in ['multilayout','multiLayout','multi_layout']: layout=kwargs.get(k,layout)
    for k in ['file','filename']: output=kwargs.get(k,output)
    if not isinstance(output,str): output=''
    filetype = '.pdf'
    if output!='': filetype=output.split('.')[-1]

    if comparison != None and isinstance(comparison,splat.Spectrum) == False and isinstance(comparison,list) == False: 
        print('plotSpectrum() Warning: comparison spectrum should be a splat Spectrum object, you passed {}'.format(comparison))
        comparison = None

# some plotting constants
    xlabel_default = 'Wavelength'
    ylabel_deafult = 'Flux'

# telluric bands in micron
    telluric_bands = [[1.1,1.2]*u.micron,[1.3,1.5]*u.micron,[1.75,2.0]*u.micron]

# assign features by group
    if not isinstance(features,list): features = [features]
    if ldwarf==True or mdwarf==True: features.extend(['k','na','feh','tio','co','h2o','h2'])
    if tdwarf==True: features.extend(['k','ch4','h2o','h2'])
    if young==True: features.extend(['vo'])
    if binary==True: features.extend(['sb'])

# clean repeats in features while maintaining order - set does not do this
    if len(features)>0:
        fea = []
        for i in features:
            if i not in fea: fea.append(i)
        features = fea


# if a list is passed, use this list
    splist = copy.deepcopy(inp)
    if isinstance(splist,list) == False: splist = [splist]
    
# set up for multiplot
    if len(splist) == 1: multiplot = False
    
# array of lists => force multiplot
    elif len(splist) > 1 and isinstance(splist[0],list) == True: multiplot = True
    else: pass

# reformat array of spectra of multiplot is used (i.e., user forgot to set)
    if multiplot == True and isinstance(splist[0],splat.Spectrum):
        splist = [[s] for s in splist]

    elif multiplot == False and isinstance(splist[0],splat.Spectrum):
        splist = [splist]
        
# flatten array if multiplot is not set
    elif multiplot == False and isinstance(splist[0],list) and len(splist) > 1:
        splist = [[item for sublist in splist for item in sublist]]       # flatten
    else: pass

# total number of spectra - use to assign default legends
    allsps = [item for sublist in splist for item in sublist]   # Total number of spectra
    if len(legend) == 0: legend=[sp.name for sp in allsps]
    if len(legend) < len(allsps):
        legend.extend([allsps[i].name for i in range(len(legend),len(allsps)-len(legend))])
    

# now run a loop through the input subarrays
    plt.close('all')

# set up here for multiple file output
    nplot = 1
    if multipage == True or multiplot == True:
        nplot = layout[0]*layout[1]
        numpages = int(len(splist) / nplot) + 1
        if (len(splist) % nplot == 0):
               numpages -= 1
        fig = []
        
    if multipage == True and filetype.lower() == 'pdf':
        pdf_pages = PdfPages(output)
        
    if multipage == False:
        if len(splist) > 1:
            filebase = output.replace('.{}'.format(filetype),'')
            files = [filebase+'{}.'.format(i+1)+filetype for i in numpy.arange(len(splist))]
        else:
            files = [output]

    pg_n = 0        # page counter
    plt_n = 0       # plot per page counter
    lg_n = 0        # legend per plot counter

    for plts,sp in enumerate(splist):
# set specific plot parameters
        if not isinstance(sp[0],splat.Spectrum):
            raise ValueError('\nInput to plotSpectrum has wrong format:\n\n{}\n\n'.format(sp[0]))

# set up plotting defaults for the list of spectra - REPLACE THIS
        if not isinstance(zeropoint,list): zeropoint = [zeropoint]*len(sp)
        if len(zeropoint) < len(sp): zeropoint.extend([zeropoint[-1] for x in range(len(sp)-len(zeropoint))])
        if not isinstance(color,list): color = [color]*len(sp)
        if len(color) < len(sp): color.extend([color[-1] for x in range(len(sp)-len(color))])
        if not isinstance(linestyle,list): linestyle = [linestyle]*len(sp)
        if len(linestyle) < len(sp): linestyle.extend([linestyle[-1] for x in range(len(sp)-len(linestyle))])
        if not isinstance(linewidth,list): linewidth = [linewidth]*len(sp)
        if len(linewidth) < len(sp): linewidth.extend([linewidth[-1] for x in range(len(sp)-len(linewidth))])
        if not isinstance(alpha,list): alpha = [alpha]*len(sp)
        if len(alpha) < len(sp): alpha.extend([alpha[-1] for x in range(len(sp)-len(alpha))])
        if not isinstance(color_noise,list): color_noise = [color_noise]*len(sp)
        if len(color_noise) < len(sp): color_noise.extend([color_noise[-1] for x in range(len(sp)-len(color_noise))])
        if not isinstance(linestyle_noise,list): linestyle_noise = [linestyle_noise]*len(sp)
        if len(linestyle_noise) < len(sp): linestyle_noise.extend([linestyle_noise[-1] for x in range(len(sp)-len(linestyle_noise))])
        if not isinstance(linewidth_noise,list): linewidth_noise = [linewidth_noise]*len(sp)
        if len(linewidth_noise) < len(sp): linewidth_noise.extend([linewidth_noise[-1] for x in range(len(sp)-len(linewidth_noise))])
        if not isinstance(alpha_noise,list): alpha_noise = [alpha_noise]*len(sp)
        if len(alpha_noise) < len(sp): alpha_noise.extend([alpha_noise[-1] for x in range(len(sp)-len(color_noise))])
        if not isinstance(color_comparison,list): color_comparison = [color_comparison]*len(sp)
        if len(color_comparison) < len(sp): color_comparison.extend([color_comparison[-1] for x in range(len(sp)-len(color_comparison))])
        if not isinstance(linestyle_comparison,list): linestyle_comparison = [linestyle_comparison]*len(sp)
        if len(linestyle_comparison) < len(sp): linestyle_comparison.extend([linestyle_comparison[-1] for x in range(len(sp)-len(linestyle_comparison))])
        if not isinstance(linewidth_comparison,list): linewidth_comparison = [linewidth_comparison]*len(sp)
        if len(linewidth_comparison) < len(sp): linewidth_comparison.extend([linewidth_comparison[-1] for x in range(len(sp)-len(linewidth_comparison))])
        if not isinstance(alpha_comparison,list): alpha_comparison = [alpha_comparison]*len(sp)
        if len(alpha_comparison) < len(sp): alpha_comparison.extend([alpha_comparison[-1] for x in range(len(sp)-len(alpha_comparison))])

# settings that work if the spectrum was read in as legitmate Spectrum object
        try:
            xlabel = kwargs.get('xlabel','{} ({})'.format(sp[0].wave_label,sp[0].wave.unit))
            ylabel = kwargs.get('ylabel','{} ({})'.format(sp[0].flux_label,sp[0].flux.unit))
        except:
            xlabel = kwargs.get('xlabel',xlabel_default)
            ylabel = kwargs.get('ylabel',ylabel_default)
# initial plot range
        bound = [numpy.nanmin(sp[0].wave.value),numpy.nanmax(sp[0].wave.value)]
        ymax = [numpy.nanquantile(s.flux.value,0.98) for s in sp]
        bound.extend(numpy.array([-0.02,1.3])*numpy.nanmax(ymax)+\
            numpy.array([numpy.nanmin(zeropoint),numpy.nanmax(zeropoint)+stack*(len(sp)-1)]))

# set colormap if provided
        if colormap != None:
            values = numpy.arange(len(sp))
            color_map = plt.get_cmap(colormap)
            norm  = colmap.Normalize(vmin=0, vmax=1.0*values[-1])
            scalarMap = cm.ScalarMappable(norm=norm, cmap=color_map)
            for i in range(len(sp)): color[i] = scalarMap.to_rgba(values[i])

# GENERATE PLOTS
        if multiplot == True or multipage == True:
            plt_n = plts % nplot
            if (plt_n == 0):
                fig.append(plt.figure())
                pg_n += 1
            ax = fig[pg_n-1].add_subplot(layout[0], layout[1], plt_n+1)
            
# plotting a single plot with all spectra
        else:
            plt.close('all')
            plt_n = 0
            fig = []
            if len(figsize)>0: fig.append(plt.figure(figsize=figsize))
            else: fig.append(plt.figure())
            ax = fig[0].add_subplot(111)
        
        for ii, a in enumerate(sp):
# zeropoint and stack
            flx = [i+zeropoint[ii] for i in a.flux.value]
            if stack > 0: flx = [f + (len(sp)-ii-1)*stack for f in flx]
            ax.plot(a.wave.value,flx,color=color[ii],linestyle=linestyle[ii], lw=linewidth[ii], alpha=alpha[ii], zorder = 10, label = legend[lg_n])  

# add comparison
            if comparison != None:
# zeropoint and stack
                cflx = [i+zeropoint[ii] for i in comparison.flux.value]
                if stack > 0: cflx = [f + (len(sp)-ii-1)*stack for f in cflx]
                ax.plot(comparison.wave.value,cflx,color=color_comparison[ii],linestyle=linestyle_comparison[ii], lw=linewidth_comparison[ii], alpha=alpha_comparison[ii], zorder = 10)
    
# add residual
            if residual == True and len(sp) == 2:
                # Save flux values from first spectrum
                if ii == 0:
                    flx0 = [f - (len(sp)-ii-1)*stack for f in flx]
                    
                # Subtract fluxes and plot
                elif ii == 1:
                    res = [flx0[f_n] - f for f_n, f in enumerate(flx)]
                    ax.plot(a.wave.value, res, alpha = alpha_residual[ii], color = color_residual[ii], linsetyle=linestyle_residual[ii], lw=linewidth_residual[ii])
                    
                    # Fix bound[2] if residual goes below 0
                    if numpy.nanmin(res) < bound[2]:
                        b0 = numpy.argmin(a.wave.value[a.wave.value > bound[0]])
                        b1 = numpy.argmax(a.wave.value[a.wave.value < bound[1]])
                        bound[2] = numpy.nanmin(res[b0:b1])

# noise
            if show_noise == True:
                ns = [i+zeropoint[ii] for i in a.noise.value]
                ax.plot(a.wave.value,ns,color=color_noise[ii],linestyle=linestyle_noise[ii],alpha=alpha_noise[ii], lw=linewidth_noise[ii], zorder = 10)

# zeropoint
            if show_zero == True:
                ze = numpy.ones(len(a.flux))*zeropoint[ii]
                ax.plot(a.wave.value,ze,color=color[ii],linestyle=linestyle_zero,alpha=alpha_zero,lw=linewidth_zero, zorder = 10)

# save maximum flux among all spectra for plotting
# THIS IS VERY SLOW AND IT WOULD BE BETTER TO FIND AN ALTERNATE APPROACH
            if len(features)>0:
                f = interp1d(a.wave,flx,bounds_error=False,fill_value=0.)
                if ii == 0: 
                    wvmax = numpy.linspace(bound[0],bound[1],nsamples)
                    flxmax = numpy.array(f(wvmax))
                else: flxmax = numpy.maximum(flxmax,numpy.array(f(wvmax)))

# legend counter
            lg_n = lg_n + 1 # Increment legend


# label features
# THIS NEEDS TO BE FIXED WITH GRETEL'S STUFF
        if len(features) > 0:
            yoff = 0.02*(bound[3]-bound[2]) # label offset
            fontsize = int((10-numpy.nanmin([(layout[0]*layout[1]-1),6]))*fontscale)
            for ftr in features:
                ftr = ftr.lower()
                if ftr in FEATURE_LABELS:
                    ftrc = checkDict(ftr,FEATURE_LABELS)
                    if ftrc != False:
                        for ii,waveRng in enumerate(FEATURE_LABELS[ftrc]['wavelengths']):
                            wRng = waveRng.to(sp[0].wave.unit).value
# features must be contained in plot range (may change this)
                            if numpy.nanmin(wRng) > bound[0] and numpy.nanmax(wRng) < bound[1]:
                                wfeature = numpy.where(numpy.logical_and(wvmax >= numpy.nanmin(wRng),wvmax <= numpy.nanmax(wRng)))
                                if len(wvmax[wfeature]) == 0: wfeature = numpy.argmax(numpy.absolute(wvmax-numpy.nanmedian(wRng)))
                                y = numpy.nanmax(flxmax[wfeature])+yoff
                                flxmax[wfeature] = flxmax[wfeature]+3.*yoff

                                if FEATURE_LABELS[ftrc]['type'] == 'band':
                                    ax.plot(wRng,[y+yoff]*2,color='k',linestyle='-')
                                    ax.plot([wRng[0]]*2,[y,y+yoff],color='k',linestyle='-')
                                    ax.text(numpy.mean(wRng),y+1.5*yoff,FEATURE_LABELS[ftrc]['label'],horizontalalignment='center',fontsize=fontsize)
                                else:
                                    for w in wRng: ax.plot([w]*2,[y,y+yoff],color='k',linestyle='-')
                                    ax.text(numpy.mean(wRng),y+1.5*yoff,FEATURE_LABELS[ftrc]['label'],horizontalalignment='center',fontsize=fontsize)
            bound[3] = numpy.nanmax([numpy.nanmax(flxmax)+2.*yoff,bound[3]])

# add grid
        if grid == True: ax.grid()            

# axis labels 
        fontsize = (numpy.round(numpy.max([13./((layout[0]*layout[1])**0.33),5]))) * fontscale
        legend_fontsize = (13-numpy.min([(layout[0]*layout[1]-1),8])) * legend_fontscale
        ax.set_xlabel(xlabel, fontsize = fontsize)
        ax.set_ylabel(ylabel, fontsize = fontsize)
        ax.tick_params(axis='x', labelsize=fontsize)
        ax.tick_params(axis='y', labelsize=fontsize)

# add title
        if title!='': ax.set_title(title)

# log scale?
        if kwargs.get('xlog',False): ax.set_xscale('log',nonposx='clip')
        if kwargs.get('ylog',False): ax.set_yscale('log',nonposy='clip')

# place legend
        if len(legend) > 0:
            if legend_location == 'outside':
                box = ax.get_position()
                ax.set_position([box.x0, box.y0 + box.height * 0.15, box.width * 0.7, box.height * 0.7])
                ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), prop={'size':legend_fontsize})
            else:
                ax.legend(loc=legend_location, prop={'size':legend_fontsize})
                bound[3] = bound[3]+0.1*(bound[3]-bound[2])     # extend axis for in-plot legends

# overplot telluric absorption
        if telluric == True:
            yoff = 0.02*(bound[3]-bound[2]) # label offset
            for waveRng in telluric_bands:
                wR = waveRng.to(sp[0].wave.unit).value
                rect = patches.Rectangle((wR[0],bound[2]),wR[1]-wR[0],bound[3]-bound[2],facecolor=color_telluric,alpha=alpha_telluric,color=color_telluric)
                ax.add_patch(rect)
                ax.text(numpy.mean(wR),bound[2]+3*yoff,r'$\oplus$',horizontalalignment='center',fontsize=fontsize)

# overplot color swaths for pre-specified bands
        if len(band) > 0:
            for i,b in enumerate(band):
                if not isinstance(b,list): 
                    try: b = [float(b)-0.5*band_width,float(b)+0.5*band_width]
                    except:
                        print('\nWarning: plotSpectrum bands variables should be array of 2-element arrays; you passed {}'.format(band))
                        b = [0.,0.]
                rect = patches.Rectangle((b[0],bound[2]),b[1]-b[0],bound[3]-bound[2],facecolor=band_color[i],color=band_color[i],alpha=band_alpha[i])
                ax.add_patch(rect)
                if band_label_position[i].lower() == 'top':
                    ax.text(numpy.mean(b),bound[3]-3*yoff,band_label[i],horizontalalignment='center',fontsize=fontsize)
                elif band_label_position[i].lower() == 'middle':
                    ax.text(numpy.mean(b),0.5*(bound[2]+bound[3]),band_label[i],horizontalalignment='center',fontsize=fontsize)
                else:
                    ax.text(numpy.mean(b),bound[2]+3*yoff,band_label[i],horizontalalignment='center',fontsize=fontsize)

# place inset - RIGHT NOW ONLY SETTING LIMITS WITH FIRST SPECTRUM IN LIST
        if inset == True and len(inset_xrange) == 2:
            ax_inset = fig[pg_n-1].add_axes(inset_position) #, axisbg='white')
            bound2 = inset_xrange
            if len(inset_yrange) == 0:
                b0 = numpy.argmax(sp[0].wave.value > bound2[0])
                b1 = numpy.argmin(sp[0].wave.value < bound2[1])
                inset_yrange = [numpy.nanmin(sp[0].flux.value[b0:b1]),numpy.nanmax(sp[0].flux.value[b0:b1])]
            bound2.extend(inset_yrange)
            db = (bound2[3]-bound2[2])
            bound2[2] = bound2[2]-0.05*db
            bound2[3] = bound2[3]+0.05*db
            ax_inset.axis(bound2)
            inset_fontsize = fontsize*0.7

            for ii,a in enumerate(sp):
                flx = [i+zeropoint[ii] for i in a.flux.value]
                ax_inset.plot(a.wave.value,flx,color=colors[ii],linestyle=linestyle[ii],linewidth=linewidth[ii],alpha=alpha[ii])            
                ax_inset.set_xlabel('')
                ax_inset.set_ylabel('')
                ax_inset.tick_params(axis='x', labelsize=inset_fontsize)
                ax_inset.tick_params(axis='y', labelsize=inset_fontsize)
#                ax_inset.legend()

# inset feature labels
            if len(inset_features) > 0:
                yoff = 0.05*(bound2[3]-bound2[2])
                for ftr in inset_features:
                    ftrc = checkDict(ftr,FEATURE_LABELS)
                    if ftrc != False:
                        for ii,waveRng in enumerate(FEATURE_LABELS[ftrc]['wavelengths']):
                            wRng = waveRng.to(sp[0].wave.unit).value
                            if (numpy.min(wRng) > bound2[0] and numpy.max(wRng) < bound2[1]):
                                wfeature = numpy.where(numpy.logical_and(wvmax >= numpy.nanmin(wRng),wvmax <= numpy.nanmax(wRng)))
                                if len(wvmax[wfeature]) == 0: wfeature = numpy.argmax(numpy.absolute(wvmax-numpy.nanmedian(wRng)))
                                y = numpy.nanmax(flxmax[wfeature])+yoff
                                flxmax[wfeature] = flxmax[wfeature]+3.*yoff
        
                                if FEATURE_LABELS[ftrc]['type'] == 'band':
                                    ax_inset.plot(wR,[y+yoff]*2,color='k',linestyle='-')
                                    ax_inset.plot([wR[0]]*2,[y,y+yoff],color='k',linestyle='-')
                                    ax_inset.text(numpy.mean(wR),y+2*yoff,FEATURE_LABELS[ftrc]['label'],horizontalalignment='center',fontsize=inset_fontsize)
                                else:
                                    for w in waveRng:
                                        ax_inset.plot([w]*2,[y,y+yoff],color='k',linestyle='-')
                                    ax_inset.text(numpy.mean(wR),y+2*yoff,FEATURE_LABELS[ftrc]['label'],horizontalalignment='center',fontsize=inset_fontsize)
                                    waveRng = [wR[0]-0.02,wR[1]+0.02]   # for overlap
        
# update offset
                if len(inset_features) > 0: bound2[3] = numpy.nanmax([bound2[3],numpy.nanmax(flxmax)+5.*yoff])
                ax_inset.axis(bound2)

# finalize bounding
        if len(xrng) > 0: bound[0:2] = xrng
        if len(yrng) > 0: bound[2:4] = yrng
        if isUnit(bound[0]): bound = [x.value for x in bound]
        ax.axis(bound)
    
# save to file or display
# ERROR HERE - CHECK WHAT FILES
        if multipage == False:
            if files[plts] != '' and (plts % nplot == 3 or plts == len(splist)-1):
                if kwargs.get('tight',True) == True: 
                    plt.savefig(files[plts], bbox_inches='tight')
                else:
                    plt.savefig(files[plts])
    if output == '' and not kwargs.get('web',False):
        plt.show()
        if (kwargs.get('interactive',False) != False): plt.ion()
        else: plt.ioff()


# save figures in multipage format and write off pdf file
    if multipage == True:    
        for pg_n in numpy.arange(numpages):
#            fig[pg_n].text(0.5, 0.04, xlabel, ha = 'center', va = 'center')
#            fig[pg_n].text(0.06, 0.5, ylabel, ha = 'center', va = 'center', rotation = 'vertical')
            fig[pg_n].tight_layout
            fig[pg_n].suptitle(title, fontsize = int(14*fontsize), fontweight = 'bold')
            pdf_pages.savefig(fig[pg_n])
        if filetype.lower() == 'pdf':
            pdf_pages.close()

    plt.clf()
    return fig


def plotBatch(inp,output='spectra_plot.pdf',comparisons=None,classify=False,normalize=False,normrange=[0.9,1.4],layout=[2,2],basecolors=['k','m'],legend=[],fontscale=0.7,classify_kwargs={},plot_kwargs={},**kwargs):
    """
    Plots a batch of spectra into a 2x2 set of PDF files, with options of overplotting 
    comparison spectra, including best-match spectral standards inferred from `splat.classifyByStandard()`_.
    This routine is essentially a wrapper for a specific use case of `splat.plot.plotSpectrum()`_.

    Parameters 
    ----------
    inp : Spectrum or string or array of these
        A single or list of Spectrum objects or filenames, or the glob search string 
        for a set of files (e.g., '/Data/myspectra/*.fits'). 

    output : string, default = 'spectra_plot.pdf'
        Filename for PDF file output; full path should be included if not saving to local directory

    comparisons : Spectrum or string or array of these, default = None
        list of Spectrum objects or filenames for comparison spectra. If comparisons list is shorter 
        than source list, then the last comparison source will be repeated. If the comparisons list is 
        longer, the list will be truncated. 
    
    classify : boolean, default = False
        Set to True to classify sources based on comparison to spectral standards using `splat.classifyByStandard()`_.
        Specific parameters for classification can be specifed with the `classify_kwargs` keyword

    normalize : boolean, default = False
        Set to True to normalize source and (if passed) comparison spectra.

    normrange : list, default = [0.9,1.4]
        Range in micron to compute the normalization

    layout : list, default = [2,2]
        Layout of plots on page (# per row x # per column)

    basecolors : list, default = ['k','m']
        Baseline color pair for plotted spectrum and its (optional) comparison template

    legend : list, default = []
        Set to list of legends for plots. The number of legends should equal the number of 
        sources and comparisons (if passed) in an alternating sequence. The default is to display
        the file name for each source and object name for comparison (if passed)

    fontscale : float, default = 0.7
        Scale factor to change font size; default value is optimal for 2x2 layout 

    classify_kwargs : dict, default = {}
        Dictionary of keywords for `splat.classifyByStandard()`_ 

    plot_kwargs : dict, default = {}
        Dictionary of additional keywords for `splat.plot.plotSpectrum()`_ 

    **kwargs : dict, optional
        Additional keyword parameters

    Returns 
    -------
    matplotlib figure object
        Figure object containing all of the generated plots

    Examples
    --------
    **Case 1:** These commands will grab a set of high S/N, optically-classified L5 dwarfs, and
    plot them into a single multiple-page PDF document in a 2x2 grid with the best fit 
    classification standard overplotted:

        >>> import glob, splat
        >>> import splat.plot as splot
        >>> splist = splat.getSpectrum(opt_spt='L5',snr=100)
        >>> splot.plotBatch(splist,classify=True,normalize=True,yrange=[0,1.2],output='comparison.pdf')
              
    **Case 2:** These commands will grab a list of spectral files, plot them into a single multiple-page PDF
    document in a 2x2 grid, with the best fit classification standard overplotted. The three 
    calls of `plotBatch()` produce the same outcome.

        >>> import glob, splat
        >>> import splat.plot as splot
        >>> files = glob.glob('/home/mydata/*.fits')
        >>> sp = splot.plotBatch(files,classify=True,output='comparison.pdf')
        >>> sp = splot.plotBatch('/home/mydata/*.fits',classify=True,output='comparison.pdf')
        >>> sp = splot.plotBatch([Spectrum(file=f) for f in files],classify=True,output='comparison.pdf')

    .. _`splat.plot.plotSpectrum()` : api.html#splat.plot.plotSpectrum
    .. _`splat.classifyByStandard()` : api.html#splat.classifyByStandard
    .. _`Spectrum()` : api.html#Spectrum



    See Also
    --------
    plotSpectrum : primary Spectrum plotting routine

    """

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
        raise ValueError('\nInput should be list of Spectrum objects or filenames')

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
    """
    :Purpose: ``Plot index-index plots.``

    indices should be a dictionary of {'index_name': {'ranges': [[],[]], 'value': #, 'unc': #}, ...}
    Not currently implemented
    """

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
    """
    :Purpose: ``Plot SED photometry with SpeX spectrum.``

    Not currently implemented
    """
    pass
    return



def plotSequence(spec,type_range=2, std_class='dwarf', spt='', output='', verbose=False, **kwargs):
    """
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

        A matplotlib figure object containing the plot of the spectrum compared to sequence of standards on screen or saved to file

    :Example:
       >>> import splat
       >>> import splat.plot as splot
       >>> sp = splat.getSpectrum(lucky=True)[0]
       >>> fig = splat.plotSequence(sp,output='classify.pdf')

.. _`plotSpectrum()` : api.html#splat.plot.plotSpectrum
.. _`classifyByStandard()` : api.html#splat.classifyByStandard

    """

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


    
