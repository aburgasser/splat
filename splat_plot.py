"""
.. note::
         These are the plotting functinos for the SPLAT code 
"""


# Related third party imports.
from scipy.interpolate import interp1d 
import numpy
import matplotlib.pyplot as plt


# To do:
#     masking telluric regions
#    labeling features
#     add legend
def plotSpectrum(*args, **kwargs):

# error check - make sure you're plotting something
    if (len(args) < 1):
        print 'plotSpectrum needs at least on Spectrum object to plot'
        return
#    if isinstance(args[0],Spectrum):
#        print 'plotSpectrum needs at least one Spectrum object to plot'
#        return

# this loop solves issues if a list of spectrum objects is provided accidentally
    sp = []
    for i,a in enumerate(args):
        if (type(a) == list):
            sp.append(a[0])
        else:
            sp.append(a)

# absorption features
    feature_labels = { \
        'h2o': {'label': r'H$_2$O', 'type': 'band', 'wavelengths': [[0.92,0.95],[1.08,1.20],[1.325,1.450],[1.72,2.14]]}, \
        'ch4': {'label': r'CH$_4$', 'type': 'band', 'wavelengths': [[1.1,1.24],[1.28,1.44],[1.6,1.76],[2.05,2.35]]}, \
        'co': {'label': r'CO$$', 'type': 'band', 'wavelengths': [[2.28,2.42]]}, \
        'tio': {'label': r'TiO$$', 'type': 'band', 'wavelengths': [[0.76,0.80],[0.825,0.831]]}, \
        'vo': {'label': r'VO$$', 'type': 'band', 'wavelengths': [[1.04,1.08]]}, \
        'feh': {'label': r'FeH$$', 'type': 'band', 'wavelengths': [[0.86,0.90],[0.98,1.03],[1.19,1.31],[1.57,1.64]]}, \
        'h2': {'label': r'H$_2$', 'type': 'band', 'wavelengths': [[2.05,2.6]]}, \
        'oh': {'label': r'tell$$', 'type': 'band', 'wavelengths': [[1.37,1.45],[1.82,2.0]]}, \
        'sb': {'label': r'*$$', 'type': 'band', 'wavelengths': [[1.6,1.64]]}, \
        'h': {'label': r'H I$$', 'type': 'line', 'wavelengths': [[1.004,1.005],[1.093,1.094],[1.281,1.282],[1.944,1.945],[2.166,2.166]]},\
        'na': {'label': r'Na I$$', 'type': 'line', 'wavelengths': [[0.8186,0.8195],[1.136,1.137],[2.206,2.209]]}, \
        'k': {'label': r'K I$$', 'type': 'line', 'wavelengths': [[0.7699,0.7665],[1.169,1.177],[1.244,1.252]]}}


# keyword parameters
    nsamples = kwargs.get('nsamples',1000)
    title = kwargs.get('title','')
    zeropoint = kwargs.get('zeropoint',[0. for x in range(len(sp))])
    xlabel = kwargs.get('xlabel','{} ({})'.format(sp[0].wlabel,sp[0].wunit))
    ylabel = kwargs.get('ylabel','{} {} ({})'.format(sp[0].fscale,sp[0].flabel,sp[0].funit))
#    xrange = kwargs.get('xrange',[x.value for x in sp[0].waveRange()])
    xrange = kwargs.get('xrange',[0.85,2.4])
    bound = xrange
    ymax = [s.fluxMax().value for s in sp]
    yrange = kwargs.get('yrange',[0,1.2*(numpy.nanmax(ymax)+numpy.nanmax(zeropoint))])
    bound.extend(yrange)
    grid = kwargs.get('grid',False)
    colors = kwargs.get('colors',['k' for x in range(len(sp))])
    if (len(colors) < len(args)):
        colors.extend(['k' for x in range(len(sp)-len(colors))])
    colorsUnc = kwargs.get('colors',['k' for x in range(len(sp))])
    if (len(colorsUnc) < len(sp)):
        colorsUnc.extend(['k' for x in range(len(sp)-len(colorsUnc))])
    linestyle = kwargs.get('linestyle',['steps' for x in range(len(sp))])
    if (len(linestyle) < len(sp)):
        linestyle.extend(['steps' for x in range(len(sp)-len(linestyle))])
    filename = kwargs.get('filename','')
    format = kwargs.get('format',filename.split('.')[-1])
    showNoise = kwargs.get('showNoise',[False for x in range(len(sp))])
    if not isinstance(showNoise, list):
        showNoise = [showNoise]
    if (len(showNoise) < len(sp)):
        showNoise.extend([True for x in range(len(sp)-len(showNoise))])
    showZero = kwargs.get('showZero',[False for x in numpy.arange(len(sp))])
    if not isinstance(showZero, list):
        showZero = [showZero]
    if (len(showZero) < len(sp)):
        showZero.extend([True for x in range(len(sp)-len(showZero))])
    features = kwargs.get('features',[])
    if (kwargs.get('ldwarf',False) or kwargs.get('mdwarf',False)):
        features = list(set(features)|set(['H2O','TiO','CO','K','NA','FEH','H2']))
    if (kwargs.get('tdwarf',False)):
        features = list(set(features)|set(['H2O','CH4','K','H2']))
    if (kwargs.get('young',False)):
        features = list(set(features)|set(['VO']))
    if (kwargs.get('telluric',False)):
        features = list(set(features)|set(['OH']))
    if (kwargs.get('binary',False)):
        features = list(set(features)|set(['SB']))
#    mask = kwargs.get('mask',False)                # not yet implemented
#    labels = kwargs.get('labels','')            # not yet implemented


#    plt.clf()
# loop through sources
    plt.subplots(1)
    for ii,a in enumerate(sp):
        flx = [i+zeropoint[ii] for i in a.flux.value]
        plt.plot(a.wave.value,flx,color=colors[ii],linestyle=linestyle[ii])
# show noise
        if (showNoise[ii]):
            ns = [i+zeropoint[ii] for i in a.noise.value]
            plt.plot(a.wave.value,ns,color=colorsUnc[ii],linestyle=linestyle[ii],alpha=0.3)
# zeropoint
        if (showZero[ii]):
            ze = numpy.ones(len(a.flux))*zeropoint[ii]
            plt.plot(a.wave.value,ze,color=colors[ii],linestyle=':',alpha=0.3)
# determine maximum flux for all spectra
        f = interp1d(a.wave,flx,bounds_error=False,fill_value=0.)
        if (ii == 0):
            flxmax = flx
            wvmax = a.wave
        else:
            flxmax = numpy.maximum(flxmax,f(wvmax))
	         
# grid
    if (grid):
        plt.grid()

# label features
    f = interp1d(wvmax,flxmax,bounds_error=False,fill_value=0.)
    yoff = 0.02*numpy.nanmax(yrange)
    xset = [99.,99.]
    oflg = False
    for ftr in features:
        ftr = ftr.lower()
        if ftr in feature_labels:
            for ii,waveRng in enumerate(feature_labels[ftr]['wavelengths']):
                if (max(waveRng) > min(xrange)):
                    x = (numpy.arange(0,nsamples+1.0)/nsamples)* \
                        (numpy.nanmax(waveRng)-numpy.nanmin(waveRng)+0.1)+numpy.nanmin(waveRng)-0.05
                    y = numpy.nanmax(f(x))+0.5*yoff
#  THIS PART IS SUPPOSED TO SOLVE OVERLAP PROBLEM BUT NOT YET WORKING
                    if oflg:
#                        flg = [1 if ((numpy.nanmin(waveRng) < numpy.nanmin(w) and numpy.max(waveRng) > numpy.nanmax(w)) or (numpy.nanmin(waveRng) > numpy.nanmin(w) and numpy.nanmax(waveRng) < numpy.nanmax(w))) else 0 for w in xset]
#                        flg = [1 if ((waveRng[0] < w[0] and waveRng[1] > w[1]) or (waveRng[0] > w[0] and waveRng[1] < w[1])) else 0 for w in xset]
                        flg1 = [1 if (waveRng[0] < w[0] and waveRng[1] > w[1]) else 0 for w in xset]
                        flg2 = [1 if (min(waveRng) > min(w) and max(waveRng) < max(w)) else 0 for w in xset]
                        flg = flg1+flg2
#                        flg = [1 if (abs(numpy.mean(waveRng)-numpy.mean(w))<0.1) else 0 for w in xset]
                        y = y+4*yoff*numpy.sum(flg)
                        print ftr, numpy.sum(flg1), numpy.sum(flg2) #, [numpy.mean(a)-numpy.mean(waveRng) for a in xset]
                    else:
                        oflg = True
                    xset = [xset,waveRng]
#                w = numpy.where(numpy.nanmin(numpy.abs(numpy.subtract(xset,[numpy.mean(waveRng)]*len(xset)))) < 0.1)
#                if (len(w[0]) > 0):
#                    print ftr, waveRng, w, len(w[0]), xset[w[0]], numpy.abs(numpy.subtract(xset,[numpy.mean(waveRng)]*len(xset)))
#                    y = y+4*yoff*len(w)
#                if numpy.nanmin(numpy.abs(numpy.subtract(xset,[numpy.mean(waveRng)]*len(xset)))) < 0.2:
#                    y = y+4*yoff
#                xset.append(numpy.mean(waveRng))
                    if feature_labels[ftr]['type'] == 'band':
                        plt.plot(waveRng,[y+yoff]*2,color='k',linestyle='-')
                        plt.plot([waveRng[0]]*2,[y,y+yoff],color='k',linestyle='-')
                        plt.text(numpy.mean(waveRng),y+1.5*yoff,feature_labels[ftr]['label'],horizontalalignment='center')
                    else:
                        for w in waveRng:
                            plt.plot([w]*2,[y,y+yoff],color='k',linestyle='-')
                        plt.text(numpy.mean(waveRng),y+1.5*yoff,feature_labels[ftr]['label'],horizontalalignment='center')


# axis labels
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.axis(bound)
    plt.title(title)
    
# save to filen or display
    if (len(filename) > 0): 
        plt.savefig(filename, format=format)
    
    else:
        plt.show()
        if (kwargs.get('interactive',False) != False):
            plt.ion()        # make window interactive by default
        else:
            plt.ioff()
    return

