"""
.. note::
         These are the plotting functinos for the SPLAT code 
"""


# Related third party imports.
import matplotlib.pyplot as plt
import matplotlib.colors as colmap
import matplotlib.cm as cm
from matplotlib.backends.backend_pdf import PdfPages
import numpy
from scipy.interpolate import interp1d 
import splat

# To do:
#     masking telluric regions
#    labeling features
#     add legend
def plotSpectrum(*args, **kwargs):

# keyword parameters
    nsamples = kwargs.get('nsamples',1000)
    multiplot = kwargs.get('multiplot',False)       # create multiple plots
    multipage = kwargs.get('multipage',False)       # create a multiple page sequence of plots
    multilayout = kwargs.get('multilayout',[1,1])       # layout of multiple plots, [# horizontal, # vertical]
    stack = kwargs.get('stack',0)                   # stack spectra on top of each other
    legend = kwargs.get('legend','')                # legend text
    grid = kwargs.get('grid',False)                 # plot internal grid lines
    filename = kwargs.get('filename','')            # output filename
    filename = kwargs.get('file',filename)
    title = kwargs.get('title','')
    filebase = filename.split('.')[0]               # filebase for multiple files
    filetype = kwargs.get('format',filename.split('.')[-1])
    filetype.lower()
    if filetype == '':
        filetype = 'eps'
    comparison = kwargs.get('comparison',False)
    if comparison.__class__.__name__ != 'Spectrum':
        comparison = False
        
#    mask = kwargs.get('mask',False)                # not yet implemented

# features to label on spectra
    feature_labels = { \
        'h2o': {'label': r'H$_2$O', 'type': 'band', 'wavelengths': [[0.92,0.95],[1.08,1.20],[1.325,1.550],[1.72,2.14]]}, \
        'ch4': {'label': r'CH$_4$', 'type': 'band', 'wavelengths': [[1.1,1.24],[1.28,1.44],[1.6,1.76],[2.05,2.35]]}, \
        'co': {'label': r'CO', 'type': 'band', 'wavelengths': [[2.28,2.42]]}, \
        'tio': {'label': r'TiO', 'type': 'band', 'wavelengths': [[0.76,0.80],[0.825,0.831]]}, \
        'vo': {'label': r'VO', 'type': 'band', 'wavelengths': [[1.04,1.08]]}, \
        'feh': {'label': r'FeH', 'type': 'band', 'wavelengths': [[0.86,0.90],[0.98,1.03],[1.19,1.31],[1.57,1.64]]}, \
        'h2': {'label': r'H$_2$', 'type': 'band', 'wavelengths': [[2.05,2.6]]}, \
        'oh': {'label': r'tell', 'type': 'band', 'wavelengths': [[1.37,1.45],[1.82,2.0]]}, \
        'sb': {'label': r'*', 'type': 'band', 'wavelengths': [[1.6,1.64]]}, \
        'h': {'label': r'H I', 'type': 'line', 'wavelengths': [[1.004,1.005],[1.093,1.094],[1.281,1.282],[1.944,1.945],[2.166,2.166]]},\
        'na': {'label': r'Na I', 'type': 'line', 'wavelengths': [[0.8186,0.8195],[1.136,1.137],[2.206,2.209]]}, \
        'k': {'label': r'K I', 'type': 'line', 'wavelengths': [[0.7699,0.7665],[1.169,1.177],[1.244,1.252]]}, \
        'telluric': {'label': r'O$_2\oplus$', 'type': 'band', 'wavelengths': [[1.1,1.2],[1.3,1.5],[1.75,2.0],[2.4,2.5]]}}

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

# PROCESS INPUTS
# FORMAT FOR PLOTTING WILL BE ARRAY OF SPECTRUM OBJECTS = ONE PLOT

# error check - make sure you're plotting something
    if (len(args) < 1):
        print 'plotSpectrum needs at least one Spectrum object to plot'
        return

    splist = []
    cnt=0
    for a in args:
        if a.__class__.__name__ == 'Spectrum':      # a spectrum object
            splist.append([a])
        elif isinstance(a,list):
            splist.append(a)
            cnt+=1
        else:
            print '\nplotSpectrum: Ignoring input object {} as it is neither a Spectrum object nor a list\n\n'.format(a)

# force multiplot if more than one list of objects passed (i.e., user forgot to set)
    if cnt > 1:
        multiplot = True
        
# rearrange according to setting of multiplot
    if multiplot == False:
        splist = [[item for sublist in splist for item in sublist]]       # flatten


    
# now run a loop through the input subarrays
    plt.close('all')

# set up here for multiple file output
    if multipage == True:
        nplot = multilayout[0]*multilayout[1]
        if filetype == 'pdf':
            pdf_pages = PdfPages(filename)
        numpages = int(len(splist) / nplot) + 1
        if (len(splist) % nplot == 0):
                numpages -= 1
        fig = range(numpages)
    else:
        if len(splist) > 1:
            files = [filebase+'{}.'.format(i+1)+filetype for i in range(len(splist))]

    pg_n = 0
    plt_n = 0
    plt.close('all')
    for plts,sp in enumerate(splist):

# set specific plot parameters
        zeropoint = kwargs.get('zeropoint',[0. for x in range(len(sp))])
        xlabel = kwargs.get('xlabel','{} ({})'.format(sp[0].wlabel,sp[0].wunit))
        ylabel = kwargs.get('ylabel','{} {} ({})'.format(sp[0].fscale,sp[0].flabel,sp[0].funit))
        xrange = kwargs.get('xrange',[0.85,2.4])
        bound = xrange
        ymax = [s.fluxMax().value for s in sp]
        yrng = kwargs.get('yrange',map(lambda x: x*(numpy.nanmax(ymax)+numpy.nanmax(zeropoint)),[-0.02,1.2]))
        bound.extend(yrng)
        linestyle = kwargs.get('linestyle',['steps' for x in range(len(sp))])
        if (len(linestyle) < len(sp)):
            linestyle.extend(['steps' for x in range(len(sp)-len(linestyle))])
# colors
        colors = kwargs.get('colors',['k' for x in range(len(sp))])
        if (len(colors) < len(args)):
            colors.extend(['k' for x in range(len(sp)-len(colors))])
        colorsUnc = kwargs.get('colorsUnc',['k' for x in range(len(sp))])
        if (len(colorsUnc) < len(sp)):
            colorsUnc.extend(['k' for x in range(len(sp)-len(colorsUnc))])
        if (kwargs.get('colorScheme') != None):
            colorScheme = kwargs.get('colorScheme')
            values = range(len(sp))
            color_map = plt.get_cmap(colorScheme)
            norm  = colmap.Normalize(vmin=0, vmax=values[-1])
            scalarMap = cm.ScalarMappable(norm=norm, cmap=color_map)
            for i in range(len(sp)):
                if colors[i] == '':
                    colors[i] = scalarMap.to_rgba(values[i])
        else:
            for i in range(len(sp)):
                if colors[i] == '':
                    colors[i] = 'k'
        labels = kwargs.get('labels',['' for x in range(len(sp))])
        if(len(labels) < len(args)):
            labels.extend(['' for x in range(len(sp)-len(labels))])
#        for x in range(len(sp)):
#            if labels[x] == '':
#                labels[x] = sp[x].shortname


        showNoise = kwargs.get('showNoise',[False for x in range(len(sp))])
        if not isinstance(showNoise, list):
            showNoise = [showNoise]
        if (len(showNoise) < len(sp)):
            showNoise.extend([True for x in range(len(sp)-len(showNoise))])
        showZero = kwargs.get('showZero',[True for x in numpy.arange(len(sp))])
        if not isinstance(showZero, list):
            showZero = [showZero]
        if (len(showZero) < len(sp)):
            showZero.extend([True for x in range(len(sp)-len(showZero))])

# GENERATE PLOTS
        if (multipage == True):
            plt_n = plts % nplot
            if (plt_n == 0):
#                ax = range(nplot)
#                t = tuple([tuple([i+b*multilayout[1] for i in range(multilayout[1])]) for b in range(multilayout[0])])
                fig[pg_n], ax = plt.subplots(multilayout[0], multilayout[1], sharex = True, sharey = True)
                ax = [item for sublist in ax for item in sublist]       # flatten
                print ax
                pg_n += 1
            
# plotting a single plot with all spectra
        else:
            plt.close('all')
            ax = range(1)
            plt_n = 0
            if (kwargs.get('figsize') != None):
                fig = plt.figure(figsize = kwargs.get('figsize'))
            else:
                fig = plt.figure()
            ax[0] = fig.add_subplot(111)
        
        for ii, a in enumerate(sp):
            flx = [i+zeropoint[ii] for i in a.flux.value]
#stack
            if stack > 0:
                flx = [f + ii*stack for f in flx]
                if kwargs.get('yrange') == None:
                    bound[3] += bound[3] * stack * ii

            ax[plt_n].plot(a.wave.value,flx,color=colors[ii],linestyle=linestyle[ii], zorder = 10, label = labels[ii])  

# noise
            if (showNoise[ii]):
                ns = [i+zeropoint[ii] for i in a.noise.value]
                ax[plt_n].plot(a.wave.value,ns,color=colorsUnc[ii],linestyle=linestyle[ii],alpha=0.3, zorder = 10)

# zeropoint
            if (showZero[ii]):
                ze = numpy.ones(len(a.flux))*zeropoint[ii]
                ax[plt_n].plot(a.wave.value,ze,color=colors[ii],linestyle=':',alpha=0.3, zorder = 10)

# determine maximum flux for all spectra
            f = interp1d(a.wave,flx,bounds_error=False,fill_value=0.)
            if (ii == 0):
                flxmax = flx
                wvmax = a.wave
            else:
                flxmax = numpy.maximum(flxmax,f(wvmax))

# label features
# THIS NEEDS TO BE FIXED WITH GRETEL'S STUFF
        f = interp1d(wvmax,flxmax,bounds_error=False,fill_value=0.)
        yoff = 0.02*numpy.nanmax(yrng)
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
#                           flg = [1 if ((numpy.nanmin(waveRng) < numpy.nanmin(w) and numpy.max(waveRng) > numpy.nanmax(w)) or (numpy.nanmin(waveRng) > numpy.nanmin(w) and numpy.nanmax(waveRng) < numpy.nanmax(w))) else 0 for w in xset]
#                           flg = [1 if ((waveRng[0] < w[0] and waveRng[1] > w[1]) or (waveRng[0] > w[0] and waveRng[1] < w[1])) else 0 for w in xset]
                            flg1 = [1 if (waveRng[0] < w[0] and waveRng[1] > w[1]) else 0 for w in xset]
                            flg2 = [1 if (min(waveRng) > min(w) and max(waveRng) < max(w)) else 0 for w in xset]
                            flg = flg1+flg2
#                           flg = [1 if (abs(numpy.mean(waveRng)-numpy.mean(w))<0.1) else 0 for w in xset]
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
                            ax[plt_n].plot(waveRng,[y+yoff]*2,color='k',linestyle='-')
                            ax[plt_n].plot([waveRng[0]]*2,[y,y+yoff],color='k',linestyle='-')
                            ax[plt_n].text(numpy.mean(waveRng),y+1.5*yoff,feature_labels[ftr]['label'],horizontalalignment='center')
                        else:
                            for w in waveRng:
                                ax[plt_n].plot([w]*2,[y,y+yoff],color='k',linestyle='-')
                            ax[plt_n].text(numpy.mean(waveRng),y+1.5*yoff,feature_labels[ftr]['label'],horizontalalignment='center')


# grid
        if (grid):
            ax[plt_n].grid()
            
        ax[plt_n].axis(bound)
#        ax[plt_n].set_title(a.shortname, fontsize = 12)
        print bound

# axis labels - this may not work
        ax[plt_n].set_xlabel(xlabel)
        ax[plt_n].set_ylabel(ylabel)
#        plt.title(title)
    
# save to file or display
        if multipage == False:
            if filebase != '':
                plt.savefig(files[plts], format=filetype)
            else:
                plt.show()
                if (kwargs.get('interactive',False) != False):
                    plt.ion()        # make window interactive by default
                else:
                    plt.ioff()



# save figures in multipage format and write off pdf file
    if (multipage == True):    
        for pg_n in range(numpages):
#            fig[pg_n].text(0.5, 0.04, xlabel, ha = 'center', va = 'center')
#            fig[pg_n].text(0.06, 0.5, ylabel, ha = 'center', va = 'center', rotation = 'vertical')
            fig[pg_n].suptitle(title, fontsize = 14, fontweight = 'bold')
            pdf_pages.savefig(fig[pg_n])
        if filetype == 'pdf':
            pdf_pages.close()



    return

