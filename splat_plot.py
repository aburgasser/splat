"""
.. note::
         These are the plotting functions for the SPLAT code 
"""


# Related third party imports.
import matplotlib.cm as cm
import matplotlib.colors as colmap
import matplotlib.patches as patches
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy
from scipy.interpolate import interp1d 
import splat

# To do:
#     add legend
#   
def plotSpectrum(*args, **kwargs):

# keyword parameters
    nsamples = kwargs.get('nsamples',1000)
    multiplot = kwargs.get('multiplot',False)       # create multiple plots
    multipage = kwargs.get('multipage',False)       # create a multiple page sequence of plots
    multilayout = kwargs.get('multilayout',[1,1])       # layout of multiple plots, [# horizontal, # vertical]
    multilayout = kwargs.get('layout',multilayout)       # layout of multiple plots, [# horizontal, # vertical]
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
        'ch4': {'label': r'CH$_4$', 'type': 'band', 'wavelengths': [[1.1,1.24],[1.28,1.44],[1.6,1.76],[2.2,2.35]]}, \
        'co': {'label': r'CO', 'type': 'band', 'wavelengths': [[2.28,2.39]]}, \
        'tio': {'label': r'TiO', 'type': 'band', 'wavelengths': [[0.76,0.80],[0.825,0.831]]}, \
        'vo': {'label': r'VO', 'type': 'band', 'wavelengths': [[1.04,1.08]]}, \
        'feh': {'label': r'FeH', 'type': 'band', 'wavelengths': [[0.86,0.90],[0.98,1.03],[1.19,1.25],[1.57,1.64]]}, \
        'h2': {'label': r'H$_2$', 'type': 'band', 'wavelengths': [[2.05,2.6]]}, \
        'sb': {'label': r'*', 'type': 'band', 'wavelengths': [[1.6,1.64]]}, \
        'h': {'label': r'H I', 'type': 'line', 'wavelengths': [[1.004,1.005],[1.093,1.094],[1.281,1.282],[1.944,1.945],[2.166,2.166]]},\
        'na': {'label': r'Na I', 'type': 'line', 'wavelengths': [[0.8186,0.8195],[1.136,1.137],[2.206,2.209]]}, \
        'nai': {'label': r'Na I', 'type': 'line', 'wavelengths': [[0.8186,0.8195],[1.136,1.137],[2.206,2.209]]}, \
        'na1': {'label': r'Na I', 'type': 'line', 'wavelengths': [[0.8186,0.8195],[1.136,1.137],[2.206,2.209]]}, \
        'k': {'label': r'K I', 'type': 'line', 'wavelengths': [[0.7699,0.7665],[1.169,1.177],[1.244,1.252]]}, \
        'ki': {'label': r'K I', 'type': 'line', 'wavelengths': [[0.7699,0.7665],[1.169,1.177],[1.244,1.252]]}, \
        'k1': {'label': r'K I', 'type': 'line', 'wavelengths': [[0.7699,0.7665],[1.169,1.177],[1.244,1.252]]}}

    features = kwargs.get('features',[])
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
        print 'plotSpectrum needs at least one Spectrum object to plot'
        return

# if a list is passed, use this list
    elif (len(args) == 1 and isinstance(args[0],list)):
        splist = args[0]
    
# if a set of objects is passed, turn into a list
    else:
        splist = []
        for a in args:
            if a.__class__.__name__ == 'Spectrum':      # a spectrum object
                splist.append([a])
            elif isinstance(a,list):
                splist.append(a)
            else:
                print '\nplotSpectrum: Ignoring input object {} as it is neither a Spectrum object nor a list\n\n'.format(a)

# set up for multiplot
    if (len(splist) == 1):
        multiplot = False
    
# array of lists => force multiplot
    elif (len(splist) > 1 and isinstance(splist[0],list)):
        multiplot = True

# reformat array of spectra of multiplot is used (i.e., user forgot to set)
    if (multiplot == True and splist[0].__class__.__name__ == 'Spectrum'):
        splist = [[s] for s in splist]

    elif (multiplot == False and splist[0].__class__.__name__ == 'Spectrum'):
        splist = [splist]
        
# flatten array if multiplot is not set
    elif (multiplot == False and isinstance(splist[0],list) and len(splist) > 1):
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
        if (sp[0].__class__.__name__ != 'Spectrum'):
            raise ValueError('\nInput to plotSpectrum has wrong format:\n\n{}\n\n'.format(args[0]))
        zeropoint = kwargs.get('zeropoint',[0. for x in range(len(sp))])

# settings that work if the spectrum was read in as legitmate Spectrum object
        try:
            xlabel = kwargs.get('xlabel','{} ({})'.format(sp[0].wlabel,sp[0].wunit))
            ylabel = kwargs.get('ylabel','{} {} ({})'.format(sp[0].fscale,sp[0].flabel,sp[0].funit))
        except:
            xlabel = kwargs.get('xlabel','Wavelength (unknown units)')
            ylabel = kwargs.get('ylabel','Flux (unknown units)')
        xrange = kwargs.get('xrange',[0.85,2.42])
        bound = xrange
        ymax = [s.fluxMax().value for s in sp]
        yrng = kwargs.get('yrange',map(lambda x: x*(numpy.nanmax(ymax)+numpy.nanmax(zeropoint)),[-0.02,1.2]))
        bound.extend(yrng)
        linestyle = kwargs.get('linestyle',['steps' for x in range(len(sp))])
        if (len(linestyle) < len(sp)):
            linestyle.extend(['steps' for x in range(len(sp)-len(linestyle))])

# colors
        colors = kwargs.get('colors',['k' for x in range(len(sp))])
        if (len(colors) < len(sp)):
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
                colors[i] = scalarMap.to_rgba(values[i])
        labels = kwargs.get('labels',['' for x in range(len(sp))])
        if(len(labels) < len(args)):
            labels.extend(['' for x in range(len(sp)-len(labels))])


# show uncertainties
        showNoise = kwargs.get('showNoise',[False for x in range(len(sp))])
        showNoise = kwargs.get('noise',showNoise)
        showNoise = kwargs.get('uncertainty',showNoise)
        if not isinstance(showNoise, list):
            showNoise = [showNoise]
        if (len(showNoise) < len(sp)):
            showNoise.extend([True for x in range(len(sp)-len(showNoise))])

# zero points - by default true
        showZero = kwargs.get('showZero',[True for x in numpy.arange(len(sp))])
        if not isinstance(showZero, list):
            showZero = [showZero]
        if (len(showZero) < len(sp)):
            showZero.extend([True for x in range(len(sp)-len(showZero))])

# GENERATE PLOTS
        if (multipage == True):
            plt_n = plts % nplot
            if (plt_n == 0):# and plts != len(splist)):
#                ax = range(nplot)
#                t = tuple([tuple([i+b*multilayout[1] for i in range(multilayout[1])]) for b in range(multilayout[0])])
#                fig[pg_n], ax = plt.subplots(multilayout[0], multilayout[1], sharex = True, sharey = True)
                fig[pg_n] = plt.figure()
                pg_n += 1
            ax = fig[pg_n-1].add_subplot(multilayout[0], multilayout[1], plt_n+1)
            
# plotting a single plot with all spectra
        else:
            plt.close('all')
#            ax = range(1)
            plt_n = 0
            if (kwargs.get('figsize') != None):
                fig = plt.figure(figsize = kwargs.get('figsize'))
            else:
                fig = plt.figure()
            ax = fig.add_subplot(111)
        
        for ii, a in enumerate(sp):
            flx = [i+zeropoint[ii] for i in a.flux.value]
#stack
            if stack > 0:
                flx = [f + ii*stack for f in flx]
                if kwargs.get('yrange') == None:
                    bound[3] = 1.2 + stack * ii
            

            ax.plot(a.wave.value,flx,color=colors[ii],linestyle=linestyle[ii], zorder = 10, label = labels[ii])  

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
                flxmax = flx
                wvmax = a.wave
            else:
                flxmax = numpy.maximum(flxmax,f(wvmax))

# label features
# THIS NEEDS TO BE FIXED WITH GRETEL'S STUFF
        yoff = 0.02*numpy.nanmax(yrng)
        fontsize = 10-numpy.min([(multilayout[0]*multilayout[1]-1),6])
        for ftr in features:
            ftr = ftr.lower()
            if ftr in feature_labels:
                for ii,waveRng in enumerate(feature_labels[ftr]['wavelengths']):
                    if (numpy.min(waveRng) > bound[0] and numpy.max(waveRng) < bound[1]):
                        x = (numpy.arange(0,nsamples+1.0)/nsamples)* \
                            (numpy.nanmax(waveRng)-numpy.nanmin(waveRng)+0.1)+numpy.nanmin(waveRng)-0.05
                        f = interp1d(wvmax,flxmax,bounds_error=False,fill_value=0.)
                        y = numpy.nanmax(f(x))+0.5*yoff

                        if feature_labels[ftr]['type'] == 'band':
                            ax.plot(waveRng,[y+yoff]*2,color='k',linestyle='-')
                            ax.plot([waveRng[0]]*2,[y,y+yoff],color='k',linestyle='-')
                            ax.text(numpy.mean(waveRng),y+1.5*yoff,feature_labels[ftr]['label'],horizontalalignment='center',fontsize=fontsize)
                        else:
                            for w in waveRng:
                                ax.plot([w]*2,[y,y+yoff],color='k',linestyle='-')
                            ax.text(numpy.mean(waveRng),y+1.5*yoff,feature_labels[ftr]['label'],horizontalalignment='center',fontsize=fontsize)
                            waveRng = [waveRng[0]-0.05,waveRng[1]+0.05]   # for overlap

# update offset
                        foff = [y+3*yoff if (w >= waveRng[0] and w <= waveRng[1]) else 0 for w in wvmax.value]
                        flxmax = [numpy.max([x,y]) for x, y in zip(flxmax, foff)]

# overplot telluric absorption
        bound[3] = numpy.max(flxmax)+2.*yoff
        if (kwargs.get('telluric',False) == True):
            twv = [[1.1,1.2],[1.3,1.5],[1.75,2.0]]
            for waveRng in twv:
                rect = patches.Rectangle((waveRng[0],bound[2]),waveRng[1]-waveRng[0],bound[3]-bound[2],facecolor='grey', alpha=0.1)
                ax.add_patch(rect)
                ax.text(numpy.mean(waveRng),bound[2]+3*yoff,r'$\oplus$',horizontalalignment='center',fontsize=fontsize)


# grid
        if (grid):
            ax.grid()            
        ax.axis(bound)

# axis labels 
        fontsize = 13-numpy.min([(multilayout[0]*multilayout[1]-1),8])
        ax.set_xlabel(xlabel, fontsize = fontsize)
        ax.set_ylabel(ylabel, fontsize = fontsize)
        ax.tick_params(axis='x', labelsize=fontsize)
        ax.tick_params(axis='y', labelsize=fontsize)
    
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

