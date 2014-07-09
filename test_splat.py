#!/usr/bin/python -tt

# code to test SPLAT functionality
import splat
import sys

sys.stderr.write('\n\n>>>>>>>>>>>> TESTING SPLAT CODE <<<<<<<<<<<<\n\n')
sp = splat.loadSpectrum('spex_prism_0415-0935_030917.txt')
sp.info()
sys.stderr.write('...loadSpectrum successful\n\n')
sp = splat.getSpectrum(shortname='0415-0935')[0]
sp.info()
sys.stderr.write('...getSpectrum successful\n\n')
ind = splat.measureIndexSet(sp,set='burgasser')
sys.stderr.write('Spectral indices:\n')
for k in ind.keys():
	print '{:s}: {:4.3f}+/-{:4.3f}'.format(k, ind[k][0], ind[k][1])
sys.stderr.write('...index measurement successful\n')
spt, spt_e = splat.classifyByIndex(sp,set='burgasser')
sys.stderr.write('\n...index classification of 2MASS J0415-0935 = {:s}+/-{:2.1f}; successful\n'.format(splat.typeToNum(spt),spt_e))
sp.fluxCalibrate('2MASS J',15.0)
mag = splat.filterMag(sp,'MKO J')
sys.stderr.write('\n...apparent magnitude MKO J = {:3.2f} from 2MASS J = 15.0; filter calibration successful\n'.format(mag))
mdl = splat.loadModel(teff=750,logg=4.5,set='btsettl08')
sys.stderr.write('\n...subsampled model generation successful\n')
mdl.normalize()
sp.normalize()
sys.stderr.write('\n...normalization successful\n')
splat.plotSpectrum(sp,mdl,colors=['k','r'],title='If this appears everything is OK: close window')
sys.stderr.write('\n...plotting successful\n')
sys.stderr.write('\n>>>>>>>>>>>> SPLAT TEST SUCCESSFUL; HAVE FUN! <<<<<<<<<<<<\n\n')
