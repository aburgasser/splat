# -*- coding: utf-8 -*-
from __future__ import print_function, division

"""
.. note::
         These are the web functions for SPLAT 
"""

# imports: internal
import os
import sys

# imports: external
from bokeh import mpl
from bokeh.models import HBox, VBox, Paragraph, Range1d
from bokeh.models.widgets import Panel, Tabs
from bokeh.embed import components
from flask import Flask, render_template, request

import splat
from splat.initialize import *
from splat.utilities import *
import splat.plot as splot

app = Flask(__name__)

app.config['UPLOAD_FOLDER'] = '/tmp/'
# These are the extension that we are accepting to be uploaded
app.config['ALLOWED_EXTENSIONS'] = set(['txt', 'pdf', 'png', 'jpg', 'jpeg', 'gif', 'fits'])


@app.route('/')
def load_html():
	return render_template('input.html', error='')
	

@app.route('/values', methods=['GET','POST'])
def load_spectra():
	#'/Users/guillaumeshippee/Desktop/splat-master/reference/Spectra/10001_10443.fits' #change
	if request.method == 'GET':
		return render_template('input.html', error='')
	else:
#		for k in list(request.form.keys()):
#			print(k,request.form[k])

# search by file "upload"	
		if request.form['submit'] == 'Load File':
			try:
				path = request.form['path']
				sp = splat.Spectrum(file=str(path))
				sp = [sp]
			except:
				return render_template('input.html',  error = "\n\nProblem with file upload button")
		
# search by file path specification	
		if request.form['submit'] == 'Load File ':
			try:
				path = request.form['path']
				sp = splat.Spectrum(file=str(path))
				sp = [sp]
			except:
				return render_template('input.html',  error = "\n\nProblem with file path specification")
		
# search by spectrum key	
		if request.form['submit'] == 'Load by ID':
			try:
				sp = splat.Spectrum(int(str(request.form['key'])))
				sp = [sp]
			except:
				return render_template('input.html',  error = "\n\nProblem with key specification")
		
# search by date observed	
		if request.form['submit'] == 'Load by Date':
			try:
				sp = splat.getSpectrum(date = str(request.form['date']))
			except:
				return render_template('input.html',  error = "\n\nProblem with key specification")
		
# search by shortname	
		elif request.form['submit'] == 'Load by Shortname':
			try:
				sp = splat.getSpectrum(shortname = str(request.form['shortname']))
			except:
				return render_template('input.html',  error = "\n\nProblem with specifying file by shortname")
	
# search by name	
		elif request.form['submit'] == 'Load by Name':
			try:
				sp = splat.getSpectrum(name = str(request.form['name']))
			except:
				return render_template('input.html',  error = "\n\nProblem with specifying file by name")
	
			
# search by options
		elif request.form['submit'] == 'Load by Options':
			sp1 = request.form['sp1']
			sp2 = request.form['sp2']
			mag1 = request.form['mag1']
			mag2 =  request.form['mag2']


			if sp2 == '':
				sp = sp1
			elif sp1 == '':
				sp2 = sp1
			elif sp1 == '' and sp2 == '':
				sp = ''
			else:
				sp = [sp1, sp2]

			if mag2 == '':
				mag = mag1
			elif sp1 == '':
				mag2 = mag1
			elif sp1 == '' and sp2 == '':
				mag = ''
			else:
				mag = [mag1, mag2]

			kwargs = {'spt': sp, 'jmag' : mag, 'snr' : request.form['snr'], 'date' : request.form['date'] }
			kwargs = {k: v for k, v in kwargs.items() if v}

			try:
				sp = splat.getSpectrum(**kwargs)
			except:
				return render_template('input.html', error = "\n\nProblem with option search")
			
# lucky pull	
		elif request.form['submit'] == 'Get Lucky!':
			sp = splat.getSpectrum(lucky=True)

		if len(sp) == 0:
			return render_template('input.html', error = "\n\nNo spectra matched search constratins")

		try:
			tab = []

			for s in sp:
				spectral_type = splat.classifyByStandard(s)[0]
				mpl_fig = splot.plotSpectrum(s, web=True, uncertainty = True)[0]
				bokehfig = mpl.to_bokeh(fig=mpl_fig)
				bokehfig.set(x_range=Range1d(.8,2.4),y_range=Range1d(0,s.fluxMax().value*1.2))
#				sys.stdout = open("out1.txt", "w") 
#				sys.stdout = sys.__stdout__
#				with open("out1.txt", "r") as f:
#					content = f.read() 
				content = s.info(printout=False)
#				print(content)
				p = Paragraph(text=content)
				widget = VBox(bokehfig, p)
				tab.append(Panel(child=widget, title=str(s.name)))
#				tab.append(Panel(child=widget, title=str(spectral_type)+ " Star"))
				
			plottabs = Tabs(tabs= tab)
			script, div_dict = components({"plot" : plottabs})
		except:
				return render_template('input.html', error = "\n\nProblem Plotting Spectra")
		
		return render_template('output.html', star_type = spectral_type, script=script,  div=div_dict)	

if __name__ == '__main__':
	port = int(os.environ.get('PORT', 5000))
	app.run(host='0.0.0.0', port=port, debug=False)

