import os
import splat
import sys
from bokeh import mpl
from bokeh.models import HBox, VBox, Paragraph, Range1d
from bokeh.models.widgets import Panel, Tabs
from bokeh.embed import components
from flask import Flask, render_template, request

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
		if request.form['submit'] == 'Load Spectra':
			try:
				path = request.form['path']
				sp = splat.Spectrum(path)
			except:
				return render_template('input.html',  error = "Could not load spectra")
		

		elif request.form['submit'] == 'Load Spectra ':
			try:
				name = request.form['name']
				sp = splat.getSpectrum(shortname = name)

			except:
				return render_template('input.html',  error = "Could not load spectra")
	
			

		elif request.form['submit'] == 'Load Spectra  ':
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
				if len(sp) == 0:
					return render_template('input.html', error = "Could not load spectra")
			except:
				return render_template('input.html', error = "Could not load spectra")
			
		
		elif request.form['submit'] == 'Get Lucky!':
			sp = splat.getSpectrum(lucky=True)

		if len(sp) == 0:
			return render_template('input.html', error = "Could not load spectra")

		try:
			tab = []

			for i in range(len(sp)):
				spectral_type = splat.classifyByStandard(sp[i])[0]
				mpl_fig = splat.plotSpectrum(sp[i], web=True, uncertainty = True, mdwarf=True)[0]
				bokehfig = mpl.to_bokeh(fig=mpl_fig)
				bokehfig.set(x_range=Range1d(.8,2.4))
				sys.stdout = open("out1.txt", "w") 
				sp[i].info()
				sys.stdout = sys.__stdout__
				with open("out1.txt", "r") as f:
					content = f.read() 
				p = Paragraph(text=content, width=200, height=100)
				widget = VBox(bokehfig, p)
				tab.append(Panel(child=widget, title=str(spectral_type)+ " Star"))
				
			tabs = Tabs(tabs= tab)
			script, div_dict = components({"plot" : tabs})
		except:
				return render_template('input.html', error = "Error Plotting Spectra")
		
		return render_template('out.html', star_type = spectral_type, script=script,  div=div_dict)	

if __name__ == '__main__':
	port = int(os.environ.get('PORT', 5000))
	app.run(host='0.0.0.0', port=port, debug=False)
		