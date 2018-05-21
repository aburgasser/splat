Installation and Dependencies
===============================================

SPLAT is best forked from the github site https://github.com/aburgasser/splat/, 
which is updated on a regular basis. There is a possibility of ``pip`` installation, but this is not working quite yet.

There are two additional steps to complete to get full functionality:

- Copy the file ``.splat_access`` into your home directory - this is your access key if you have priveleged access to unpublished data. 

- Set up your environment variables to find the SPLAT code; this can be done in three ways (in the order of which the code looks for this directory):

	- Best: create an environment variable called SPLAT_PATH and set it to the root directory of the SPLAT code (in bsh environment, add the line ``export SPLAT_PATH=/Users/adam/projects/splat`` to your .bashrc or .bash_profile)

	- If you use the PYTHONPATH environment variable, add the root directory of the SPLAT code to it (in bsh environment add the line ``export PYTHONPATH=/Users/adam/projects/splat:${PYTHONPATH}`` to your .bashrc or .bash_profile)

	- add the root directory of the SPLAT code to your system PATH variable (in bsh environment add the line ``export PATH=/Users/adam/projects/splat:${PATH}`` to your .bashrc or .bash_profile)


The SPLAT code uses the following external packages that are contained in the `Ananconda <https://docs.continuum.io/>`_ installation:

* `astropy <http://www.astropy.org/>`_
* `bokeh <http://bokeh.pydata.org/en/latest/>`_
* `flask <http://flask.pocoo.org/>`_
* `matplotlib <http://matplotlib.org/>`_
* `numpy <http://www.numpy.org/>`_
* `pandas <http://pandas.pydata.org/>`_
* `requests <http://docs.python-requests.org/en/master/>`_
* `scipy <https://www.scipy.org/>`_

The following external packages must be installed separately (using pip):

* `astroquery <https://astroquery.readthedocs.io/en/latest/>`_
* `corner <http://corner.readthedocs.io/en/latest/>`_
* `emcee <http://dan.iel.fm/emcee/current/>`_
* `triangle <https://pypi.python.org/pypi/triangle_plot>`_

It also use the following internal (python) packages (no installation should be needed):

* base64
* copy
* csv
* datetime
* division (for Python v2.X)
* glob
* math
* os
* print_function (for Python v2.X)
* random
* re
* shutil
* sys
* string (for Python v2.X)
* time
* warnings


SPLAT has been tested in Python 2.7.X and 3.5.X. 

SPLAT has not yet reached v1.0, so bugs are common. Please help us squish them by 
sending bug reports to aburgasser@ucsd.edu or start an issue on the github site.



* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

