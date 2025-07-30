Installation and Dependencies
===============================================

SPLAT is best forked from the github site https://github.com/aburgasser/splat/, 
which is updated on a regular basis. The command is:

>> git clone https://github.com/aburgasser/splat.git

Note that ``pip`` installation IS NOT CURRENTLY WORKING, so best to install from github.

In addition, you will need to set up your environment variables to find the SPLAT code

For MAC/UNIX:
-------------
	- Best: create an environment variable called SPLAT_PATH and set it to the root directory of the SPLAT code; in bash environment, add the line ``export SPLAT_PATH=[path to splat]`` to your .bashrc or .bash_profile, where [path to splat] is the full path to your SPLAT python code (e.g., /Users/adam/python_codes/splat/). If you have a newwer Mac OS, which uses the zsh environment, add the same line to your .zshrc or .zsh_profile

	- If you use the PYTHONPATH environment variable, add the root directory of the SPLAT code to it by adding the line ``export PYTHONPATH=[path to splat]:${PYTHONPATH}`` to your .bashrc / .bash_profile or .zshrc / .zsh_profile files

	- add the root directory of the SPLAT code to your system PATH variable by adding the line ``export PATH=[path to splat]:${PATH}`` to your .bashrc / .bash_profile or .zshrc / .zsh_profile files

For Windows:
------------
	- SPLAT has been found to work fine in the Anaconda environment, and you can set environmental variables by entering the following lines on your Anaconda terminal:

		* conda env config vars set SPLAT_PATH=[path to splat]
		* conda env config vars set PYTHONPATH=[path to splat]

	where again [path to splat] is the full path to your SPLAT python code (e.g., /Users/adam/python_codes/splat/).


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



*Search*


* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

