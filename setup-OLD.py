
#!/usr/bin/env python
# Licensed under a 3-clause BSD style license - see LICENSE.rst

from __future__ import unicode_literals
from setuptools import setup, find_packages
from splat.initialize import VERSION
#from sphinx.setup_command import BuildDoc
#cmdclass = {'build_sphinx': BuildDoc}


name = "SPLAT"
version = VERSION

setup(
  name=name,
  version=version,
  packages = find_packages(),
#  cmdclass = cmdclass,
#  packages = find_packages(exclude=['docs','tests']),

  # Project uses reStructuredText, so ensure that the docutils get
  # installed or upgraded on the target machine
  install_requires = [
    'asdf-astropy',
#    'sphinx==6.2.1',    
    'astropy',
    'astroquery',
    'bokeh',
    'corner',
    'emcee',
    'flask',
    'matplotlib',
    'numpy',
    'pandas',
    'requests',
    'scipy',
  ],

  package_dir = {'splat': 'splat'},    
  package_data = {'splat': ['db/*','docs/*','resources/*','tutorials/*','build/*']},
  include_package_data=True,

  zip_safe = True,
  use_2to3 = False,
  classifiers=[
      'Development Status :: 3 - Alpha',
      'Intended Audience :: Science/Research',
      'License :: OSI Approved :: MIT License',
      'Operating System :: OS Independent',
      'Programming Language :: Python :: 2',
      'Programming Language :: Python :: 2.6',
      'Programming Language :: Python :: 2.7',
      'Programming Language :: Python :: 3',
      'Programming Language :: Python :: 3.3',
      'Programming Language :: Python :: 3.4',
      'Programming Language :: Python :: 3.5',
      'Programming Language :: Python :: 3.6',
      'Topic :: Scientific/Engineering :: Astronomy',
      'Topic :: Scientific/Engineering :: Physics'
  ],

  # metadata for upload to PyPI
  author = "Adam Burgasser",
  author_email = "aburgasser@ucsd.edu",
  description = "SpeX Prism Library Analysis Toolkit (SPLAT)",
  long_description = "Code base for the SpeX Prism Library Analysis Toolkit (SPLAT), which provides access to near-infrared spectral data for over 3000 sources obtained with the NASA/SpeX instrument primarily of low-mass stars and brown dwarfs; relevant atmosphere and evolutionary models; and analysis functions for spectral measurements and emprical calibrations.",
  license = "MIT",
#    download_url='%s/astropy-%s.tar.gz' % (DOWNLOAD_BASE_URL, VERSION),
  keywords = ['splat','spectroscopy', 'spectral analysis', 'astronomy','astrophysics',\
              'ultracool dwarfs','low mass stars', 'brown dwarfs', 'spex','prism', 'classification'],
  url = "http://splat.physics.ucsd.edu/splat/",   # project home page, if any

  # command_options={
  #       'build_sphinx': {
  #           'project': ('setup.py', name),
  #           'version': ('setup.py', version),
  #           'source_dir': ('setup.py', './docs'),
  #           'build_dir': ('setup.py', './docs')}},
)
