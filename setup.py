import os
from distutils.core import setup

# Utility function to read the README file.
# Used for the long_description.  It's nice, because now 1) we have a top level
# README file and 2) it's easier to type in the README file than to put a raw
# string in below ...
def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

setup(
    name = "splat",
    version = "0.3.3",
    author = "Adam Burgasser",
    author_email = "aburgasser@gmail.com",
    download_url = "http://pono.ucsd.edu/splat/code/splat.py",
    keywords = ["astronomy", "astrophysics", "spectroscopy", "near infrared", "spex", "brown dwarfs"],
    description = ("SpeX Prism Library Analysis Toolkit"),
    license = "BSD",
    url = "http://pono.ucsd.edu/splat",
#    packages=['splat'],
    long_description=read('README'),
    classifiers=[
        "Development Status :: 2 - Pre-Alpha",
        "Intended Audience :: Science/Research",
        "License :: Public Domain",
        "Natural Language :: English",
        "Programming Language :: Python :: 2.7",
        "Topic :: Scientific/Engineering :: Astronomy",
    ],
)

