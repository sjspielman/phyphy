#!/usr/bin/env python

##############################################################################
##  phyhy: Python HyPhy: Facilitating the execution and parsing of standard HyPhy analyses.
##
##  Written by Stephanie J. Spielman (stephanie.spielman@temple.edu) 
##############################################################################

'''
Setup.py script (uses setuptools) for building, testing, and installing phyphy.

To build and install the package as root (globally), enter (from this directory!) - 
    sudo python setup.py build
    sudo python setup.py test   # OPTIONAL BUT RECOMMENDED. Please file an issue for any failed tests! 
    sudo python setup.py install
    

To install for a particular user (locally), enter - 
    python setup.py build
    python setup.py test   # OPTIONAL BUT RECOMMENDED. Please file an issue for any failed tests! 
    python setup.py build --user # where user is the current account
'''

_VERSION="0.1"

from setuptools import setup
setup(name = 'phyphy', 
    version = _VERSION, 
    description = 'Facilitating the execution and parsing of standard HyPhy (>=2.3.7) analyses',
    author = 'Stephanie J. Spielman', 
    author_email = 'stephanie.spielman@temple.edu', 
    url = 'https://github.com/sjspielman/phyphy',
    download_url = 'https://github.com/sjspielman/phyphy/tarball/' + _VERSION,
    platforms = 'Tested on Mac OS X.',
    package_dir = {'phyphy':'src'},
    packages = ['phyphy'],
    package_data = {'tests': ['test_jsons/*']},
    install_requires=['Biopython', 'ete3>=3.1'],
    test_suite = "tests"
)
