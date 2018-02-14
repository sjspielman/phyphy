.. pyvolve documentation master file, created by
   sphinx-quickstart on Mon Jan 19 10:26:47 2015.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

``phyphy`` documentation 
====================================================================

`phyphy <http://github.com/sjspielman/phyphy>`_ is a **P**\ ython package for facilitating  **HyPhy** `(http://hyphy.org/) <http://hyphy.org/>`_ analysis execution and parsing. ``phyphy`` is compatible with **HyPhy version >= 2.3.7**. 

Overview
-----------

``phyphy`` can automate the execution and output parsing for the following standard analyses:

+ `FEL (Fixed Effects Likelihood) <http://hyphy.org/methods/selection-methods/#fel>`_ 
+ `MEME (Mixed Effects Model of Evolution) <http://hyphy.org/methods/selection-methods/#meme>`_ 
+ `SLAC (Single-Likelihood Ancestor Counting) <http://hyphy.org/methods/selection-methods/#slac>`_ 
+ `FUBAR (Fast, Unconstrained Bayesian AppRoximation) <http://hyphy.org/methods/selection-methods/#fubar>`_ 
+ `aBSREL (adaptive Branch-Site Random Effects Likelihood) <http://hyphy.org/methods/selection-methods/#absrel>`_ 
+ `BUSTED (Branch-Site Unrestricted Statistical Test for Episodic Diversification) <http://hyphy.org/methods/selection-methods/#busted>`_ 
+ `RELAX (Test for selection RELAXation) <http://hyphy.org/methods/selection-methods/#relax>`_ 
+ `LEISR (Likelihood Estimation of Individual Site Rates) <https://doi.org/10.1101/206011>`_


``phyphy`` usage relies on three primary modules:

+ ``Hyphy``
    + This module can be **optionally** used to specify that a *non-canonical* (i.e. not installed to ``/usr/local/``) installation or a local build of HyPhy. This module can additionally be used to specify that HyPhy be run quietly and/or without outputting its standard log files, ``messages.log`` and ``errors.log``.
+ ``Analysis``
    + This module contains the analysis methods to execute, named according to the analysis. For example, the ``FEL`` class within the ``Analysis`` module would be used to execute an FEL analysis.
+ ``Extractor``
    + This module contains functionality to parse a given analysis output. ``Extractor`` makes it simple to extract information from a given analysis JSON output, including:
        + Fitted model parameters
        + Fitted phylogenies
        + `Newick-extended format <https://home.cc.umanitoba.ca/~psgendb/doc/atv/NHX.pdf>`_ phylogenies with branch features for downstream visualization
        + CSV files, for methods FEL, MEME, SLAC, FUBAR, LEISR, aBSREL

API Documentation
---------------------

.. toctree::
   :maxdepth: 2

   modules


Download/Installation
--------------------------

``phyphy`` is freely available under a FreeBSD license. The easiest way to install ``phyphy`` is using ``pip`` (or ``pip3``, for use with Python3 instead of Python2.7). 

.. code-block:: bash

   sudo pip install phyphy
   
Alternatively, the most recent release of Phyphy is available for download from `https://github.com/sjspielman/phyphy/releases <https://github.com/sjspielman/pyvolve/releases>`_. Once Pyvolve has been downloaded, navigate to the directory in the terminal. To install for all users, enter this command:


.. code-block:: bash

   sudo python setup.py install


Alternatively, to install locally for a specific user (or if you do not have root privileges), enter this command:

.. code-block:: python
   
   python setup.py install --user 


Optional tests (Python2.7 compatible only!) may be run with the command (which may or may not require ``sudo``, depending on your install choice), 

.. code-block:: python
   
   python setup.py test

Further note that ``phyphy`` has the following dependencies which must be installed (``pip`` should take care of these for you): 

1. `BioPython <http://biopython.org/wiki/Main_Page>`_
2. `ete3 <http://etetoolkit.org/>`_


Issues and Questions
---------------------

Please file all bugs, issues, feature requests, and/or questions in the Issues section of the ``phyphy`` github repository: `https://github.com/sjspielman/phyphy/issues <https://github.com/sjspielman/phyphy/issues>`_. Note that you will need a github account to file.


Citation
---------
Forthcoming
       


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

