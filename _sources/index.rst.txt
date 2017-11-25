.. pyvolve documentation master file, created by
   sphinx-quickstart on Mon Jan 19 10:26:47 2015.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to phyphy's documentation!
===================================

TEXT TEXT TEXT 

--------------------------

phyphy is freely available under a TBD license. The easiest way to install phyphy is using `pip` or `easy_install`:

.. code-block:: bash

   sudo pip install phyphy
   #or
   sudo easy_install install phyphy
   
Alternatively, the most recent release of Phyphy is available for download from `https://github.com/sjspielman/phyphy/releases <https://github.com/sjspielman/pyvolve/releases>`_. Once Pyvolve has been downloaded, navigate to the directory in the terminal. To install for all users, enter this command:


.. code-block:: bash

   sudo python setup.py install


Alternatively, to install locally for a specific user (or if you do not have root privileges), enter this command:

.. code-block:: python
   
   python setup.py install --user 


Optional tests may be run with the command (which may or may not require ``sudo``, depending on your install choice), 

.. code-block:: python
   
   python setup.py test


Dependencies
-------------

Phyphy has several dependencies which must be installed: 

1. `BioPython <http://biopython.org/wiki/Main_Page>`_
2. `Dendropy <https://www.dendropy.org/>`_

Note that installing Phyphy with pip and/or easy_install will give you these dependencies if they are missing.

Issues and Questions
---------------------

Please file all bugs, issues, feature requests, and/or questions in the Issues section of Phyphy's github repository: `https://github.com/sjspielman/phyphy/issues <https://github.com/sjspielman/phyphy/issues>`_. Note that you will need a github account to file.


Citation
---------
Forthcoming


Contents
-----------

.. toctree::
   :maxdepth: 2

   modules

        


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

