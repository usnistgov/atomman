AtomMan
=======

Atomistic Manipulation Toolkit
 
AtomMan: the Atomistic Manipulation Toolkit is a Python library for 
creating, representing, manipulating, and analyzing large-scale atomic 
systems of atoms. The focus of the package is to facilitate the rapid design 
and development of simulations that are fully documented and easily adaptable 
to new potentials, configurations, etc.  The code has no requirements that 
limit which systems it can be used on, i.e. it should work on Linux, Mac and 
Windows computers.

Features:

1. Allows for efficient and fast calculations on millions of atoms, each with many freely defined per-atom properties.
2. Create dislocation monopoles and evaluate them with differential displacement and Nye tensor plots.
3. Generate point defects.
4. Call LAMMPS directly from Python and instantly retrieve the resulting data or LAMMPS error statement.
5. Easily convert systems to/from the other Python atomic environments of ASE and PyMatGen.
6. Can create systems based on CIF crystal structure files, and LAMMPS atom and dump files.
7. Built-in unit conversions.

Installation
------------

As of version 1.2, the atomman package is Python 2/3 compatible. It makes heavy use of numpy, so
it's easiest to download a Python environment like Anaconda.

The latest release can be installed using pip::

    pip install atomman

This pip command should install atomman and any other required packages, but
occasionally a requirement may have to be installed separately. The list of required packages are given below.

Alternatively, all code and documentation can be downloaded from GitHub. 
    
    - The stable releases are available at `https://github.com/usnistgov/atomman`_.
    
    - The working development versions are at `https://github.com/lmhale99/atomman`_.
    
Documentation
-------------

Tutorial Jupyter Notebooks can be found in the doc/tutorial directory.  The tutorials come in two flavors:
    
    - The tutorials starting with ##. provide a general overview/example of the various capabilities.
    
    - The tutorials starting with ##.#. give more detailed descriptions and list options available to the tools mentioned in the overview tutorials.
    
Documentation for the code itself can be found in the doc/html directory.

To view the code documentation, use `browser-sync <https://browsersync.io/>`_ which needs `yarn <https://yarnpkg.com/lang/en/docs/install/>`_::

    yarn
    yarn browser-sync doc/html -w


Required packages
-----------------

This is a list of the required Python packages

    - `xmltodict`_ converts XML files to Python dictionaries. Used by 
      DataModelDict.
    
    - `DataModelDict`_ class allowing for easy transformations between 
      XML/JSON/Python representations of structured data models.
      
    - `numericalunits`_ forms the basis for unit conversions.  
      
    - `numpy`_, `scipy`_,  `pandas_`, and `matplotlib`_ Python scientific tools
      for representing, manipulating and plotting data.
    
Optional packages
-----------------

This is a list of additional Python packages that can add functionality

    - `diffpy.Structure`_ - CIF reader. Required for loading systems from
      CIF files.
    
    - `ase`_ - The Atomic Simulation Environment for interacting with small 
      systems and DFT calculations. Required for converting to/from ase.Atoms 
      objects.
    
    - `pymatgen`_ - The Python Materials Genomics package used by the Materials
      Project for DFT calculations. Required for converting to/from 
      pymatgen.Structure objects.
      
    - `cython`_ - Allows for construction of c/Python hybrid code for faster calculations.  Alternate cython versions of some of the calculation heavy functions can be built if cython is installed.
       
.. _https://github.com/usnistgov/atomman: https://github.com/usnistgov/atomman
.. _https://github.com/lmhale99/atomman: https://github.com/lmhale99/atomman
.. _Introduction to atomman: https://github.com/usnistgov/atomman/blob/master/docs/tutorial/1%20Basics.ipynb
.. _Interacting with LAMMPS: https://github.com/usnistgov/atomman/blob/master/docs/tutorial/2%20LAMMPS%20Functionality.ipynb
.. _xmltodict: https://github.com/martinblech/xmltodict
.. _DataModelDict: https://github.com/usnistgov/DataModelDict
.. _numericalunits: https://pypi.python.org/pypi/numericalunits
.. _numpy: http://www.numpy.org/
.. _scipy: https://www.scipy.org/
.. _pandas: http://pandas.pydata.org/
.. _matplotlib: http://matplotlib.org/
.. _diffpy.Structure: http://www.diffpy.org/diffpy.Structure/
.. _ase: https://wiki.fysik.dtu.dk/ase/
.. _pymatgen: http://pymatgen.org/
