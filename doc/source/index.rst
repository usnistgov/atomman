=======================================
atomman: Atomistic Manipulation Toolkit
=======================================

Description
===========

atomman is a Python library for creating, representing, manipulating, and analyzing large-scale atomic 
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
============

As of version 1.2, the atomman package is Python 2/3 compatible. It makes heavy use of numpy, so
it's easiest to download a Python environment like Anaconda.

The latest release can be installed using pip::

    pip install atomman

This pip command should install atomman and any other required packages, but
occasionally a requirement may have to be installed separately. The list of required packages are given below.

Alternatively, all code and documentation can be downloaded from GitHub. 
    
    - The stable releases are available at `https://github.com/usnistgov/atomman`_.
    
    - The working development versions are at `https://github.com/lmhale99/atomman`_.

Tutorials
=========

.. toctree::
   :maxdepth: 2
   
   tutorial/index
    
    
Code Documentation
==================

.. toctree::
   :maxdepth: 3
   
   atomman

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`

.. _https://github.com/usnistgov/atomman: https://github.com/usnistgov/atomman
.. _https://github.com/lmhale99/atomman: https://github.com/lmhale99/atomman