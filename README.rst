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

1. Allows for efficient and fast calculations on millions of atoms, each with
   many freely defined per-atom properties.

2. Create dislocation monopoles and evaluate them with differential
   displacement and Nye tensor plots.

3. Generate point defects.

4. Call LAMMPS directly from Python and instantly retrieve the resulting data
   or LAMMPS error statement.

5. Easily convert systems to/from the other Python atomic environments of ASE
   and PyMatGen.

6. Can create systems based on CIF crystal structure files, and LAMMPS atom and
   dump files.

7. Built-in unit conversions.

Installation
------------

As of version 1.2, the atomman package is Python 2/3 compatible. It makes heavy
use of numpy, so it's easiest to download a Python environment like Anaconda.

The latest release can be installed using pip::

    pip install atomman

This pip command should install atomman and any other required packages, but
occasionally a requirement may have to be installed separately. The list of
required packages are given below.

Alternatively, all code and documentation can be downloaded from GitHub.

- The stable releases are available at
  `https://github.com/usnistgov/atomman <https://github.com/usnistgov/atomman>`__.

- The working development versions are at
  `https://github.com/lmhale99/atomman <https://github.com/lmhale99/atomman>`__.

Documentation
-------------

Web-based documentation for the atomman package is available at
`https://www.ctcms.nist.gov/potentials/atomman <https://www.ctcms.nist.gov/potentials/atomman>`__.

Source code for the documentation can be found in the
`github doc directory <https://github.com/usnistgov/atomman/tree/master/doc/>`__.
The doc directory contains the information both as the source RestructuredText
files and as unformatted HTML. If you download a copy, you can view the HTML
version offline by

    cd {atomman_path}/doc/html
    python -m http.server (for python 3)
    python -m SimpleHttpServer (for python 2)

Then, opening localhost:8000 in a web browser.

The documentation consists of two main components:

1. **Tutorial Jupyter Notebooks:**
   `Online html version <https://www.ctcms.nist.gov/potentials/atomman/tutorial/index.html>`__,
   `Downloadable Notebook version <https://github.com/usnistgov/atomman/tree/master/doc/tutorial>`__.
   The tutorials starting with ##. provide a general overview/example of the
   various capabilities.  The tutorials starting with ##.#. give more detailed
   descriptions and list options available to the tools mentioned in the
   overview tutorials.

2. **Code Documentation:**
   `Online html version <https://www.ctcms.nist.gov/potentials/atomman/atomman.html>`__.
   This provides a rendering of the Python docstrings for the included
   functions and classes.

Required packages
-----------------

This is a list of the required Python packages

- `xmltodict <https://github.com/martinblech/xmltodict>`__

- `DataModelDict <https://github.com/usnistgov/DataModelDict>`__

- `numericalunits <https://pypi.python.org/pypi/numericalunits>`__

- `numpy <http://www.numpy.org/>`__

- `scipy <https://www.scipy.org/>`__

- `pandas <http://pandas.pydata.org/>`__

- `matplotlib <http://matplotlib.org/>`__

Optional packages
-----------------

This is a list of additional Python packages that can add functionality

- `cython <http://cython.org/>`__:
  Allows for construction of c/Python hybrid code for faster calculations.
  Alternate cython versions of some of the calculation heavy
  functions can be built if cython is installed.

- `diffpy.Structure <http://www.diffpy.org/diffpy.Structure/>`__:
  CIF reader. Required for loading systems from CIF files.

- `ase <https://wiki.fysik.dtu.dk/ase/>`__:
  The Atomic Simulation Environment for interacting with small systems
  and DFT calculations. Required for converting to/from ase.Atoms objects.

- `pymatgen <http://pymatgen.org/>`__:
  The Python Materials Genomics package used by the Materials
  Project for DFT calculations. Required for converting to/from
  pymatgen.Structure objects.

- `spglib <https://atztogo.github.io/spglib/python-spglib.html>`__:
  A Python interface to the spglib spacegroup analysis code.  spglib
  can be used to analyze and determine the spacegroup for an atomic system.
  Required for converting to/from spglib.cell objects.

- `lammps <https://lammps.sandia.gov/doc/Python_library.html>`__:
  The Python library interface for the LAMMPS molecular dynamics simulation
  code that allows for a LAMMPS session to be controlled directly from Python.
  Required for direct loading/dumping of a system to/from a LAMMPS session.
