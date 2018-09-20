"""
Attributes
----------
rootdir : str
    The absolute path to the iprPy package's root directory used to locate
    contained data files.
"""

# Standard Python libraries
from __future__ import absolute_import, print_function, division
import os

# Define rootdir
rootdir = os.path.dirname(os.path.abspath(__file__))

# Read version from VERSION file
with open(os.path.join(rootdir, 'VERSION')) as version_file:
    __version__ = version_file.read().strip()

# atomman imports
from . import compatibility
from . import unitconvert
from . import tools
from . import lammps
from .dump import *
from .dump import __all__ as dump_all
from .core import *
from .core import __all__ as core_all
from .load import *
from .load import __all__ as load_all
from . import plot
from . import defect

__all__ = ['__version__', 'rootdir'] + dump_all + core_all + load_all
__all__.sort()

# Define default working units
unitconvert.reset_units(length = 'angstrom', mass = 'amu', energy='eV', charge='e')
