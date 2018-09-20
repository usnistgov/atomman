# Standard Python libraries
from __future__ import (absolute_import, print_function,
                        division, unicode_literals)

# http://www.numpy.org/
import numpy as np

# atomman imports
from . import dvect

def displacement(system_0, system_1, box_reference='final', code=None):
    """
    Compute the displacement vectors between all matching atoms for two systems.
    
    Parameters
    ----------
    system_0 : atomman.System
        The initial system to calculate displacements from.
    system_1 : atomman.System
        The final system to calculate displacements to.
    box_reference : str or None
        Specifies which system's boundary conditions to use.  'initial' uses
        system_0's box and pbc.  'final' uses system_1's box and pbc (Default)
        None computes the straight difference between the positions without
        accounting for periodic boundaries.
    code : str, optional
        Option for specifying which code version of dvect to use (see dvect's
        documentation for values).
    
    Returns
    -------
    numpy.ndarray
        The displacement vectors for all atoms.
        
    Raises
    ------
    ValueError
        If the systems have different numbers of atoms or for invalid
        box_reference values.
    """
    if system_0.natoms != system_1.natoms:
        raise ValueError('systems have different number of atoms')
    
    if box_reference == 'final':
        disp = dvect(system_0.atoms.pos, system_1.atoms.pos, system_1.box, system_1.pbc, code=code)
    elif box_reference == 'initial':
        disp = dvect(system_0.atoms.pos, system_1.atoms.pos, system_0.box, system_0.pbc, code=code)
    elif box_reference is None:
        disp = system_1.atoms.pos - system_0.atoms.pos
    else:
        raise ValueError("box_reference must be 'final', 'initial', or None")
    
    return disp
    