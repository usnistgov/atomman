# Standard Python libraries
from __future__ import (absolute_import, print_function,
                        division, unicode_literals)

# http://www.numpy.org/
import numpy as np

# atomman imports
from ... import Atoms, Box, System

def load(cell, symbols=()):
    """
    Convert an spglib cell tuple into an atomman.System
    
    Parameters
    ----------
    cell : tuple
        A tuple containing 3x3 lattice vectors, 3XN box relative positions,
        and N numeric atomic types.
    symbols : list
        The elemental symbols to pair with the unique atom types/numbers.
    
    Returns
    -------
    system : atomman.System
        A atomman representation of a system.
    """
    
    # Separate out cell tuple
    lattice = cell[0]
    positions = cell[1]
    numbers = cell[2]
    
    # Define box
    box = Box(vects = lattice)
    
    # Define atoms
    atoms = Atoms(atype=numbers, pos=positions)
    
    # Build system
    system = System(atoms=atoms, box=box, symbols=symbols, scale=True)
    
    return system