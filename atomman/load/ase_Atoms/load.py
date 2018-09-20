# Standard Python libraries
from __future__ import (absolute_import, print_function,
                        division, unicode_literals)

# http://www.numpy.org/
import numpy as np

# atomman imports
from ... import Atoms, Box, System

# https://wiki.fysik.dtu.dk/ase/
try:
    import ase
    has_ase = True
except:
    has_ase = False

def load(aseatoms, symbols=None, prop={}):
    """
    Convert an ase.Atoms into an atomman.System.
    
    Parameters
    ----------
    aseatoms : ase.Atoms
        An ase representation of a collection of atoms.
    symbols : tuple, optional
        Allows the list of element symbols to be assigned during loading.
        Useful if the symbols for the model differ from the standard element
        tags.
    prop : dict, optional
        Dictionary containing any extra per-atom properties to include.
        
    Returns
    -------
    system : atomman.System
        A atomman representation of a system.
    """
    
    assert has_ase, 'ase not imported'
    
    # Get box/cell information
    box = Box(vects = aseatoms.get_cell())
    pbc = aseatoms.get_pbc()
    
    # Get element information
    all_elements = np.array(aseatoms.get_chemical_symbols())
    elements, atype = np.unique(all_elements, return_inverse = True)
    if symbols is None:
        symbols = elements

    # Get atomic information
    prop['atype'] = atype + 1
    prop['pos'] = aseatoms.get_positions()
    atoms = Atoms(prop=prop)
    
    # Build system
    system = System(atoms=atoms, box=box, pbc=pbc, symbols=symbols)
    
    return system