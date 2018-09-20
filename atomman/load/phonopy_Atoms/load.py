# Standard Python libraries
from __future__ import (absolute_import, print_function,
                        division, unicode_literals)

# http://www.numpy.org/
import numpy as np

# atomman imports
from ... import Atoms, Box, System

# https://atztogo.github.io/phonopy/
try:
    import phonopy
    has_phonopy = True
except:
    has_phonopy = False

def load(phonopyatoms, symbols=None, prop={}):
    """
    Convert a phonopy.structure.atoms.PhonopyAtoms into an atomman.System.
    
    Parameters
    ----------
    phonopyatoms : phonopy.structure.atoms.PhonopyAtoms
        A phonopy representation of a collection of atoms, based on ase.Atoms.
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
    
    assert has_phonopy, 'phonopy not imported'
    
    # Get box/cell information
    box = Box(vects = phonopyatoms.get_cell())
    pbc = phonopyatoms.get_pbc()
    
    # Get element information
    all_elements = np.array(phonopyatoms.get_chemical_symbols())
    elements, atype = np.unique(all_elements, return_inverse = True)
    if symbols is None:
        symbols = elements

    # Get atomic information
    prop['atype'] = atype + 1
    prop['pos'] = phonopyatoms.get_positions()
    atoms = Atoms(prop=prop)
    
    # Build system
    system = System(atoms=atoms, box=box, pbc=pbc, symbols=symbols)
    
    return system