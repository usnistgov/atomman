# Standard Python libraries
from __future__ import (absolute_import, print_function,
                        division, unicode_literals)

# http://www.numpy.org/
import numpy as np

# atomman imports
from ... import Atoms, Box, System

#http://pymatgen.org
try:
    import pymatgen as pmg
    has_pmg = True
except:
    has_pmg = False

def load(structure, symbols=None):
    """
    Convert a pymatgen.Structure into an atomman.System.
    
    Parameters
    ----------
    structure : pymatgen.Structure
        A pymatgen representation of a structure.
    symbols : tuple, optional
        Allows the list of element symbols to be assigned during loading.
        Useful if the symbols for the model differ from the standard element
        tags.
    
    Returns
    -------
    system : atomman.System
        A atomman representation of a system.
    """
    
    assert has_pmg, 'pymatgen not imported'
    
    # Get box/lattice information
    box = Box(vects = structure.lattice.matrix)
    
    # Get element information
    all_elements =  np.array([str(symbol) for symbol in structure.species])
    elements, atype = np.unique(all_elements, return_inverse = True)
    if symbols is None:
        symbols = elements

    # Get atomic information
    prop = structure.site_properties
    prop['atype'] = atype + 1
    prop['pos'] = structure.cart_coords
    atoms = Atoms(prop=prop)
    
    # Build system
    system = System(atoms=atoms, box=box, symbols=symbols)
    
    return system