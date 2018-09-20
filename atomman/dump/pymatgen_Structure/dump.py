# Standard Python libraries
from __future__ import (absolute_import, print_function,
                        division, unicode_literals)

# http://www.numpy.org/
import numpy as np

# http://pymatgen.org
try:
    import pymatgen as pmg
    has_pmg = True
except:
    has_pmg = False

def dump(system, symbols=None):
    """
    Convert an atomman.System into a pymatgen.Structure.
    
    Parameters
    ----------
    system : atomman.System
        A atomman representation of a system.
    symbols : tuple, optional
        List of the element symbols that correspond to the atom types.  If not
        given, will use system.symbols if set, otherwise no element content
        will be included.
    
    Returns
    -------
    structure : pymatgen.Structure
        A pymatgen representation of a structure.
    """
    
    assert has_pmg, 'pymatgen not imported'
    
    # Get box/lattice information
    lattice = pmg.Lattice(system.box.vects)
    
    # Get symbols information
    if symbols is None:
        symbols = system.symbols
    symbols = np.asarray(symbols)
    if None in symbols:
        raise ValueError('Symbols needed for all atypes')
    
    # Convert short symbols list to full species list
    atype = system.atoms.atype
    species = symbols[atype-1]
    
    # Get atomic information
    sites = system.atoms_prop(key='pos', scale=True)
    prop = {}
    for p in system.atoms_prop():
        if p != 'atype' and p != 'pos':
            prop[p] = system.atoms_prop(key=p)
    
    # Build structure
    structure = pmg.Structure(lattice, species, sites, site_properties=prop)
    
    return structure