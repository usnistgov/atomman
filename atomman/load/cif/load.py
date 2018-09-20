# Standard Python libraries
from __future__ import (absolute_import, print_function,
                        division, unicode_literals)
import os

# http://www.numpy.org/
import numpy as np

# http://www.diffpy.org/diffpy.structure/
try:
    import diffpy.Structure
    has_diffpy = True
except:
    has_diffpy = False

# atomman imports
from ... import Atoms, Box, System
from ...tools import uber_open_rmode

def load(cif, symbols=None):
    """
    Reads in a CIF crystal file.  Requires diffpy.Structure.
    
    Parameters
    ----------
    cif : str or file-like object
        The cif content to read.
    symbols : tuple, optional
        Allows the list of element symbols to be assigned during loading.
        Useful if the symbols for the model differ from the standard element
        tags.
    
    Returns
    -------
    system : atomman.System
        An atomman representation of a system.
    """
    
    assert has_diffpy, 'diffpy.Structure not imported'
    
    # Read in cif file to diffpy Structure
    dps = diffpy.Structure.structure.Structure()
    with uber_open_rmode(cif) as f:
        dps.readStr(f.read())
    
    all_elements = dps.element
    elements, all_atype = np.unique(all_elements, return_inverse=True)
    all_atype += 1
    all_pos = dps.xyz
    
    atype = []
    pos = []
    for a1, p1 in zip(all_atype, all_pos):
        noMatch = True
        for a2, p2 in zip(atype, pos):
            if a1 == a2 and np.allclose(p1, p2):
                noMatch = False
                break
        if noMatch:
            atype.append(a1)
            pos.append(p1)
    
    atoms = Atoms(atype=atype, pos=pos)
    
    box = Box(a=dps.lattice.a,
              b=dps.lattice.b,
              c=dps.lattice.c,
              alpha=dps.lattice.alpha,
              beta=dps.lattice.beta,
              gamma=dps.lattice.gamma)
    
    return System(atoms=atoms, box=box, scale=True, symbols = elements)