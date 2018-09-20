# Standard Python libraries
from __future__ import (absolute_import, print_function,
                        division, unicode_literals)

# http://www.numpy.org/
import numpy as np

# atomman imports
from ... import Atoms, Box, System
from ...compatibility import range
from ...tools import uber_open_rmode

def load(poscar, symbols=None, prop={}):
    """
    Reads a poscar-style coordination file for a system.
    
    Parameters
    ----------
    poscar : str or file-like object
        The POSCAR content to read
    symbols : tuple, optional
        Allows the list of element symbols to be assigned during loading.
        Useful if the symbols for the model differ from the standard element
        tags or if the poscar file has no elemental information.
    prop : dict, optional
        Dictionary containing any extra per-atom properties to include.
    
    Returns
    -------
    system : atomman.System
        The system object associated with the poscar file.
    """
    
    # Read in all lines of the file
    with uber_open_rmode(poscar) as f:
        lines = f.readlines()
    
    # Interpret box information
    box_scale = float(lines[1])
    avect = np.array(lines[2].split(), dtype='float64') * box_scale
    bvect = np.array(lines[3].split(), dtype='float64') * box_scale
    cvect = np.array(lines[4].split(), dtype='float64') * box_scale
    box = Box(avect=avect, bvect=bvect, cvect=cvect)
    
    # Read in elements, number of types, and style info
    try:
        typenums = np.array(lines[5].split(), dtype='int64')
        elements = [None for n in range(len(typenums))]
        style = lines[6]
        start_i = 7
    except:
        elements = lines[5].split()
        typenums = np.array(lines[6].split(), dtype='int64')
        style = lines[7]
        start_i = 8
    
    # Decode elements if needed
    for i in range(len(elements)):
        try:
            elements[i] = elements[i].decode('UTF-8')
        except:
            pass
    
    # Handle elements and symbols
    if symbols is None:
        symbols = elements
    
    # Build atype list
    atype = np.array([], dtype='int64')
    for i in range(len(typenums)):
        atype = np.hstack((atype, np.full(typenums[i], i+1, dtype='int64')))
    
    # Check which coordinate style to use
    try:
        style = style.decode('UTF-8')
    except:
        pass
    if style[0] in 'cCkK':
        scale = False
    else:
        scale = True
    
    # Read in positions
    natoms = np.sum(typenums)
    pos = np.empty((natoms, 3), dtype='float64')
    count = 0
    for i in range(start_i, len(lines)):
        terms = lines[i].split()
        if len(terms) > 0:
            pos[count,:] = np.array(terms, dtype='float64')
        count += 1
    
    # Add atype, pos to prop dictionary
    prop['atype'] = atype
    prop['pos'] = pos
    
    atoms = Atoms(prop=prop)
    system = System(atoms=atoms, box=box, scale=scale, symbols=symbols)
    
    return system