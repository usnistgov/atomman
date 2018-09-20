# Standard Python libraries
from __future__ import (absolute_import, print_function,
                        division, unicode_literals)
from collections import OrderedDict

# http://www.numpy.org/
import numpy as np

# https://github.com/usnistgov/DataModelDict
from DataModelDict import DataModelDict as DM

# atomman imports
import atomman.unitconvert as uc
from ... import Atoms, Box, System
from ...compatibility import range, iteritems

def load(model, symbols=None, key='atomic-system', index=0):
    """
    Read in a data model containing a crystal structure.
    
    Parameters
    ----------
    model : str, file-like object or DataModelDict
        The data model to read.
    symbols : tuple, optional
        Allows the list of element symbols to be assigned during loading.
    key : str, optional
        The key identifying the root element for the system definition.
        Default value is 'atomic-system'.
    index : int, optional.
        If the full model has multiple key entries, the index specifies which
        to access.  Default value is 0 (first, or only entry).
    
    Returns
    -------
    system : atomman.System
        The system object associated with the data model.
    """
    
    # Pull system model out of data model using key and index
    a_sys = DM(model).finds(key)
    if len(a_sys) == 0:
        raise KeyError(key + ' not found in model')
    try:
        a_sys = a_sys[index]
    except:
        raise IndexError('Invalid index ' + str(index) + ' for model key ' + key)
    
    # Extract crystal system and box values
    c_system = list(a_sys['cell'].keys())[0]
    cell = a_sys['cell'][c_system]
    
    if c_system == 'cubic': 
        a = b = c = uc.value_unit(cell['a'])
        alpha = beta = gamma = 90.0
    
    elif c_system == 'hexagonal':
        a = b = uc.value_unit(cell['a'])
        c = uc.value_unit(cell['c'])
        alpha = beta = 90.0
        gamma = 120.0
    
    elif c_system == 'tetragonal':
        a = b = uc.value_unit(cell['a'])
        c = uc.value_unit(cell['c'])
        alpha = beta = gamma = 90.0
    
    elif c_system == 'trigonal' or c_system == 'rhombohedral':
        a = b = c = uc.value_unit(cell['a'])
        alpha = beta = gamma = cell['alpha']
    
    elif c_system == 'orthorhombic':
        a = uc.value_unit(cell['a'])
        b = uc.value_unit(cell['b'])
        c = uc.value_unit(cell['c'])
        alpha = beta = gamma = 90.0
    
    elif c_system == 'monoclinic':
        a = uc.value_unit(cell['a'])
        b = uc.value_unit(cell['b'])
        c = uc.value_unit(cell['c'])
        alpha = gamma = 90.0
        beta = cell['beta']
    
    elif c_system == 'triclinic':
        a = uc.value_unit(cell['a'])
        b = uc.value_unit(cell['b'])
        c = uc.value_unit(cell['c'])
        alpha = cell['alpha']
        beta = cell['beta']
        gamma = cell['gamma']
    
    box = Box(a=a, b=b, c=c, alpha=alpha, beta=beta, gamma=gamma)
    
    # Count atypes and generate list of symbols if given
    atoms = []
    scale = None
    
    all_atypes = np.array(a_sys.finds('component'))
    all_symbols = np.array(a_sys.finds('symbol'))
    all_elements = np.array(a_sys.finds('element'))
    
    if len(all_atypes) == 0:
        if len(all_symbols) != 0:
            lsymbols, atypes = np.unique(all_symbols, return_inverse=True)
        elif len(all_elements) != 0:
            lsymbols, atypes = np.unique(all_elements, return_inverse=True)
        else:
            raise ValueError('No atom components, symbols or elements listed')
    
    else:
        atypes = all_atypes
        lsymbols = [None for i in range(max(all_atypes))]
        
        if len(all_elements) != 0 and len(all_symbols) == 0:
            all_symbols = all_elements
        
        if len(all_symbols) != 0:
            assert len(all_symbols) == len(atypes)
            sym_dict = {}
            for atype, symbol in zip(atypes, all_symbols):
                if atype not in sym_dict:
                    sym_dict[atype] = symbol
                else:
                    assert sym_dict[atype] == symbol
            
            for atype, symbol in iteritems(sym_dict):
                lsymbols[atype-1] = symbol
    
    # Use lsymbols if symbols parameter is not given.
    if symbols is None:
        symbols = lsymbols
    
    # Read per-atom properties
    natoms = len(atypes)
    prop = OrderedDict()
    prop['atype'] = atypes
    prop['pos'] = np.zeros((natoms, 3), dtype='float64')
    count = 0
    
    pos_units = []
    for atom in a_sys.iteraslist('atom'):
        
        # Read in pos for atom and unit info
        prop['pos'][count] = uc.value_unit(atom['position'])
        pos_units.append(atom['position'].get('unit', None))
        
        # Add other per-atom properties
        for property in atom.iteraslist('property'):
            if property['name'] not in prop:
                value = uc.value_unit(property)
                shape = (natoms, ) + value.shape
                prop[property['name']] = np.zeros(shape, dtype=value.dtype)
                
            prop[property['name']][count] = uc.value_unit(property)
        count += 1
    
    pos_unit = np.unique(pos_units)
    assert len(pos_unit) == 1, 'Mixed units for positions'
    if pos_unit[0] == 'scaled':
        scale=True
    else:
        scale=False
    
    atoms = Atoms(**prop)
    system = System(box=box, atoms=atoms, scale=scale, symbols=symbols)
    
    return system