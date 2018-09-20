# Standard Python libraries
from __future__ import (absolute_import, print_function,
                        division, unicode_literals)
import os

# https://github.com/usnistgov/DataModelDict
from DataModelDict import DataModelDict as DM

# atomman imports
import atomman.unitconvert as uc
from ...compatibility import range, iteritems, stringtype
from ...tools import identifyfamily

def dump(system, **kwargs):
    """
    Return a DataModelDict 'cell' representation of the system.
    
    Parameters
    ----------
    system : atomman.System
        The system to generate the data model for.
    f : str or file-like object, optional
        File path or file-like object to write the content to.  If not given,
        then the content is returned as a DataModelDict.
    format : str, optional
        File format 'xml' or 'json' to save if f is given.  If not given, will
        be inferred from f if f is a filename, or taken as 'json'.
    indent : int, optional
        Indentation option to use for XML/JSON content if f is given.
    box_unit : str, optional
        Length unit to use for the box. Default value is 'angstrom'.
    symbols : list, optional
        list of atom-model symbols corresponding to the atom types.  If not
        given, will use system.symbols.
    elements : list, optional
        list of element tags corresponding to the atom types.
    prop_units : dict, optional
        dictionary where the keys are the property keys to include, and
        the values are units to use. If not given, only the positions in
        scaled units are included.
    a_std : float, optional
        Standard deviation of a lattice constant to include if available.
    b_std : float, optional
        Standard deviation of b lattice constant to include if available.
    c_std : float, optional
        Standard deviation of c lattice constant to include if available.
    
    Returns
    -------
    DataModelDict
        A 'cell' data model of the system.
    """
    
    # Set default values
    box_unit = kwargs.get('box_unit', 'angstrom')
    indent = kwargs.get('indent', None)
    symbols = kwargs.get('symbols', system.symbols)
    if isinstance(symbols, stringtype):
        symbols = [symbols]
    assert len(symbols) == system.natypes, 'Number of symbols does not match number of atom types'
    
    elements = kwargs.get('elements', [None for i in range(system.natypes)])
    if not isinstance(elements, list):
        elements = [elements]
    assert len(elements) == system.natypes, 'Number of elements does not match number of atom types'
    
    prop_units = kwargs.get('prop_units', {})
    if 'pos' not in prop_units:
        prop_units['pos'] = 'scaled'
    
    # Extract system values
    a = system.box.a
    b = system.box.b
    c = system.box.c
    alpha = system.box.alpha
    beta =  system.box.beta
    gamma = system.box.gamma
    
    # Check for box standard deviations
    if 'a_std' in kwargs and 'b_std' in kwargs and 'c_std' in kwargs:
        errors = True
        a_std = kwargs['a_std']
        b_std = kwargs['b_std']
        c_std = kwargs['c_std']
    else:
        errors = False
        a_std = None
        b_std = None
        c_std = None
    
    # Initialize DataModelDict
    model = DM()
    model['cell'] = cell = DM()
    
    # Test crystal family
    c_family = identifyfamily(system.box)
    if c_family is None:
        c_family = 'triclinic'
    cell[c_family] = DM()
    
    if c_family == 'cubic':
        a_ave = (a + b + c) / 3
        if errors is True:
            a_std_ave = (a_std + b_std + c_std) / 3
        else:
            a_std_ave = None
        
        cell[c_family]['a'] = uc.model(a_ave, box_unit, error=a_std_ave)
    
    elif c_family == 'tetragonal':
        a_ave = (a + b) / 2
        if errors is True:
            a_std_ave = (a_std + b_std) / 2
        else:
            a_std_ave = None
        
        cell[c_family]['a'] = uc.model(a_ave, box_unit, error=a_std_ave)
        cell[c_family]['c'] = uc.model(c, box_unit, error=c_std)
    
    elif c_family == 'orthorhombic':
        cell[c_family]['a'] = uc.model(a, box_unit, error=a_std)
        cell[c_family]['b'] = uc.model(b, box_unit, error=b_std)
        cell[c_family]['c'] = uc.model(c, box_unit, error=c_std)
    
    elif c_family == 'hexagonal':
        a_ave = (a + b) / 2
        if errors is True:
            a_std_ave = (a_std + b_std) / 2
        else:
            a_std_ave = None
        cell[c_family]['a'] = uc.model(a_ave, box_unit, error=a_std_ave)
        cell[c_family]['c'] = uc.model(c, box_unit, error=c_std)
    
    elif c_family == 'rhombohedral':
        a_ave = (a + b + c) / 3
        alpha_ave = (alpha + beta + gamma) / 3
        if errors is True:
            a_std_ave = (a_std + b_std + c_std) / 3
        else:
            a_std_ave = None
        
        cell[c_family]['a'] = uc.model(a_ave, box_unit, error=a_std_ave)
        cell[c_family]['alpha'] = alpha_ave
    
    elif c_family == 'monoclinic':
        cell[c_family]['a'] = uc.model(a, box_unit, error=a_std)
        cell[c_family]['b'] = uc.model(b, box_unit, error=b_std)
        cell[c_family]['c'] = uc.model(c, box_unit, error=c_std)
        cell[c_family]['beta'] = beta
    
    elif c_family == 'triclinic':
        cell[c_family]['a'] = uc.model(a, box_unit, error=a_std)
        cell[c_family]['b'] = uc.model(b, box_unit, error=b_std)
        cell[c_family]['c'] = uc.model(c, box_unit, error=c_std)
        cell[c_family]['alpha'] = alpha
        cell[c_family]['beta'] = beta
        cell[c_family]['gamma'] = gamma
    else:
        raise ValueError('Unknown crystal family')
    
    atype = system.atoms.atype
    aindex = atype - 1
    
    # Build list of atoms and per-atom properties
    for i in range(system.natoms):
        atom = DM()
        
        atom['component'] = int(atype[i])
        
        symbol = symbols[aindex[i]]
        if symbol is not None:
            atom['symbol'] = symbol
        
        element = elements[aindex[i]]
        if element is not None:
            atom['element'] = element
        
        atom['position'] = DM()
        if prop_units['pos'] == 'scaled':
            atom['position']['value'] = list(system.atoms_prop(a_id=i, key='pos', scale=True))
        else:
            atom['position']['value'] = list(uc.get_in_units(system.atoms.pos[i], prop_units['pos']))
        atom['position']['unit'] = prop_units['pos']
        
        for key, unit in iteritems(prop_units):
            if key != 'pos' and key != 'atype':
                value = system.atoms.view[key][i]
                prop = DM()
                prop['name'] = key 
                prop.update(uc.model(value, unit))
                atom.append('property', prop)
        
        model.append('atom', atom)
    
    model = DM([('atomic-system', model)])
    
    # Return DataModelDict or str
    if 'f' not in kwargs:
        if 'format' not in kwargs:
            return model
        elif kwargs['format'].lower() == 'xml':
            return model.xml(indent=indent)
        elif kwargs['format'].lower() == 'json':
            return model.json(indent=indent)
    
    # Write to file
    else:
        f = kwargs['f']
        if 'format' not in kwargs:
            try:
                format = os.path.splitext(f)[1][1:]
            except:
                format = 'json'
        else:
            format = kwargs['format']
        
        if hasattr(f, 'write'):
            if format.lower() == 'xml':
                return model.xml(fp=f, indent=indent)
            elif format.lower() == 'json':
                return model.json(fp=f, indent=indent)
        
        else:
            with open(f, 'w') as fp:
                if format.lower() == 'xml':
                    return model.xml(fp=fp, indent=indent)
                elif format.lower() == 'json':
                    return model.json(fp=fp, indent=indent)
        
    
    return 