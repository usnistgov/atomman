# Standard Python libraries
from __future__ import (absolute_import, print_function,
                        division, unicode_literals)
from copy import deepcopy

# atomman imports
from ...compatibility import stringtype, range
from ...tools import indexstr

def process_prop_info(prop_name=None, table_name=None, shape=None, unit=None, dtype=None, prop_info=None):
    """Handles common setting of prop_info terms"""
    
    if prop_info is not None:
        # Check that competing parameters are not given
        try:
            assert prop_name is None
            assert table_name is None
            assert shape is None
            assert unit is None
            assert dtype is None
        except:
            raise ValueError('prop_info cannot be given with prop_name, table_name, shape, unit, or dtype')
        prop_info = deepcopy(prop_info)
    
    # Combine list parameters into prop_info
    else:
        
        # Check that prop_name is given
        if prop_name is None:
            raise ValueError('prop_info or prop_name is required.')
        numprops = len(prop_name)
        
        # Check that sizes match
        try:
            assert table_name is None or len(table_name) == numprops
            assert shape is None or len(shape) == numprops
            assert unit is None or len(unit) == numprops
            assert dtype is None or len(dtype) == numprops
        except:
            raise ValueError('any of prop_name, table_name, shape, unit, and dtype given must be the same length')
        
        # Build prop_info
        prop_info = []
        for i in range(numprops):
            pinfo = {}
            pinfo['prop_name'] = prop_name[i]
            if table_name is not None:
                pinfo['table_name'] = table_name[i]
            if shape is not None:
                pinfo['shape'] = shape[i]
            if unit is not None:
                pinfo['unit'] = unit[i]
            if dtype is not None:
                pinfo['dtype'] = dtype[i]
            prop_info.append(pinfo)
    
    # Fill in missing default values
    for prop in prop_info:
        
        # prop_name is required
        if 'prop_name' not in prop:
            raise KeyError('prop_name required for each entry in prop_info')
        
        # Generate table_name from prop_name and shape
        if 'table_name' not in prop:
            prop['shape'] = prop.get('shape', ())
            prop['table_name'] = []
            for index, istr in indexstr(prop['shape']):
                prop['table_name'].append(prop['prop_name'] + istr)
        # Make certain table_name is a list
        else:
            try:
                assert not isinstance(prop['table_name'], stringtype)
                prop['table_name'] = list(prop['table_name'])
            except:
                prop['table_name'] = [prop['table_name']]
        
        # Take shape as length of table_name
        if 'shape' not in prop:
            numtnames = len(prop['table_name'])
            if numtnames == 1:
                prop['shape'] = ()
            else:
                prop['shape'] = (numtnames, )
        # Make certain shape is a tuple
        else:
            try:
                prop['shape'] = tuple(prop['shape'])
            except:
                prop['shape'] = (prop['shape'],)
        
        # Default unit and dtype values are None
        prop['unit'] = prop.get('unit', None)
        prop['dtype'] = prop.get('dtype', None)
    
    return prop_info