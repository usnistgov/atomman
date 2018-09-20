# Standard Python libraries
from __future__ import (absolute_import, print_function,
                        division, unicode_literals)
from copy import deepcopy

# atomman imports
from ...lammps import style
from ...compatibility import stringtype, range
from ...tools import indexstr

def process_prop_info(prop_name=None, table_name=None, shape=None, unit=None,
                      dtype=None, prop_info=None, lammps_units='metal'):
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
        for i in range(len(prop_name)):
            
            # Search for prop_name match 
            for sinfo in standard_conversions(lammps_units):
                match = False
                if sinfo['prop_name'] == prop_name[i]:
                    match = True
                    break
            
            pinfo = {}
            pinfo['prop_name'] = prop_name[i]
            
            if table_name is not None:
                pinfo['table_name'] = table_name[i]
            elif match is True and 'table_name' in sinfo:
                pinfo['table_name'] = sinfo['table_name']
            
            if shape is not None:
                pinfo['shape'] = shape[i]
            elif match is True and 'shape' in sinfo:
                pinfo['shape'] = sinfo['shape']
                
            if unit is not None:
                pinfo['unit'] = unit[i]
            elif match is True and 'unit' in sinfo:
                pinfo['unit'] = sinfo['unit']
                
            if dtype is not None:
                pinfo['dtype'] = dtype[i]
            elif match is True and 'dtype' in sinfo:
                pinfo['dtype'] = sinfo['dtype']
            
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

def standard_conversions(lammps_units='metal'):
    """
    Gives pre-defined conversion information for standard LAMMPS variables.
    """
    
    lammps_unit = style.unit(lammps_units)
    
    return [{"prop_name": "a_id",
             "table_name": "id"},
            
            {"prop_name": "m_id",
             "table_name": "mol"},
            
            {"prop_name": "p_id",
             "table_name": "proc"},
            
            {"prop_name": "p_id_plus1",
             "table_name": "procp1"},
            
            {"prop_name": "atype",
             "table_name": "type"},
            
            {"prop_name": "element",
             "table_name": "element"},
            
            {"prop_name": "mass",
             "table_name": "mass",
             "unit": lammps_unit['mass']},
            
            {"prop_name": "pos",
             "table_name": ["x", "y", "z"],
             "unit": lammps_unit['length']},
            
            {"prop_name": "spos",
             "table_name": ["xs", "ys", "zs"],
             "unit": 'scaled'},
            
            {"prop_name": "upos",
             "table_name": ["xu", "yu", "zu"],
             "unit": lammps_unit['length']},
            
            {"prop_name": "supos",
             "table_name": ["xsu", "ysu", "zsu"],
             "unit": 'scaled'},
            
            {"prop_name": "boximage",
             "table_name": ["ix", "iy", "iz"],
             "unit": 'scaled'},
            
            {"prop_name": "velocity",
             "table_name": ["vx", "vy", "vz"],
             "unit": lammps_unit['velocity']},
            
            {"prop_name": "force",
             "table_name": ["fx", "fy", "fz"],
             "unit": lammps_unit['force']},
            
            {"prop_name": "charge",
             "table_name": "q",
             "unit": lammps_unit['charge']},
            
            {"prop_name": "mu",
             "table_name": ["mux", "muy", "muz"],
             "unit": lammps_unit['dipole']},
            
            {"prop_name": "mu_mag",
             "table_name": "mu",
             "unit": lammps_unit['dipole']},
            
            {"prop_name": "radius",
             "table_name": "radius",
             "unit": lammps_unit['length']},
            
            {"prop_name": "diameter",
             "table_name": "diameter",
             "unit": lammps_unit['length']},
            
            {"prop_name": "ang_velocity",
             "table_name": ["omegax", "omegay", "omegaz"],
             "unit": lammps_unit['ang-vel']},
            
            {"prop_name": "ang_momentum",
             "table_name": ["angmomx", "angmomy", "angmomz"],
             "unit": lammps_unit['ang-mom']},
            
            {"prop_name": "torque",
             "table_name": ["tqx", "tqy", "tqz"],
             "unit": lammps_unit['force'] + '*' + lammps_unit['length']}]