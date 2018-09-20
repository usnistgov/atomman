# Standard Python libraries
from __future__ import (absolute_import, print_function,
                        division, unicode_literals)

# atomman imports
from ...lammps import style

def velocities_prop_info(atom_style='atomic', units='metal'):
    """
    Returns the prop_info associated with a Velocity section in a LAMMPS data file
    
    Parameters
    ----------
    atom_style : str, optional
        The LAMMPS atom_style option.  Default value is 'atomic'.
    units : str, optional
        The LAMMPS units option.  Default value is 'metal'.
        
    Returns
    -------
    prop_info : list of dict
        The proper conversion metadata as used by atomman.convert.atoms_table
    """
    lammps_unit = style.unit(units)
    
    if atom_style in ('angle', 'atomic', 'body', 'bond', 'charge', 'dipole',
                      'full', 'line', 'meso', 'molecular', 'peri', 'smd',
                      'template', 'tri', 'wavepacket'):
        prop_info = [{"prop_name": "a_id",
                      "table_name": "id"},
                     
                     {"prop_name": "velocity",
                      "table_name": ["vx", "vy", "vz"],
                      "unit": lammps_unit['velocity']}]
    
    elif atom_style == 'electron':
        prop_info = [{"prop_name": "a_id",
                      "table_name": "id"},
                     
                     {"prop_name": "velocity",
                      "table_name": ["vx", "vy", "vz"],
                      "unit": lammps_unit['velocity']},
                     
                     {"prop_name": "eradial_velocity",
                      "table_name": "ervel",
                      "unit": lammps_unit['velocity']}]
    
    elif atom_style == 'ellipsoid':
        prop_info = [{"prop_name": "a_id",
                      "table_name": "id"},
                     
                     {"prop_name": "velocity",
                      "table_name": ["vx", "vy", "vz"],
                      "unit": lammps_unit['velocity']},
                     
                     {"prop_name": "ang_momentum",
                      "table_name": ["lx", "ly", "lz"],
                      "unit": lammps_unit['ang-mom']}]
    
    elif atom_style == 'sphere':
        prop_info = [{"prop_name": "a_id",
                      "table_name": "id"},
                     
                     {"prop_name": "velocity",
                      "table_name": ["vx", "vy", "vz"],
                      "unit": lammps_unit['velocity']},
                     
                     {"prop_name": "ang_velocity",
                      "table_name": ["wx", "wy", "wz"],
                      "unit": lammps_unit['ang-vel']}]
    
    elif atom_style[:6] == 'hybrid':
        substyles = atom_style.split()
        prop_info = velocities_prop_info('atomic')
        for substyle in substyles[1:]:
            prop_names = []
            for prop in prop_info:
                prop_names.append(prop['prop_name'])
            subprop_info = velocities_prop_info(substyle)
            for prop in subprop_info:
                if prop['prop_name'] not in prop_names:
                    prop_info.append(prop)
    else:
        raise ValueError('atom_style ' + atom_style + ' not supported')
    
    return prop_info