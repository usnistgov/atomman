# Standard Python libraries
from __future__ import (absolute_import, print_function,
                        division, unicode_literals)

# atomman imports
from ...lammps import style

def atoms_prop_info(atom_style='atomic', units='metal'):
    """
    Returns the prop_info associated with an Atoms section in a LAMMPS data file
    
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
    
    if atom_style == 'angle':
        prop_info = [{"prop_name": "a_id",
                      "table_name": "id"},
                     
                     {"prop_name": "m_id",
                      "table_name": "mol"},
                     
                     {"prop_name": "atype",
                      "table_name": "type"},
                     
                     {"prop_name": "pos",
                      "table_name": ["x", "y", "z"],
                      "unit": lammps_unit['length']}]
    
    elif atom_style == 'atomic':
        prop_info = [{"prop_name": "a_id",
                      "table_name": "id"},
                     
                     {"prop_name": "atype",
                      "table_name": "type"},
                     
                     {"prop_name": "pos",
                      "table_name": ["x", "y", "z"],
                      "unit": lammps_unit['length']}]
    
    elif atom_style == 'body':
        prop_info = [{"prop_name": "a_id",
                      "table_name": "id"},
                     
                     {"prop_name": "atype",
                      "table_name": "type"},
                     
                     {"prop_name": "bflag",
                      "table_name": "bodyflag"},
                     
                     {"prop_name": "mass",
                      "table_name": "mass",
                      "unit": lammps_unit['mass']},
                     
                     {"prop_name": "pos",
                      "table_name": ["x", "y", "z"],
                      "unit": lammps_unit['length']}]
    
    elif atom_style == 'bond':
        prop_info = [{"prop_name": "a_id",
                      "table_name": "id"},
                     
                     {"prop_name": "m_id",
                      "table_name": "mol"},
                     
                     {"prop_name": "atype",
                      "table_name": "type"},
                     
                     {"prop_name": "pos",
                      "table_name": ["x", "y", "z"],
                      "unit": lammps_unit['length']}]
    
    elif atom_style == 'charge':
        prop_info = [{"prop_name": "a_id",
                      "table_name": "id"},
                     
                     {"prop_name": "atype",
                      "table_name": "type"},
                     
                     {"prop_name": "charge",
                      "table_name": "q",
                      "unit": lammps_unit['charge']},
                     
                     {"prop_name": "pos",
                      "table_name": ["x", "y", "z"], 
                      "unit": lammps_unit['length']}]
    
    elif atom_style == 'dipole':
        prop_info = [{"prop_name": "a_id",
                      "table_name": "id"},
                     
                     {"prop_name": "atype",
                      "table_name": "type"},
                     
                     {"prop_name": "charge",
                      "table_name": "q",
                      "unit": lammps_unit['charge']},
                     
                     {"prop_name": "pos",
                      "table_name": ["x", "y", "z"],
                      "unit": lammps_unit['length']},
                     
                     {"prop_name": "mu",
                      "table_name": ["mux", "muy", "muz"],
                      "unit": lammps_unit['dipole']}]
    
    elif atom_style == 'electron':
        prop_info = [{"prop_name": "a_id",
                      "table_name": "id"},
                     
                     {"prop_name": "atype",
                      "table_name": "type"},
                     
                     {"prop_name": "charge",
                      "table_name": "q",
                      "unit": lammps_unit['charge']},
                     
                     {"prop_name": "espin",
                      "table_name": "spin"},
                     
                     {"prop_name": "eradius",
                      "table_name": "eradius", 
                      "unit": lammps_unit['length']},
                     
                     {"prop_name": "pos",
                      "table_name": ["x", "y", "z"],
                      "unit": lammps_unit['length']}]
    
    elif atom_style == 'ellipsoid':
        prop_info = [{"prop_name": "a_id",
                      "table_name": "id"},
                     
                     {"prop_name": "atype",
                      "table_name": "type"},
                     
                     {"prop_name": "eflag",
                      "table_name": "ellipsoidflag"},
                     
                     {"prop_name": "density",
                      "table_name": "density",
                      "unit": lammps_unit['density']},
                     
                     {"prop_name": "pos",
                      "table_name": ["x", "y", "z"],
                      "unit": lammps_unit['length']}]
    
    elif atom_style == 'full':
        prop_info = [{"prop_name": "a_id",
                      "table_name": "id"},
                     
                     {"prop_name": "m_id",
                      "table_name": "mol"},
                     
                     {"prop_name": "atype",
                      "table_name": "type"},
                     
                     {"prop_name": "charge",
                      "table_name": "q",
                      "unit": lammps_unit['charge']},
                     
                     {"prop_name": "pos",
                      "table_name": ["x", "y", "z"],
                      "unit": lammps_unit['length']}]
    
    elif atom_style == 'line':
        prop_info = [{"prop_name": "a_id",
                      "table_name": "id"},
                     
                     {"prop_name": "m_id",
                      "table_name": "mol"},
                     
                     {"prop_name": "atype",
                      "table_name": "type"},
                     
                     {"prop_name": "lflag",
                      "table_name": "lineflag"},
                     
                     {"prop_name": "density",
                      "table_name": "density",
                      "unit": lammps_unit['density']},
                     
                     {"prop_name": "pos",
                      "table_name": ["x", "y", "z"],
                      "unit": lammps_unit['length']}]
    
    elif atom_style == 'meso':
        prop_info = [{"prop_name": "a_id",
                      "table_name": "id"},
                     
                     {"prop_name": "atype",
                      "table_name": "type"},
                     
                     {"prop_name": "rho",
                      "table_name": "rho"},
                     
                     {"prop_name": "e",
                      "table_name": "e"},
                     
                     {"prop_name": "cv",
                      "table_name": "cv"},
                     
                     {"prop_name": "pos",
                      "table_name": ["x", "y", "z"],
                      "unit": lammps_unit['length']}]
    
    elif atom_style == 'molecular':
        prop_info = [{"prop_name": "a_id",
                      "table_name": "id"},
                     
                     {"prop_name": "m_id",
                      "table_name": "mol"},
                     
                     {"prop_name": "atype",
                      "table_name": "type"},
                     
                     {"prop_name": "pos",
                      "table_name": ["x", "y", "z"],
                      "unit": lammps_unit['length']}]
    
    elif atom_style == 'peri':
        prop_info = [{"prop_name": "a_id",
                      "table_name": "id"},
                     
                     {"prop_name": "atype",
                      "table_name": "type"},
                     
                     {"prop_name": "volume",
                      "table_name": "volume",
                      "unit": lammps_unit['volume']},
                     
                     {"prop_name": "density",
                      "table_name": "density",
                      "unit": lammps_unit['density']},
                     
                     {"prop_name": "pos",
                      "table_name": ["x", "y", "z"],
                      "unit": lammps_unit['length']}]
    
    elif atom_style == 'smd':
        prop_info = [{"prop_name": "a_id",
                      "table_name": "id"},
                     
                     {"prop_name": "atype",
                      "table_name": "type"},
                     
                     {"prop_name": "m_id",
                      "table_name": "mol"},
                     
                     {"prop_name": "volume",
                      "table_name": "volume",
                      "unit": lammps_unit['volume']},
                     
                     {"prop_name": "mass",
                      "table_name": "mass",
                      "unit": lammps_unit['mass']},
                     
                     {"prop_name": "kradius",
                      "table_name": "kernalradius",
                      "unit": lammps_unit['length']},
                     
                     {"prop_name": "cradius",
                      "table_name": "contactradius",
                      "unit": lammps_unit['length']},
                     
                     {"prop_name": "pos",
                      "table_name": ["x", "y", "z"],
                      "unit": lammps_unit['length']}]
    
    elif atom_style == 'sphere':
        prop_info = [{"prop_name": "a_id",
                      "table_name": "id"},
                     
                     {"prop_name": "atype",
                      "table_name": "type"},
                     
                     {"prop_name": "diameter",
                      "table_name": "diameter",
                      "unit": lammps_unit['length']},
                     
                     {"prop_name": "density",
                      "table_name": "density",
                      "unit": lammps_unit['density']},
                     
                     {"prop_name": "pos",
                      "table_name": ["x", "y", "z"],
                      "unit": lammps_unit['length']}]
    
    elif atom_style == 'template':
        prop_info = [{"prop_name": "a_id",
                      "table_name": "id"},
                     
                     {"prop_name": "m_id",
                      "table_name": "mol"},
                     
                     {"prop_name": "m_template",
                      "table_name": "templateindex"},
                     
                     {"prop_name": "a_template",
                      "table_name": "templateatom"},
                     
                     {"prop_name": "atype",
                      "table_name": "type"},
                     
                     {"prop_name": "pos",
                      "table_name": ["x", "y", "z"],
                      "unit": lammps_unit['length']}]
    
    elif atom_style == 'tri':
        prop_info = [{"prop_name": "a_id",
                      "table_name": "id"},
                     
                     {"prop_name": "m_id",
                      "table_name": "mol"},
                     
                     {"prop_name": "atype",
                      "table_name": "type"},
                     
                     {"prop_name": "tflag",
                      "table_name": "triangleflag"},
                     
                     {"prop_name": "density",
                      "table_name": "density",
                      "unit": lammps_unit['density']},
                     
                     {"prop_name": "pos",
                      "table_name": ["x", "y", "z"],
                      "unit": lammps_unit['length']}]
    
    elif atom_style == 'wavepacket':
        prop_info = [{"prop_name": "a_id",
                      "table_name": "id"},
                     
                     {"prop_name": "atype",
                      "table_name": "type"},
                     
                     {"prop_name": "charge",
                      "table_name": "q",
                      "unit": lammps_unit['charge']},
                     
                     {"prop_name": "espin",
                      "table_name": "spin"},
                     
                     {"prop_name": "eradius",
                      "table_name": "eradius", 
                      "unit": lammps_unit['length']},
                     
                     {"prop_name": "e_id",
                      "table_name": "etag"},
                     
                     {"prop_name": "cs_re",
                      "table_name": "cs_re"},
                     
                     {"prop_name": "cs_im",
                      "table_name": "cs_im"},
                     
                     {"prop_name": "pos",
                      "table_name": ["x", "y", "z"],
                      "unit": lammps_unit['length']}]
    
    elif atom_style[:6] == 'hybrid':
        substyles = atom_style.split()
        prop_info = atoms_prop_info('atomic')
        for substyle in substyles[1:]:
            prop_names = []
            for prop in prop_info:
                prop_names.append(prop['prop_name'])
            subprop_info = atoms_prop_info(substyle)
            for prop in subprop_info:
                if prop['prop_name'] not in prop_names:
                    prop_info.append(prop)
    else:
        raise ValueError('atom_style ' + atom_style + ' not supported')
    
    return prop_info