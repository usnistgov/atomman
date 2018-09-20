# Standard Python libraries
from __future__ import (absolute_import, print_function,
                        division, unicode_literals)
from io import open
from collections import OrderedDict

# atomman imports
import atomman.unitconvert as uc
from ...lammps import style
from .process_prop_info import process_prop_info
from ...compatibility import range, ispython2
from ...tools import indexstr

def dump(system, f=None, lammps_units='metal', scale=False, prop_name=None,
         table_name=None, shape=None, unit=None, dtype=None,
         prop_info=None, float_format ='%.13f', return_prop_info=False):
    """
    Write a LAMMPS-style atom data file from a System.
    
    Parameters
    ----------
    system : atomman.System 
        The system to write to the atom data file.
    f : str or file-like object, optional
        File path or file-like object to write the content to.  If not given,
        then the content is returned as a str.
    lammps_units : str, optional
        The LAMMPS units option associated with the table values.  This is used
        for the box dimensions and default units for standard dump properties 
        (not compute/fix definitions). 
    scale : bool, optional
        Flag indicating if atom positions are to be scaled relative to the box
        (True) or given in absolute Cartesian values (False, default).
    prop_name : list, optional
        The Atoms properties to include.  If neither prop_name or prop_info are
        given, all system properties will be included.
    table_name : list, optional
        The dump table column name(s) that correspond to each prop_name.  If not
        given, the table_name values will be based on the prop_name and shape
        values.
    shape : list, optional
        The shape of each per-atom property.  If not given, will be inferred
        from the length of each table_name value.
    unit : list, optional
        Lists the units for each prop_name as stored in the table.  For a
        value of None, no conversion will be performed for that property.  For
        a value of 'scaled', the corresponding table values will be taken in
        box-scaled units.  If not given, unit values will be taken based on
        lammps_units if prop_name corresponds to a standard LAMMPS property, 
        otherwise will be set to None (no conversion).
    dtype : list, optional
        Allows for the data type of each property to be explicitly given.
        Values of None will infer the data type from the corresponding
        property values.  If not given, all values will be None.
    prop_info : list of dict, optional
        Structured form of property conversion parameters, in which each
        dictionary in the list corresponds to a single atoms property.  Each
        dictionary must have a 'prop_name' field, and can optionally have
        'table_name', 'shape', 'unit', and 'dtype' fields.
    float_format : str, optional
        c-style formatting string for floating point values.  Default value is
        '%.13f'.
    return_prop_info : bool, optional
        Flag indicating if the filled-in prop_info is to be returned.  Having
        this allows for 1:1 load/dump conversions.  Default value is False
        (prop_info is not returned).
        
    Returns
    -------
    prop_info : list of dict
        The filled-in prop_info structure (if return_prop_info is True).
    """
    lammps_unit = style.unit(lammps_units)
    
    # Set default values
    if prop_info is None:
        if prop_name is None:
            prop_name = ['a_id'] + system.atoms_prop()
        
        if shape is None and table_name is None:
            shape = []
            for name in prop_name:
                if name == 'a_id':
                    shape.append(())
                elif name in ['spos', 'upos', 'supos']:
                    shape.append((3,))
                else:
                    shape.append(system.atoms.view[name].shape[1:])
    
    # Process conversion parameters
    prop_info = process_prop_info(prop_name=prop_name, table_name=table_name,
                                  shape=shape, unit=unit, dtype=dtype,
                                  prop_info=prop_info,
                                  lammps_units=lammps_units)
    
    # Write timestep info
    content = 'ITEM: TIMESTEP\n'
    try:
        content += '%i\n' % system.timestep
    except:
        content += '0\n'
    
    # Write number of atoms
    content += 'ITEM: NUMBER OF ATOMS\n'
    content += '%i\n' % (system.natoms)
    
    # Extract and convert box values
    xlo = uc.get_in_units(system.box.xlo, lammps_unit['length'])
    xhi = uc.get_in_units(system.box.xhi, lammps_unit['length'])
    ylo = uc.get_in_units(system.box.ylo, lammps_unit['length'])
    yhi = uc.get_in_units(system.box.yhi, lammps_unit['length'])
    zlo = uc.get_in_units(system.box.zlo, lammps_unit['length'])
    zhi = uc.get_in_units(system.box.zhi, lammps_unit['length'])
    xy = uc.get_in_units(system.box.xy, lammps_unit['length'])
    xz = uc.get_in_units(system.box.xz, lammps_unit['length'])
    yz = uc.get_in_units(system.box.yz, lammps_unit['length'])
    
    # Compute absolute box bounds
    xlo_bound = xlo + min((0.0, xy, xz, xy + xz))
    xhi_bound = xhi + max((0.0, xy, xz, xy + xz))
    ylo_bound = ylo + min((0.0, yz))
    yhi_bound = yhi + max((0.0, yz))
    zlo_bound = zlo
    zhi_bound = zhi
    
    is_orthogonal = (xy == 0.0 and xz == 0.0 and yz == 0.0)
    
    # Write system boundary info
    if is_orthogonal:
        content += 'ITEM: BOX BOUNDS'
    else:
        content += 'ITEM: BOX BOUNDS xy xz yz'
        
    # Write pbc info
    for i in range(3):
        if system.pbc[i]:
            content += ' pp'
        else:
            content += ' fm'
    content += '\n'
        
    # Write system boundary info
    if is_orthogonal:
        xf2 = float_format + ' ' + float_format + '\n'
        content += xf2 % (xlo_bound, xhi_bound)
        content += xf2 % (ylo_bound, yhi_bound)
        content += xf2 % (zlo_bound, zhi_bound)
    else:
        xf3 = float_format + ' ' + float_format + ' ' + float_format + '\n'
        content += xf3 % (xlo_bound, xhi_bound, xy)
        content += xf3 % (ylo_bound, yhi_bound, xz)
        content += xf3 % (zlo_bound, zhi_bound, yz)
    
    # Write atomic header info and prepare outarray for writing
    header = 'ITEM: ATOMS'
    for prop in prop_info:
        header += ' ' + ' '.join(prop['table_name'])
    header += '\n'
    content += header
    
    content += table_dump(system, prop_info=prop_info, float_format=float_format)
    
    returns = []
    
    # Save to the file-like object
    if hasattr(f, 'write'):
        f.write(content)
    
    # Save to the file name
    elif f is not None:
        with open(f, 'w') as fp:
            fp.write(content)
    
    # Return as a string
    else:
        returns.append(content)
    
    if return_prop_info is True:
        returns.append(prop_info)
        
    if len(returns) == 1:
        return returns[0]
    elif len(returns) > 1:
        return tuple(returns)
        
def table_dump(system, f=None, prop_info=None, float_format ='%.13f'):
    """
    Converts a system's atoms' values to a string table.  Modified from
    table.dump to handle alternate pos fields that dump files can use.
    
    Parameters
    ----------
    system : atomman.System
        An atomman representation of a system.
    f : str or file-like object, optional
        File path or file-like object to write the content to.  If not given,
        then the content is returned as a str.
    prop_info : list of dict, optional
        Structured form of property conversion parameters, in which each
        dictionary in the list corresponds to a single atoms property.  Each
        dictionary must have a 'prop_name' field, and can optionally have
        'table_name', 'shape', 'unit', and 'dtype' fields.
    float_format : str, optional
        c-style formatting string for floating point values.  Default value is
        '%.13f'.
    
    Returns
    -------
    str
        The generated data table.  Only returned if fp is None.
    """
    # Set parameters
    natoms = system.natoms
    key_rename = OrderedDict()
    
    # Build list of properties to scale and alternate pos
    scale = []
    altpos = []
    for prop in prop_info:
        if prop['prop_name'] in ['spos', 'upos', 'supos']:
            altpos.append(prop['prop_name'])
        if prop['unit'] == 'scaled':
            scale.append(prop['prop_name'])
            prop['unit'] = None
    
    # Transform to dataframe
    df = system.atoms_df(scale)
    
    # Add a_id values
    df['a_id'] = range(1, natoms+1)
    
    # Add alternate pos terms
    if 'upos' in altpos:
        df['upos[0]'] = df['pos[0]']
        df['upos[1]'] = df['pos[1]']
        df['upos[2]'] = df['pos[2]']
    if 'spos' in altpos or 'supos' in altpos:
        spos = system.atoms_prop(key='pos', scale=True)
        if 'spos' in altpos:
            df['spos[0]'] = spos[:,0]
            df['spos[1]'] = spos[:,1]
            df['spos[2]'] = spos[:,2]
        if 'supos' in altpos:
            df['supos[0]'] = spos[:,0]
            df['supos[1]'] = spos[:,1]
            df['supos[2]'] = spos[:,2]
    
    # Loop over all properties
    for prop in prop_info:
        pname = prop['prop_name']
        
        # loop over all table names and property indexes
        for tname, (index, istr) in zip(prop['table_name'],
                                        indexstr(prop['shape'])):
            
            # Build name change dict
            key_rename[pname + istr] = tname
            
            # Convert units if needed
            if prop['unit'] is not None:
                df[pname + istr] = uc.get_in_units(df[pname + istr], prop['unit'])
    
    # Rename and reorganize
    df = df.rename(columns=key_rename)[list(key_rename.values())]
  
    # Generate table
    sep = ' '
    if ispython2:
        sep = sep.encode('utf-8')
    return df.to_csv(path_or_buf=f, sep=sep, index=None, header=False,
                     float_format=float_format)