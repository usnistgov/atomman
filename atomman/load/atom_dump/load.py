# Standard Python libraries
from __future__ import (absolute_import, print_function,
                        division, unicode_literals)
from collections import OrderedDict
                        
# atomman imports
import atomman.unitconvert as uc
from ... import Atoms, Box, System
from .process_prop_info import process_prop_info, standard_conversions
from ...lammps import style
from .. import load_table
from ...tools import uber_open_rmode

def load(data, symbols=None, lammps_units='metal', prop_name=None,
         table_name=None, shape=None, unit=None, dtype=None,
         prop_info=None, return_prop_info=False):
    """
    Reads in a LAMMPS atomic dump file into a System.
    
    Parameters
    ----------
    data : str or file-like object
        The content, file path or file-like object containing the content to
        read.
    symbols : tuple, optional
        Allows the list of element symbols to be assigned during loading.
    lammps_units : str
        The LAMMPS units option associated with the parameters.  Default value
        is 'metal'.
    prop_name : list, optional
         The Atoms properties to generate.
    table_name : list, optional
        The table column name(s) that correspond to each prop_name.  If
        prop_name, table_name and prop_info are not given, prop_name and
        table_name will be read in from data.
    shape : list, optional
        The shape of each per-atom property.  If not given, will be taken from
        standard LAMMPS parameter names, or left at () for direct 
        property-table conversion.
    unit : list, optional
        Lists the units for each prop_name as stored in the table.  For a
        value of None, no conversion will be performed for that property.  For
        a value of 'scaled', the corresponding table values will be taken in
        box-scaled units.  If not given, all unit values will be set to None
        (i.e. no conversions).
    dtype : list, optional
        Allows for the data type of each property to be explicitly given.
        Values of None will infer the data type from the corresponding
        property values.  If not given, all values will be None.
    prop_info : list of dict, optional
        Structured form of property conversion parameters, in which each
        dictionary in the list corresponds to a single atoms property.  Each
        dictionary must have a 'prop_name' field, and can optionally have
        'table_name', 'shape', 'unit', and 'dtype' fields.
    return_prop_info : bool, optional
        Flag indicating if the full prop_info is to be returned.  Default value
        is False.
        
    Returns
    -------
    system : atomman.System
        The generated system.
    prop_info : list of dict
        The full prop_info detailing the property-table conversion. Returned
        if return_prop_info is True.
    """
    
    lammps_unit = style.unit(lammps_units)
    
    # Initialize parameter values
    pbc = None
    box = None
    natoms = None
    atomsstart = None
    xy = 0.0
    xz = 0.0
    yz = 0.0
    
    readnatoms = False
    readtimestep = False
    bcount = 3
    
    # Read str and files in the same way
    with uber_open_rmode(data) as fp:
        
        # Loop over all lines in fp
        for i, line in enumerate(fp):
            terms = line.decode('UTF-8').split()
            
            # Skip blank lines
            if len(terms) > 0:
                
                # Read number of atoms if time to do so
                if readnatoms:
                    natoms = int(terms[0])
                    readnatoms = False
                
                # Read timestep if time to do so
                elif readtimestep:
                    #timestep = int(terms[0])
                    readtimestep = False
                
                # Read x boundary condition values if time to do so
                elif bcount == 0:
                    xlo = uc.set_in_units(float(terms[0]), lammps_unit['length'])
                    xhi = uc.set_in_units(float(terms[1]), lammps_unit['length'])
                    if len(terms) == 3:
                        xy = uc.set_in_units(float(terms[2]),
                                             lammps_unit['length'])
                    bcount += 1
                
                # Read y boundary condition values if time to do so
                elif bcount == 1:
                    ylo = uc.set_in_units(float(terms[0]), lammps_unit['length'])
                    yhi = uc.set_in_units(float(terms[1]), lammps_unit['length'])
                    if len(terms) == 3:
                        xz = uc.set_in_units(float(terms[2]),
                                             lammps_unit['length'])
                    bcount += 1
                
                # Read z boundary condition values if time to do so
                elif bcount == 2:
                    zlo = uc.set_in_units(float(terms[0]), lammps_unit['length'])
                    zhi = uc.set_in_units(float(terms[1]), lammps_unit['length'])
                    if len(terms) == 3:
                        yz = uc.set_in_units(float(terms[2]),
                                             lammps_unit['length'])
                        
                        # Convert from max, min to hi, lo
                        xlo = xlo - min((0.0, xy, xz, xy + xz))
                        xhi = xhi - max((0.0, xy, xz, xy + xz))
                        ylo = ylo - min((0.0, yz))
                        yhi = yhi - max((0.0, yz))
                    bcount += 1
                
                # Otherwise, only check lines starting with ITEM
                elif terms[0] == 'ITEM:':
                    
                    # ITEM: TIMESTEP indicates it is time to read the timestep
                    if terms[1] == 'TIMESTEP':
                        readtimestep = True
                    
                    # ITEM: NUMBER indicates it is time to read natoms
                    elif terms[1] == 'NUMBER':
                        readnatoms = True
                    
                    # ITEM: BOX gives pbc and indicates it is time to read box parameters
                    elif terms[1] == 'BOX':
                        pbc = [True, True, True]
                        for i in range(3):
                            if terms[i + len(terms) - 3] != 'pp':
                                pbc[i] = False
                        bcount = 0
                        
                    # ITEM: ATOMS gives list of property names and indicates it is time to read atomic values
                    elif terms[1] == 'ATOMS':
                        
                        # Read list of property names
                        name_list = terms[2:]
                        
                        # Create default prop_name and table_name if needed
                        if prop_info is None and prop_name is None:
                            assert table_name is None, 'table_name cannot be given without prop_name'
                            prop_name, table_name = matchprops(name_list)
                            
                        atomsstart = i + 1
    
    # Create system
    box = Box(xlo=xlo, xhi=xhi,
              ylo=ylo, yhi=yhi,
              zlo=zlo, zhi=zhi,
              xy=xy, xz=xz, yz=yz)
    atoms = Atoms(natoms=natoms)
    system = System(box=box, atoms=atoms, pbc=pbc)
    
    # Generate prop_info
    prop_info = process_prop_info(prop_name=prop_name,
                                  table_name=table_name,
                                  shape=shape, unit=unit, dtype=dtype,
                                  prop_info=prop_info,
                                  lammps_units=lammps_units)
    
    # Remove duplicate pos fields
    firstpos = True
    short_prop_info = []
    for pinfo in prop_info:
        if pinfo['prop_name'] in ['pos', 'spos', 'upos', 'supos']:
            if firstpos:
                pinfo['prop_name'] = 'pos'
            else:
                continue
        short_prop_info.append(pinfo)
    
    # Read atoms into system
    system = load_table(data, box=system.box, symbols=symbols, system=system,
                        prop_info=short_prop_info, skiprows=atomsstart,
                        nrows=natoms)
    
    return system
    
def matchprops(items):
    """
    Takes a list of table_names, pairs them up and matches them to prop_names.
    
    Parameters
    ----------
    items : list
        One dimensional list of all table names.
    
    Returns
    -------
    prop_name : list
        The list of system property names corresponding to items.
    table_name : list
        The list of items paired up and corresponding to the prop_name list.
    """
    
    prop2table = OrderedDict()
    for item in items:
        # Search for prop_name match 
        for sinfo in standard_conversions():
            match = False
            table_names = sinfo['table_name']
            if not isinstance(table_names, list):
                table_names = [table_names]
            if item in table_names:
                match = True
                break
        
        if match is True:
            name = sinfo['prop_name']
        else:
            name = item
        
        if name not in prop2table:
            if match is True:
                for table_name in table_names:
                    assert table_name in items, 'Incomplete propery ' + str(name)
            prop2table[name] = []
        prop2table[name].append(item)
    
    prop_name = list(prop2table.keys())
    table_name = list(prop2table.values())
    for i in range(len(table_name)):
        if len(table_name[i]) == 1:
            table_name[i] = table_name[i][0]
    
    return prop_name, table_name