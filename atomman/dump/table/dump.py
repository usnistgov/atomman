# Standard Python libraries
from __future__ import (absolute_import, print_function,
                        division, unicode_literals)
from collections import OrderedDict

# atomman imports
import atomman.unitconvert as uc
from .process_prop_info import process_prop_info
from ...compatibility import range, ispython2
from ...tools import indexstr

def dump(system, f=None, prop_name=None, table_name=None, shape=None,
         unit=None, dtype=None, prop_info=None, header=False,
         float_format ='%.13f', return_prop_info=False):
    """
    Converts a system's atoms' values to a string table.
    
    Parameters
    ----------
    system : atomman.System
        An atomman representation of a system.
    f : str or file-like object, optional
        File path or file-like object to write the content to.  If not given,
        then the content is returned as a str.
    prop_name : list, optional
        The Atoms properties to include.  Must be given if prop_info is not.
    table_name : list, optional
        The table column name(s) that correspond to each prop_name.  If not
        given, the table_name values will be based on the prop_name values.
    shape : list, optional
        The shape of each per-atom property.  If not given, will be inferred
        from the length of each table_name value.
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
    header : bool, optional
        Flag indicating whether to include the column names in the outputted
        table.  Default value is False (no column names).
    float_format : str, optional
        c-style formatting string for floating point values.  Default value is
        '%.13f'.
    return_prop_info : bool, optional
        Flag indicating if the filled-in prop_info is to be returned.  Having
        this allows for 1:1 load/dump conversions.  Default value is False
        (prop_info is not returned).
        
    Returns
    -------
    str
        The generated data table.  Only returned if fp is None.
    """
    # Set parameters
    natoms = system.natoms
    key_rename = OrderedDict()
    
    # Set default values
    if prop_info is None:
        if prop_name is None:
            prop_name = system.atoms_prop()
        
        if shape is None and table_name is None:
            shape = []
            for name in prop_name:
                shape.append(system.atoms.view[name].shape[1:])
    
    # Process conversion parameters
    prop_info = process_prop_info(prop_name=prop_name, table_name=table_name,
                                  shape=shape, unit=unit, dtype=dtype,
                                  prop_info=prop_info)
    
    
    # Build list of properties to scale
    scale = []
    for prop in prop_info:
        if prop['unit'] == 'scaled':
            scale.append(prop['prop_name'])
            prop['unit'] = None
    
    # Transform to dataframe
    df = system.atoms_df(scale)
    
    # Add a_id values
    df['a_id'] = range(1, natoms+1)
    
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
    table = df.to_csv(path_or_buf=f,
                      sep=sep,
                      index=None,
                      header=header,
                      float_format=float_format,
                      encoding='ascii'
                      )
    
    returns = []
    
    if table is not None:
        returns.append(table)
    
    if return_prop_info is True:
        returns.append(prop_info)
        
    if len(returns) == 1:
        return returns[0]
    elif len(returns) > 1:
        return tuple(returns)
    
    return