# Standard Python libraries
from __future__ import (absolute_import, print_function,
                        division)

from .ase_Atoms import load as load_ase_Atoms
from .cif import load as load_cif
from .pymatgen_Structure import load as load_pymatgen_Structure
from .table import load as load_table
from .system_model import load as load_system_model
from .poscar import load as load_poscar
from .atom_data import load as load_atom_data
from .atom_dump import load as load_atom_dump
from .spglib_cell import load as load_spglib_cell
from .phonopy_Atoms import load as load_phonopy_Atoms

__all__ = ['load', 'load_ase_Atoms', 'load_pymatgen_Structure', 'load_table',
           'load_system_model', 'load_poscar', 'load_atom_data',
           'load_atom_dump', 'load_cif', 'load_spglib_cell', 'load_phonopy_Atoms']

def load(style, input, **kwargs):
    """
    Load a System.
    
    Parameters
    ----------
    style : str
        Indicates the format of the content to load as an atomman.System
    input : str, file-like object or object
        The content to load.
    kwargs
        Any extra keyword arguments to pass to the underlying load methods.
        
    Returns
    -------
    system : atomman.System
        The system object associated with the data model.
    """
    
    if style == 'system_model':
        return load_system_model(input, **kwargs)
    
    elif style == 'cif':
        return load_cif(input, **kwargs)
    
    elif style == 'atom_data':
        return load_atom_data(input, **kwargs)
    
    elif style == 'atom_dump':
        return load_atom_dump(input, **kwargs)
    
    elif style == 'table':
        return load_table(input, **kwargs)
    
    elif style == 'ase_Atoms':
        return load_ase_Atoms(input, **kwargs)
    
    elif style == 'phonopy_Atoms':
        return load_phonopy_Atoms(input, **kwargs)
    
    elif style == 'pymatgen_Structure':
        return load_pymatgen_Structure(input, **kwargs)
    
    elif style == 'poscar':
        return load_poscar(input, **kwargs)
    
    elif style == 'spglib_cell':
        return load_spglib_cell(input, **kwargs)
    
    else:
        raise ValueError('Unsupported style')