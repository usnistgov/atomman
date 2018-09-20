# Standard Python libraries
from __future__ import (absolute_import, print_function,
                        division)

from .ase_Atoms import dump as dump_ase_Atoms
from .pymatgen_Structure import dump as dump_pymatgen_Structure
from .table import dump as dump_table
from .system_model import dump as dump_system_model
from .poscar import dump as dump_poscar
from .atom_data import dump as dump_atom_data
from .atom_dump import dump as dump_atom_dump
from .spglib_cell import dump as dump_spglib_cell
from .phonopy_Atoms import dump as dump_phonopy_Atoms

__all__ = ['dump', 'dump_ase_Atoms', 'dump_pymatgen_Structure', 'dump_table',
           'dump_system_model', 'dump_poscar', 'dump_atom_data',
           'dump_atom_dump', 'dump_spglib_cell', 'dump_phonopy_Atoms']

def dump(style, system, **kwargs):
    """
    Convert a System to another format.
    
    Parameters
    ----------
    style : str
        Indicates the format of the content to dump the atomman.System as.
    system : atomman.System
        The system to convert.
    kwargs : any, optional
        Any extra keyword arguments to pass to the underlying dump methods.
        
    Returns
    -------
    str, object or tuple
        Any content returned by the underlying dump methods.
    """
    
    if style == 'system_model':
        return dump_system_model(system, **kwargs)
    
    elif style == 'atom_data':
        return dump_atom_data(system, **kwargs)
    
    elif style == 'atom_dump':
        return dump_atom_dump(system, **kwargs)
    
    elif style == 'table':
        return dump_table(system, **kwargs)
    
    elif style == 'ase_Atoms':
        return dump_ase_Atoms(system, **kwargs)

    elif style == 'phonopy_Atoms':
        return dump_phonopy_Atoms(system, **kwargs)
    
    elif style == 'pymatgen_Structure':
        return dump_pymatgen_Structure(system, **kwargs)
    
    elif style == 'poscar':
        return dump_poscar(system, **kwargs)
    
    elif style == 'spglib_cell':
        return dump_spglib_cell(system, **kwargs)
        
    else:
        raise ValueError('Unsupported style')