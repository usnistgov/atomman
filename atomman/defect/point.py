# Standard Python libraries
from __future__ import (absolute_import, print_function,
                        division, unicode_literals)
from copy import deepcopy

# http://www.numpy.org/
import numpy as np

# atomman imports
from ..core import Atoms, System
import atomman.unitconvert as uc

__all__ = ['point', 'vacancy', 'interstitial', 'substitutional', 'dumbbell']

def point(system, ptd_type='v', pos=None, ptd_id=None,
          db_vect=None, scale=False, atol=None, **kwargs):
    """
    Generates a new System where a point defect has been inserted.
    
    Parameters
    ----------
    system : atomman.System
        The base System to add the defect to.
    ptd_type : str, optional
        Key indicating which type of defect to add. Default value is 'v':
        - 'v' : vacancy.
        - 'i' : positional interstitial.
        - 's' : substitutional.
        - 'db' : dumbbell interstitial.
    pos : array-like object, optional
        Position for adding the defect atom (all styles).
    ptd_id : int, optional
        Atom id where defect is added.  Alternative to using pos
        ('v', 's', 'db' styles).
    db_vect : array-like object, optional
        Vector associated with the dumbbell interstitial ('db' style).
    scale : bool, optional
        Indicates if pos and db_vect are Cartesian (False) or box-relative
        (True). Default value is False.
    atol : float, optional
        Absolute tolerance for position-based searching. Default value is 0.01
        angstroms.
    **kwargs : any, optional
        Keyword arguments corresponding to per-atom property values for the
        new atom ('i', 's', 'db' styles).
    
    Raises
    ------
    AssertionError
        If parameters are given for styles that don't allow them.
    ValueError
        If an invalid ptd_type is given.
    
    Returns
    -------
    atomman.System
        A new system containing the defect.
    """
    
    # If vacancy
    if ptd_type == 'v':
        assert db_vect is None, "db_vect not allowed with ptd_type=='v'"
        assert len(kwargs) == 0, "per-atom properties not allowed with ptd_type=='v'"
        return vacancy(system, pos=pos, ptd_id=ptd_id, scale=scale, atol=atol)
    
    # If interstitial
    elif ptd_type == 'i':
        assert ptd_id is None, "ptd_id not allowed with ptd_type=='i'"
        assert db_vect is None, "db_vect not allowed with ptd_type=='i'"
        return interstitial(system, pos=pos, scale=scale, atol=atol, **kwargs)
    
    # If substitutional
    elif ptd_type == 's':
        assert db_vect is None, "db_vect not allowed with ptd_type=='s'"
        return substitutional(system, pos=pos, ptd_id=ptd_id, scale=scale,
                              atol=atol, **kwargs)
    
    # If dumbbell
    elif ptd_type == 'db':
        return dumbbell(system, pos=pos, ptd_id=ptd_id, db_vect=db_vect,
                        scale=scale, atol=atol, **kwargs)
    
    else:
        raise ValueError('Invalid ptd_type. Options are: v, i, s, or db')

def vacancy(system, pos=None, ptd_id=None, scale=False, atol=None):
    """
    Generates a new system by adding a vacancy point defect.
    1. Removes the indicated atom from the system
    2. Adds per-atom property old_id if it doesn't exist corresponding to the
    atom ids in the original system.
    
    Parameters
    ----------
    system : atomman.System
        The base System to add the defect to.
    pos : array-like object, optional
        Position of the atom to be removed.  Either pos or ptd_id must be
        given.
    ptd_id : int, optional
        Id of the atom to be removed.  Either pos or ptd_id must be given.
    scale : bool, optional
        Indicates if pos is Cartesian (False) or box-relative (True). Default
        value is False.
    atol : float, optional
        Absolute tolerance for position-based searching. Default value is 0.01
        angstroms.
    
    Returns
    -------
    atomman.System
        A new system with the vacancy added.
    """
    # Set default atol
    if atol is None:
        atol = uc.set_in_units(0.01, 'angstrom')
    
    # Identify the id of the atom at pos
    if pos is not None:
        if ptd_id is not None:
            raise ValueError('pos and ptd_id cannot both be supplied')
        
        if scale:
            pos = system.unscale(pos)
        
        dist = np.linalg.norm(system.dvect(pos, system.atoms.pos), axis=1)
        ptd_id = np.where(np.isclose(dist, 0.0, atol=atol))
        if len(ptd_id) == 1 and len(ptd_id[0]) == 1:
            ptd_id = ptd_id[0][0]
        else:
            raise ValueError('Unique atom at pos not identified')
    
    elif ptd_id is not None:
        if ptd_id < 0:
            ptd_id += system.natoms
        if ptd_id < 0 or ptd_id >= system.natoms:
            raise ValueError('invalid ptd_id')
    
    else:
        raise ValueError('Either pos or ptd_id required')
    
    # Generate atomic index list for defect
    index = list(range(system.natoms))
    try:
        index.pop(ptd_id)
    except:
        raise TypeError('ptd_id must be an integer type')
    
    # Build new system
    d_system = System(box=deepcopy(system.box), pbc=deepcopy(system.pbc),
                      atoms=deepcopy(system.atoms[index]),
                      symbols=system.symbols)
    
    # Add property old_id with each atom's original id
    if 'old_id' not in d_system.atoms_prop():
        d_system.atoms.old_id = index
    
    return d_system

def interstitial(system, pos, scale=False, atol=None, **kwargs):
    """
    Generates a new system by adding an interstitial point defect.
    1. Adds a new atom to the end of the Atoms list.
    2. Adds per-atom property old_id if it doesn't exist corresponding to the
    atom ids in the original system.
    3. Sets any of the new atom's per-atom properties to values given as
    kwargs.  Any undefined properties are given zero values except atype,
    which is set to 1.
    
    Parameters
    ----------
    system : atomman.System
        The base System to add the defect to.
    pos : array-like object
        Position of the atom being added.
    scale : bool, optional
        Indicates if pos is Cartesian (False) or box-relative (True).  Default
        value is False.
    atol : float, optional
        Absolute tolerance for position-based searching. Default value is 0.01
        angstroms.
    **kwargs : any, optional
        Keyword arguments corresponding to per-atom property values for the
        new atom.  By default, atype==1 and all other properties are set to
        be all zeros for the property's shape.
    
    Returns
    -------
    atomman.System
        A new system with the interstitial added.
    """
    # Set default atol
    if atol is None:
        atol = uc.set_in_units(0.01, 'angstrom')
    
    if scale:
        pos = system.unscale(pos)
    
    # Check that no atoms are already at pos
    dist = np.linalg.norm(system.dvect(pos, system.atoms.pos), axis=1)
    ptd_id = np.where(np.isclose(dist, 0.0, atol=atol))
    if not (len(ptd_id) == 1 and len(ptd_id[0]) == 0):
        raise ValueError('atom already at pos')
    
    # Generate atomic index list for defect
    index = list(range(system.natoms))
    index.append(0)
    
    # Build new system with atom 0 copied
    d_system = System(box=deepcopy(system.box), pbc=deepcopy(system.pbc),
                      atoms=deepcopy(system.atoms[index]),
                      symbols=system.symbols)
    
    # Add property old_id with each atom's original id
    if 'old_id' not in d_system.atoms_prop():
        d_system.atoms.old_id = index
    
    # Set values for the new atom
    for prop in d_system.atoms_prop():
        if prop == 'atype':
            d_system.atoms.atype[-1] = kwargs.pop('atype', 1)
        elif prop == 'pos':
            d_system.atoms.pos[-1] = pos
        elif prop == 'old_id':
            d_system.atoms.old_id[-1] = kwargs.pop('old_id',
                                                   d_system.atoms.old_id.max()+1)
        else:
            d_system.atoms.view[prop][-1] = kwargs.pop(prop,
                                                       np.zeros_like(d_system.atoms.view[prop][-1]))
    
    return d_system

def substitutional(system, pos=None, ptd_id=None, atype=1,scale=False,
                   atol=None, **kwargs):
    """
    Generates a new system by adding a substitutional point defect.
    1. Moves the indicated atom to the end of the list and changes its atype
    to the value given.
    2. Adds per-atom property old_id if it doesn't exist corresponding to the
    atom ids in the original system.
    3. Sets any of the moved atom's per-atom properties to values given as
    kwargs.  Any undefined properties are left unchanged.
    
    Parameters
    ----------
    system : atomman.System
        The base System to add the defect to.
    pos : array-like object, optional
        Position of the atom being modified.  Either pos or ptd_id must be
        given.
    ptd_id : int, optional
        Id of the atom to be modified.  Either pos or ptd_id must be given.
    atype : int, optional
        Integer atomic type to change the identified atom to.  Must be
        different than the atom's current id.  Default value is 1.
    scale : bool, optional
        Indicates if pos is Cartesian (False) or box-relative (True).  Default
        value is False.
    atol : float, optional
        Absolute tolerance for position-based searching. Default value is 0.01
        angstroms.
    **kwargs : any, optional
        Keyword arguments corresponding to per-atom property values for the
        modified atom.  By default, all properties (except atype) are left
        unchanged.
    
    Returns
    -------
    atomman.System
        A new system with the substitutional added.
    """
    # Set default atol
    if atol is None:
        atol = uc.set_in_units(0.01, 'angstrom')
    
    # Identify the id of the atom at pos
    if pos is not None:
        if ptd_id is not None:
            raise ValueError('pos and ptd_id cannot both be supplied')
        
        if scale:
            pos = system.unscale(pos)
        
        dist = np.linalg.norm(system.dvect(pos, system.atoms.pos), axis=1)
        ptd_id = np.where(np.isclose(dist, 0.0, atol=atol))
        if len(ptd_id) == 1 and len(ptd_id[0]) == 1:
            ptd_id = ptd_id[0][0]
        else:
            raise ValueError('Unique atom at pos not identified')
    
    elif ptd_id is not None:
        if ptd_id < 0:
            ptd_id += system.natoms
        if ptd_id < 0 or ptd_id >= system.natoms:
            raise ValueError('invalid ptd_id')
    
    else:
        raise ValueError('Either pos or ptd_id required')
    
    if system.atoms.atype[ptd_id] == atype:
        raise ValueError('identified atom is already of the specified atype')
    
    # Generate atomic index list for defect
    index = list(range(system.natoms))
    index.pop(ptd_id)
    index.append(ptd_id)
    
    # Build new system
    d_system = System(box=deepcopy(system.box), pbc=deepcopy(system.pbc),
                      atoms=deepcopy(system.atoms[index]),
                      symbols=system.symbols)
    
    # Add property old_id with each atom's original id
    if 'old_id' not in d_system.atoms_prop():
        d_system.atoms.old_id = index
    
    # Set values for the new atom
    for prop in d_system.atoms_prop():
        if prop == 'atype':

            d_system.atoms.atype[-1] = atype
        elif prop == 'pos':
            pass
        else:
            d_system.atoms.view[prop][-1] = kwargs.pop(prop,
                                                       d_system.atoms.view[prop][-1])
    
    return d_system

def dumbbell(system, pos=None, ptd_id=None, db_vect=None, scale=False,
             atol=None, **kwargs):
    """
    Generates a new system by adding a dumbbell interstitial point defect.
    
    1. Copies the indicated atom and moves both the original and copy to the
    end of the Atoms list.
    2. Displaces the dumbbell atoms position's by +-db_vect.
    3. Adds per-atom property old_id if it doesn't exist corresponding to the
    atom ids in the original system.
    4. Sets any of the new atom's per-atom properties to values given as
       kwargs.  Any undefined properties are left unchanged.
    
    Parameters
    ----------
    system : atomman.System
        The base System to add the defect to.
    pos : array-like object, optional
        Position of the atom being modified.  Either pos or ptd_id must be
        given.
    ptd_id : int, optional
        Id of the atom to be modified.  Either pos or ptd_id must be given.
    db_vect : array-like object
        Vector shift to apply to the atoms in the dumbbell.
    scale : bool, optional
        Indicates if pos and db_vect are Cartesian (False) or box-relative
        (True).  Default value is False.
    atol : float, optional
        Absolute tolerance for position-based searching. Default value is 0.01
        angstroms.
    \*\*kwargs : any, optional
        Keyword arguments corresponding to per-atom property values for the
        new atom in the dumbbell.  By default, all properties are left
        unchanged (i.e. same as atom that was copied).
    
    Returns
    -------
    atomman.System
        A new system with the dumbbell added.
    """
    # Set default atol
    if atol is None:
        atol = uc.set_in_units(0.01, 'angstrom')
    
    # Identify the id of the atom at pos
    if pos is not None:
        if ptd_id is not None:
            raise ValueError('pos and ptd_id cannot both be supplied')
        
        if scale:
            pos = system.unscale(pos)
        
        dist = np.linalg.norm(system.dvect(pos, system.atoms.pos), axis=1)
        ptd_id = np.where(np.isclose(dist, 0.0, atol=atol))
        if len(ptd_id) == 1 and len(ptd_id[0]) == 1:
            ptd_id = ptd_id[0][0]
        else:
            raise ValueError('Unique atom at pos not identified')
    
    elif ptd_id is not None:
        if ptd_id < 0:
            ptd_id += system.natoms
        if ptd_id < 0 or ptd_id >= system.natoms:
            raise ValueError('invalid ptd_id')
    
    else:
        raise ValueError('Either pos or ptd_id required')
        
    # Unscale db_vect
    if scale:
        db_vect = system.unscale(db_vect)
    
    # Generate atomic index list for defect
    index = list(range(system.natoms))
    index.pop(ptd_id)
    index.append(ptd_id)
    index.append(ptd_id)
    
    # Build new system
    d_system = System(box=deepcopy(system.box), pbc=deepcopy(system.pbc),
                      atoms=deepcopy(system.atoms[index]),
                      symbols=system.symbols)
    
    # Add property old_id with each atom's original id
    if 'old_id' not in d_system.atoms_prop():
        d_system.atoms.old_id = index
    
    # Set values for the new atom
    for prop in d_system.atoms_prop():
        if prop == 'pos':
            d_system.atoms.pos[-2] -= db_vect
            d_system.atoms.pos[-1] += db_vect
        elif prop == 'old_id':
            d_system.atoms.old_id[-1] = kwargs.pop('old_id',
                                                   d_system.atoms.old_id.max()+1)
        else:
            d_system.atoms.view[prop][-1] = kwargs.pop(prop,
                                                       d_system.atoms.view[prop][-1])
    
    return d_system