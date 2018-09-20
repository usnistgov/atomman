# Standard Python libraries
from __future__ import (absolute_import, print_function,
                        division, unicode_literals)

# atomman imports
import atomman.unitconvert as uc
from .atoms_prop_info import atoms_prop_info
from .velocities_prop_info import velocities_prop_info
from ... import Atoms, Box, System
from ...lammps import style
from .. import load_table
from ...tools import uber_open_rmode

def load(data, pbc=(True, True, True), symbols=None, atom_style='atomic', units='metal'):
    """
    Read a LAMMPS-style atom data file.
    
    Parameters
    ----------
    data : str or file-like object
        The atom data content to read.  Can be str content, path name, or open
        file-like object.
    pbc : list of bool
        Three boolean values indicating which System directions are periodic.
        Default value is (True, True, True).
    symbols : tuple, optional
        Allows the list of element symbols to be assigned during loading.
    atom_style :str
        The LAMMPS atom_style option associated with the data file.  Default
        value is 'atomic'.
    units : str 
        The LAMMPS units option associated with the data file. Default value
        is 'metal'.
    
    Returns
    -------
    atomman.System
        The corresponding system.  Note all property values will be
        automatically converted to atomman.unitconvert's set working units.
    """
    
    # Get units information
    units_dict = style.unit(units)
    
    # Initialize parameter values
    atomsstart = None
    velocitiesstart = None
    xy = 0.0
    xz = 0.0
    yz = 0.0
    
    # Read str and files in the same way
    with uber_open_rmode(data) as fp:
        
        # Loop over all lines in fp
        for i, line in enumerate(fp):
            terms = line.decode('UTF-8').split()

            # Skip blank lines
            if len(terms)>0:
                
                # Read number of atoms 
                if len(terms) == 2 and terms[1] == 'atoms':
                    natoms = int(terms[0])
                
                # Read number of atom types
                elif len(terms) == 3 and terms[1] == 'atom' and terms[2] == 'types': 
                    #natypes = int(terms[0])
                    pass
                
                # Read boundary info
                elif len(terms) == 4 and terms[2] == 'xlo' and terms[3] == 'xhi':
                    xlo = uc.set_in_units(float(terms[0]), units_dict['length'])
                    xhi = uc.set_in_units(float(terms[1]), units_dict['length'])
                elif len(terms) == 4 and terms[2] == 'ylo' and terms[3] == 'yhi':
                    ylo = uc.set_in_units(float(terms[0]), units_dict['length'])
                    yhi = uc.set_in_units(float(terms[1]), units_dict['length'])
                elif len(terms) == 4 and terms[2] == 'zlo' and terms[3] == 'zhi':
                    zlo = uc.set_in_units(float(terms[0]), units_dict['length'])
                    zhi = uc.set_in_units(float(terms[1]), units_dict['length'])
                elif len(terms) == 6 and terms[3] == 'xy' and terms[4] == 'xz' and terms[5] == 'yz':
                    xy = uc.set_in_units(float(terms[0]), units_dict['length'])
                    xz = uc.set_in_units(float(terms[1]), units_dict['length'])
                    yz = uc.set_in_units(float(terms[2]), units_dict['length'])
                
                # Identify starting line number for Atoms data
                elif len(terms) == 1 and terms[0] == 'Atoms':
                    atomsstart = i + 1
                
                # Identify starting line number for Velocity data
                elif len(terms) == 1 and terms[0] == 'Velocities':
                    velocitiesstart = i + 1
    
    # Create system
    box = Box(xlo=xlo, xhi=xhi,
              ylo=ylo, yhi=yhi,
              zlo=zlo, zhi=zhi,
              xy=xy, xz=xz, yz=yz)
    atoms = Atoms(natoms=natoms)
    system = System(box=box, atoms=atoms, pbc=pbc)
    
    # Read in Atoms info
    if atomsstart is not None:
        prop_info = atoms_prop_info(atom_style, units)
        system = load_table(data, box=system.box, system=system, symbols=symbols,
                            prop_info=prop_info, skiprows=atomsstart, nrows=natoms)
    else:
        raise ValueError('No Atoms section found!')
    
    # Read in Velocities info
    if velocitiesstart is not None:
        prop_info = velocities_prop_info(atom_style, units)
        system = load_table(data, box=system.box, system=system,
                            prop_info=prop_info,
                            skiprows=velocitiesstart, nrows=natoms)
    
    return system