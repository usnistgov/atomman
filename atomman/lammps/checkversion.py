# Standard Python libraries
from __future__ import (absolute_import, print_function,
                        division, unicode_literals)
import os
import shutil
import tempfile

# atomman imports
from . import run

def checkversion(lammps_command):
    """
    Gets the LAMMPS version and date information by testing lammps_command.
    
    Parameters
    ----------
    lammps_command: str
        A LAMMPS executable.
        
    Returns
    -------
    version_info: dict
        Dictionary containing 'version', the str LAMMPS version, and
        'date', the corresponding  datetime.date for the LAMMPS 
        version.
        
    Raises
    ------
    ValueError
        If lammps fails to run
    """
    # Define emptyscript and logfile paths
    tempdir = tempfile.mkdtemp()
    emptyscript = os.path.join(tempdir, 'empty.in')
    logfile = os.path.join(tempdir, 'empty.lammps')
    
    # Create emptyscript
    with open(emptyscript, 'w') as f:
        f.write('')
    
    # Run lammps_command with emptyscript
    try:
        log = run(lammps_command, emptyscript, logfile=logfile)
    except:
        raise ValueError('Failed to run simulation with lammps_command '+lammps_command)
    finally:
        shutil.rmtree(tempdir)
    
    # Extract lammps version and date info
    version_info = {}
    version_info['version'] = log.lammps_version
    version_info['date'] = log.lammps_date
    return version_info