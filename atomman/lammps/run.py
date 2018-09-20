# Standard Python libraries
from __future__ import (absolute_import, print_function,
                        division, unicode_literals)
import os
import glob
import shutil
import subprocess as sp

# atomman imports
from .Log import Log
from ..compatibility import range, stringtype

def run(lammps_command, script_name, mpi_command=None,
        restart_script_name=None, return_style='object', logfile='log.lammps',
        flatten=None):
    """
    Calls LAMMPS to run. Returns data model containing LAMMPS output.
    
    Parameters
    ----------
    lammps_command : str
        The LAMMPS inline run command (sans -in script_name).
    script_name : str
        Path of the LAMMPS input script to use.
    mpi_command : str or None, optional
        The MPI inline command to run LAMMPS in parallel. Default value is 
        None (no mpi).
    restart_script_name : str or None, optional
        Alternative script to use for restarting if logfile already exists.
        Default value is None (no restarting).
    return_style :str, optional
        The format for the returned data. Default value is 'object'.
        'object' -- returns an atomman.lammps.Log object.
        'model' -- returns a DataModelDict.
    logfile : str, optional
        Specifies the path to the logfile to write to.  Default value is
        'log.lammps'.
    flatten : str or None, optional
        Specifies if the simulations are to be flattened, and which flatten
        style to use:
        - None does not flatten the simulations (default).
        - 'first' uses the values from the earliest simulation.
        - 'last' uses the values from the latest simulation.
    
    Returns
    -------
    atomman.lammps.Log or DataModelDict
        The content either as a Log object or in data model format.
    """
    
    assert return_style in ['object', 'model'], 'Invalid return_style value'
    if flatten is not None:
        assert flatten in ['first', 'last'], 'Invalid flatten value'
    
    # Check if restart_script_name is given
    if restart_script_name is not None:
        
        # Check if simulation was previously started by looking for log.lammps
        if os.path.isfile(logfile):
            logname, logext = os.path.splitext(logfile)
            # Replace script_name with restart_script_name
            script_name = restart_script_name
            
            # Search for any earlier log files with the name log-*.lammps
            logids = []
            for oldlog in glob.iglob(logname + '-*' + logext):
                logids.append(int(os.path.splitext(os.path.basename(oldlog))[0][len(logname)+1:]))
            
            # Rename old logfile to keep it from being overwritten
            if len(logids) == 0:
                lognum = 1
            else:
                lognum = max(logids)+1
            shutil.move(logfile, logname + '-' + str(lognum) + logext)
        else:
            lognum = 0
    else:
        lognum = 0
    
    # Convert lammps_command into list of terms
    if isinstance(lammps_command, stringtype):
        lammps_command = lammps_command.split(' ')
    elif not isinstance(lammps_command, list):
        lammps_command = [lammps_command]
    
    # Convert script_name into list of terms
    if isinstance(script_name, stringtype):
        script_name = script_name.split(' ')
    elif not isinstance(script_name, list):
        script_name = [script_name]
    
    # Convert mpi_command into list of terms
    if mpi_command is None:
        mpi_command = []
    elif isinstance(mpi_command, stringtype):
        mpi_command = mpi_command.split(' ')
    elif not isinstance(mpi_command, list):
        mpi_command = [mpi_command]
    
    # Extra terms
    extra = []
    
    # Check if logfile is not log.lammps
    if logfile != 'log.lammps':
        # Convert logfile into list of terms
        extra += ['-log'] + logfile.split(' ')
    
    # Try to run lammps as a subprocess
    try:
        output = sp.check_output(mpi_command + lammps_command + extra + ['-in'] + script_name)
    
    # Pass LAMMPS error to a Python error if failed
    except sp.CalledProcessError as e:
        if e.output != '':
            lines = e.output.decode("utf-8").split('\n')
            raise ValueError('Invalid LAMMPS input: \n%s' % lines[-2])
        else:
            raise OSError('Failed to run LAMMPS')
    
    # Initialize Log object
    log = Log()
    
    # Read in from all old log files
    for i in range(1, lognum+1):
        log.read(logname + '-' + str(i) + logext)
    
    # Read in from current logfile
    log.read(logfile)
    
    if flatten is not None:
        log.flatten(flatten)
    
    if return_style == 'model':
        return log.model()
    else:
        return log