# Standard Python libraries
from __future__ import (absolute_import, print_function,
                        division, unicode_literals)
import datetime

# http://www.numpy.org/
import numpy as np

# https://pandas.pydata.org/
import pandas as pd

# https://github.com/usnistgov/DataModelDict
from DataModelDict import DataModelDict as DM

# atomman imports
from ..tools import uber_open_rmode
from ..compatibility import range

class Log(object):
    """Object for representing a LAMMPS log output"""
    
    def __init__(self, log_info=None):
        """
        Initializes a Log object.
        
        Parameters
        ----------
        log_info : str or file-like object, optional
            The LAMMPS log content to read in.  If None (default), then the
            Log object is created but empty.
        """
        
        # Initialize simulation properties
        self.__simulations = []
        self.__lammps_version = None
        self.__lammps_date = None
        
        # Read log data if supplied
        if log_info is not None:
            self.read(log_info)
    
    def read(self, log_info, append=True):
        """
        Parses a LAMMPS screen output/log file.
        
        Parameters
        ----------
        log_info : str or file-like object
            The LAMMPS log content to read in.
        append : bool, optional
            Flag indicating if the content being read in is appended to the
            current data (True, default), or if it overwrites any existing
            saved content (False).
        """
        
        # Reset properties and values if append is False
        if append is False:
            self.__simulations = []
            self.__lammps_version = None
            self.__lammps_date = None
        
        # Strings found directly before and after run and minimize simulations
        sim_start_trigger = ['Memory usage per processor =',
                             'Per MPI rank memory allocation (min/avg/max) =']
        sim_end_trigger = ['Loop time of']
        
        # Handle file names, strings and open file-like objects equivalently
        with uber_open_rmode(log_info) as log_info:
            
            # Initialize parameters
            headers = []
            footers = []
            i = 0
            
            # For all lines in file/output
            for line in log_info:
                line = line.decode('UTF-8')
                # Skip blank lines
                if len(line.split()) == 0:
                    continue
                
                # Save the LAMMPS version information
                if line[:8] == 'LAMMPS (' and self.lammps_version is None:
                    month = {'Jan': 1, 'Feb': 2, 'Mar': 3, 'Apr': 4, 
                             'May': 5, 'Jun': 6, 'Jul': 7, 'Aug': 8,
                             'Sep': 9, 'Oct': 10,'Nov': 11,'Dec': 12}
                    self.__lammps_version = line.strip()[8:-1]
                    d = self.lammps_version.split('-')[0].split()
                    self.__lammps_date = datetime.date(int(d[2]), month[d[1]], int(d[0]))
                
                # Check for strings listed prior to run and minimize simulations
                if any([trigger in line for trigger in sim_start_trigger]):
                    headers.append(i+1)
                
                # Check for strings listed after run and minimize simulations
                elif any([trigger in line for trigger in sim_end_trigger]):
                    footers.append(i-1)
                
                i += 1
            
            # Add last line to footers for incomplete logs
            footers.append(i)
            
            # Reset file pointer
            log_info.seek(0)
            
            # For all lines in file/output
            for header, footer in zip(headers, footers):
                
                # Initialize simulation data dictionary
                sim = {}
                
                # Use pandas to read all thermo data at once
                sim['thermo'] = pd.read_csv(log_info, header=header,
                                            nrows=footer-header,
                                            delim_whitespace=True,
                                            skip_blank_lines=True)
                
                # Reset file pointer
                log_info.seek(0)
                
                # Append simulation results
                self.__simulations.append(sim)
    
    @property
    def simulations(self):
        """list of dict: parsed data for each separate LAMMPS run/minimize action"""
        return self.__simulations
    
    @property
    def lammps_version(self):
        """str : The LAMMPS version used."""
        return self.__lammps_version
    
    @property
    def lammps_date(self):
        """datetime.date : The date associated with the LAMMPS version used."""
        return self.__lammps_date
    
    def flatten(self, style='last'):
        """
        Combines all simulations into one.
        
        Parameters
        ----------
        style : str, optional
            Specifies which values to use for duplicate time steps:
            - 'first' uses the values from the earliest simulation.
            - 'last' uses the values from the latest simulation (default).
        """
        # Check that all simulations with thermo data have step values
        for sim in self.simulations:
            if 'thermo' in sim and len(sim['thermo']) > 0:
                assert 'Step' in sim['thermo'], 'All simulation thermos must have Step key in order to flatten'
        
        # Combine the data into merged_df
        merged_df = self.simulations[0]['thermo']
        for i in range(1, len(self.simulations)):
            if 'thermo' in self.simulations[i]:
                thermo = self.simulations[i]['thermo']
                if style == 'first':
                    merged_df = pd.concat([merged_df, thermo[thermo.Step > merged_df.Step.max()]], ignore_index=True)
                elif style == 'last':
                    merged_df = pd.concat([merged_df[merged_df.Step < thermo.Step.min()], thermo], ignore_index=True)
                else:
                    raise ValueError('Unsupported style')
        
        self.__simulations = [{'thermo':merged_df}]
    
    def model(self, flatten=None):
        """
        Returns an XML/JSON equivalent data mode of the information
        
        Parameters
        ----------
        flatten : str or None, optional
            Specifies if the simulations are to be flattened, and which
            flatten style to use:
            - None does not flatten the simulations (default).
            - 'first' uses the values from the earliest simulation.
            - 'last' uses the values from the latest simulation.
            
        Returns
        -------
        DataModelDict
            The Log content in data model form.
        """
        if flatten is not None:
            self.flatten(flatten)
        
        # Create DataModelDict root
        log_model = DM()
        log_model['LAMMPS-log-thermo-data'] = DM()
        
        # Loop over all simulations
        for sim in self.simulations:
            sim_model = DM()
            
            # Convert to DataModelDict
            sim_model['thermo'] = DM()
            thermo = sim['thermo']
            for j in thermo:
                sim_model['thermo'][str(j)] = thermo[j].values.tolist()
            
            # Append simulation results to DataModelDict root
            log_model['LAMMPS-log-thermo-data'].append('simulation', sim_model)
            
        return log_model