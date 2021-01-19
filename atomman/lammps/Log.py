# coding: utf-8
# Standard Python libraries
import datetime

# http://www.numpy.org/
import numpy as np

# https://pandas.pydata.org/
import pandas as pd

# https://github.com/usnistgov/DataModelDict
from DataModelDict import DataModelDict as DM

# atomman imports
from ..tools import uber_open_rmode



class Simulation():
    """
    Object for representing the LAMMPS log output for a single MS/MD run
    """

    def __init__(self, thermo=None, perfo=None):
        """
        Initializes a Simulation object with performance information
        """
        self.__thermo = None
        self.__performance = None
        self.__keys = []
        if thermo is not None:
            self.thermo = thermo
        if perfo is not None:
            self.performance = perfo

    def keys(self):
        """List of attribute keys that have been set"""
        return tuple(self.__keys)

    def __getitem__(self, key):
        """Accesses class attributes as dict keys""" 

        if key not in self.keys():
            raise KeyError(key)

        return getattr(self, key)

    def __iter__(self):
        """Iterates over set attribute keys"""
        for key in self.keys():
            yield(key)

    @property
    def thermo(self):
        """pandas.DataFrame: The Simulation's thermo data"""
        return self.__thermo

    @property
    def performance(self):
        """pandas.DataFrame: The Simulation's performance data"""
        return self.__performance

    @thermo.setter
    def thermo(self, value):
        if isinstance(value, pd.DataFrame):
            self.__thermo = value
        else:
            self.__thermo = pd.DataFrame(value)
        if 'thermo' not in self.keys():
            self.__keys.append('thermo')
            
    @performance.setter
    def performance(self, value):
        if isinstance(value, pd.DataFrame):
            self.__performance = value
        else:
            self.__performance = pd.DataFrame(value)
        if 'performance' not in self.keys():
            self.__keys.append('performance')


class Log():
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
        self.__performance_data = []
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
            self.__performance_data = []
            self.__lammps_version = None
            self.__lammps_date = None
        
        # Strings found directly before and after run and minimize simulations
        sim_start_trigger = ['Memory usage per processor =',
                             'Per MPI rank memory allocation (min/avg/max) =']
        sim_end_trigger = ['Loop time of']
        
        perf_start_trigger = ['MPI task timing breakdown']
        perf_start_trigger_old_version = ['Pair  time (%)']
        perf_end_trigger = ['Nlocal:']
        is_old_version = False

        # Handle file names, strings and open file-like objects equivalently
        with uber_open_rmode(log_info) as log_info:
            
            # Initialize parameters
            headers = []
            footers = []
            headers_perf = []
            footers_perf = []
            i = 0
            
            # For all lines in file/output
            for line in log_info:
                line = line.decode('UTF-8')
                
                # Skip blank lines
                if len(line.split()) == 0:
                    continue
                
                # Save the LAMMPS version information
                if line[:8] == 'LAMMPS (' and self.lammps_version is None:
                    self.__read_lammps_version(line)
                
                # Check for strings listed prior to run and minimize simulations
                if any([trigger in line for trigger in sim_start_trigger]):
                    headers.append(i+1)
                
                # Check for strings listed after run and minimize simulations
                elif any([trigger in line for trigger in sim_end_trigger]):
                    footers.append(i-1)

                # Check for strings listed prior to  perfomance data
                if any([trigger in line for trigger in perf_start_trigger]):
                    headers_perf.append(i+1)
                if any([trigger in line for trigger in perf_start_trigger_old_version]):
                    headers_perf.append(i)
                    is_old_version = True
                
                # Check for strings listed after performance data
                elif any([trigger in line for trigger in perf_end_trigger]):
                    footers_perf.append(i-1)
                
                i += 1
            
            # Add last line to footers for incomplete logs
            footers.append(i)
            
            # Reset file pointer
            log_info.seek(0)
            
            # Read data from all simulations
            for header, footer, header_perf, footer_perf in zip(headers, footers, headers_perf, footers_perf):
                self.__read_simulation(log_info, header, footer, header_perf, footer_perf, is_old_version)
           
    
    def __read_lammps_version(self, line):
        """
        Subfunction for reading the LAMMPS version from the log file
        """
        month = {'Jan': 1, 'Feb': 2, 'Mar': 3, 'Apr': 4, 
                    'May': 5, 'Jun': 6, 'Jul': 7, 'Aug': 8,
                    'Sep': 9, 'Oct': 10,'Nov': 11,'Dec': 12}
        self.__lammps_version = line.strip()[8:-1]
        d = self.lammps_version.split('-')[0].split()
        self.__lammps_date = datetime.date(int(d[2]), month[d[1]], int(d[0]))

    def __read_simulation(self, log_info, header, footer, header_perf, footer_perf, is_old_version):
        """
        Subfunction for reading the log data associated with a simulation run.
        """
        # Use pandas to read all thermo data at once
        
        thermo = pd.read_csv(log_info, header=header,
                                nrows=footer-header,
                                delim_whitespace=True,
                                skip_blank_lines=True)
                          
        log_info.seek(0)

        # Use pandas to read all performance data at once        
        if not is_old_version:
            performance = pd.read_csv(log_info, header=header_perf,
                                nrows=footer_perf-header_perf,
                                sep = '|',
                                skip_blank_lines=True)  
            performance = performance.drop([0])
            performance.rename(columns=lambda x: x.strip(), inplace=True)
            performance = performance.set_index('Section')
            performance = performance.replace(r'^\s*$',0.0,regex=True)
            performance=performance.astype(float) 
        else: 
            performance = pd.read_csv(log_info, header=header_perf,
                                nrows=footer_perf-header_perf,
                                sep = '=',
                                skip_blank_lines=True)
            performance = performance.columns.to_frame().T.append(performance,ignore_index=True)
            performance.rename(columns=lambda x: x.strip(), inplace=True)
            performance = performance.replace(')','')
            performance.columns = ['Section','time']
            performance = performance.set_index('Section')
            #performance['time'].replace(')','',inplace=True)
            performance[['avg. Time','percentage']] = performance.time.str.split('(',expand=True)
            performance[['%','symbol']] = performance.percentage.str.split(')',expand=True)
            del performance['time']
            del performance['symbol']
            del performance['percentage']



                           

        # Extract any other simulation information
        # ...
        
        # Reset file pointer
        log_info.seek(0)
        
        # Append simulation results
        self.__simulations.append(Simulation(thermo=thermo, perfo=performance))



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
            if sim.thermo is not None and len(sim.thermo) > 0:
                assert 'Step' in sim.thermo, 'All simulation thermos must have Step key in order to flatten'
        
        # Combine the data into merged_df
        merged_df = self.simulations[0].thermo
        
        for i in range(1, len(self.simulations)):
            if self.simulations[i].thermo is not None:
                thermo = self.simulations[i].thermo
                if style == 'first':
                    merged_df = pd.concat([merged_df, thermo[thermo.Step > merged_df.Step.max()]], ignore_index=True)
                elif style == 'last':
                    merged_df = pd.concat([merged_df[merged_df.Step < thermo.Step.min()], thermo], ignore_index=True)
                else:
                    raise ValueError('Unsupported style')
        
        #This would add the flattened data as a new simulation
        #self.__simulations.append(Simulation(thermo = merged_df))
        
        return Simulation(thermo=merged_df)
    
    def model(self, flatten=None):
        """
        Returns an XML/JSON equivalent data mode of the thermo information
        
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
            simulations=[self.flatten(flatten)]
        else:
            simulations=self.simulations
        
        # Create DataModelDict root
        log_model = DM()
        log_model['LAMMPS-log-thermo-data'] = DM()
        
        # Loop over all simulations
        for sim in simulations:
            sim_model = DM()
            
            # Convert to DataModelDict
            sim_model['thermo'] = DM()
            thermo = sim.thermo
            for j in thermo:
                sim_model['thermo'][str(j)] = thermo[j].values.tolist()
            
            # Append simulation results to DataModelDict root
            log_model['LAMMPS-log-thermo-data'].append('simulation', sim_model)
            
        return log_model