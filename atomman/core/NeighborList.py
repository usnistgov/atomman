# Standard Python libraries
from __future__ import (absolute_import, print_function,
                        division, unicode_literals)
from copy import deepcopy

# http://www.numpy.org/
import numpy as np

# atomman imports
from .nlist import nlist
from ..tools import uber_open_rmode
from ..compatibility import range

class NeighborList(object):
    """Class that finds and stores the neighbor atoms for a system."""
    
    def __init__(self, **kwargs):
        """
        Class initializer.  Calls NeighborList.load() if
        'model' is given, otherwise calls NeighborList.build().
        
        Parameters
        ----------
        system : atomman.System, optional
            The system to calculate the neighbor list for. Must be given if
            model is not given.
        cutoff : float, optional
            Radial cutoff distance for identifying neighbors.  Must be given if
            model is not given.
        model : str or file-like object, optional
            Gives the file path or content to load.  If given, no other
            parameters are allowed.
        cmult : int, optional
            Parameter associated with the binning routine.  Default value is most
            likely the fastest.
        code: str, optional
            Option for specifying which underlying code function of nlist to use:
            - 'cython' uses the version of the function built in cython (faster).
            - 'python' uses the purely python version.
            Default is 'cython' if the code can be imported, otherwise 'python'.
        initialsize : int, optional
            The number of neighbor positions to initially assign to each atom.
            Default value is 20.
        deltasize : int, optional
            Specifies the number of extra neighbor positions to allow each atom
            when the number of neighbors exceeds the underlying array size.
            Default value is 10.
        """
        if 'model' in kwargs:
            model = kwargs.pop('model')
            self.load(model, **kwargs)
        else:
            system = kwargs.pop('system')
            cutoff = kwargs.pop('cutoff')
            self.build(system, cutoff, **kwargs)
    
    @property
    def coord(self):
        """The atomic coordination numbers"""
        return deepcopy(self.__coord)
    
    def __len__(self):
        """len returns the number of atoms"""
        return len(self.__coord)
    
    def __getitem__(self, key):
        """Get returns the list of neighbors for the specified atom."""
        return self.__neighbors[key, :self.coord[key]]
    
    def build(self, system, cutoff, cmult=1, code=None, initialsize=20,
              deltasize=10):
        """
        Builds the neighbor list for a system.
        
        Parameters
        ----------
        system : atomman.System 
            The system to calculate the neighbor list for.
        cutoff : float
            Radial cutoff distance for identifying neighbors.
        cmult : int, optional
            Parameter associated with the binning routine. Default value is most
            likely the fastest.
        code: str, optional
            Option for specifying which underlying code function of nlist to use:
            - 'cython' uses the version of the function built in cython (faster).
            - 'python' uses the purely python version.
            Default is 'cython' if the code can be imported, otherwise 'python'.
        initialsize : int, optional
            The number of neighbor positions to initially assign to each atom.
            Default value is 20.
        deltasize : int, optional
            Specifies the number of extra neighbor positions to allow each atom
            when the number of neighbors exceeds the underlying array size.
            Default value is 10.
        """
        # Call nlist
        neighbors = nlist(system, cutoff, cmult=cmult, code=code,
                          initialsize=initialsize, deltasize=deltasize)
        
        # Split coord and neighbors
        self.__coord = neighbors[:, 0]
        self.__neighbors = neighbors[:, 1:]
        
    def load(self, model):
        """
        Read in a neighbor list from a file.
        
        Parameters
        ----------
        model : str or file-like object
            Gives the file path or content to load.
        """
        # First pass determines number of atoms and max number of neighbors
        nterms = 0
        natoms = 0
        with uber_open_rmode(model) as fin:
            for line in fin:
                line = line.decode('UTF-8')
                terms = line.split()
                n_n = len(terms)
                if terms[0][0] != '#' and n_n > 0:
                    natoms += 1
                    if n_n > nterms:
                        nterms = n_n
            
            self.__coord = np.zeros(natoms, dtype=int)
            self.__neighbors = np.empty((natoms, nterms), dtype=int)
            
            # Second pass gets values
            fin.seek(0)
            for line in fin:
                line = line.decode('UTF-8')
                terms = line.split()
                if len(terms) > 0 and terms[0][0] != '#':
                    i = int(terms[0])
                    self.__coord[i] = len(terms) - 1
                    for j in range(1, len(terms)):
                        self.__neighbors[i, j-1] = terms[j]
    
    def dump(self, fname):
        """
        Saves the neighbor list to a file.
        
        Parameters
        ----------
        fname : str
            The file name to save the content to.
        """
        with open(fname, 'w') as fp:
            fp.write('# Neighbor list:\n')
            fp.write('# The first column gives an atom index.\n')
            fp.write('# The rest of the columns are the indexes of the identified neighbors.\n')
            for i in range(len(self)):
                fp.write('%i' % i)
                for j in self[i]:
                    fp.write(' %i' % j)
                fp.write('\n')