# Standard Python libraries
from __future__ import (absolute_import, print_function,
                        division, unicode_literals)

# atomman imports
from ..compatibility import range

def indexstr(shape):
    """
    Iterates through all unique indicies of an array with a given shape.
    
    Parameters
    ----------
    shape : tuple of int
        The array shape to iterate through.
        
    Yields
    ------
    index : tuple of int
        A unique index set of the array.
    istr : str
        A string representation of index with numbers in [].
    """
    
    if tuple(shape) == ():
        # Yield for empty shape
        yield (), ''
    else:
        # Loop over all values of the first index of shape
        for i in range(shape[0]):
            index1 = (i,)
            istr1 = '['+str(i)+']'
            
            # Recursively go through other indicies of shape
            for index2, istr2 in indexstr(shape[1:]):
                
                # Combine index components
                index = index1 + index2
                istr = istr1 + istr2
                
                # Yield composite
                yield index, istr