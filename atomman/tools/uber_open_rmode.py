# Standard Python libraries
from __future__ import (absolute_import, print_function,
                        division, unicode_literals)
from io import BytesIO, open
import os

# atomman imports
from ..compatibility import unicode

class uber_open_rmode():
    """
    Context manager for reading data from file-like objects, file names,
    and data strings in the same manner.
    """
    
    def __init__(self, data):
        """
        Initialize context manager.
        
        Parameters
        ----------
        data : file-like object or str
            Specifies what content to read.
        """
        self.data = data
        
    def __enter__(self):
        """Define different open actions."""
        def isfile(data):
            try:
                return os.path.isfile(self.data)
            except:
                return False
        
        # Any open file-like objects in 'r' mode will have a read attribute
        if hasattr(self.data, 'read'):
            self.open_file = self.data
            self.to_close = False
        
        # If data is path to a file, open the file
        elif isfile(self.data):
            self.open_file = open(self.data, 'rb')
            self.to_close = True
        
        # If data is not unicode, read using BytesIO
        elif not isinstance(self.data, unicode):
            self.open_file = BytesIO(self.data)
            self.to_close = True
            
        # If data is unicode, encode and read using BytesIO
        else:
            self.open_file = BytesIO(self.data.encode('utf-8'))
            self.to_close = True
        
        return self.open_file
    
    def __exit__(self, *args):
        """Close file if one was opened."""
        if self.to_close:
            self.open_file.close()