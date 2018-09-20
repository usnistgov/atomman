# Standard Python libraries
from __future__ import (absolute_import, print_function,
                        division, unicode_literals)

def filltemplate(template, variable, s_delimiter, e_delimiter):
    """
    Takes a template and fills in values for delimited template variables.
    
    Parameters
    ----------
    template : string or file-like object
        The template file or file content to fill in.
    variable : dict 
        Dictionary with keys defining the delimited template variable terms,
        and values the values to replace the variable terms with.
    s_delimiter : str
        The leading delimiter for identifying the template variable terms.
    e_delimiter : str 
        The trailing delimiter for identifying the template variable terms.
    
    Returns
    -------
    str
        The template with all delimited variable terms replaced with their
        corresponding defined values from variable.
        
    Raises
    ------
    KeyError
        If delimited term found in template that has no value in variable.
    ValueError
        If parsing of s_delimiter, e_delimiter pairs fails.
    """
    
    # Convert to string if a file-like object
    try:
        template = template.read()
    except AttributeError:
        pass
    
    # Loop until done
    while True:
        
        # Search for starting delimiter
        try:
            s = template.index(s_delimiter)
        except ValueError: 
            s = None
        else:
            s = s + len(s_delimiter)
        
        # Search for ending delimiter
        try:
            e = template.index(e_delimiter)
        except ValueError:
            e = None        
        
        # Replace delimited string with value
        if s is not None and e is not None and s < e:
            name = template[s: e]
            var = s_delimiter + name + e_delimiter
            try:
                value = str(variable[name])
            except KeyError:
                raise KeyError(name + ' not found in variable dictionary')
            template = template.replace(var, value)
        
        # Finish if no delimiters remain
        elif s is None and e is None:
            break
            
        # Issue errors
        elif s is None:
            raise ValueError('ending delimiter found without starting delimiter')
        elif e is None:
            raise ValueError('starting delimiter found without ending delimiter')
        else:
            raise ValueError('ending delimiter found before starting delimiter')
            
    return template