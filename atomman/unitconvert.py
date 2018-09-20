# Standard Python libraries
from __future__ import (absolute_import, print_function,
                        division, unicode_literals)
import ast

# https://pypi.python.org/pypi/numericalunits
import numericalunits as nu

# https://github.com/usnistgov/DataModelDict
from DataModelDict import DataModelDict as DM

# http://www.numpy.org/
import numpy as np

# atomman imports
from .compatibility import iteritems, range, stringtype

def build_unit():
    """
    Saves numericalunits attributes to global dictionary unit so the values
    can be retrieved by their string names.
    """
    
    # Define global dictionary
    global unit
    unit = {}
    
    # Copy all float attributes of numericalunits to unit
    for key, value in iteritems(nu.__dict__):
        try:
            key = key.decode('UTF-8')
        except:
            pass
        if key[:1] != '_' and isinstance(value, float):
            unit[key] = value

def reset_units(seed=None, **kwargs):
    """
    Extends numericalunits.reset_units() by allowing for working units to be
    defined.  If no working units are specified, then random working units are
    used just like the default numericalunits behavior.  Otherwise, use the
    specified working units and SI.
    
    Parameters
    ----------
    seed : int, optional
        random number seed to use in generating random working units.
        seed='SI' will use SI units.  Cannot be given with the other
        parameters.
    length : str, optional
        Unit of length to use for the working units.
    mass : str, optional
        Unit of mass to use for the working units.
    time : str, optional
        Unit of time to use for the working units.
    energy : str, optional
        Unit of energy to use for the working units.
    charge : str, optional
        Unit of charge to use for the working units.
        
    Raises
    ------
    ValueError
        If seed is given with any other parameters, or if more than four of
        the working unit parameters are given.
    """
    
    # Generate random base working units
    if (len(kwargs) == 0):
        
        nu.reset_units(seed)
        build_unit()
    
    # Generate SI + defined working units
    elif seed is None:
        
        # Check that no more than 4 working units are defined
        if len(kwargs) > 4:
            raise ValueError('Only four working units can be defined')
        
        # Set base units to 1 (working units to SI)
        nu.reset_units('SI')
        build_unit()
        
        # Scale base units by working units
        if 'length' in kwargs:
            nu.m = unit['m'] / unit[kwargs['length']]
        
        if 'mass' in kwargs:
            nu.kg = unit['kg'] / unit[kwargs['mass']]
        
        if 'time' in kwargs:
            nu.s = unit['s'] / unit[kwargs['time']]
        
        if 'charge' in kwargs:
            nu.C = unit['C'] / unit[kwargs['charge']]
        
        # Scale derived units by working units
        if 'energy' in kwargs:
            J = unit['J'] / unit[kwargs['energy']]
            
            # Scale base units by derived units
            if 'mass' not in kwargs:
                nu.kg = J * nu.s**2 / nu.m**2
            elif 'time' not in kwargs:
                nu.s = (nu.kg * nu.m**2 / J)**0.5
            elif 'length' not in kwargs:
                nu.m = (J * nu.s**2 / nu.kg)
        
        # Rebuild derived units and unit dictionary
        nu.set_derived_units_and_constants()
        build_unit()
        
    else:
        raise ValueError('seed cannot be given with any other parameters')

def set_literal(term):
    """
    Convert string 'value unit' to numbers in working units.
    
    Parameters
    ----------
    term : str
        String containing value and associated unit. If unit is not given,
        then the value is converted to a float and assumed to be in working
        units.
        
    Returns
    -------
    numpy.ndarray
        The numerical value of term in working units.
        
    Raises
    ------
    ValueError
        If no valid float value can be parsed.
    """
    
    # Set splitting point j to end of term (i.e. assume no units given)
    j = len(term)
    
    # Loop until done
    while True:
        
        # Split term into value, unit terms
        value = term[:j].strip()
        unit = term[j:].strip()
        if len(unit) == 0:
            unit = None

        # Return number if value, unit pair is valid
        try: 
            return set_in_units(ast.literal_eval(value), unit)
        except: 
            # Find the next splitting point
            try: 
                j = term[:j].rindex(' ')
            except: 
                raise ValueError('Failed to parse term')

def set_in_units(value, units):
    """
    Convert value from specified units to working units.
    
    Parameters
    ----------
    value : int, float, numpy.ndarray, etc.
        A numerical value or list/array of values.
    units : str
        The units that value is in.
        
    Returns
    -------
    numpy.ndarray
        The given value converted from the specified units to working units.
    """
    units = parse(units)
    return np.asarray(value) * units

def get_in_units(value, units):
    """
    Convert value from working units to specified units.
    
    Parameters
    ----------
    value : int, float, numpy.ndarray, etc.
        A numerical value or list/array of values.
    units : str
        The units to convert value to (from working units).
        
    Returns
    -------
    numpy.ndarray
        The given value converted to the specified units from working units.
    """
    units = parse(units)
    return np.asarray(value) / units

def value_unit(term):
    """
    Reads numerical value from dictionary containing 'value' and 'unit' keys.
    
    Parameters
    ----------
    term : dict
        Dictionary containing 'value' and 'unit' keys.
        
    Returns
    -------
    numpy.ndarray
        The result of calling set_in_units() by passing the dictionary keys 
        'value' and 'unit' as parameters.
    
    """
    unit = term.get('unit', None)
    value = set_in_units(term['value'], unit)
    
    if 'shape' in term:
        shape = tuple(term['shape'])
        value = value.reshape(shape)
    
    return value
    
def error_unit(term):
    """
    Reads numerical error from dictionary containing 'error' and 'unit' keys.
    
    Parameters
    ----------
    term : dict
        Dictionary containing 'error' and 'unit' keys.
        
    Returns
    -------
    numpy.ndarray
        The result of calling set_in_units() by passing the dictionary keys 
        'error' and 'unit' as parameters.
    
    """
    unit = term.get('unit', None)
    error = set_in_units(term['error'], unit)
    
    if 'shape' in term:
        shape = tuple(term['shape'])
        error = error.reshape(shape)
    
    return error
    
def model(value, units, error=None):
    """
    Generates DataModelDict representation of data.
    
    Parameters
    ----------
    value : float, numpy.ndarray, etc.
        A numerical value or list/array of values.
    units : str
        The units to convert value to (from working units).
    error : float, numpy.ndarray, etc., optional
        A value error to include.  If given, must be the same
        size/shape as value.
    
    Returns
    -------
    DataModelDict
        Model representation of the value(s).
    """
    
    datamodel = DM()
    
    value = get_in_units(value, units)
    if error is not None:
        error = get_in_units(error, units)
    
    # Single value
    if value.ndim == 0:
        datamodel['value'] = value
        if error is not None:
            datamodel['error'] = error
    
    # 1D array
    elif value.ndim == 1:
        datamodel['value'] = list(value)
        if error is not None:
            datamodel['error'] = list(error)
            
    # Higher-order array requires shape
    else:
        shape = value.shape
        datamodel['value'] = list(value.flatten())
        if error is not None:
            datamodel['error'] = list(error.flatten())
        datamodel['shape'] = list(shape)
    
    datamodel['unit'] = units
    return datamodel

def parse(units):
    """
    Convert units as strings (or None) into scaling numbers.  This function
    allows for complex unit definitions with operators:
    
    - '()' for defining order of operations
    - '*' for multiplication.
    - '/' for division.
    - '^' for powers.
    
    Parameters
    ----------
    units : str or None
        String consisting of defined unit names, operators, and numerical 
        values to interpret.
        
    Returns
    -------
    float
        The scaling factor for converting numbers in the given units to
        working units. If units is None or == 'scaled', then this value is
        1.0.
    """
    
    # Units of None does no scaling
    if units is None or units == 'scaled':
        return 1

    # Parse string and return number value
    elif isinstance(units, stringtype):
        i = 0
        terms = []
        
        # Break into terms
        while i < len(units):
            
            # parse terms in parentheses first
            if units[i] == '(':
                j = i+1
                pcount = 0
                while True:
                    if j == len(units):
                        raise ValueError('Invalid () terms.')
                    elif units[j] == ')':
                        if pcount == 0:
                            break
                        else:
                            pcount -= 1
                    elif units[j] == '(':
                        pcount += 1
                    j += 1
                    
                terms.append(parse(units[i+1:j]))
                i = j+1
            
            # append string terms
            elif units[i].isalpha():
                term = ''
                while i < len(units) and units[i] not in ' */^\n\r\t':
                    term += units[i]
                    i += 1
                terms.append(unit[term])
            
            # append numeric terms
            elif units[i].isdigit() or units[i] == '-' or units[i] == '.':
                term = ''
                while i < len(units) and units[i] not in ' */^\n\r\t':
                    term += units[i]
                    i += 1
                terms.append(float(term))
            
            # append operators
            elif units[i] in '*/^':
                terms.append(units[i])
                i += 1
            
            # ignore excess white characters
            elif units[i] in ' \n\r\t':
                i += 1
            
            # issue error for unmatched ) parentheses
            elif units[i] == ')':
                raise ValueError('Invalid () terms.')
            
            else:
                raise ValueError('Unknown character: %s' % units[i])
        
        # Compute powers
        while '^' in terms:
            c = terms.index('^')
            value = [terms[c-1] ** terms[c+1]]
            terms = terms[:c-1] + value + terms[c+2:]
        
        # Compute multiplication and division
        while len(terms) > 1:
            if terms[1] == '*':
                value = [terms[0] * terms[2]]
                terms = value + terms[3:]
            elif terms[1] == '/':
                value = [terms[0] / terms[2]]
                terms = value + terms[3:]
            else:
                raise ValueError('Invalid string format')
        
        return terms[0]

    # Else assume units is already a number
    else:
        return units

# Initial reset.  Only called first time module is loaded.
reset_units()