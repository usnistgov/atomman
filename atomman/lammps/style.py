# Standard Python libraries
from __future__ import (absolute_import, print_function,
                        division, unicode_literals)
from collections import OrderedDict

def unit(units='metal'):
    """
    Returns the units information associated with a LAMMPS units option.
    
    Parameters
    ----------
    units : str, optional
        The LAMMPS unit option.  Default value is 'metal'
    
    Returns
    -------
    dict
        The keys are physical quantity types and the values are the associated
        units to use.
    """
    params = OrderedDict()
    
    if units == 'lj':
        params['mass'] =                None
        params['length'] =              None
        params['time'] =                None
        params['energy'] =              None
        params['velocity'] =            None
        params['force'] =               None
        params['torque'] =              None
        params['temperature'] =         None
        params['pressure'] =            None
        params['dynamic viscosity'] =   None
        params['charge'] =              None
        params['dipole'] =              None
        params['electric field'] =      None
        params['density'] =             None
    
    elif units == 'real':
        params['mass'] =                'g/mol'
        params['length'] =              'angstrom'
        params['time'] =                'fs'
        params['energy'] =              'kcal/mol'
        params['velocity'] =            'angstrom/fs'
        params['force'] =               'kcal/(mol*angstrom)'
        params['torque'] =              'kcal/mol'
        params['temperature'] =         'K'
        params['pressure'] =            'atm'
        params['dynamic viscosity'] =   'Pa*s/10'
        params['charge'] =              'e'
        params['dipole'] =              'e*angstrom'
        params['electric field'] =      'V/angstrom'
        params['density'] =             'g/cm^3'
    
    elif units == 'metal':
        params['mass'] =                'g/mol'
        params['length'] =              'angstrom'
        params['time'] =                'ps'
        params['energy'] =              'eV'
        params['velocity'] =            'angstrom/ps'
        params['force'] =               'eV/angstrom'
        params['torque'] =              'eV'
        params['temperature'] =         'K'
        params['pressure'] =            'bar'
        params['dynamic viscosity'] =   'Pa*s/10'
        params['charge'] =              'e'
        params['dipole'] =              'e*angstrom'
        params['electric field'] =      'V/angstrom'
        params['density'] =             'g/cm^3'
    
    elif units == 'si':
        params['mass'] =                'kg'
        params['length'] =              'm'
        params['time'] =                's'
        params['energy'] =              'J'
        params['velocity'] =            'm/s'
        params['force'] =               'N'
        params['torque'] =              'N*m'
        params['temperature'] =         'K'
        params['pressure'] =            'Pa'
        params['dynamic viscosity'] =   'Pa*s'
        params['charge'] =              'C'
        params['dipole'] =              'C*m'
        params['electric field'] =      'V/m'
        params['density'] =             'kg/m^3'
    
    elif units == 'cgs':
        params['mass'] =                'g'
        params['length'] =              'cm'
        params['time'] =                's'
        params['energy'] =              'erg'
        params['velocity'] =            'cm/s'
        params['force'] =               'dyn'
        params['torque'] =              'dyn*cm'
        params['temperature'] =         'K'
        params['pressure'] =            '0.1*Pa'
        params['dynamic viscosity'] =   'Pa*s/10'
        params['charge'] =              '10*c0*C'
        params['dipole'] =              '10*c0*C*cm'
        params['electric field'] =      'c0*uV/cm'
        params['density'] =             'g/cm^3'
    
    elif units == 'electron':
        params['mass'] =                'amu'
        params['length'] =              'aBohr'
        params['time'] =                'fs'
        params['energy'] =              '2*Ry'
        params['velocity'] =            '2*Ry*aBohr/hbar'
        params['force'] =               '2*Ry/aBohr'
        params['torque'] =              'dyn*cm'
        params['temperature'] =         'K'
        params['pressure'] =            'Pa'
        params['charge'] =              'e'
        params['dipole'] =              '1e-21/c0*C*m'
        params['electric field'] =      'V/cm'
    
    elif units == 'micro':
        params['mass'] =                'pg'
        params['length'] =              'um'
        params['time'] =                'us'
        params['energy'] =              'pg*um^2/us^2'
        params['velocity'] =            'um/us'
        params['force'] =               'pg*um/us^2'
        params['torque'] =              'pg*um^2/us^2'
        params['temperature'] =         'K'
        params['pressure'] =            'pg/(um*us^2)'
        params['dynamic viscosity'] =   'pg/(um*us)'
        params['charge'] =              '1e-12*C'
        params['dipole'] =              '1e-12*C*um'
        params['electric field'] =      'V/um'
        params['density'] =             'pg/um^3'
    
    elif units == 'nano':
        params['mass'] =                '1e-18*g'
        params['length'] =              'nm'
        params['time'] =                'ns'
        params['energy'] =              '1e-18*g*nm^2/ns^2'
        params['velocity'] =            'nm/ns'
        params['force'] =               '1e-18*g*nm/ns^2'
        params['torque'] =              '1e-18*g*nm^2/ns^2'
        params['temperature'] =         'K'
        params['pressure'] =            '1e-18*g/(nm*ns^2)'
        params['dynamic viscosity'] =   '1e-18*g/(nm*ns)'
        params['charge'] =              'e'
        params['dipole'] =              'e*nm'
        params['electric field'] =      'V/nm'
        params['density'] =             '1e-18*g/nm^3'
    
    else:
        raise ValueError('units ' + units + ' not supported')
    
    params['ang-mom'] = params['length'] + '*' + params['velocity'] + '*' + params['mass']
    params['ang-vel'] = '1/' + params['time']
    
    return params

def timestep(units='metal'):
    """
    Returns the default timestep value for a LAMMPS unit style option.
    
    Parameters
    ----------
    units : str, optional
        The LAMMPS unit option.  Default value is 'metal'
    
    Returns
    -------
    float
        The default timestep value.
    """
    if units == 'lj':
        return 0.005
    elif units == 'real':
        return 1.0
    elif units == 'metal':
        return 0.001
    elif units == 'si':
        return 1.0e-8
    elif units == 'cgs':
        return 1.0e-8
    elif units == 'electron':
        return 0.001
    elif units == 'micro':
        return 2.0
    elif units == 'nano':
        return 0.00045
    else:
        raise ValueError('units ' + units + ' not supported')