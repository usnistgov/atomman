import os
from setuptools import setup, find_packages
from setuptools.extension import Extension

try:
    import numpy
except:
    raise ImportError('Install numpy first!')

try:
    from Cython.Build import cythonize
except:
    # List cython extensions
    extensions = [Extension("atomman.core.cythonized",
                           ["atomman/core/cythonized.c"],
                           include_dirs=[numpy.get_include()])]
else:
    # List cython extensions
    extensions = [Extension("atomman.core.cythonized",
                           ["atomman/core/cythonized.pyx"],
                           include_dirs=[numpy.get_include()])]
    extensions = cythonize(extensions)

def getversion():
    """Fetches version information from VERSION file"""
    with open(os.path.join('atomman', 'VERSION')) as version_file:
        version = version_file.read().strip()
    return version

def getreadme():
    """Fetches description from README.rst file"""
    with open('README.rst') as f:
        return f.read()

setup(name = 'atomman',
      version = getversion(),
      description = 'Atomistic Manipulation Toolkit',
      long_description = getreadme(),
      ext_modules = extensions,
      classifiers = [
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'Natural Language :: English',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3.6',
        'Topic :: Scientific/Engineering :: Physics'
      ],
      keywords = [
        'atom', 
        'atomic', 
        'atomistic', 
        'molecular dynamics'
      ], 
      url = 'https://github.com/usnistgov/atomman/',
      author = 'Lucas Hale',
      author_email = 'lucas.hale@nist.gov',
      packages = find_packages(),
      install_requires = [
        'xmltodict',
        'DataModelDict',
        'numpy', 
        'matplotlib',
        'scipy',
        'pandas',
        'cython',
        'numericalunits'
      ],
      include_package_data = True,
      zip_safe = False)