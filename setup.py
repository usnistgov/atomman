import os
import sys
from setuptools import setup, find_packages
from setuptools.extension import Extension

try:
    from Cython.Build import cythonize
except:
    USE_CYTHON = False
    EXT = '.c'
else:
    USE_CYTHON = True
    EXT = '.pyx'

extensions = [Extension("*", ["atomman/core/*" + EXT]),
              Extension("*", ["atomman/defect/*" + EXT])]

if USE_CYTHON:
    extensions = cythonize(extensions, language_level=sys.version_info[0], force=True)

def getversion():
    """Fetches version information from VERSION file"""
    with open(os.path.join('atomman', 'VERSION'), encoding='UTF-8') as version_file:
        version = version_file.read().strip()
    return version

def getreadme():
    """Fetches description from README.rst file"""
    with open('README.rst', encoding='UTF-8') as f:
        return f.read()

setup(
    name = 'atomman',
    version = getversion(),
    description = 'Atomistic Manipulation Toolkit',
    long_description = getreadme(),
    ext_modules = extensions,
    classifiers = [
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Science/Research',
        'Natural Language :: English',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
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
        'numericalunits',
        'numpy>=1.15',
        'matplotlib',
        'scipy',
        'pandas',
        'cython',
        'requests',
        'toolz',
        'potentials==0.3.4'
    ],
    include_package_data = True,
    zip_safe = False
)
