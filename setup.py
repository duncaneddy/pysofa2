#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from setuptools import setup, Extension, find_packages
import glob
import pathlib

long_desc = """*pysofa2* is a module for using `International Astronomical
 Union <http://www.iau.org/>`_'s `SOFA library <http://www.iausofa.org/>`_ in 
python. It is inspired by original `pysofa <https://pypi.org/project/pysofa/>`_ 
implementation of Frederic Grollier. It extends the original pysofa module by
distributing and building the SOFA C library with the package, eliminiting the
need for the end-user to build the libraries separaetly from the package 
installtion.

Like *pysofa*, *pysofa2* is not a port of SOFA routines but a wrapper around the 
SOFA C library. Thus, no calculations are made into the pysofa software, they are
all delegated to the underlying SOFA_C library.

*pysofa2* not endorsed by the International Astronomical Union. 
In addition to *pysofa2*'s MIT license, any use of this module should comply 
with `SOFA's license and terms of use <http://www.iausofa.org/copyr.pdf>`_. 
Especially, but not exclusively, any published work or commercial products 
which includes results achieved by using *pysofa* shall acknowledge that the 
SOFA software was used in obtaining those results."""

sofa_lib = Extension("pysofa2._sofa_c",
                       glob.glob('./src/*.c'),
                       depends=["./src/sofa.h", "./src/sofam.h"],
                       include_dirs=["./src"])

# Setup information
setup(
    name='pysofa2',
    version='18.01.30.0',
    packages = find_packages(),
    install_requires = [
        'numpy>=1.14.0'
    ],
    ext_modules = [sofa_lib],
    include_package_data = True,
    author               = 'Duncan Eddy',
    author_email         = 'duncan.eddy@gmail.com',
    maintainer           = 'Duncan Eddy',
    maintainer_email     = 'duncan.eddy@gmail.com',
    description          = 'A wrapper of the International Astronomical Union\'s SOFA lbrary',
    long_description     = long_desc,
    url                  = ''
)