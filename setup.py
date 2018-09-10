#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from setuptools import setup, Extension, find_packages
import glob
import pathlib

# Load description
with open('README.md') as f:
    long_desc = f.read()

sofa_lib = Extension("pysofa2._sofa_c",
                       glob.glob('./src/*.c'),
                       depends=["./src/sofa.h", "./src/sofam.h"],
                       include_dirs=["./src"])

# Setup information
setup(
    name='pysofa2',
    version='18.01.30.4',
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
    long_description_content_type='text/markdown',
    url                  = 'https://gitlab.com/deddy/pysofa2'
)