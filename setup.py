#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from setuptools import setup, Extension, find_packages
import pathlib
import glob

# Build Library
sofa_lib = Extension("pysofa2._sofa_c",
                       glob.glob('./src/*.c'),
                       depends=["./src/sofa.h", "./src/sofam.h"],
                       include_dirs=["./src"])


# Load description
about = {}
with open(pathlib.Path.cwd() / 'pysofa2' / '__about__.py') as fobj:
    exec(fobj.read(), about)

# Read Requirements
requirements = []
with open(pathlib.Path.cwd() / 'requirements.txt') as fp:
    for line in fp:
        requirements.append(line.strip())

# Setup information
setup(
    name='pysofa2',
    version=about['__version__'],
    packages = find_packages(),
    install_requires = requirements,
    ext_modules = [sofa_lib],
    include_package_data = True,
    author               = about['__author__'],
    author_email         = about['__author_email__'],
    maintainer           = about['__author__'],
    maintainer_email     = about['__author_email__'],
    description          = about['__description__'],
    long_description     = about['__long_description__'],
    long_description_content_type = about['__long_description_content_type__'],
    url                  = about['__url__']
)