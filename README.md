| Testing | Coverage | PyPi |
| :-----: | :------: | :--: |
| [![Build Status](https://travis-ci.org/duncaneddy/pysofa2.svg?branch=master)](https://travis-ci.org/duncaneddy/pysofa2) | [![Coverage Status](https://coveralls.io/repos/github/duncaneddy/pysofa2/badge.svg?branch=master)](https://coveralls.io/github/duncaneddy/pysofa2?branch=master) | [![PyPI version](https://badge.fury.io/py/pysofa2.svg)](https://badge.fury.io/py/pysofa2) |

# pysofa2

*pysofa2* is a module for using International Astronomical
 Union's (<http://www.iau.org/>) SOFA library (<http://www.iausofa.org/>) in 
python. It is inspired by original [pysofa](https://pypi.org/project/pysofa/)
implementation of Frederic Grollier. It extends the original pysofa module by
distributing and building the SOFA C library with the package, eliminiting the
need for the end-user to build the libraries separaetly from the package 
installtion.

Like *pysofa*, *pysofa2* is not a port of SOFA routines but a wrapper around the 
SOFA C library. Thus, no calculations are made in the pysofa software, this package
simply calls the underlying SOFA_C library.

## Documentation

For package documentation please refer to the underlying SOFA documentation at:
[http://www.iausofa.org/](http://www.iausofa.org/)

## Versioning

This package diverges from semantic versioning, and instead is versioned based on
the underlying SOFA release. SOFA versions based on YEAR-MONTH-DAY of the release.
For pysofa2

## Licensing

*pysofa2* not endorsed by the International Astronomical Union. 
In addition to *pysofa2*'s MIT license, any use of this module should comply 
with SOFA's license and [terms of use](http://www.iausofa.org/copyr.pdf). 
Especially, but not exclusively, any published work or commercial products 
which includes results achieved by using *pysofa* shall acknowledge that the 
SOFA software was used in obtaining those results.