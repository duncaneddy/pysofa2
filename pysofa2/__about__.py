__version__      = "23.10.11.1"
__author__       = "Duncan Eddy"
__author_email__ = "duncan.eddy@gmail.com"
__description__  = "A wrapper of the International Astronomical Union\'s SOFA lbrary."
__url__          = "https://gitlab.com/deddy/pysofa2"
__long_description_content_type__ = 'text/markdown'
__long_description__ = '''*pysofa2* is a module for using International Astronomical
 Union's (<http://www.iau.org/>) SOFA library (<http://www.iausofa.org/>) in 
python. It is inspired by original [pysofa](https://pypi.org/project/pysofa/)
implementation of Frederic Grollier. It extends the original pysofa module by
distributing and building the SOFA C library with the package, eliminiting the
need for the end-user to build the libraries separaetly from the package 
installtion.

Like *pysofa*, *pysofa2* is not a port of SOFA routines but a wrapper around the 
SOFA C library. Thus, no calculations are made into the pysofa software, they are
all delegated to the underlying SOFA_C library.

*pysofa2* not endorsed by the International Astronomical Union. 
In addition to *pysofa2*'s MIT license, any use of this module should comply 
with SOFA's license and [terms of use](http://www.iausofa.org/copyr.pdf). 
Especially, but not exclusively, any published work or commercial products 
which includes results achieved by using *pysofa* shall acknowledge that the 
SOFA software was used in obtaining those results.'''
