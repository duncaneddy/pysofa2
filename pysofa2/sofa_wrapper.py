# -*- coding: utf-8 -*-

import warnings as _warnings
import os as _os
import numpy as _np
import pathlib as _pathlib
import glob as _glob

# Import CTypes for interfacing
import ctypes as _ct
# from   ctypes.util import find_library as _find_library
from   numpy.ctypeslib import ndpointer as _ndpointer


# Attempt to find sofa library to load
_sofalib_filename = _glob.glob(str(_pathlib.Path(__file__).parent) + '/_sofa_c*')

if len(_sofalib_filename) == 0:
    raise ImportError('Unable to find the shared C library "pysofa2._sofa_c".')

_sofa = _ct.CDLL(_sofalib_filename[0])

# Require C Shape
def _req_shape_c(a, dtype=None, shape=None, req=None):
    return _np.require(a, dtype=dtype, requirements=req).reshape(shape, order='C')


##########################
# SOFA Function Wrappers #
##########################

# These functions are based of of the headers in sofa.h

# /* Astronomy/Calendars */
_sofa.iauCal2jd.argtypes = [_ct.c_int, #iy
                            _ct.c_int, #im
                            _ct.c_int, #id
                            _ct.POINTER(_ct.c_double), #djm0
                            _ct.POINTER(_ct.c_double)] #djm
_sofa.iauCal2jd.restype = _ct.c_int
_cal2jd_msg = {
            -1: 'minimum year allowed is -4799',
            -2: 'month must be in 1..12',
            -3: 'day is out of range for this month'}

def Cal2jd(iy, im, id):
    '''Convert Gregorian Calendar to Julian Date

    Args:
        iy (int): year in Gregorian calendar
        im (int): month in Gregorian calendar
        id (int): day in Gregorian calendar

    Returns:
        djm0 (double) MJD zero-point: always 2400000.5
        djm (double) Modified Julian Date for 0 hrs.

    SOFA Notes:

        1. The algorithm used is valid from -4800 March 1, but the implementation
        rejects dates before -4799 January 1.

        2. The Julian Date is returned in two pices, in the usual SOFA manner,
        which is designed to preserved time resolution. The Julian Date is 
        available as a single number by adding djm0 and djm.

        3. In early earas the conversion is from the "Proleptic Gregorian
        Calendar"; no account is taken of the dates(s) of adoption of the
        Gregorian Calendar, nor is the AD/BC numbering convention observed.
    '''

    # Initialize returned values
    djm0 = _ct.c_double()
    djm  = _ct.c_double()

    s = _sofa.iauCal2jd(iy, im, id, _ct.byref(djm0), _ct.byref(djm))
    
    # Raise status code error
    if s != 0:
        raise ValueError(_cal2jd_msg[s])
    
    return djm0.value, djm.value

_sofa.iauEpb.argtypes = [_ct.c_double, #dj1
                         _ct.c_double] #dj2
_sofa.iauEpb.restype = _ct.c_double

def Epb(dj1, dj2):
    '''Convert Julian Date to Bessellian Epoch

    Args:
        dj1 (double): Julian Date (see note)
        dj2 (double): Julian Date (see note)

    Returns:
        epb (double): Besselian Epoch

    SOFA Notes:
        The Julian Date is supplied in two pieces, in the usual SOFA manner,
        which is designed to preserve time resolution. The Julian Date is 
        available as a single number by adding dj1 and dj2. The maximum
        resolution is achieved if dj1 is 2451545.0 (J2000.0).
    '''

    epb = _sofa.iauEpb(dj1, dj2)

    return epb

_sofa.iauEpb2jd.argtypes = [_ct.c_double, # epb
                            _ct.POINTER(_ct.c_double), # djm0
                            _ct.POINTER(_ct.c_double)] # djm
_sofa.iauEpb2jd.restype = None

def Epb2jd(epb):
    '''Convert Besselian Epoch to Julian Date

    Args:
        epb (double): Besselian Epoch (e.g. 1957.3)

    Returns:
        djm0 (double): MJD zero-point: always 2400000.5
        djm (double): Modified Julian Date

    SOFA Notes:
        The Julian Date is returned in two pieces, in the usual SOFA manner, 
        which is designed to preserve time resolution. The Julian Date is 
        available as a signle number by adding djm0 and djm.
    '''

    # Initialize pointer returns
    djm0 = _ct.c_double()
    djm  = _ct.c_double()

    # Main function call
    _sofa.iauEpb2jd(epb, _ct.byref(djm0), _ct.byref(djm))

    return djm0.value, djm.value


_sofa.iauEpj.argtypes = [
    _ct.c_double, # dj1
    _ct.c_double  # dj2
]
_sofa.iauEpj.restype = _ct.c_double

def Epj(dj1, dj2):
    '''Julian Date to Julian Epoch

    Args:
        dj1 (double): Julian Date (see note)
        dj2 (double): Julian Date (see note)

    Returns
        epj (double): Julian Epoch

    SOFA Notes:
        The Julian Date is supplied in two pieces, in the usual SOFA manner, 
        which is designed to preserve time resolution.  The Julian Date is 
        available as a single number by adding dj1 and dj2.  The maximum 
        resolution is achieved if dj1 is 2451545.0 (J2000.0).
    '''

    # Main funciton call
    epj = _sofa.iauEpj(dj1, dj2)

    return epj

_sofa.iauEpj2jd.argtypes = [
    _ct.c_double, # epj
    _ct.POINTER(_ct.c_double), # djm0
    _ct.POINTER(_ct.c_double)  # djm
]
_sofa.iauEpj2jd.restype = None

def Epj2jd(epj):
    '''Julian Epoch to Julian Date

    Args:
        epj (double): Julian Epoch (e.g. 1996.8)

    Returns:
        djm0 (double): MJD zero-point always 2400000.5
        djm (double): Modified Julian Date

    Note:
        The Julian Date is returned in two pieces, in the usual SOFA manner, 
        which is designed to preserve time resolution. The Julian Date is 
        available as a single number by adding djm0 and djm.
    '''

    # Intialize pointer variables
    djm0 = _ct.c_double()
    djm  = _ct.c_double()

    # Main function call
    _sofa.iauEpj2jd(epj, _ct.byref(djm0), _ct.byref(djm))

    return djm0.value, djm.value

_sofa.iauJd2cal.argtypes = [
    _ct.c_double, # dj1
    _ct.c_double, # dj2
    _ct.POINTER(_ct.c_int), # iy
    _ct.POINTER(_ct.c_int), # im
    _ct.POINTER(_ct.c_int), # id
    _ct.POINTER(_ct.c_double)  # fd
]
_sofa.iauJd2cal.restype = _ct.c_int
_jd2cal_msg = {
    -1: 'invalid Julian Date. Valid range -68569.5 to 1e9'
}

def Jd2cal(dj1, dj2):
    '''Julian Date to Gregorian year, month, day, and fraction of a day

    Args:
        dj1 (double): Julian Date (see notes 1, 2)
        dj2 (double): Julian Date (see notes 1, 2)

    Returns:
        iy (int): Year
        id (int): Month
        id (int): Day
        fd (double): fraction of day

    Notes:

        1. The earliest valid date is -68569.5 (-4900 March 1). The largest
        value accepted is 1e9.

        2. The Julian Date is approtioned in any convienient way between the
        arguments dj1 and dj2. For example, JD=2450123.7 could be expressed in
        any of these ways, among others:

           dj1             dj2

        2450123.7           0.0       (JD method)
        2451545.0       -1421.3       (J2000 method)
        2400000.5       50123.2       (MJD method)
        2450123.5           0.2       (date & time method)

        3. In the early eras the conversion is from the "proleptic Gregorian
        calendar"; no account is taken of the date(s) of adoption of the
        Gregorian calendar, or is the AD/BC numbering convention observed.
    '''

    # Initialized returned values
    iy = _ct.c_int()
    im = _ct.c_int()
    id = _ct.c_int()
    fd = _ct.c_double()

    # Main function call
    status = _sofa.iauJd2cal(dj1, dj2, _ct.byref(iy),  _ct.byref(im),
                                 _ct.byref(id),  _ct.byref(fd))

    # Raise any error code warnings
    if status != 0:
        raise ValueError(_jd2cal_msg[status])

    return iy.value, im.value, id.value, fd.value

_sofa.iauJdcalf.argtypes = [
    _ct.c_int, # ndp
    _ct.c_double, # dj1
    _ct.c_double, # dj2
    _ct.c_int * 4 # iymdf[4]
]
_sofa.iauJdcalf.restype = _ct.c_int
_jd2calf_msg = {
    -1: 'invalid Julian Date. Valid range -68569.5 to 1e9',
    +1: 'Number of decimal places not 0-9 (interpreted as 0)'
}

def Jd2calf(ndp, dj1, dj2):
    '''Julian Date to Gregorian Calendar conversion rounded to a 
    specified precision.

    Args:
        ndp (int): Number of decimal places of days in fraction
        dj1 (double): Julian date (see note 1)
        dj2 (double): Julian date (see note 1)

    Returns:
        iy (int): Year
        id (int): Month
        id (int): Day
        fd (double): fraction of day

    Notes:

        1. The Julian Date is approtioned in any convienient way between the
        arguments dj1 and dj2. For example, JD=2450123.7 could be expressed in
        any of these ways, among others:

           dj1             dj2

        2450123.7           0.0       (JD method)
        2451545.0       -1421.3       (J2000 method)
        2400000.5       50123.2       (MJD method)
        2450123.5           0.2       (date & time method)

        2. In the early eras the conversion is from the "proleptic Gregorian
        calendar"; no account is taken of the date(s) of adoption of the
        Gregorian calendar, or is the AD/BC numbering convention observed.

        3. Refer to the function Jd2cal

        4. ndp should be 4 or less if internal overflows are to be avoided on
        machines which use 16-bit integers.
    '''

    # Initialize internal variables
    iymdf = (_ct.c_int * 4)()

    # Main funciton call
    status = _sofa.iauJdcalf(ndp, dj1, dj2, iymdf)

    # Print status warning
    if status < 0:
        raise ValueError(_jd2calf_msg[status])
    elif status > 0:
        raise _warnings.warn(_jd2calf_msg[status], UserWarning)

    # Return output
    return iymdf[0], iymdf[1], iymdf[2], iymdf[3]

# /* Astronomy/Astrometry */
# void iauAb(double pnat[3], double v[3], double s, double bm1,
#            double ppr[3]);
# void iauApcg(double date1, double date2,
#              double ebpv[2][3], double ehp[3],
#              iauASTROM *astrom);
# void iauApcg13(double date1, double date2, iauASTROM *astrom);
# void iauApci(double date1, double date2,
#              double ebpv[2][3], double ehp[3],
#              double x, double y, double s,
#              iauASTROM *astrom);
# void iauApci13(double date1, double date2,
#                iauASTROM *astrom, double *eo);
# void iauApco(double date1, double date2,
#              double ebpv[2][3], double ehp[3],
#              double x, double y, double s, double theta,
#              double elong, double phi, double hm,
#              double xp, double yp, double sp,
#              double refa, double refb,
#              iauASTROM *astrom);
# int iauApco13(double utc1, double utc2, double dut1,
#               double elong, double phi, double hm, double xp, double yp,
#               double phpa, double tc, double rh, double wl,
#               iauASTROM *astrom, double *eo);
# void iauApcs(double date1, double date2, double pv[2][3],
#              double ebpv[2][3], double ehp[3],
#              iauASTROM *astrom);
# void iauApcs13(double date1, double date2, double pv[2][3],
#                iauASTROM *astrom);
# void iauAper(double theta, iauASTROM *astrom);
# void iauAper13(double ut11, double ut12, iauASTROM *astrom);
# void iauApio(double sp, double theta,
#              double elong, double phi, double hm, double xp, double yp,
#              double refa, double refb,
#              iauASTROM *astrom);
# int iauApio13(double utc1, double utc2, double dut1,
#               double elong, double phi, double hm, double xp, double yp,
#               double phpa, double tc, double rh, double wl,
#               iauASTROM *astrom);
# void iauAtci13(double rc, double dc,
#                double pr, double pd, double px, double rv,
#                double date1, double date2,
#                double *ri, double *di, double *eo);
# void iauAtciq(double rc, double dc, double pr, double pd,
#               double px, double rv, iauASTROM *astrom,
#               double *ri, double *di);
# void iauAtciqn(double rc, double dc, double pr, double pd,
#                double px, double rv, iauASTROM *astrom,
#                int n, iauLDBODY b[], double *ri, double *di);
# void iauAtciqz(double rc, double dc, iauASTROM *astrom,
#                double *ri, double *di);
# int iauAtco13(double rc, double dc,
#               double pr, double pd, double px, double rv,
#               double utc1, double utc2, double dut1,
#               double elong, double phi, double hm, double xp, double yp,
#               double phpa, double tc, double rh, double wl,
#               double *aob, double *zob, double *hob,
#               double *dob, double *rob, double *eo);
# void iauAtic13(double ri, double di,
#                double date1, double date2,
#                double *rc, double *dc, double *eo);
# void iauAticq(double ri, double di, iauASTROM *astrom,
#               double *rc, double *dc);
# void iauAticqn(double ri, double di, iauASTROM *astrom,
#                int n, iauLDBODY b[], double *rc, double *dc);
# int iauAtio13(double ri, double di,
#               double utc1, double utc2, double dut1,
#               double elong, double phi, double hm, double xp, double yp,
#               double phpa, double tc, double rh, double wl,
#               double *aob, double *zob, double *hob,
#               double *dob, double *rob);
# void iauAtioq(double ri, double di, iauASTROM *astrom,
#               double *aob, double *zob,
#               double *hob, double *dob, double *rob);
# int iauAtoc13(const char *type, double ob1, double ob2,
#               double utc1, double utc2, double dut1,
#               double elong, double phi, double hm, double xp, double yp,
#               double phpa, double tc, double rh, double wl,
#               double *rc, double *dc);
# int iauAtoi13(const char *type, double ob1, double ob2,
#               double utc1, double utc2, double dut1,
#               double elong, double phi, double hm, double xp, double yp,
#               double phpa, double tc, double rh, double wl,
#               double *ri, double *di);
# void iauAtoiq(const char *type,
#               double ob1, double ob2, iauASTROM *astrom,
#               double *ri, double *di);
# void iauLd(double bm, double p[3], double q[3], double e[3],
#            double em, double dlim, double p1[3]);
# void iauLdn(int n, iauLDBODY b[], double ob[3], double sc[3],
#             double sn[3]);
# void iauLdsun(double p[3], double e[3], double em, double p1[3]);
# void iauPmpx(double rc, double dc, double pr, double pd,
#              double px, double rv, double pmt, double pob[3],
#              double pco[3]);
# int iauPmsafe(double ra1, double dec1, double pmr1, double pmd1,
#               double px1, double rv1,
#               double ep1a, double ep1b, double ep2a, double ep2b,
#               double *ra2, double *dec2, double *pmr2, double *pmd2,
#               double *px2, double *rv2);
# void iauPvtob(double elong, double phi, double height, double xp,
#               double yp, double sp, double theta, double pv[2][3]);
# void iauRefco(double phpa, double tc, double rh, double wl,
#               double *refa, double *refb);

# /* Astronomy/Ephemerides */
# int iauEpv00(double date1, double date2,
#              double pvh[2][3], double pvb[2][3]);
# int iauPlan94(double date1, double date2, int np, double pv[2][3]);

# /* Astronomy/FundamentalArgs */

_sofa.iauFad03.argtypes = [
    _ct.c_double # t
]
_sofa.iauFad03.restype = _ct.c_double

def Fad03(t):
    '''Fundamental argument, mean elongation of the Moon from the Sun

    Args:
        t (double): TDB, Julian centuries since J2000.0 (Note 1)

    Returns:
        D (double): radians (Note 2)

    Notes:

        1. Though t is trictly TDB, it is usually more convenient to use TT,
        which makes no significant difference

        2. The experession used is as adopted in IERS Conventions (2003) and is
        from Simon et al. (1994).
    '''

    return _sofa.iauFad03(t)

_sofa.iauFae03.argtypes = [
    _ct.c_double # t
]
_sofa.iauFae03.restype = _ct.c_double

def Fae03(t):
    '''Fundamental argument, mean longitude of Earth

    Args:
        t (double): TDB, Julian centuries since J2000.0 (Note 1)

    Returns:
        x (double): mean longitude of Earth, radians (Note 2)

    Notes:

        1. Though t is trictly TDB, it is usually more convenient to use TT,
        which makes no significant difference

        2. The experession used is as adopted in IERS Conventions (2003) and 
        comes from Souchay et al. (1999) after Simon et al. (1994).
    '''

    return _sofa.iauFae03(t)


_sofa.iauFaf03.argtypes = [
    _ct.c_double # t
]
_sofa.iauFaf03.restype = _ct.c_double

def Faf03(t):
    '''Fundamental argument, mean longitude of the Moon minus the mean longitude
    of the ascending node

    Args:
        t (double): TDB, Julian centuries since J2000.0 (Note 1)

    Returns:
        F (double): radians (Note 2)

    Notes:

        1. Though t is trictly TDB, it is usually more convenient to use TT,
        which makes no significant difference

        2. The experession used is as adopted in IERS Conventions (2003) and 
        is from Simon et al. (1994).
    '''

    return _sofa.iauFaf03(t)

_sofa.iauFaju03.argtypes = [
    _ct.c_double # t
]
_sofa.iauFaju03.restype = _ct.c_double

def Faju03(t):
    '''Fundamental argument, mean longitude of Jupiter

    Args:
        t (double): TDB, Julian centuries since J2000.0 (Note 1)

    Returns:
        x (double): mean longitude of Jupiter, radians (Note 2)

    Notes:

        1. Though t is trictly TDB, it is usually more convenient to use TT,
        which makes no significant difference

        2. The experession used is as adopted in IERS Conventions (2003) and 
        comes from Souchay et al. (1999) after Simon et al. (1994).
    '''

    return _sofa.iauFaju03(t)

_sofa.iauFal03.argtypes = [
    _ct.c_double # t
]
_sofa.iauFal03.restype = _ct.c_double

def Fal03(t):
    '''Fundamental argument, mean anomaly of the Moon

    Args:
        t (double): TDB, Julian centuries since J2000.0 (Note 1)

    Returns:
        l (double): radians (Note 2)

    Notes:

        1. Though t is trictly TDB, it is usually more convenient to use TT,
        which makes no significant difference

        2. The experession used is as adopted in IERS Conventions (2003) and 
        is from Simon et al. (1994).
    '''

    return _sofa.iauFal03(t)

_sofa.iauFalp03.argtypes = [
    _ct.c_double # t
]
_sofa.iauFalp03.restype = _ct.c_double

def Falp03(t):
    '''Fundamental argument, mean anomaly of the Sun

    Args:
        t (double): TDB, Julian centuries since J2000.0 (Note 1)

    Returns:
        l' (double): radians (Note 2)

    Notes:

        1. Though t is trictly TDB, it is usually more convenient to use TT,
        which makes no significant difference

        2. The experession used is as adopted in IERS Conventions (2003) and 
        is from Simon et al. (1994).
    '''

    return _sofa.iauFalp03(t)

_sofa.iauFama03.argtypes = [
    _ct.c_double # t
]
_sofa.iauFama03.restype = _ct.c_double

def Fama03(t):
    '''Fundamental argument, mean longitude of Mars

    Args:
        t (double): TDB, Julian centuries since J2000.0 (Note 1)

    Returns:
        x (double): mean longitude of Mars, radians (Note 2)

    Notes:

        1. Though t is trictly TDB, it is usually more convenient to use TT,
        which makes no significant difference

        2. The experession used is as adopted in IERS Conventions (2003) and 
        comes from Souchay et al. (1999) after Simon et al. (1994).
    '''

    return _sofa.iauFama03(t)

_sofa.iauFame03.argtypes = [
    _ct.c_double # t
]
_sofa.iauFame03.restype = _ct.c_double

def Fame03(t):
    '''Fundamental argument, mean longitude of Mercury

    Args:
        t (double): TDB, Julian centuries since J2000.0 (Note 1)

    Returns:
        x (double): mean longitude of Mercury, radians (Note 2)

    Notes:

        1. Though t is trictly TDB, it is usually more convenient to use TT,
        which makes no significant difference

        2. The experession used is as adopted in IERS Conventions (2003) and 
        comes from Souchay et al. (1999) after Simon et al. (1994).
    '''

    return _sofa.iauFame03(t)

_sofa.iauFane03.argtypes = [
    _ct.c_double # t
]
_sofa.iauFane03.restype = _ct.c_double

def Fane03(t):
    '''Fundamental argument, mean longitude of Neptune

    Args:
        t (double): TDB, Julian centuries since J2000.0 (Note 1)

    Returns:
        x (double): mean longitude of Neptune, radians (Note 2)

    Notes:

        1. Though t is trictly TDB, it is usually more convenient to use TT,
        which makes no significant difference

        2. The experession used is as adopted in IERS Conventions (2003) and 
        is adapted from Simon et al. (1994).
    '''

    return _sofa.iauFane03(t)

_sofa.iauFaom03.argtypes = [
    _ct.c_double # t
]
_sofa.iauFaom03.restype = _ct.c_double

def Faom03(t):
    '''Fundamental argument, mean longitude of of the Moon's ascending node

    Args:
        t (double): TDB, Julian centuries since J2000.0 (Note 1)

    Returns:
        Omega (double): radians (Note 2)

    Notes:

        1. Though t is trictly TDB, it is usually more convenient to use TT,
        which makes no significant difference

        2. The experession used is as adopted in IERS Conventions (2003) and 
        is adapted from Simon et al. (1994).
    '''

    return _sofa.iauFaom03(t)

_sofa.iauFapa03.argtypes = [
    _ct.c_double # t
]
_sofa.iauFapa03.restype = _ct.c_double

def Fapa03(t):
    '''Fundamental argument, general accumulated precession in longitude

    Args:
        t (double): TDB, Julian centuries since J2000.0 (Note 1)

    Returns:
        x (double): general precession in longitude, radians (Note 2)

    Notes:

        1. Though t is trictly TDB, it is usually more convenient to use TT,
        which makes no significant difference

        2. The experession used is as adopted in IERS Conventions (2003). 
        It is taken from Kinoshita & Souchay (1990) and comes originally from 
        Lieske et al. (1977).
    '''

    return _sofa.iauFapa03(t)

_sofa.iauFasa03.argtypes = [
    _ct.c_double # t
]
_sofa.iauFasa03.restype = _ct.c_double

def Fasa03(t):
    '''Fundamental argument, mean longitude of Saturn

    Args:
        t (double): TDB, Julian centuries since J2000.0 (Note 1)

    Returns:
        x (double): mean longitude of Saturn, radians (Note 2)

    Notes:

        1. Though t is trictly TDB, it is usually more convenient to use TT,
        which makes no significant difference

        2. The experession used is as adopted in IERS Conventions (2003) and 
        comes from Souchay et al. (1999) after Simon et al. (1994).
    '''

    return _sofa.iauFasa03(t)

_sofa.iauFaur03.argtypes = [
    _ct.c_double # t
]
_sofa.iauFaur03.restype = _ct.c_double

def Faur03(t):
    '''Fundamental argument, mean longitude of Uranus

    Args:
        t (double): TDB, Julian centuries since J2000.0 (Note 1)

    Returns:
        x (double): mean longitude of Uranus, radians (Note 2)

    Notes:

        1. Though t is trictly TDB, it is usually more convenient to use TT,
        which makes no significant difference

        2. The experession used is as adopted in IERS Conventions (2003) and 
        is adapted from Simon et al. (1994).
    '''

    return _sofa.iauFaur03(t)


_sofa.iauFave03.argtypes = [
    _ct.c_double # t
]
_sofa.iauFave03.restype = _ct.c_double

def Fave03(t):
    '''Fundamental argument, mean longitude of Venus

    Args:
        t (double): TDB, Julian centuries since J2000.0 (Note 1)

    Returns:
        x (double): mean longitude of Venus, radians (Note 2)

    Notes:

        1. Though t is trictly TDB, it is usually more convenient to use TT,
        which makes no significant difference

        2. The experession used is as adopted in IERS Conventions (2003) and 
        comes from Souchay et al. (1999) after Simon et al. (1994).
    '''

    return _sofa.iauFave03(t)

# /* Astronomy/PrecNutPolar */
# void iauBi00(double *dpsibi, double *depsbi, double *dra);
# void iauBp00(double date1, double date2,
#              double rb[3][3], double rp[3][3], double rbp[3][3]);
# void iauBp06(double date1, double date2,
#              double rb[3][3], double rp[3][3], double rbp[3][3]);
# void iauBpn2xy(double rbpn[3][3], double *x, double *y);
# void iauC2i00a(double date1, double date2, double rc2i[3][3]);
# void iauC2i00b(double date1, double date2, double rc2i[3][3]);
# void iauC2i06a(double date1, double date2, double rc2i[3][3]);
# void iauC2ibpn(double date1, double date2, double rbpn[3][3],
#                double rc2i[3][3]);
# void iauC2ixy(double date1, double date2, double x, double y,
#               double rc2i[3][3]);

_sofa.iauC2ixys.argtypes = [
    _ct.c_double, #x
    _ct.c_double, #y
    _ct.c_double, #s
    _ndpointer(shape=(3,3), dtype=float, flags='C') #rc2i
] 

def C2ixys(x, y, s):
    '''Form the celestial to intermediate-frame-of-date matrix given the CIP
    X,Y and the CIO locator s.

    Args:
        x,y (double): Celestial Intermediate Pole (Note 1)
        s (double): the CIO locator s (Note 2)

    Returns:
        rc2i (:obj:`np.ndarray`): celestial-to-intermediate matrix (Note 3)

    Notes:
        1. The Celestial Intermediate Pole coordinates are the x,y
        components of the unit vector in the Geocentric Celestial
        Reference System.
   
        2. The CIO locator s (in radians) positions the Celestial
        Intermediate Origin on the equator of the CIP.
   
        3. The matrix rc2i is the first stage in the transformation from
        celestial to terrestrial coordinates:
   
           [TRS] = RPOM * R_3(ERA) * rc2i * [CRS]
   
                 = RC2T * [CRS]
   
        where [CRS] is a vector in the Geocentric Celestial Reference
        System and [TRS] is a vector in the International Terrestrial
        Reference System (see IERS Conventions 2003), ERA is the Earth
        Rotation Angle and RPOM is the polar motion matrix.
    '''

    # Initialize return value
    rc2i = _np.zeros(shape=(3,3), dtype=float, order='C')

    # Main function call
    _sofa.iauC2ixys(x, y, s, rc2i)

    return _np.array(rc2i, dtype=float)

# void iauC2t00a(double tta, double ttb, double uta, double utb,
#                double xp, double yp, double rc2t[3][3]);
# void iauC2t00b(double tta, double ttb, double uta, double utb,
#                double xp, double yp, double rc2t[3][3]);
# void iauC2t06a(double tta, double ttb, double uta, double utb,
#                double xp, double yp, double rc2t[3][3]);
# void iauC2tcio(double rc2i[3][3], double era, double rpom[3][3],
#                double rc2t[3][3]);
# void iauC2teqx(double rbpn[3][3], double gst, double rpom[3][3],
#                double rc2t[3][3]);
# void iauC2tpe(double tta, double ttb, double uta, double utb,
#               double dpsi, double deps, double xp, double yp,
#               double rc2t[3][3]);
# void iauC2txy(double tta, double ttb, double uta, double utb,
#               double x, double y, double xp, double yp,
#               double rc2t[3][3]);
# double iauEo06a(double date1, double date2);
# double iauEors(double rnpb[3][3], double s);
# void iauFw2m(double gamb, double phib, double psi, double eps,
#              double r[3][3]);
# void iauFw2xy(double gamb, double phib, double psi, double eps,
#               double *x, double *y);
# void iauLtp(double epj, double rp[3][3]);
# void iauLtpb(double epj, double rpb[3][3]);
# void iauLtpecl(double epj, double vec[3]);
# void iauLtpequ(double epj, double veq[3]);
# void iauNum00a(double date1, double date2, double rmatn[3][3]);
# void iauNum00b(double date1, double date2, double rmatn[3][3]);
# void iauNum06a(double date1, double date2, double rmatn[3][3]);
# void iauNumat(double epsa, double dpsi, double deps, double rmatn[3][3]);
# void iauNut00a(double date1, double date2, double *dpsi, double *deps);
# void iauNut00b(double date1, double date2, double *dpsi, double *deps);
# void iauNut06a(double date1, double date2, double *dpsi, double *deps);
# void iauNut80(double date1, double date2, double *dpsi, double *deps);
# void iauNutm80(double date1, double date2, double rmatn[3][3]);
# double iauObl06(double date1, double date2);
# double iauObl80(double date1, double date2);
# void iauP06e(double date1, double date2,
#              double *eps0, double *psia, double *oma, double *bpa,
#              double *bqa, double *pia, double *bpia,
#              double *epsa, double *chia, double *za, double *zetaa,
#              double *thetaa, double *pa,
#              double *gam, double *phi, double *psi);
# void iauPb06(double date1, double date2,
#              double *bzeta, double *bz, double *btheta);
# void iauPfw06(double date1, double date2,
#               double *gamb, double *phib, double *psib, double *epsa);
# void iauPmat00(double date1, double date2, double rbp[3][3]);
# void iauPmat06(double date1, double date2, double rbp[3][3]);
# void iauPmat76(double date1, double date2, double rmatp[3][3]);
# void iauPn00(double date1, double date2, double dpsi, double deps,
#              double *epsa,
#              double rb[3][3], double rp[3][3], double rbp[3][3],
#              double rn[3][3], double rbpn[3][3]);
# void iauPn00a(double date1, double date2,
#               double *dpsi, double *deps, double *epsa,
#               double rb[3][3], double rp[3][3], double rbp[3][3],
#               double rn[3][3], double rbpn[3][3]);
# void iauPn00b(double date1, double date2,
#               double *dpsi, double *deps, double *epsa,
#               double rb[3][3], double rp[3][3], double rbp[3][3],
#               double rn[3][3], double rbpn[3][3]);
# void iauPn06(double date1, double date2, double dpsi, double deps,
#              double *epsa,
#              double rb[3][3], double rp[3][3], double rbp[3][3],
#              double rn[3][3], double rbpn[3][3]);
# void iauPn06a(double date1, double date2,
#               double *dpsi, double *deps, double *epsa,
#               double rb[3][3], double rp[3][3], double rbp[3][3],
#               double rn[3][3], double rbpn[3][3]);
# void iauPnm00a(double date1, double date2, double rbpn[3][3]);
# void iauPnm00b(double date1, double date2, double rbpn[3][3]);
# void iauPnm06a(double date1, double date2, double rnpb[3][3]);
# void iauPnm80(double date1, double date2, double rmatpn[3][3]);

_sofa.iauPom00.argtypes = [
    _ct.c_double, #xp
    _ct.c_double, #yp
    _ct.c_double, #sp
    _ndpointer(shape=(3,3), dtype=float, flags='C') #rpom
]

def Pom00(xp, yp, sp):
    '''Form the matrix of polar motion for a given date, IAU 2000.

    Args:
        xp,yp (double): coordinates of the pole (radians, Note 1)
        sp (double): the TIO locator s' (radians, Note 2)

    Returns:
        rpom (:obj:np.matrix): polar-motion matrix (Note 3)

    Notes:
        1. The arguments xp and yp are the coordinates (in radians) of the
        Celestial Intermediate Pole with respect to the International
        Terrestrial Reference System (see IERS Conventions 2003),
        measured along the meridians to 0 and 90 deg west respectively.
   
        2. The argument sp is the TIO locator s', in radians, which
        positions the Terrestrial Intermediate Origin on the equator.  It
        is obtained from polar motion observations by numerical
        integration, and so is in essence unpredictable.  However, it is
        dominated by a secular drift of about 47 microarcseconds per
        century, and so can be taken into account by using s' = -47*t,
        where t is centuries since J2000.0.  The function iauSp00
        implements this approximation.
   
        3. The matrix operates in the sense V(TRS) = rpom * V(CIP), meaning
        that it is the final rotation when computing the pointing
        direction to a celestial source.        
    '''

    # Initialize return matrix
    rpom = _np.zeros(shape=(3,3), dtype=float, order='C')

    # Main funciton call
    _sofa.iauPom00(float(xp), float(yp), float(sp), rpom)

    return _np.array(rpom, dtype=float)

# void iauPr00(double date1, double date2,
#              double *dpsipr, double *depspr);
# void iauPrec76(double date01, double date02,
#                double date11, double date12,
#                double *zeta, double *z, double *theta);

_sofa.iauS00.argtypes = [
    _ct.c_double, # date1
    _ct.c_double, # date2
    _ct.c_double, # x
    _ct.c_double  # y
]
_sofa.iauS00.restype = _ct.c_double

def S00(date1, date2, x, y):
    '''Compute CIO locator s, positioning the Celestial Intermediate Origin on
    the equator of the Celestial Intermediate Pole, given the CIP's X,Y
    coordinates.  Compatible with IAU 2000A precession-nutation.

    Args:
        date1,date2 (double): TT as a 2-part Julian Date (Note 1)
        x,y (double): CIP coordinates (Note 3)

    Returns:
        s (double): CIO locator s in radians (Note 2)

    Notes:

        1. The TT date date1+date2 is a Julian Date, apportioned in any
        convenient way between the two arguments.  For example,
        JD(TT)=2450123.7 could be expressed in any of these ways,
        among others:
   
               date1          date2
   
            2450123.7           0.0       (JD method)
            2451545.0       -1421.3       (J2000 method)
            2400000.5       50123.2       (MJD method)
            2450123.5           0.2       (date & time method)
   
        The JD method is the most natural and convenient to use in
        cases where the loss of several decimal digits of resolution
        is acceptable.  The J2000 method is best matched to the way
        the argument is handled internally and will deliver the
        optimum resolution.  The MJD method and the date & time methods
        are both good compromises between resolution and convenience.
   
        2. The CIO locator s is the difference between the right ascensions
        of the same point in two systems:  the two systems are the GCRS
        and the CIP,CIO, and the point is the ascending node of the
        CIP equator.  The quantity s remains below 0.1 arcsecond
        throughout 1900-2100.
   
        3. The series used to compute s is in fact for s+XY/2, where X and Y
        are the x and y components of the CIP unit vector;  this series
        is more compact than a direct series for s would be.  This
        function requires X,Y to be supplied by the caller, who is
        responsible for providing values that are consistent with the
        supplied date.
   
        4. The model is consistent with the IAU 2000A precession-nutation.
    '''

    return _sofa.iauS00(date1, date2, x, y)

_sofa.iauS00a.argtypes = [
    _ct.c_double, # date1
    _ct.c_double  # date2
]
_sofa.iauS00a.restype = _ct.c_double

def S00a(date1, date2):
    '''Compute CIO locator s, positioning the Celestial Intermediate Origin on
    the equator of the Celestial Intermediate Pole, using the IAU 2000A
    precession-nutation model.

    Args:
        date1,date2 (double): TT as a 2-part Julian Date (Note 1)

    Returns:
        s (double): CIO locator s in radians (Note 2)

    Notes:

        1. The TT date date1+date2 is a Julian Date, apportioned in any
        convenient way between the two arguments.  For example,
        JD(TT)=2450123.7 could be expressed in any of these ways,
        among others:
   
               date1          date2
   
            2450123.7           0.0       (JD method)
            2451545.0       -1421.3       (J2000 method)
            2400000.5       50123.2       (MJD method)
            2450123.5           0.2       (date & time method)
   
        The JD method is the most natural and convenient to use in
        cases where the loss of several decimal digits of resolution
        is acceptable.  The J2000 method is best matched to the way
        the argument is handled internally and will deliver the
        optimum resolution.  The MJD method and the date & time methods
        are both good compromises between resolution and convenience.
   
        2. The CIO locator s is the difference between the right ascensions
        of the same point in two systems.  The two systems are the GCRS
        and the CIP,CIO, and the point is the ascending node of the
        CIP equator.  The CIO locator s remains a small fraction of
        1 arcsecond throughout 1900-2100.
   
        3. The series used to compute s is in fact for s+XY/2, where X and Y
        are the x and y components of the CIP unit vector;  this series
        is more compact than a direct series for s would be.  The present
        function uses the full IAU 2000A nutation model when predicting
        the CIP position.  Faster results, with no significant loss of
        accuracy, can be obtained via the function iauS00b, which uses
        instead the IAU 2000B truncated model.
    '''

    return _sofa.iauS00a(date1, date2)

_sofa.iauS00b.argtypes = [
    _ct.c_double, # date1
    _ct.c_double  # date2
]
_sofa.iauS00b.restype = _ct.c_double

def S00b(date1, date2):
    '''Compute CIO locator s, positioning the Celestial Intermediate Origin on
    the equator of the Celestial Intermediate Pole, using the IAU 2000A
    precession-nutation model.

    Args:
        date1,date2 (double): TT as a 2-part Julian Date (Note 1)

    Returns:
        s (double): CIO locator s in radians (Note 2)

    Notes:

        1. The TT date date1+date2 is a Julian Date, apportioned in any
        convenient way between the two arguments.  For example,
        JD(TT)=2450123.7 could be expressed in any of these ways,
        among others:
   
               date1          date2
   
            2450123.7           0.0       (JD method)
            2451545.0       -1421.3       (J2000 method)
            2400000.5       50123.2       (MJD method)
            2450123.5           0.2       (date & time method)
   
        The JD method is the most natural and convenient to use in
        cases where the loss of several decimal digits of resolution
        is acceptable.  The J2000 method is best matched to the way
        the argument is handled internally and will deliver the
        optimum resolution.  The MJD method and the date & time methods
        are both good compromises between resolution and convenience.
   
        2. The CIO locator s is the difference between the right ascensions
        of the same point in two systems.  The two systems are the GCRS
        and the CIP,CIO, and the point is the ascending node of the
        CIP equator.  The CIO locator s remains a small fraction of
        1 arcsecond throughout 1900-2100.
   
        3. The series used to compute s is in fact for s+XY/2, where X and Y
        are the x and y components of the CIP unit vector;  this series
        is more compact than a direct series for s would be.  The present
        function uses the IAU 2000B truncated nutation model when
        predicting the CIP position.  The function iauS00a uses instead
        the full IAU 2000A model, but with no significant increase in
        accuracy and at some cost in speed.
    '''

    return _sofa.iauS00b(date1, date2)
    
_sofa.iauS06.argtypes = [
    _ct.c_double, # date1
    _ct.c_double, # date2
    _ct.c_double, # x
    _ct.c_double  # y
]
_sofa.iauS06.restype = _ct.c_double

def S06(date1, date2, x, y):
    '''The CIO locator s, positioning the Celestial Intermediate Origin on
    the equator of the Celestial Intermediate Pole, given the CIP's X,Y
    coordinates.  Compatible with IAU 2006/2000A precession-nutation. Cannonical
    model.

    Args:
        date1,date2 (double): TT as a 2-part Julian Date (Note 1)
        x,y (double): CIP coordinates (Note 3)

    Returns:
        s (double): CIO locator s in radians (Note 2)

    Notes:

        1. The TT date date1+date2 is a Julian Date, apportioned in any
        convenient way between the two arguments.  For example,
        JD(TT)=2450123.7 could be expressed in any of these ways,
        among others:
   
               date1          date2
   
            2450123.7           0.0       (JD method)
            2451545.0       -1421.3       (J2000 method)
            2400000.5       50123.2       (MJD method)
            2450123.5           0.2       (date & time method)
   
        The JD method is the most natural and convenient to use in
        cases where the loss of several decimal digits of resolution
        is acceptable.  The J2000 method is best matched to the way
        the argument is handled internally and will deliver the
        optimum resolution.  The MJD method and the date & time methods
        are both good compromises between resolution and convenience.
   
        2. The CIO locator s is the difference between the right ascensions
        of the same point in two systems:  the two systems are the GCRS
        and the CIP,CIO, and the point is the ascending node of the
        CIP equator.  The quantity s remains below 0.1 arcsecond
        throughout 1900-2100.
   
        3. The series used to compute s is in fact for s+XY/2, where X and Y
        are the x and y components of the CIP unit vector;  this series
        is more compact than a direct series for s would be.  This
        function requires X,Y to be supplied by the caller, who is
        responsible for providing values that are consistent with the
        supplied date.
   
        4. The model is consistent with the "P03" precession (Capitaine et
        al. 2003), adopted by IAU 2006 Resolution 1, 2006, and the
        IAU 2000A nutation (with P03 adjustments).
    '''

    return _sofa.iauS06(date1, date2, x, y)

_sofa.iauS06a.argtypes = [
    _ct.c_double, # date1
    _ct.c_double  # date2
]
_sofa.iauS06a.restype = _ct.c_double

def S06a(date1, date2):
    '''The CIO locator s, positioning the Celestial Intermediate Origin on
    the equator of the Celestial Intermediate Pole, using the IAU 2006
    precession and IAU 2000A nutation models.

    Args:
        date1,date2 (double): TT as a 2-part Julian Date (Note 1)

    Returns:
        s (double): CIO locator s in radians (Note 2)

    Notes:

        1. The TT date date1+date2 is a Julian Date, apportioned in any
        convenient way between the two arguments.  For example,
        JD(TT)=2450123.7 could be expressed in any of these ways,
        among others:
   
               date1          date2
   
            2450123.7           0.0       (JD method)
            2451545.0       -1421.3       (J2000 method)
            2400000.5       50123.2       (MJD method)
            2450123.5           0.2       (date & time method)
   
        The JD method is the most natural and convenient to use in
        cases where the loss of several decimal digits of resolution
        is acceptable.  The J2000 method is best matched to the way
        the argument is handled internally and will deliver the
        optimum resolution.  The MJD method and the date & time methods
        are both good compromises between resolution and convenience.
   
        2. The CIO locator s is the difference between the right ascensions
        of the same point in two systems.  The two systems are the GCRS
        and the CIP,CIO, and the point is the ascending node of the
        CIP equator.  The CIO locator s remains a small fraction of
        1 arcsecond throughout 1900-2100.
   
        3. The series used to compute s is in fact for s+XY/2, where X and Y
        are the x and y components of the CIP unit vector;  this series is
        more compact than a direct series for s would be.  The present
        function uses the full IAU 2000A nutation model when predicting
        the CIP position.
    '''

    return _sofa.iauS06a(date1, date2)

_sofa.iauSp00.argtypes = [
    _ct.c_double, # date1
    _ct.c_double  # date2
]
_sofa.iauSp00.restype = _ct.c_double

def Sp00(date1, date2):
    '''The TIO locator s', positioning the Terrestrial Intermediate Origin
    on the equator of the Celestial Intermediate Pole.

    Args:
        date1,date2 (double): TT as a 2-part Julian Date (Note 1)

    Returns:
        s (double): CIO locator s in radians (Note 2)

    Notes:

        1. The TT date date1+date2 is a Julian Date, apportioned in any
        convenient way between the two arguments.  For example,
        JD(TT)=2450123.7 could be expressed in any of these ways,
        among others:
   
               date1          date2
   
            2450123.7           0.0       (JD method)
            2451545.0       -1421.3       (J2000 method)
            2400000.5       50123.2       (MJD method)
            2450123.5           0.2       (date & time method)
   
        The JD method is the most natural and convenient to use in
        cases where the loss of several decimal digits of resolution
        is acceptable.  The J2000 method is best matched to the way
        the argument is handled internally and will deliver the
        optimum resolution.  The MJD method and the date & time methods
        are both good compromises between resolution and convenience.
   
        2. The TIO locator s' is obtained from polar motion observations by
        numerical integration, and so is in essence unpredictable.
        However, it is dominated by a secular drift of about
        47 microarcseconds per century, which is the approximation
        evaluated by the present function.
    '''

    return _sofa.iauSp00(date1, date2)

_sofa.iauXy06.argtypes = [
    _ct.c_double, # date1
    _ct.c_double  # date2
]
_sofa.iauXy06.restype = None

def Xy06(date1, date2):
    '''X,Y coordinates of celestial intermediate pole from series based
    on IAU 2006 precession and IAU 2000A nutation.

    Args:
        date1,date2 (double): TT as a 2-part Julian Date (Note 1)

    Returns:
        x,y (double): CIP X,Y coordinates (Note 2)

    Notes:

        1. The TT date date1+date2 is a Julian Date, apportioned in any
        convenient way between the two arguments.  For example,
        JD(TT)=2450123.7 could be expressed in any of these ways,
        among others:
   
               date1          date2
   
            2450123.7           0.0       (JD method)
            2451545.0       -1421.3       (J2000 method)
            2400000.5       50123.2       (MJD method)
            2450123.5           0.2       (date & time method)
   
        The JD method is the most natural and convenient to use in
        cases where the loss of several decimal digits of resolution
        is acceptable.  The J2000 method is best matched to the way
        the argument is handled internally and will deliver the
        optimum resolution.  The MJD method and the date & time methods
        are both good compromises between resolution and convenience.
   
        2. The X,Y coordinates are those of the unit vector towards the
        celestial intermediate pole.  They represent the combined effects
        of frame bias, precession and nutation.
   
        3. The fundamental arguments used are as adopted in IERS Conventions
        (2003) and are from Simon et al. (1994) and Souchay et al. (1999).

        4. This is an alternative to the angles-based method, via the SOFA
        function iauFw2xy and as used in iauXys06a for example.  The two
        methods agree at the 1 microarcsecond level (at present), a
        negligible amount compared with the intrinsic accuracy of the
        models.  However, it would be unwise to mix the two methods
        (angles-based and series-based) in a single application. 
    '''

    # Initialized returned values
    x = _ct.c_double()
    y = _ct.c_double()

    # Main function call
    _sofa.iauXy06(date1, date2, _ct.byref(x),  _ct.byref(y))

    return x.value, y.value

_sofa.iauXys00a.argtypes = [
    _ct.c_double, # date1
    _ct.c_double  # date2
]
_sofa.iauXys00a.restype = None

def Xys00a(date1, date2):
    '''For a given TT date, compute the X,Y coordinates of the Celestial
    Intermediate Pole and the CIO locator s, using the IAU 2000A
    precession-nutation model.

    Args:
        date1,date2 (double): TT as a 2-part Julian Date (Note 1)

    Returns:
        x,y (double): Celestial Intermediate Pole (Note 2)
        s (double): CIO locator s in radians (Note 2)

    Notes:

        1. The TT date date1+date2 is a Julian Date, apportioned in any
        convenient way between the two arguments.  For example,
        JD(TT)=2450123.7 could be expressed in any of these ways,
        among others:
   
               date1          date2
   
            2450123.7           0.0       (JD method)
            2451545.0       -1421.3       (J2000 method)
            2400000.5       50123.2       (MJD method)
            2450123.5           0.2       (date & time method)
   
        The JD method is the most natural and convenient to use in
        cases where the loss of several decimal digits of resolution
        is acceptable.  The J2000 method is best matched to the way
        the argument is handled internally and will deliver the
        optimum resolution.  The MJD method and the date & time methods
        are both good compromises between resolution and convenience.
   
        2. The Celestial Intermediate Pole coordinates are the x,y
        components of the unit vector in the Geocentric Celestial
        Reference System.

        3. The CIO locator s (in radians) positions the Celestial
        Intermediate Origin on the equator of the CIP.

        4. A faster, but slightly less accurate result (about 1 mas for
        X,Y), can be obtained by using instead the iauXys00b function.
    '''

    # Initialized returned values
    x = _ct.c_double()
    y = _ct.c_double()
    s = _ct.c_double()

    # Main function call
    _sofa.iauXys00a(date1, date2, _ct.byref(x), _ct.byref(y), _ct.byref(s))

    return x.value, y.value, s.value

_sofa.iauXys00b.argtypes = [
    _ct.c_double, # date1
    _ct.c_double  # date2
]
_sofa.iauXys00b.restype = None

def Xys00b(date1, date2):
    '''For a given TT date, compute the X,Y coordinates of the Celestial
    Intermediate Pole and the CIO locator s, using the IAU 2000B
    precession-nutation model.

    Args:
        date1,date2 (double): TT as a 2-part Julian Date (Note 1)

    Returns:
        x,y (double): Celestial Intermediate Pole (Note 2)
        s (double): CIO locator s in radians (Note 2)

    Notes:

        1. The TT date date1+date2 is a Julian Date, apportioned in any
        convenient way between the two arguments.  For example,
        JD(TT)=2450123.7 could be expressed in any of these ways,
        among others:
   
               date1          date2
   
            2450123.7           0.0       (JD method)
            2451545.0       -1421.3       (J2000 method)
            2400000.5       50123.2       (MJD method)
            2450123.5           0.2       (date & time method)
   
        The JD method is the most natural and convenient to use in
        cases where the loss of several decimal digits of resolution
        is acceptable.  The J2000 method is best matched to the way
        the argument is handled internally and will deliver the
        optimum resolution.  The MJD method and the date & time methods
        are both good compromises between resolution and convenience.
   
        2. The Celestial Intermediate Pole coordinates are the x,y
        components of the unit vector in the Geocentric Celestial
        Reference System.

        3. The CIO locator s (in radians) positions the Celestial
        Intermediate Origin on the equator of the CIP.

        4. The present function is faster, but slightly less accurate (about
        1 mas in X,Y), than the iauXys00a function.
    '''

    # Initialized returned values
    x = _ct.c_double()
    y = _ct.c_double()
    s = _ct.c_double()

    # Main function call
    _sofa.iauXys00b(date1, date2, _ct.byref(x), _ct.byref(y), _ct.byref(s))

    return x.value, y.value, s.value

_sofa.iauXys06a.argtypes = [
    _ct.c_double, # date1
    _ct.c_double  # date2
]
_sofa.iauXys06a.restype = None

def Xys06a(date1, date2):
    '''For a given TT date, compute the X,Y coordinates of the Celestial
    Intermediate Pole and the CIO locator s, using the IAU 2006
    precession and IAU 2000A nutation models.

    Args:
        date1,date2 (double): TT as a 2-part Julian Date (Note 1)

    Returns:
        x,y (double): Celestial Intermediate Pole (Note 2)
        s (double): CIO locator s in radians (Note 2)

    Notes:

        1. The TT date date1+date2 is a Julian Date, apportioned in any
        convenient way between the two arguments.  For example,
        JD(TT)=2450123.7 could be expressed in any of these ways,
        among others:
   
               date1          date2
   
            2450123.7           0.0       (JD method)
            2451545.0       -1421.3       (J2000 method)
            2400000.5       50123.2       (MJD method)
            2450123.5           0.2       (date & time method)
   
        The JD method is the most natural and convenient to use in
        cases where the loss of several decimal digits of resolution
        is acceptable.  The J2000 method is best matched to the way
        the argument is handled internally and will deliver the
        optimum resolution.  The MJD method and the date & time methods
        are both good compromises between resolution and convenience.
   
        2. The Celestial Intermediate Pole coordinates are the x,y
        components of the unit vector in the Geocentric Celestial
        Reference System.

        3. The CIO locator s (in radians) positions the Celestial
        Intermediate Origin on the equator of the CIP.

        4. Series-based solutions for generating X and Y are also available:
        see Capitaine & Wallace (2006) and iauXy06.
    '''

    # Initialized returned values
    x = _ct.c_double()
    y = _ct.c_double()
    s = _ct.c_double()

    # Main function call
    _sofa.iauXys06a(date1, date2, _ct.byref(x), _ct.byref(y), _ct.byref(s))

    return x.value, y.value, s.value

# /* Astronomy/RotationAndTime */

_sofa.iauEe00.argtypes = [
    _ct.c_double, # date1
    _ct.c_double, # date2
    _ct.c_double, # epsa
    _ct.c_double  # dpsi
]
_sofa.iauEe00.restype = _ct.c_double

def Ee00(date1, date2, epsa, dpsi):
    '''The equation of the equinoxes, compatible with IAU 2000 resolutions,
    given the nutation in longitude and the mean obliquity.

    Args:
        date1,date2 (double): TT as a 2-part Julian Date (Note 1)
        epsa (double): mean obliquity (Note 2)
        dpsi (double): nutation in longitude (Note 3)

    Returns:
        ee (double): Equation of the equinoxes (Note 4)

    Notes:

        1. The TT date date1+date2 is a Julian Date, apportioned in any
        convenient way between the two arguments.  For example,
        JD(TT)=2450123.7 could be expressed in any of these ways,
        among others:
   
               date1          date2
   
            2450123.7           0.0       (JD method)
            2451545.0       -1421.3       (J2000 method)
            2400000.5       50123.2       (MJD method)
            2450123.5           0.2       (date & time method)
   
        The JD method is the most natural and convenient to use in
        cases where the loss of several decimal digits of resolution
        is acceptable.  The J2000 method is best matched to the way
        the argument is handled internally and will deliver the
        optimum resolution.  The MJD method and the date & time methods
        are both good compromises between resolution and convenience.

        2. The obliquity, in radians, is mean of date.

        3. The result, which is in radians, operates in the following sense:
   
           Greenwich apparent ST = GMST + equation of the equinoxes

        4. The result is compatible with the IAU 2000 resolutions.  For
        further details, see IERS Conventions 2003 and Capitaine et al.
        (2002).
    '''

    return _sofa.iauEe00(date1, date2, epsa, dpsi)

_sofa.iauEe00a.argtypes = [
    _ct.c_double, # date1
    _ct.c_double  # date2
]
_sofa.iauEe00a.restype = _ct.c_double

def Ee00a(date1, date2):
    '''Equation of the equinoxes, compatible with IAU 2000 resolutions.

    Args:
        date1,date2 (double): TT as a 2-part Julian Date (Note 1)

    Returns:
        ee (double): Equation of the equinoxes (Note 2)

    Notes:

        1. The TT date date1+date2 is a Julian Date, apportioned in any
        convenient way between the two arguments.  For example,
        JD(TT)=2450123.7 could be expressed in any of these ways,
        among others:
   
               date1          date2
   
            2450123.7           0.0       (JD method)
            2451545.0       -1421.3       (J2000 method)
            2400000.5       50123.2       (MJD method)
            2450123.5           0.2       (date & time method)
   
        The JD method is the most natural and convenient to use in
        cases where the loss of several decimal digits of resolution
        is acceptable.  The J2000 method is best matched to the way
        the argument is handled internally and will deliver the
        optimum resolution.  The MJD method and the date & time methods
        are both good compromises between resolution and convenience.

        2. The result, which is in radians, operates in the following sense:
   
           Greenwich apparent ST = GMST + equation of the equinoxes

        3. The result is compatible with the IAU 2000 resolutions.  For
        further details, see IERS Conventions 2003 and Capitaine et al.
        (2002).
    '''

    return _sofa.iauEe00a(date1, date2)

_sofa.iauEe00b.argtypes = [
    _ct.c_double, # date1
    _ct.c_double  # date2
]
_sofa.iauEe00b.restype = _ct.c_double

def Ee00b(date1, date2):
    '''Equation of the equinoxes, compatible with IAU 2000 resolutions but
    using the truncated nutation model IAU 2000B.

    Args:
        date1,date2 (double): TT as a 2-part Julian Date (Note 1)

    Returns:
        ee (double): Equation of the equinoxes (Note 2)

    Notes:

        1. The TT date date1+date2 is a Julian Date, apportioned in any
        convenient way between the two arguments.  For example,
        JD(TT)=2450123.7 could be expressed in any of these ways,
        among others:
   
               date1          date2
   
            2450123.7           0.0       (JD method)
            2451545.0       -1421.3       (J2000 method)
            2400000.5       50123.2       (MJD method)
            2450123.5           0.2       (date & time method)
   
        The JD method is the most natural and convenient to use in
        cases where the loss of several decimal digits of resolution
        is acceptable.  The J2000 method is best matched to the way
        the argument is handled internally and will deliver the
        optimum resolution.  The MJD method and the date & time methods
        are both good compromises between resolution and convenience.

        2. The result, which is in radians, operates in the following sense:
   
           Greenwich apparent ST = GMST + equation of the equinoxes

        3. The result is compatible with the IAU 2000 resolutions.  For
        further details, see IERS Conventions 2003 and Capitaine et al.
        (2002).
    '''

    return _sofa.iauEe00b(date1, date2)

_sofa.iauEe06a.argtypes = [
    _ct.c_double, # date1
    _ct.c_double  # date2
]
_sofa.iauEe06a.restype = _ct.c_double

def Ee06a(date1, date2):
    '''Equation of the equinoxes, compatible with IAU 2000 resolutions and
    IAU 2006/2000A precession-nutation.

    Args:
        date1,date2 (double): TT as a 2-part Julian Date (Note 1)

    Returns:
        ee (double): Equation of the equinoxes (Note 2)

    Notes:

        1. The TT date date1+date2 is a Julian Date, apportioned in any
        convenient way between the two arguments.  For example,
        JD(TT)=2450123.7 could be expressed in any of these ways,
        among others:
   
               date1          date2
   
            2450123.7           0.0       (JD method)
            2451545.0       -1421.3       (J2000 method)
            2400000.5       50123.2       (MJD method)
            2450123.5           0.2       (date & time method)
   
        The JD method is the most natural and convenient to use in
        cases where the loss of several decimal digits of resolution
        is acceptable.  The J2000 method is best matched to the way
        the argument is handled internally and will deliver the
        optimum resolution.  The MJD method and the date & time methods
        are both good compromises between resolution and convenience.

        2. The result, which is in radians, operates in the following sense:
   
           Greenwich apparent ST = GMST + equation of the equinoxes
    '''

    return _sofa.iauEe06a(date1, date2)

_sofa.iauEect00.argtypes = [
    _ct.c_double, # date1
    _ct.c_double  # date2
]
_sofa.iauEect00.restype = _ct.c_double

def Eect00(date1, date2):
    '''Equation of the equinoxes complementary terms, consistent with
    IAU 2000 resolutions.

    Args:
        date1,date2 (double): TT as a 2-part Julian Date (Note 1)

    Returns:
        ee (double): complementary terms (Note 2)

    Notes:

        1. The TT date date1+date2 is a Julian Date, apportioned in any
        convenient way between the two arguments.  For example,
        JD(TT)=2450123.7 could be expressed in any of these ways,
        among others:
   
               date1          date2
   
            2450123.7           0.0       (JD method)
            2451545.0       -1421.3       (J2000 method)
            2400000.5       50123.2       (MJD method)
            2450123.5           0.2       (date & time method)
   
        The JD method is the most natural and convenient to use in
        cases where the loss of several decimal digits of resolution
        is acceptable.  The J2000 method is best matched to the way
        the argument is handled internally and will deliver the
        optimum resolution.  The MJD method and the date & time methods
        are both good compromises between resolution and convenience.

        2. The "complementary terms" are part of the equation of the
        equinoxes (EE), classically the difference between apparent and
        mean Sidereal Time:
   
           GAST = GMST + EE
   
        with:
   
           EE = dpsi * cos(eps)
   
        where dpsi is the nutation in longitude and eps is the obliquity
        of date.  However, if the rotation of the Earth were constant in
        an inertial frame the classical formulation would lead to
        apparent irregularities in the UT1 timescale traceable to side-
        effects of precession-nutation.  In order to eliminate these
        effects from UT1, "complementary terms" were introduced in 1994
        (IAU, 1994) and took effect from 1997 (Capitaine and Gontier,
        1993):
   
           GAST = GMST + CT + EE
   
        By convention, the complementary terms are included as part of
        the equation of the equinoxes rather than as part of the mean
        Sidereal Time.  This slightly compromises the "geometrical"
        interpretation of mean sidereal time but is otherwise
        inconsequential.
   
        The present function computes CT in the above expression,
        compatible with IAU 2000 resolutions (Capitaine et al., 2002, and
        IERS Conventions 2003).
    '''

    return _sofa.iauEect00(date1, date2)

_sofa.iauEqeq94.argtypes = [
    _ct.c_double, # date1
    _ct.c_double  # date2
]
_sofa.iauEqeq94.restype = _ct.c_double

def Eqeq94(date1, date2):
    '''Equation of the equinoxes, IAU 1994 model.

    Args:
        date1,date2 (double): TDB Date (Note 1)

    Returns:
        ee (double): equation of the equinoxes (Note 2)

    Notes:

        1. The TT date date1+date2 is a Julian Date, apportioned in any
        convenient way between the two arguments.  For example,
        JD(TT)=2450123.7 could be expressed in any of these ways,
        among others:
   
               date1          date2
   
            2450123.7           0.0       (JD method)
            2451545.0       -1421.3       (J2000 method)
            2400000.5       50123.2       (MJD method)
            2450123.5           0.2       (date & time method)
   
        The JD method is the most natural and convenient to use in
        cases where the loss of several decimal digits of resolution
        is acceptable.  The J2000 method is best matched to the way
        the argument is handled internally and will deliver the
        optimum resolution.  The MJD method and the date & time methods
        are both good compromises between resolution and convenience.

        2. The result, which is in radians, operates in the following sense:
   
           Greenwich apparent ST = GMST + equation of the equinoxes
    '''

    return _sofa.iauEqeq94(date1, date2)

_sofa.iauEra00.argtypes = [
    _ct.c_double, # dj1
    _ct.c_double  # dj2
]
_sofa.iauEra00.restype = _ct.c_double

def Era00(dj1, dj2):
    '''Earth rotation angle (IAU 2000 model).

    Args:
        dj1,dj2 (double): UT1 as a 2-part Julian Date (see note)

    Returns:
        theta (double): Earth rotation angle (radians), range 0-2pi

    Notes:

        1. The UT1 date dj1+dj2 is a Julian Date, apportioned in any
        convenient way between the arguments dj1 and dj2.  For example,
        JD(UT1)=2450123.7 could be expressed in any of these ways,
        among others:
   
                dj1            dj2
   
            2450123.7           0.0       (JD method)
            2451545.0       -1421.3       (J2000 method)
            2400000.5       50123.2       (MJD method)
            2450123.5           0.2       (date & time method)
   
        The JD method is the most natural and convenient to use in
        cases where the loss of several decimal digits of resolution
        is acceptable.  The J2000 and MJD methods are good compromises
        between resolution and convenience.  The date & time method is
        best matched to the algorithm used:  maximum precision is
        delivered when the dj1 argument is for 0hrs UT1 on the day in
        question and the dj2 argument lies in the range 0 to 1, or vice
        versa.

        2. The algorithm is adapted from Expression 22 of Capitaine et al.
        2000.  The time argument has been expressed in days directly,
        and, to retain precision, integer contributions have been
        eliminated.  The same formulation is given in IERS Conventions
        (2003), Chap. 5, Eq. 14.
    '''

    return _sofa.iauEra00(dj1, dj2)

_sofa.iauGmst00.argtypes = [
    _ct.c_double, # uta
    _ct.c_double, # utb
    _ct.c_double, # tta
    _ct.c_double  # ttb
]

_sofa.iauGmst00.restype = _ct.c_double

def Gmst00(uta, utb, tta, ttb):
    '''Greenwich mean sidereal time (model consistent with IAU 2000
    resolutions).

    Args:
        uta,utb (double): UT1 as a 2-part Julian Date (Notes 1,2)
        tta,ttb (double): TT as a 2-part Julian Date (Notes 1,2)

    Returns:
        theta (double): Greenwich mean sidereal time (radians)

    Notes:

        1. The UT1 and TT dates uta+utb and tta+ttb respectively, are both
        Julian Dates, apportioned in any convenient way between the
        argument pairs.  For example, JD=2450123.7 could be expressed in
        any of these ways, among others:
   
               Part A         Part B
   
            2450123.7           0.0       (JD method)
            2451545.0       -1421.3       (J2000 method)
            2400000.5       50123.2       (MJD method)
            2450123.5           0.2       (date & time method)
   
        The JD method is the most natural and convenient to use in
        cases where the loss of several decimal digits of resolution
        is acceptable (in the case of UT;  the TT is not at all critical
        in this respect).  The J2000 and MJD methods are good compromises
        between resolution and convenience.  For UT, the date & time
        method is best matched to the algorithm that is used by the Earth
        Rotation Angle function, called internally:  maximum precision is
        delivered when the uta argument is for 0hrs UT1 on the day in
        question and the utb argument lies in the range 0 to 1, or vice
        versa.

        2. Both UT1 and TT are required, UT1 to predict the Earth rotation
        and TT to predict the effects of precession.  If UT1 is used for
        both purposes, errors of order 100 microarcseconds result.

        3. This GMST is compatible with the IAU 2000 resolutions and must be
        used only in conjunction with other IAU 2000 compatible
        components such as precession-nutation and equation of the
        equinoxes.

        4. The result is returned in the range 0 to 2pi.

        5. The algorithm is from Capitaine et al. (2003) and IERS Conventions 2003.
    '''

    return _sofa.iauGmst00(uta, utb, tta, ttb)

_sofa.iauGmst06.argtypes = [
    _ct.c_double, # uta
    _ct.c_double, # utb
    _ct.c_double, # tta
    _ct.c_double  # ttb
]

_sofa.iauGmst06.restype = _ct.c_double

def Gmst06(uta, utb, tta, ttb):
    '''Greenwich mean sidereal time (consistent with IAU 2006 precession).

    Args:
        uta,utb (double): UT1 as a 2-part Julian Date (Notes 1,2)
        tta,ttb (double): TT as a 2-part Julian Date (Notes 1,2)

    Returns:
        theta (double): Earth rotation angle (radians), range 0-2pi

    Notes:

        1. The UT1 and TT dates uta+utb and tta+ttb respectively, are both
        Julian Dates, apportioned in any convenient way between the
        argument pairs.  For example, JD=2450123.7 could be expressed in
        any of these ways, among others:
   
               Part A         Part B
   
            2450123.7           0.0       (JD method)
            2451545.0       -1421.3       (J2000 method)
            2400000.5       50123.2       (MJD method)
            2450123.5           0.2       (date & time method)
   
        The JD method is the most natural and convenient to use in
        cases where the loss of several decimal digits of resolution
        is acceptable (in the case of UT;  the TT is not at all critical
        in this respect).  The J2000 and MJD methods are good compromises
        between resolution and convenience.  For UT, the date & time
        method is best matched to the algorithm that is used by the Earth
        Rotation Angle function, called internally:  maximum precision is
        delivered when the uta argument is for 0hrs UT1 on the day in
        question and the utb argument lies in the range 0 to 1, or vice
        versa.

        2. Both UT1 and TT are required, UT1 to predict the Earth rotation
        and TT to predict the effects of precession.  If UT1 is used for
        both purposes, errors of order 100 microarcseconds result.

        3. This GMST is compatible with the IAU 2006 precession and must not
        be used with other precession models.

        4. The result is returned in the range 0 to 2pi.
    '''

    return _sofa.iauGmst06(uta, utb, tta, ttb)

_sofa.iauGmst82.argtypes = [
    _ct.c_double, # dj1
    _ct.c_double  # dj2
]

_sofa.iauGmst82.restype = _ct.c_double

def Gmst82(dj1, dj2):
    '''Universal Time to Greenwich mean sidereal time (IAU 1982 model).

    Args:
        dj1,dj2 (double): UT1 Julian Date (see note)

    Returns:
        theta (double): Earth rotation angle (radians), range 0-2pi

    Notes:

        1. TThe UT1 date dj1+dj2 is a Julian Date, apportioned in any
        convenient way between the arguments dj1 and dj2.  For example,
        JD(UT1)=2450123.7 could be expressed in any of these ways,
        among others:
   
                dj1            dj2
   
            2450123.7          0          (JD method)
             2451545        -1421.3       (J2000 method)
            2400000.5       50123.2       (MJD method)
            2450123.5         0.2         (date & time method)
   
        The JD method is the most natural and convenient to use in
        cases where the loss of several decimal digits of resolution
        is acceptable.  The J2000 and MJD methods are good compromises
        between resolution and convenience.  The date & time method is
        best matched to the algorithm used:  maximum accuracy (or, at
        least, minimum noise) is delivered when the dj1 argument is for
        0hrs UT1 on the day in question and the dj2 argument lies in the
        range 0 to 1, or vice versa.

        2. The algorithm is based on the IAU 1982 expression.  This is
        always described as giving the GMST at 0 hours UT1.  In fact, it
        gives the difference between the GMST and the UT, the steady
        4-minutes-per-day drawing-ahead of ST with respect to UT.  When
        whole days are ignored, the expression happens to equal the GMST
        at 0 hours UT1 each day.

        3. In this function, the entire UT1 (the sum of the two arguments
        dj1 and dj2) is used directly as the argument for the standard
        formula, the constant term of which is adjusted by 12 hours to
        take account of the noon phasing of Julian Date.  The UT1 is then
        added, but omitting whole days to conserve accuracy.
    '''

    return _sofa.iauGmst82(dj1, dj2)

_sofa.iauGst00a.argtypes = [
    _ct.c_double, # uta
    _ct.c_double, # utb
    _ct.c_double, # tta
    _ct.c_double  # ttb
]

_sofa.iauGst00a.restype = _ct.c_double

def Gst00a(uta, utb, tta, ttb):
    '''Greenwich apparent sidereal time (consistent with IAU 2000 resolutions).

    Args:
        uta,utb (double): UT1 as a 2-part Julian Date (Notes 1,2)
        tta,ttb (double): TT as a 2-part Julian Date (Notes 1,2)

    Returns:
        gast (double): Greenwich apparent sidereal time (radians)

    Notes:

        1. The UT1 and TT dates uta+utb and tta+ttb respectively, are both
        Julian Dates, apportioned in any convenient way between the
        argument pairs.  For example, JD=2450123.7 could be expressed in
        any of these ways, among others:
   
               Part A         Part B
   
            2450123.7           0.0       (JD method)
            2451545.0       -1421.3       (J2000 method)
            2400000.5       50123.2       (MJD method)
            2450123.5           0.2       (date & time method)
   
        The JD method is the most natural and convenient to use in
        cases where the loss of several decimal digits of resolution
        is acceptable (in the case of UT;  the TT is not at all critical
        in this respect).  The J2000 and MJD methods are good compromises
        between resolution and convenience.  For UT, the date & time
        method is best matched to the algorithm that is used by the Earth
        Rotation Angle function, called internally:  maximum precision is
        delivered when the uta argument is for 0hrs UT1 on the day in
        question and the utb argument lies in the range 0 to 1, or vice
        versa.

        2. Both UT1 and TT are required, UT1 to predict the Earth rotation
        and TT to predict the effects of precession-nutation.  If UT1 is
        used for both purposes, errors of order 100 microarcseconds
        result.

        3. This GAST is compatible with the IAU 2000 resolutions and must be
        used only in conjunction with other IAU 2000 compatible
        components such as precession-nutation.

        4. The result is returned in the range 0 to 2pi.

        5. The algorithm is from Capitaine et al. (2003) and IERS Conventions 2003.
    '''

    return _sofa.iauGst00a(uta, utb, tta, ttb)

_sofa.iauGst00b.argtypes = [
    _ct.c_double, # uta
    _ct.c_double  # utb
]

_sofa.iauGst00b.restype = _ct.c_double

def Gst00b(uta, utb):
    '''Greenwich apparent sidereal time (consistent with IAU 2000
    resolutions but using the truncated nutation model IAU 2000B).

    Args:
        uta,utb (double): UT1 as a 2-part Julian Date (Notes 1,2)

    Returns:
        gast (double): Greenwich apparent sidereal time (radians)

    Notes:

        1. The UT1 date uta+utb is a Julian Date, apportioned in any
        convenient way between the argument pair.  For example,
        JD=2450123.7 could be expressed in any of these ways, among
        others:
   
                uta            utb
   
            2450123.7           0.0       (JD method)
            2451545.0       -1421.3       (J2000 method)
            2400000.5       50123.2       (MJD method)
            2450123.5           0.2       (date & time method)
   
        The JD method is the most natural and convenient to use in cases
        where the loss of several decimal digits of resolution is
        acceptable.  The J2000 and MJD methods are good compromises
        between resolution and convenience.  For UT, the date & time
        method is best matched to the algorithm that is used by the Earth
        Rotation Angle function, called internally:  maximum precision is
        delivered when the uta argument is for 0hrs UT1 on the day in
        question and the utb argument lies in the range 0 to 1, or vice
        versa.

        2. The result is compatible with the IAU 2000 resolutions, except
        that accuracy has been compromised for the sake of speed and
        convenience in two respects:
   
        . UT is used instead of TDB (or TT) to compute the precession
          component of GMST and the equation of the equinoxes.  This
          results in errors of order 0.1 mas at present.
   
        .  The IAU 2000B abridged nutation model (McCarthy & Luzum, 2001)
          is used, introducing errors of up to 1 mas.

        3. This GAST is compatible with the IAU 2000 resolutions and must be
        used only in conjunction with other IAU 2000 compatible
        components such as precession-nutation.

        4. The result is returned in the range 0 to 2pi.

        5. The algorithm is from Capitaine et al. (2003) and IERS Conventions 2003.
    '''

    return _sofa.iauGst00b(uta, utb)

_sofa.iauGst06.argtypes = [
    _ct.c_double, #uta
    _ct.c_double, #utb
    _ct.c_double, #tta
    _ct.c_double, #ttb
    _ndpointer(shape=(3,3), dtype=float, flags='C')
] #rnpb
_sofa.iauGst06.restype = _ct.c_double

def Gst06(uta, utb, tta, ttb, rnpb):
    '''Greenwich apparent sidereal time, IAU 2006, given the NPB matrix.

    Args:
        uta,utb (double): UT1 as a 2-part Julian Date (Notes 1,2)
        tta,ttb (double): TT as a 2-part Julian Date (Notes 1,2)
        rnpb    (:obj:np.ndarray): nutation x precession x bias matrix

    Returns:
        gast (double): Greenwich apparent sidereal time (radians)

    Notes:

        1. The UT1 and TT dates uta+utb and tta+ttb respectively, are both
        Julian Dates, apportioned in any convenient way between the
        argument pairs.  For example, JD=2450123.7 could be expressed in
        any of these ways, among others:
   
               Part A        Part B
   
            2450123.7           0.0       (JD method)
            2451545.0       -1421.3       (J2000 method)
            2400000.5       50123.2       (MJD method)
            2450123.5           0.2       (date & time method)
   
        The JD method is the most natural and convenient to use in
        cases where the loss of several decimal digits of resolution
        is acceptable (in the case of UT;  the TT is not at all critical
        in this respect).  The J2000 and MJD methods are good compromises
        between resolution and convenience.  For UT, the date & time
        method is best matched to the algorithm that is used by the Earth
        rotation angle function, called internally:  maximum precision is
        delivered when the uta argument is for 0hrs UT1 on the day in
        question and the utb argument lies in the range 0 to 1, or vice
        versa.

        2. Both UT1 and TT are required, UT1 to predict the Earth rotation
        and TT to predict the effects of precession-nutation.  If UT1 is
        used for both purposes, errors of order 100 microarcseconds
        result.

        3. Although the function uses the IAU 2006 series for s+XY/2, it is
        otherwise independent of the precession-nutation model and can in
        practice be used with any equinox-based NPB matrix.

        4. The result is returned in the range 0 to 2pi.
    '''

    return _sofa.iauGst06(uta, utb, tta, ttb, _req_shape_c(rnpb, float, (3,3)))


_sofa.iauGst06a.argtypes = [
    _ct.c_double, # uta
    _ct.c_double, # utb
    _ct.c_double, # tta
    _ct.c_double  # ttb
]

_sofa.iauGst06a.restype = _ct.c_double

def Gst06a(uta, utb, tta, ttb):
    '''Greenwich apparent sidereal time (consistent with IAU 2000 and 2006 resolutions).

    Args:
        uta,utb (double): UT1 as a 2-part Julian Date (Notes 1,2)
        tta,ttb (double): TT as a 2-part Julian Date (Notes 1,2)

    Returns:
        gast (double): Greenwich apparent sidereal time (radians)

    Notes:

        1. The UT1 and TT dates uta+utb and tta+ttb respectively, are both
        Julian Dates, apportioned in any convenient way between the
        argument pairs.  For example, JD=2450123.7 could be expressed in
        any of these ways, among others:
   
               Part A        Part B
   
            2450123.7           0.0       (JD method)
            2451545.0       -1421.3       (J2000 method)
            2400000.5       50123.2       (MJD method)
            2450123.5           0.2       (date & time method)
   
        The JD method is the most natural and convenient to use in
        cases where the loss of several decimal digits of resolution
        is acceptable (in the case of UT;  the TT is not at all critical
        in this respect).  The J2000 and MJD methods are good compromises
        between resolution and convenience.  For UT, the date & time
        method is best matched to the algorithm that is used by the Earth
        rotation angle function, called internally:  maximum precision is
        delivered when the uta argument is for 0hrs UT1 on the day in
        question and the utb argument lies in the range 0 to 1, or vice
        versa.

        2. Both UT1 and TT are required, UT1 to predict the Earth rotation
        and TT to predict the effects of precession-nutation.  If UT1 is
        used for both purposes, errors of order 100 microarcseconds
        result.

        3. This GAST is compatible with the IAU 2000/2006 resolutions and
        must be used only in conjunction with IAU 2006 precession and
        IAU 2000A nutation.

        4. The result is returned in the range 0 to 2pi.
    '''

    return _sofa.iauGst06a(uta, utb, tta, ttb)

_sofa.iauGst94.argtypes = [
    _ct.c_double, # uta
    _ct.c_double  # utb
]

_sofa.iauGst94.restype = _ct.c_double

def Gst94(uta, utb):
    '''Greenwich apparent sidereal time (consistent with IAU 1982/94 resolutions).

    Args:
        uta,utb (double): UT1 as a 2-part Julian Date (Notes 1,2)

    Returns:
        gast (double): Greenwich apparent sidereal time (radians)

    Notes:

        1. The UT1 date uta+utb is a Julian Date, apportioned in any
        convenient way between the argument pair.  For example,
        JD=2450123.7 could be expressed in any of these ways, among
        others:
   
                uta            utb
   
            2450123.7           0.0       (JD method)
            2451545.0       -1421.3       (J2000 method)
            2400000.5       50123.2       (MJD method)
            2450123.5           0.2       (date & time method)
   
        The JD method is the most natural and convenient to use in cases
        where the loss of several decimal digits of resolution is
        acceptable.  The J2000 and MJD methods are good compromises
        between resolution and convenience.  For UT, the date & time
        method is best matched to the algorithm that is used by the Earth
        Rotation Angle function, called internally:  maximum precision is
        delivered when the uta argument is for 0hrs UT1 on the day in
        question and the utb argument lies in the range 0 to 1, or vice
        versa.

        2. The result is compatible with the IAU 1982 and 1994 resolutions,
        except that accuracy has been compromised for the sake of
        convenience in that UT is used instead of TDB (or TT) to compute
        the equation of the equinoxes.

        3. This GAST must be used only in conjunction with contemporaneous
        IAU standards such as 1976 precession, 1980 obliquity and 1982
        nutation.  It is not compatible with the IAU 2000 resolutions.

        4. The result is returned in the range 0 to 2pi.
    '''

    return _sofa.iauGst94(uta, utb)

# /* Astronomy/SpaceMotion */
# int iauPvstar(double pv[2][3], double *ra, double *dec,
#               double *pmr, double *pmd, double *px, double *rv);
# int iauStarpv(double ra, double dec,
#               double pmr, double pmd, double px, double rv,
#               double pv[2][3]);

# /* Astronomy/StarCatalogs */
# void iauFk52h(double r5, double d5,
#               double dr5, double dd5, double px5, double rv5,
#               double *rh, double *dh,
#               double *drh, double *ddh, double *pxh, double *rvh);
# void iauFk5hip(double r5h[3][3], double s5h[3]);
# void iauFk5hz(double r5, double d5, double date1, double date2,
#               double *rh, double *dh);
# void iauH2fk5(double rh, double dh,
#               double drh, double ddh, double pxh, double rvh,
#               double *r5, double *d5,
#               double *dr5, double *dd5, double *px5, double *rv5);
# void iauHfk5z(double rh, double dh, double date1, double date2,
#               double *r5, double *d5, double *dr5, double *dd5);
# int iauStarpm(double ra1, double dec1,
#               double pmr1, double pmd1, double px1, double rv1,
#               double ep1a, double ep1b, double ep2a, double ep2b,
#               double *ra2, double *dec2,
#               double *pmr2, double *pmd2, double *px2, double *rv2);

# /* Astronomy/EclipticCoordinates */
# void iauEceq06(double date1, double date2, double dl, double db,
#                double *dr, double *dd);
# void iauEcm06(double date1, double date2, double rm[3][3]);
# void iauEqec06(double date1, double date2, double dr, double dd,
#                double *dl, double *db);
# void iauLteceq(double epj, double dl, double db, double *dr, double *dd);
# void iauLtecm(double epj, double rm[3][3]);
# void iauLteqec(double epj, double dr, double dd, double *dl, double *db);

# /* Astronomy/GalacticCoordinates */
# void iauG2icrs(double dl, double db, double *dr, double *dd);
# void iauIcrs2g(double dr, double dd, double *dl, double *db);

# /* Astronomy/GeodeticGeocentric */
# int iauEform(int n, double *a, double *f);
# int iauGc2gd(int n, double xyz[3],
#              double *elong, double *phi, double *height);
# int iauGc2gde(double a, double f, double xyz[3],
#               double *elong, double *phi, double *height);
# int iauGd2gc(int n, double elong, double phi, double height,
#              double xyz[3]);
# int iauGd2gce(double a, double f,
#               double elong, double phi, double height, double xyz[3]);

# /* Astronomy/Timescales */

_sofa.iauD2dtf.argtypes = [
    _ct.c_char_p, #scale
    _ct.c_int, #ndp
    _ct.c_double, #d1
    _ct.c_double, #d2
    _ct.POINTER(_ct.c_int), #iy
    _ct.POINTER(_ct.c_int), #im
    _ct.POINTER(_ct.c_int), #id
    _ct.c_int * 4 #ihmsf
]

_sofa.iauD2dtf.restype = _ct.c_int

_d2dtf_msg = {
    1: 'D2dtf: dubious year',
    -1: 'unacceptable date',
}

def D2dtf(scale, ndp, d1, d2):
    '''Format for output a 2-part Julian Date (or in the case of UTC a
    quasi-JD form that includes special provision for leap seconds).

    Args:
        scale (str): time scale ID (Note 1)
        ndp (int): resolution (Note 2)
        d1,d2 (double): time as a 2-part Julian Date (Notes 3,4)

    Returns:
        iy (int) year in Greogorian Calendar (Note 5)
        im (int) month in Greogorian Calendar (Note 5)
        id (int) day in Greogorian Calendar (Note 5)
        ihmsf (:obj:np.ndarray): hours, minutes, seconds, fraction (Note 1)

    Notes:

        1. scale identifies the time scale.  Only the value "UTC" (in upper
        case) is significant, and enables handling of leap seconds (see
        Note 4).
   
        2. ndp is the number of decimal places in the seconds field, and can
        have negative as well as positive values, such as:
   
        ndp         resolution
        -4            1 00 00
        -3            0 10 00
        -2            0 01 00
        -1            0 00 10
         0            0 00 01
         1            0 00 00.1
         2            0 00 00.01
         3            0 00 00.001
   
        The limits are platform dependent, but a safe range is -5 to +9.
   
        3. d1+d2 is Julian Date, apportioned in any convenient way between
        the two arguments, for example where d1 is the Julian Day Number
        and d2 is the fraction of a day.  In the case of UTC, where the
        use of JD is problematical, special conventions apply:  see the
        next note.
   
        4. JD cannot unambiguously represent UTC during a leap second unless
        special measures are taken.  The SOFA internal convention is that
        the quasi-JD day represents UTC days whether the length is 86399,
        86400 or 86401 SI seconds.  In the 1960-1972 era there were
        smaller jumps (in either direction) each time the linear UTC(TAI)
        expression was changed, and these "mini-leaps" are also included
        in the SOFA convention.
   
        5. The warning status "dubious year" flags UTCs that predate the
        introduction of the time scale or that are too far in the future
        to be trusted.  See iauDat for further details.
   
        6. For calendar conventions and limitations, see iauCal2jd.
    '''

    # Initialize return pointers
    iy    = _ct.c_int()
    im    = _ct.c_int()
    id    = _ct.c_int()
    ihmsf = (_ct.c_int * 4)()

    # Main function call
    status = _sofa.iauD2dtf(bytes(scale, encoding='utf-8'), ndp, d1, d2, _ct.byref(iy), _ct.byref(im), _ct.byref(id), ihmsf)
    
    # Raise warnings
    if status < 0:
        raise ValueError(_d2dtf_msg[status])
    elif status > 0:
        _warnings.warn(_d2dtf_msg[status], UserWarning, 2)

    return iy.value, im.value, id.value, ihmsf

_sofa.iauDat.argtypes = [
    _ct.c_int, #iy
    _ct.c_int, #im
    _ct.c_int, #id
    _ct.c_double, #fd
    _ct.POINTER(_ct.c_double)  #deltat
]

_sofa.iauDat.restype = _ct.c_int

_dat_msg = {
    1: 'Dat: dubious year',
    -1: 'minimum year allowed is -4799',
    -2: 'month must be in 1..12',
    -3: 'day is out of range for this month',
    -4: 'bad fraction of day',
}

def Dat(iy, im, id, fd):
    '''For a given UTC date, calculate Delta(AT) = TAI-UTC.

    Latest leap second:  2016 December 31

    Args:
        iy (int): UTC year (Notes 1 and 2)
        im (int): UTC month (Note 2)
        id (int): UTC day (Notes 2 and 3)
        fd (double): fraction of day (Note 4)

    Returns:
        deltat (double) TAI minus UTC, seconds)

    Notes:

        1. UTC began at 1960 January 1.0 (JD 2436934.5) and it is improper
        to call the function with an earlier date.  If this is attempted,
        zero is returned together with a warning status.
   
        Because leap seconds cannot, in principle, be predicted in
        advance, a reliable check for dates beyond the valid range is
        impossible.  To guard against gross errors, a year five or more
        after the release year of the present function (see the constant
        IYV) is considered dubious.  In this case a warning status is
        returned but the result is computed in the normal way.
   
        For both too-early and too-late years, the warning status is +1.
        This is distinct from the error status -1, which signifies a year
        so early that JD could not be computed.
   
        2. If the specified date is for a day which ends with a leap second,
        the TAI-UTC value returned is for the period leading up to the
        leap second.  If the date is for a day which begins as a leap
        second ends, the TAI-UTC returned is for the period following the
        leap second.
   
        3. The day number must be in the normal calendar range, for example
        1 through 30 for April.  The "almanac" convention of allowing
        such dates as January 0 and December 32 is not supported in this
        function, in order to avoid confusion near leap seconds.
   
        4. The fraction of day is used only for dates before the
        introduction of leap seconds, the first of which occurred at the
        end of 1971.  It is tested for validity (0 to 1 is the valid
        range) even if not used;  if invalid, zero is used and status -4
        is returned.  For many applications, setting fd to zero is
        acceptable;  the resulting error is always less than 3 ms (and
        occurs only pre-1972).
   
        5. The status value returned in the case where there are multiple
        errors refers to the first error detected.  For example, if the
        month and day are 13 and 32 respectively, status -2 (bad month)
        will be returned.  The "internal error" status refers to a
        case that is impossible but causes some compilers to issue a
        warning.
   
        6. In cases where a valid result is not available, zero is returned.
    '''

    # Initialize return values
    deltat = _ct.c_double()

    # Main funcitno call
    status = _sofa.iauDat(iy, im, id, fd, _ct.byref(deltat))

    # Raise warnings
    if status < 0:
        raise ValueError(_dat_msg[status])
    elif status > 0:
        _warnings.warn(_dat_msg[status], UserWarning, 2)

    return deltat.value

_sofa.iauDtdb.argtypes = [
    _ct.c_double, # date1
    _ct.c_double, # date2
    _ct.c_double, # ut
    _ct.c_double, # elong
    _ct.c_double, # u
    _ct.c_double, # v
]

_sofa.iauDtdb.restype = _ct.c_double

def Dtdb(date1, date2, ut, elong, u, v):
    '''An approximation to TDB-TT, the difference between barycentric
    dynamical time and terrestrial time, for an observer on the Earth.

    Args:
        date1,date2 (double): date, TDB (Notes 1-3)
        ut (double): universal time (UT1, fraction of one day)
        elong (double): longitude (east positive, radians)
        u (double): distance from Earth spin axis (km)
        v (double): distance north of equatorial plane (km)

    Returns:
        tdb_tt (double)  TDB-TT (seconds)

    Notes:
        1. The date date1+date2 is a Julian Date, apportioned in any
        convenient way between the two arguments.  For example,
        JD(TT)=2450123.7 could be expressed in any of these ways,
        among others:
   
               date1          date2
   
            2450123.7           0.0       (JD method)
            2451545.0       -1421.3       (J2000 method)
            2400000.5       50123.2       (MJD method)
            2450123.5           0.2       (date & time method)
   
        The JD method is the most natural and convenient to use in
        cases where the loss of several decimal digits of resolution
        is acceptable.  The J2000 method is best matched to the way
        the argument is handled internally and will deliver the
        optimum resolution.  The MJD method and the date & time methods
        are both good compromises between resolution and convenience.
   
        Although the date is, formally, barycentric dynamical time (TDB),
        the terrestrial dynamical time (TT) can be used with no practical
        effect on the accuracy of the prediction.
   
        2. TT can be regarded as a coordinate time that is realized as an
        offset of 32.184s from International Atomic Time, TAI.  TT is a
        specific linear transformation of geocentric coordinate time TCG,
        which is the time scale for the Geocentric Celestial Reference
        System, GCRS.
   
        3. TDB is a coordinate time, and is a specific linear transformation
        of barycentric coordinate time TCB, which is the time scale for
        the Barycentric Celestial Reference System, BCRS.
   
        4. The difference TCG-TCB depends on the masses and positions of the
        bodies of the solar system and the velocity of the Earth.  It is
        dominated by a rate difference, the residual being of a periodic
        character.  The latter, which is modeled by the present function,
        comprises a main (annual) sinusoidal term of amplitude
        approximately 0.00166 seconds, plus planetary terms up to about
        20 microseconds, and lunar and diurnal terms up to 2 microseconds.
        These effects come from the changing transverse Doppler effect
        and gravitational red-shift as the observer (on the Earth's
        surface) experiences variations in speed (with respect to the
        BCRS) and gravitational potential.
   
        5. TDB can be regarded as the same as TCB but with a rate adjustment
        to keep it close to TT, which is convenient for many applications.
        The history of successive attempts to define TDB is set out in
        Resolution 3 adopted by the IAU General Assembly in 2006, which
        defines a fixed TDB(TCB) transformation that is consistent with
        contemporary solar-system ephemerides.  Future ephemerides will
        imply slightly changed transformations between TCG and TCB, which
        could introduce a linear drift between TDB and TT;  however, any
        such drift is unlikely to exceed 1 nanosecond per century.
   
        6. The geocentric TDB-TT model used in the present function is that of
        Fairhead & Bretagnon (1990), in its full form.  It was originally
        supplied by Fairhead (private communications with P.T.Wallace,
        1990) as a Fortran subroutine.  The present C function contains an
        adaptation of the Fairhead code.  The numerical results are
        essentially unaffected by the changes, the differences with
        respect to the Fairhead & Bretagnon original being at the 1e-20 s
        level.
   
        The topocentric part of the model is from Moyer (1981) and
        Murray (1983), with fundamental arguments adapted from
        Simon et al. 1994.  It is an approximation to the expression
        ( v / c ) . ( r / c ), where v is the barycentric velocity of
        the Earth, r is the geocentric position of the observer and
        c is the speed of light.
   
        By supplying zeroes for u and v, the topocentric part of the
        model can be nullified, and the function will return the Fairhead
        & Bretagnon result alone.
   
        7. During the interval 1950-2050, the absolute accuracy is better
        than +/- 3 nanoseconds relative to time ephemerides obtained by
        direct numerical integrations based on the JPL DE405 solar system
        ephemeris.
   
        8. It must be stressed that the present function is merely a model,
        and that numerical integration of solar-system ephemerides is the
        definitive method for predicting the relationship between TCG and
        TCB and hence between TT and TDB.
    '''

    return _sofa.iauDtdb(date1, date2, ut, elong, u, v)

_sofa.iauDtf2d.argtypes = [
    _ct.c_char_p, #scale
    _ct.c_int, #iy
    _ct.c_int, #im
    _ct.c_int, #id
    _ct.c_int, #ihr
    _ct.c_int, #imn
    _ct.c_double, #sec
    _ct.POINTER(_ct.c_double), #d1
    _ct.POINTER(_ct.c_double)  #d2
]

_sofa.iauDtf2d.restype = _ct.c_int

_dtf2d_msg = {
    3: 'time is after end of day (Note 5) and dubious year (Note 6)',
    2: 'time is after end of day (Note 5)',
    1: 'dubious year (Note 6)',
    -1: 'bad year',
    -2: 'bad month',
    -3: 'bad day',
    -4: 'bad hour',
    -5: 'bad minute',
    -6: 'bad second (<0)',
}

def Dtf2d(scale, iy, im, id, ihr, imn, sec):
    '''Encode date and time fields into 2-part Julian Date (or in the case
    of UTC a quasi-JD form that includes special provision for leap
    seconds).

    Args:
        scale (str): time scale ID (Note 1)
        iy,im,id (int): year, month, day in Gregorian calendar (Note 2)
        ihr,imn (int): hour, minute
        sec (double): seconds

    Returns:
        d1,d2 (double): 2-part Julian Date (Notes 3,4)

    Notes:
        1. scale identifies the time scale.  Only the value "UTC" (in upper
        case) is significant, and enables handling of leap seconds (see
        Note 4).
   
        2. For calendar conventions and limitations, see iauCal2jd.
   
        3. The sum of the results, d1+d2, is Julian Date, where normally d1
        is the Julian Day Number and d2 is the fraction of a day.  In the
        case of UTC, where the use of JD is problematical, special
        conventions apply:  see the next note.
   
        4. JD cannot unambiguously represent UTC during a leap second unless
        special measures are taken.  The SOFA internal convention is that
        the quasi-JD day represents UTC days whether the length is 86399,
        86400 or 86401 SI seconds.  In the 1960-1972 era there were
        smaller jumps (in either direction) each time the linear UTC(TAI)
        expression was changed, and these "mini-leaps" are also included
        in the SOFA convention.
   
        5. The warning status "time is after end of day" usually means that
        the sec argument is greater than 60.0.  However, in a day ending
        in a leap second the limit changes to 61.0 (or 59.0 in the case
        of a negative leap second).
   
        6. The warning status "dubious year" flags UTCs that predate the
        introduction of the time scale or that are too far in the future
        to be trusted.  See iauDat for further details.
   
        7. Only in the case of continuous and regular time scales (TAI, TT,
        TCG, TCB and TDB) is the result d1+d2 a Julian Date, strictly
        speaking.  In the other cases (UT1 and UTC) the result must be
        used with circumspection;  in particular the difference between
        two such results cannot be interpreted as a precise time
        interval.
    '''

    # Initialize return values
    d1 = _ct.c_double()
    d2 = _ct.c_double()

    # Main function call
    status = _sofa.iauDtf2d(bytes(scale, encoding='utf-8'), iy, im, id, ihr, imn, sec, _ct.byref(d1), _ct.byref(d2))

    # Raise warnings
    if status < 0:
        raise ValueError(_dtf2d_msg[status])
    elif status > 0:
        _warnings.warn(_dtf2d_msg[status], UserWarning, 2)

    return d1.value, d2.value

_sofa.iauTaitt.argtypes = [
    _ct.c_double, # tai1
    _ct.c_double, # tai2
    _ct.POINTER(_ct.c_double), #tt1
    _ct.POINTER(_ct.c_double), #tt2
]

_sofa.iauTaitt.restype = _ct.c_int

def Taitt(tai1, tai2):
    '''Time scale transformation:  International Atomic Time, TAI, to
    Terrestrial Time, TT.

    Args:
        tai1,tai2 (double): TAI as a 2-part Julian Date

    Returns:
        tt1,tt2 (double): TT as a 2-part Julian Date

    Notes:

        1. tai1+tai2 is Julian Date, apportioned in any convenient way
        between the two arguments, for example where tai1 is the Julian
        Day Number and tai2 is the fraction of a day.  The returned
        tt1,tt2 follow suit.
    '''

    # Initialize return pointers
    tt1 = _ct.c_double()
    tt2 = _ct.c_double()

    # Main function call
    status = _sofa.iauTaitt(tai1, tai2, _ct.byref(tt1), _ct.byref(tt2))

    return tt1.value, tt2.value

_sofa.iauTaiut1.argtypes = [
    _ct.c_double, # tai1
    _ct.c_double, # tai2
    _ct.c_double, # dta
    _ct.POINTER(_ct.c_double), #ut11
    _ct.POINTER(_ct.c_double), #ut12
]

_sofa.iauTaiut1.restype = _ct.c_int

def Taiut1(tai1, tai2, dta):
    '''Time scale transformation:  International Atomic Time, TAI, to
    Universal Time, UT1.

    Args:
        tai1,tai2 (double): TAI as a 2-part Julian Date
        dta (double): UT1-TAI in seconds

    Returns:
        ut11,ut12 (double): UT1 as a 2-part Julian Date

    Notes:

        1. tai1+tai2 is Julian Date, apportioned in any convenient way
        between the two arguments, for example where tai1 is the Julian
        Day Number and tai2 is the fraction of a day.  The returned
        UT11,UT12 follow suit.
   
        2. The argument dta, i.e. UT1-TAI, is an observed quantity, and is
        available from IERS tabulations.
    '''

    # Initialize return pointers
    ut11 = _ct.c_double()
    ut12 = _ct.c_double()

    # Main function call
    status = _sofa.iauTaiut1(tai1, tai2, dta, _ct.byref(ut11), _ct.byref(ut12))

    return ut11.value, ut12.value

_sofa.iauTaiutc.argtypes = [
    _ct.c_double, # tai1
    _ct.c_double, # tai2
    _ct.POINTER(_ct.c_double), #utc1
    _ct.POINTER(_ct.c_double), #utc2
]

_sofa.iauTaiutc.restype = _ct.c_int

_taiutc_msg = {
    1: 'dubious year (Note 4)',
    -1: 'unacceptable date'
}

def Taiutc(tai1, tai2):
    '''Time scale transformation:  International Atomic Time, TAI, to
    Coordinated Universal Time, UTC.

    Args:
        tai1,tai2 (double): TAI as a 2-part Julian Date (Note 1)

    Returns:
        utc1,utc2 (double): UTC as a 2-part quasi Julian Date (Notes 1-3)

    Notes:

        1. tai1+tai2 is Julian Date, apportioned in any convenient way
        between the two arguments, for example where tai1 is the Julian
        Day Number and tai2 is the fraction of a day.  The returned utc1
        and utc2 form an analogous pair, except that a special convention
        is used, to deal with the problem of leap seconds - see the next
        note.
   
        2. JD cannot unambiguously represent UTC during a leap second unless
        special measures are taken.  The convention in the present
        function is that the JD day represents UTC days whether the
        length is 86399, 86400 or 86401 SI seconds.  In the 1960-1972 era
        there were smaller jumps (in either direction) each time the
        linear UTC(TAI) expression was changed, and these "mini-leaps"
        are also included in the SOFA convention.
   
        3. The function iauD2dtf can be used to transform the UTC quasi-JD
        into calendar date and clock time, including UTC leap second
        handling.
   
        4. The warning status "dubious year" flags UTCs that predate the
        introduction of the time scale or that are too far in the future
        to be trusted.  See iauDat for further details.
    '''

    # Initialize return pointers
    utc1 = _ct.c_double()
    utc2 = _ct.c_double()

    # Main function call
    status = _sofa.iauTaiutc(tai1, tai2, _ct.byref(utc1), _ct.byref(utc2))

    # Raise warnings
    if status < 0:
        raise ValueError(_taiutc_msg[status])
    elif status > 0:
        _warnings.warn(_taiutc_msg[status], UserWarning, 2)

    return utc1.value, utc2.value

_sofa.iauTcbtdb.argtypes = [
    _ct.c_double, # tcb1
    _ct.c_double, # tcb2
    _ct.POINTER(_ct.c_double), #tdb1
    _ct.POINTER(_ct.c_double), #tdb2
]

_sofa.iauTcbtdb.restype = _ct.c_int

def Tcbtdb(tcb1, tcb2):
    '''Time scale transformation:  Barycentric Coordinate Time, TCB, to
    Barycentric Dynamical Time, TDB.

    Args:
        tcb1,tcb2 (double): TCB as a 2-part Julian Date

    Returns:
        tdb1,tdb2 (double): TDB as a 2-part Julian Date

    Notes:

        1. tcb1+tcb2 is Julian Date, apportioned in any convenient way
        between the two arguments, for example where tcb1 is the Julian
        Day Number and tcb2 is the fraction of a day.  The returned
        tdb1,tdb2 follow suit.
   
        2. The 2006 IAU General Assembly introduced a conventional linear
        transformation between TDB and TCB.  This transformation
        compensates for the drift between TCB and terrestrial time TT,
        and keeps TDB approximately centered on TT.  Because the
        relationship between TT and TCB depends on the adopted solar
        system ephemeris, the degree of alignment between TDB and TT over
        long intervals will vary according to which ephemeris is used.
        Former definitions of TDB attempted to avoid this problem by
        stipulating that TDB and TT should differ only by periodic
        effects.  This is a good description of the nature of the
        relationship but eluded precise mathematical formulation.  The
        conventional linear relationship adopted in 2006 sidestepped
        these difficulties whilst delivering a TDB that in practice was
        consistent with values before that date.
   
        3. TDB is essentially the same as Teph, the time argument for the
        JPL solar system ephemerides.
    '''

    # Initialize return pointers
    tdb1 = _ct.c_double()
    tdb2 = _ct.c_double()

    # Main function call
    status = _sofa.iauTcbtdb(tcb1, tcb2, _ct.byref(tdb1), _ct.byref(tdb2))

    return tdb1.value, tdb2.value

_sofa.iauTcgtt.argtypes = [
    _ct.c_double, # tcg1
    _ct.c_double, # tcg2
    _ct.POINTER(_ct.c_double), #tt1
    _ct.POINTER(_ct.c_double), #tt2
]

_sofa.iauTcgtt.restype = _ct.c_int

def Tcgtt(tcg1, tcg2):
    '''Time scale transformation:  Geocentric Coordinate Time, TCG, to
    Terrestrial Time, TT.

    Args:
        tcg1,tc2 (double): TCG as a 2-part Julian Date

    Returns:
        tt1,tt2 (double): TT as a 2-part Julian Date

    Notes:

        1. tcg1+tcg2 is Julian Date, apportioned in any convenient way
        between the two arguments, for example where tcg1 is the Julian
        Day Number and tcg22 is the fraction of a day.  The returned
        tt1,tt2 follow suit.
    '''

    # Initialize return pointers
    tt1 = _ct.c_double()
    tt2 = _ct.c_double()

    # Main function call
    status = _sofa.iauTcgtt(tcg1, tcg2, _ct.byref(tt1), _ct.byref(tt2))

    return tt1.value, tt2.value

_sofa.iauTdbtcb.argtypes = [
    _ct.c_double, # tdb1
    _ct.c_double, # tdb2
    _ct.POINTER(_ct.c_double), #tcb1
    _ct.POINTER(_ct.c_double), #tcb2
]

_sofa.iauTdbtcb.restype = _ct.c_int

def Tdbtcb(tdb1, tdb2):
    '''Time scale transformation:  Barycentric Dynamical Time, TDB, to
    Barycentric Coordinate Time, TCB.

    Args:
        tdb1,tdb2 (double): TDB as a 2-part Julian Date

    Returns:
        tcb1,tcb2 (double): TCB as a 2-part Julian Date

    Notes:

        1. tdb1+tdb2 is Julian Date, apportioned in any convenient way
        between the two arguments, for example where tdb1 is the Julian
        Day Number and tdb2 is the fraction of a day.  The returned
        tcb1,tcb2 follow suit.
   
        2. The 2006 IAU General Assembly introduced a conventional linear
        transformation between TDB and TCB.  This transformation
        compensates for the drift between TCB and terrestrial time TT,
        and keeps TDB approximately centered on TT.  Because the
        relationship between TT and TCB depends on the adopted solar
        system ephemeris, the degree of alignment between TDB and TT over
        long intervals will vary according to which ephemeris is used.
        Former definitions of TDB attempted to avoid this problem by
        stipulating that TDB and TT should differ only by periodic
        effects.  This is a good description of the nature of the
        relationship but eluded precise mathematical formulation.  The
        conventional linear relationship adopted in 2006 sidestepped
        these difficulties whilst delivering a TDB that in practice was
        consistent with values before that date.
   
        3. TDB is essentially the same as Teph, the time argument for the
        JPL solar system ephemerides.
    '''

    # Initialize return pointers
    tcb1 = _ct.c_double()
    tcb2 = _ct.c_double()

    # Main function call
    status = _sofa.iauTdbtcb(tdb1, tdb2, _ct.byref(tcb1), _ct.byref(tcb2))

    return tcb1.value, tcb2.value

_sofa.iauTdbtt.argtypes = [
    _ct.c_double, # tdb1
    _ct.c_double, # tdb2
    _ct.c_double, # dtr
    _ct.POINTER(_ct.c_double), #tt1
    _ct.POINTER(_ct.c_double), #tt2
]

_sofa.iauTdbtt.restype = _ct.c_int

def Tdbtt(tdb1, tdb2, dtr):
    '''Time scale transformation:  Barycentric Dynamical Time, TDB, to
    Terrestrial Time, TT.

    Args:
        tdb1,tdb2 (double): TDB as a 2-part Julian Date
        dtr (double): TDB-TT in seconds

    Returns:
        tt1,tt2 (double): TT as a 2-part Julian Date

    Notes:

        1. tdb1+tdb2 is Julian Date, apportioned in any convenient way
        between the two arguments, for example where tdb1 is the Julian
        Day Number and tdb2 is the fraction of a day.  The returned
        tt1,tt2 follow suit.
   
        2. The argument dtr represents the quasi-periodic component of the
        GR transformation between TT and TCB.  It is dependent upon the
        adopted solar-system ephemeris, and can be obtained by numerical
        integration, by interrogating a precomputed time ephemeris or by
        evaluating a model such as that implemented in the SOFA function
        iauDtdb.   The quantity is dominated by an annual term of 1.7 ms
        amplitude.
   
        3. TDB is essentially the same as Teph, the time argument for the
        JPL solar system ephemerides.
    '''

    # Initialize return pointers
    tt1 = _ct.c_double()
    tt2 = _ct.c_double()

    # Main function call
    status = _sofa.iauTdbtt(tdb1, tdb2, dtr, _ct.byref(tt1), _ct.byref(tt2))

    return tt1.value, tt2.value

_sofa.iauTttai.argtypes = [
    _ct.c_double, # tt1
    _ct.c_double, # tt2
    _ct.POINTER(_ct.c_double), #tai1
    _ct.POINTER(_ct.c_double), #tai2
]

_sofa.iauTttai.restype = _ct.c_int

def Tttai(tt1, tt2):
    '''Time scale transformation:  Terrestrial Time, TT, to International
    Atomic Time, TAI.

    Args:
        tdb1,tdb2 (double): TDB as a 2-part Julian Date

    Returns:
        tt1,tt2 (double): TT as a 2-part Julian Date

    Notes:

        1. tt1+tt2 is Julian Date, apportioned in any convenient way between
        the two arguments, for example where tt1 is the Julian Day Number
        and tt2 is the fraction of a day.  The returned tai1,tai2 follow
        suit.
    '''

    # Initialize return pointers
    tai1 = _ct.c_double()
    tai2 = _ct.c_double()

    # Main function call
    status = _sofa.iauTttai(tt1, tt2, _ct.byref(tai1), _ct.byref(tai2))

    return tai1.value, tai2.value

_sofa.iauTttcg.argtypes = [
    _ct.c_double, # tt1
    _ct.c_double, # tt2
    _ct.POINTER(_ct.c_double), #tcg1
    _ct.POINTER(_ct.c_double), #tcg2
]

_sofa.iauTttcg.restype = _ct.c_int

def Tttcg(tt1, tt2):
    '''Time scale transformation:  Terrestrial Time, TT, to Geocentric
    Coordinate Time, TCG.

    Args:
        tt1,tt2 (double): TT as a 2-part Julian Date

    Returns:
        tcg1,tcg2 (double): TCG as a 2-part Julian Date

    Notes:

        1. tt1+tt2 is Julian Date, apportioned in any convenient way between
        the two arguments, for example where tt1 is the Julian Day Number
        and tt2 is the fraction of a day.  The returned tcg1,tcg2 follow
        suit.
    '''

    # Initialize return pointers
    tcg1 = _ct.c_double()
    tcg2 = _ct.c_double()

    # Main function call
    status = _sofa.iauTttcg(tt1, tt2, _ct.byref(tcg1), _ct.byref(tcg2))

    return tcg1.value, tcg2.value

_sofa.iauTttdb.argtypes = [
    _ct.c_double, # tt1
    _ct.c_double, # tt2
    _ct.c_double, # dtr
    _ct.POINTER(_ct.c_double), #tdb1
    _ct.POINTER(_ct.c_double), #tdb2
]

_sofa.iauTttdb.restype = _ct.c_int

def Tttdb(tt1, tt2, dtr):
    '''Time scale transformation:  Terrestrial Time, TT, to Barycentric
    Dynamical Time, TDB.

    Args:
        tdb1,tdb2 (double): TT as a 2-part Julian Date
        dtr (double): TDB-TT in seconds

    Returns:
        tdb1,tdb2 (double): TDB as a 2-part Julian Date

    Notes:

        1. tt1+tt2 is Julian Date, apportioned in any convenient way between
        the two arguments, for example where tt1 is the Julian Day Number
        and tt2 is the fraction of a day.  The returned tdb1,tdb2 follow
        suit.
   
        2. The argument dtr represents the quasi-periodic component of the
        GR transformation between TT and TCB.  It is dependent upon the
        adopted solar-system ephemeris, and can be obtained by numerical
        integration, by interrogating a precomputed time ephemeris or by
        evaluating a model such as that implemented in the SOFA function
        iauDtdb.   The quantity is dominated by an annual term of 1.7 ms
        amplitude.
   
        3. TDB is essentially the same as Teph, the time argument for the JPL
        solar system ephemerides.
    '''

    # Initialize return pointers
    tdb1 = _ct.c_double()
    tdb2 = _ct.c_double()

    # Main function call
    status = _sofa.iauTttdb(tt1, tt2, dtr, _ct.byref(tdb1), _ct.byref(tdb2))

    return tdb1.value, tdb2.value

_sofa.iauTtut1.argtypes = [
    _ct.c_double, # tt1
    _ct.c_double, # tt2
    _ct.c_double, # dt
    _ct.POINTER(_ct.c_double), #ut11
    _ct.POINTER(_ct.c_double), #ut12
]

_sofa.iauTtut1.restype = _ct.c_int

def Ttut1(tt1, tt2, dt):
    '''Time scale transformation:  Terrestrial Time, TT, to Universal Time, UT1

    Args:
        tdb1,tdb2 (double): TDB as a 2-part Julian Date
        dt (double): TT-UT1 in seconds

    Returns:
        ut11,ut12 (double): TT as a 2-part Julian Date

    Notes:

        1. tt1+tt2 is Julian Date, apportioned in any convenient way between
        the two arguments, for example where tt1 is the Julian Day Number
        and tt2 is the fraction of a day.  The returned ut11,ut12 follow
        suit.
   
        2. The argument dt is classical Delta T.
    '''

    # Initialize return pointers
    ut11 = _ct.c_double()
    ut12 = _ct.c_double()

    # Main function call
    status = _sofa.iauTtut1(tt1, tt2, dt, _ct.byref(ut11), _ct.byref(ut12))

    return ut11.value, ut12.value

_sofa.iauUt1tai.argtypes = [
    _ct.c_double, # ut11
    _ct.c_double, # ut12
    _ct.c_double, # dta
    _ct.POINTER(_ct.c_double), #tai1
    _ct.POINTER(_ct.c_double), #tai2
]

_sofa.iauUt1tai.restype = _ct.c_int

def Ut1tai(ut11, ut12, dta):
    '''Time scale transformation:  Universal Time, UT1, to International Atomic Time, TAI.

    Args:
        tdb1,tdb2 (double): UT1 as a 2-part Julian Date
        dta (double): UT1-TAI in seconds

    Returns:
        tai1,tai2 (double): TAI as a 2-part Julian Date

    Notes:

        1. ut11+ut12 is Julian Date, apportioned in any convenient way
        between the two arguments, for example where ut11 is the Julian
        Day Number and ut12 is the fraction of a day.  The returned
        tai1,tai2 follow suit.
   
        2. The argument dta, i.e. UT1-TAI, is an observed quantity, and is
        available from IERS tabulations.
    '''

    # Initialize return pointers
    tai1 = _ct.c_double()
    tai2 = _ct.c_double()

    # Main function call
    status = _sofa.iauUt1tai(ut11, ut12, dta, _ct.byref(tai1), _ct.byref(tai2))

    return tai1.value, tai2.value

_sofa.iauUt1tt.argtypes = [
    _ct.c_double, # ut11
    _ct.c_double, # ut12
    _ct.c_double, # dt
    _ct.POINTER(_ct.c_double), #tt1
    _ct.POINTER(_ct.c_double), #tt2
]

_sofa.iauUt1tt.restype = _ct.c_int

def Ut1tt(ut11, ut12, dt):
    '''Time scale transformation:  Universal Time, UT1, to Terrestrial Time, TT.

    Args:
        tdb1,tdb2 (double): UT1 as a 2-part Julian Date
        dt (double): TT-UT1 in seconds

    Returns:
        tai1,tai2 (double): TT as a 2-part Julian Date

    Notes:

        1. ut11+ut12 is Julian Date, apportioned in any convenient way
        between the two arguments, for example where ut11 is the Julian
        Day Number and ut12 is the fraction of a day.  The returned
        tt1,tt2 follow suit.
   
        2. The argument dt is classical Delta T.
    '''

    # Initialize return pointers
    tt1 = _ct.c_double()
    tt2 = _ct.c_double()

    # Main function call
    status = _sofa.iauUt1tt(ut11, ut12, dt, _ct.byref(tt1), _ct.byref(tt2))

    return tt1.value, tt2.value

_sofa.iauUt1utc.argtypes = [
    _ct.c_double, # ut11
    _ct.c_double, # ut12
    _ct.c_double, # dut1
    _ct.POINTER(_ct.c_double), #utc1
    _ct.POINTER(_ct.c_double), #utc2
]

_sofa.iauUt1utc.restype = _ct.c_int

_ut1utc_msg = {
    1: 'dubious year (Note 5)',
    -1: 'unacceptable date'
}

def Ut1utc(ut11, ut12, dta):
    '''Time scale transformation:  Universal Time, UT1, to Terrestrial Time, TT.

    Args:
        ut11,ut12 (double): UT1 as a 2-part Julian Date (Note 1)
        dut1 (double): Delta UT1: UT1-UTC in seconds (Note 2)

    Returns:
        utc1,utc2 (double): UTC as a 2-part quasi Julian Date (Notes 3,4)

    Notes:

        1. ut11+ut12 is Julian Date, apportioned in any convenient way
        between the two arguments, for example where ut11 is the Julian
        Day Number and ut12 is the fraction of a day.  The returned utc1
        and utc2 form an analogous pair, except that a special convention
        is used, to deal with the problem of leap seconds - see Note 3.
   
        2. Delta UT1 can be obtained from tabulations provided by the
        International Earth Rotation and Reference Systems Service.  The
        value changes abruptly by 1s at a leap second;  however, close to
        a leap second the algorithm used here is tolerant of the "wrong"
        choice of value being made.
   
        3. JD cannot unambiguously represent UTC during a leap second unless
        special measures are taken.  The convention in the present
        function is that the returned quasi JD day UTC1+UTC2 represents
        UTC days whether the length is 86399, 86400 or 86401 SI seconds.
   
        4. The function iauD2dtf can be used to transform the UTC quasi-JD
        into calendar date and clock time, including UTC leap second
        handling.
   
        5. The warning status "dubious year" flags UTCs that predate the
        introduction of the time scale or that are too far in the future
        to be trusted.  See iauDat for further details.
    '''

    # Initialize return pointers
    utc1 = _ct.c_double()
    utc2 = _ct.c_double()

    # Main function call
    status = _sofa.iauUt1utc(ut11, ut12, dta, _ct.byref(utc1), _ct.byref(utc2))

    # Raise warnings
    if status < 0:
        raise ValueError(_ut1utc_msg[status])
    elif status > 0:
        _warnings.warn(_ut1utc_msg[status], UserWarning, 2)

    return utc1.value, utc2.value

_sofa.iauUtctai.argtypes = [
    _ct.c_double, # ut11
    _ct.c_double, # ut12
    _ct.POINTER(_ct.c_double), #utc1
    _ct.POINTER(_ct.c_double), #utc2
]

_sofa.iauUtctai.restype = _ct.c_int

_utctai_msg = {
    1: 'dubious year (Note 3)',
    -1: 'unacceptable date'
}

def Utctai(utc1, utc2):
    '''Time scale transformation:  Coordinated Universal Time, UTC, to
    International Atomic Time, TAI.

    Args:
        utc1,utc2 (double): UTC as a 2-part quasi Julian Date (Notes 1-4)

    Returns:
        tai1,tai2 (double): TAI as a 2-part Julian Date (Note 5)

    Notes:

        1. utc1+utc2 is quasi Julian Date (see Note 2), apportioned in any
        convenient way between the two arguments, for example where utc1
        is the Julian Day Number and utc2 is the fraction of a day.
   
        2. JD cannot unambiguously represent UTC during a leap second unless
        special measures are taken.  The convention in the present
        function is that the JD day represents UTC days whether the
        length is 86399, 86400 or 86401 SI seconds.  In the 1960-1972 era
        there were smaller jumps (in either direction) each time the
        linear UTC(TAI) expression was changed, and these "mini-leaps"
        are also included in the SOFA convention.
   
        3. The warning status "dubious year" flags UTCs that predate the
        introduction of the time scale or that are too far in the future
        to be trusted.  See iauDat for further details.
   
        4. The function iauDtf2d converts from calendar date and time of day
        into 2-part Julian Date, and in the case of UTC implements the
        leap-second-ambiguity convention described above.
   
        5. The returned TAI1,TAI2 are such that their sum is the TAI Julian
        Date.
    '''

    # Initialize return pointers
    tai1 = _ct.c_double()
    tai2 = _ct.c_double()

    # Main function call
    status = _sofa.iauUtctai(utc1, utc2, _ct.byref(tai1), _ct.byref(tai2))

    # Raise warnings
    if status < 0:
        raise ValueError(_utctai_msg[status])
    elif status > 0:
        _warnings.warn(_utctai_msg[status], UserWarning, 2)

    return tai1.value, tai2.value

_sofa.iauUtcut1.argtypes = [
    _ct.c_double, # ut11
    _ct.c_double, # ut12
    _ct.c_double, # dut1
    _ct.POINTER(_ct.c_double), #utc1
    _ct.POINTER(_ct.c_double), #utc2
]

_sofa.iauUtcut1.restype = _ct.c_int

_ut1utc_msg = {
    1: 'dubious year (Note 3)',
    -1: 'unacceptable date'
}

def Utcut1(utc1, utc2, dut1):
    '''Time scale transformation:  Coordinated Universal Time, UTC, to Universal Time, UT1.

    Args:
        utc1,utc2 (double): UTC as a 2-part quasi Julian Date (Notes 1-4)
        dut1 (double): Delta UT1 = UT1-UTC in seconds (Note 5)

    Returns:
        utc1,utc2 (double): UT1 as a 2-part Julian Date (Note 6)

    Notes:

        1. utc1+utc2 is quasi Julian Date (see Note 2), apportioned in any
        convenient way between the two arguments, for example where utc1
        is the Julian Day Number and utc2 is the fraction of a day.
   
        2. JD cannot unambiguously represent UTC during a leap second unless
        special measures are taken.  The convention in the present
        function is that the JD day represents UTC days whether the
        length is 86399, 86400 or 86401 SI seconds.
   
        3. The warning status "dubious year" flags UTCs that predate the
        introduction of the time scale or that are too far in the future
        to be trusted.  See iauDat for further details.
   
        4. The function iauDtf2d converts from calendar date and time of
        day into 2-part Julian Date, and in the case of UTC implements
        the leap-second-ambiguity convention described above.
   
        5. Delta UT1 can be obtained from tabulations provided by the
        International Earth Rotation and Reference Systems Service.
        It is the caller's responsibility to supply a dut1 argument
        containing the UT1-UTC value that matches the given UTC.
   
        6. The returned ut11,ut12 are such that their sum is the UT1 Julian
        Date.
    '''

    # Initialize return pointers
    ut11 = _ct.c_double()
    ut12 = _ct.c_double()

    # Main function call
    status = _sofa.iauUtcut1(utc1, utc2, dut1, _ct.byref(ut11), _ct.byref(ut12))

    # Raise warnings
    if status < 0:
        raise ValueError(_ut1utc_msg[status])
    elif status > 0:
        _warnings.warn(_ut1utc_msg[status], UserWarning, 2)

    return ut11.value, ut12.value

# /* Astronomy/HorizonEquatorial */
# void iauAe2hd(double az, double el, double phi,
#               double *ha, double *dec);
# void iauHd2ae(double ha, double dec, double phi,
#               double *az, double *el);
# double iauHd2pa(double ha, double dec, double phi);

# /* Astronomy/Gnomonic */
# int iauTpors(double xi, double eta, double a, double b,
#              double *a01, double *b01, double *a02, double *b02);
# int iauTporv(double xi, double eta, double v[3],
#              double v01[3], double v02[3]);
# void iauTpsts(double xi, double eta, double a0, double b0,
#               double *a, double *b);
# void iauTpstv(double xi, double eta, double v0[3], double v[3]);
# int iauTpxes(double a, double b, double a0, double b0,
#              double *xi, double *eta);
# int iauTpxev(double v[3], double v0[3], double *xi, double *eta);

# /* VectorMatrix/AngleOps */
# void iauA2af(int ndp, double angle, char *sign, int idmsf[4]);
# void iauA2tf(int ndp, double angle, char *sign, int ihmsf[4]);
# int iauAf2a(char s, int ideg, int iamin, double asec, double *rad);
# double iauAnp(double a);
# double iauAnpm(double a);
# void iauD2tf(int ndp, double days, char *sign, int ihmsf[4]);
# int iauTf2a(char s, int ihour, int imin, double sec, double *rad);
# int iauTf2d(char s, int ihour, int imin, double sec, double *days);

# /* VectorMatrix/BuildRotations */

_sofa.iauRx.argtypes = [
    _ct.c_double, #phi
    _ndpointer(shape=(3,3), dtype=float, flags='C') #r
] 

def Rx(phi, r):
    '''Rotate a r-matrix about the x-axis.
    '''
    r2 = _req_shape_c(r, float, (3,3)).copy()
    _sofa.iauRx(float(phi), r2)
    return _np.array(r2, dtype=float)

_sofa.iauRy.argtypes = [
    _ct.c_double, #theta
    _ndpointer(shape=(3,3), dtype=float, flags='C') #r
] 

def Ry(theta, r):
    ''' Rotate a r-matrix about the y-axis.
    '''
    r2 = _req_shape_c(r, float, (3,3)).copy()

    # Main funciton call
    _sofa.iauRy(float(theta), r2)

    return _np.array(r2, dtype=float)


_sofa.iauRz.argtypes = [
    _ct.c_double, #psi
    _ndpointer(shape=(3,3), dtype=float, flags='C') #r
]

def Rz(psi, r):
    '''Rotate a r-matrix about the z-axis.
    '''

    r2 = _req_shape_c(r, float, (3,3)).copy()
    
    # Main function call
    _sofa.iauRz(float(psi), r2)

    return _np.array(r2, dtype=float)

# /* VectorMatrix/CopyExtendExtract */
# void iauCp(double p[3], double c[3]);
# void iauCpv(double pv[2][3], double c[2][3]);
# void iauCr(double r[3][3], double c[3][3]);
# void iauP2pv(double p[3], double pv[2][3]);
# void iauPv2p(double pv[2][3], double p[3]);

# /* VectorMatrix/Initialization */
# void iauIr(double r[3][3]);
# void iauZp(double p[3]);
# void iauZpv(double pv[2][3]);
# void iauZr(double r[3][3]);

# /* VectorMatrix/MatrixOps */
# void iauRxr(double a[3][3], double b[3][3], double atb[3][3]);
# void iauTr(double r[3][3], double rt[3][3]);

# /* VectorMatrix/MatrixVectorProducts */
# void iauRxp(double r[3][3], double p[3], double rp[3]);
# void iauRxpv(double r[3][3], double pv[2][3], double rpv[2][3]);
# void iauTrxp(double r[3][3], double p[3], double trp[3]);
# void iauTrxpv(double r[3][3], double pv[2][3], double trpv[2][3]);

# /* VectorMatrix/RotationVectors */
# void iauRm2v(double r[3][3], double w[3]);
# void iauRv2m(double w[3], double r[3][3]);

# /* VectorMatrix/SeparationAndAngle */
# double iauPap(double a[3], double b[3]);
# double iauPas(double al, double ap, double bl, double bp);
# double iauSepp(double a[3], double b[3]);
# double iauSeps(double al, double ap, double bl, double bp);

# /* VectorMatrix/SphericalCartesian */
# void iauC2s(double p[3], double *theta, double *phi);
# void iauP2s(double p[3], double *theta, double *phi, double *r);
# void iauPv2s(double pv[2][3],
#              double *theta, double *phi, double *r,
#              double *td, double *pd, double *rd);
# void iauS2c(double theta, double phi, double c[3]);
# void iauS2p(double theta, double phi, double r, double p[3]);
# void iauS2pv(double theta, double phi, double r,
#              double td, double pd, double rd,
#              double pv[2][3]);

# /* VectorMatrix/VectorOps */
# double iauPdp(double a[3], double b[3]);
# double iauPm(double p[3]);
# void iauPmp(double a[3], double b[3], double amb[3]);
# void iauPn(double p[3], double *r, double u[3]);
# void iauPpp(double a[3], double b[3], double apb[3]);
# void iauPpsp(double a[3], double s, double b[3], double apsb[3]);
# void iauPvdpv(double a[2][3], double b[2][3], double adb[2]);
# void iauPvm(double pv[2][3], double *r, double *s);
# void iauPvmpv(double a[2][3], double b[2][3], double amb[2][3]);
# void iauPvppv(double a[2][3], double b[2][3], double apb[2][3]);
# void iauPvu(double dt, double pv[2][3], double upv[2][3]);
# void iauPvup(double dt, double pv[2][3], double p[3]);
# void iauPvxpv(double a[2][3], double b[2][3], double axb[2][3]);
# void iauPxp(double a[3], double b[3], double axb[3]);
# void iauS2xpv(double s1, double s2, double pv[2][3], double spv[2][3]);
# void iauSxp(double s, double p[3], double sp[3]);
# void iauSxpv(double s, double pv[2][3], double spv[2][3]);