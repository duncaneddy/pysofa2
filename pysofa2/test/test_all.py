from pytest import approx

import numpy as _np
from pysofa2 import *

def test_Cal2jd():
    djm0, djm = Cal2jd(2003, 6, 1)

    assert djm0 == 2400000.5
    assert djm  == 52791.0

def test_Epb():
    epb = Epb(2415019.8135, 30103.18648)

    assert epb == approx(1982.418424159278580, 1e-12)

def test_Epb2jd():
    epb = 1957.3

    djm0, djm = Epb2jd(epb)

    assert djm0 == approx(2400000.5, 1e-9)
    assert djm  == approx(35948.1915101513, 1e-9)

def test_Epj():
    epj = Epj(2451545, -7392.5)

    assert epj == approx(1979.760438056125941, 1e-12)

def test_Epj2jd():
    epj = 1996.8

    djm0, djm = Epj2jd(epj)

    assert djm0 == approx(2400000.5, 1e-9)
    assert djm  == approx(50375.7, 1e-9)

def test_jd2cal():
    dj1 = 2400000.5
    dj2 = 50123.9999

    iy, im, id, fd = Jd2cal(dj1, dj2)

    assert iy == 1996
    assert im == 2
    assert id == 10
    assert fd == approx(0.9999, 1e-7)

def test_jd2calf():
    dj1 = 2400000.5
    dj2 = 50123.9999

    iy, im, id, fd = Jd2calf(4, dj1, dj2)

    assert iy == 1996
    assert im == 2
    assert id == 10
    assert fd == 9999

def test_fad03():
    assert Fad03(0.80) == approx(1.946709205396925672, 1e-12) 

def test_fae03():
    assert Fae03(0.80) == approx(1.744713738913081846, 1e-12) 

def test_faf03():
    assert Faf03(0.80) == approx(0.2597711366745499518, 1e-12) 

def test_faju03():
    assert Faju03(0.80) == approx(5.275711665202481138, 1e-12) 

def test_fal03():
    assert Fal03(0.80) == approx(5.132369751108684150, 1e-12) 

def test_falp03():
    assert Falp03(0.80) == approx(6.226797973505507345, 1e-12) 

def test_fama03():
    assert Fama03(0.80) == approx(3.275506840277781492, 1e-12) 

def test_fame03():
    assert Fame03(0.80) == approx(5.417338184297289661, 1e-12) 

def test_fane03():
    assert Fane03(0.80) == approx(2.079343830860413523, 1e-12) 

def test_faom03():
    assert Faom03(0.80) == approx(-5.973618440951302183, 1e-12)

def test_fapa03():
    assert Fapa03(0.80) == approx(0.1950884762240000000e-1, 1e-12) 

def test_fasa03():
    assert Fasa03(0.80) == approx(5.371574539440827046, 1e-12) 

def test_faur03():
    assert Faur03(0.80) == approx(5.180636450180413523, 1e-12) 

def test_fave03():
    assert Fave03(0.80) == approx(3.424900460533758000, 1e-12) 

def test_c2ixys():
    x =  0.5791308486706011000e-3
    y =  0.4020579816732961219e-4
    s = -0.1220040848472271978e-7

    rc2i = C2ixys(x, y, s)

    assert rc2i[0, 0] == approx(0.9999998323037157138, 1e-12)
    assert rc2i[0, 1] == approx(0.5581984869168499149e-9, 1e-12)
    assert rc2i[0, 2] == approx(-0.5791308491611282180e-3, 1e-12)

    assert rc2i[1, 0] == approx(-0.2384261642670440317e-7, 1e-12)
    assert rc2i[1, 1] == approx(0.9999999991917468964, 1e-12)
    assert rc2i[1, 2] == approx(-0.4020579110169668931e-4, 1e-12)

    assert rc2i[2, 0] == approx(0.5791308486706011000e-3, 1e-12)
    assert rc2i[2, 1] == approx(0.4020579816732961219e-4, 1e-12)
    assert rc2i[2, 2] == approx(0.9999998314954627590, 1e-12)

def test_pom00():
    xp =  2.55060238e-7
    yp =  1.860359247e-6
    sp = -0.1367174580728891460e-10

    rpom = Pom00(xp, yp, sp)

    assert rpom[0, 0] == approx(0.9999999999999674721, 1e-12)
    assert rpom[0, 1] == approx(-0.1367174580728846989e-10, 1e-16)
    assert rpom[0, 2] == approx(0.2550602379999972345e-6, 1e-16)

    assert rpom[1, 0] == approx(0.1414624947957029801e-10, 1e-16)
    assert rpom[1, 1] == approx(0.9999999999982695317, 1e-12)
    assert rpom[1, 2] == approx(-0.1860359246998866389e-5, 1e-16)

    assert rpom[2, 0] == approx(-0.2550602379741215021e-6, 1e-16)
    assert rpom[2, 1] == approx(0.1860359247002414021e-5, 1e-16)
    assert rpom[2, 2] == approx(0.9999999999982370039, 1e-12)

def test_s00():
    x = 0.5791308486706011000e-3
    y = 0.4020579816732961219e-4

    s = S00(2400000.5, 53736.0, x, y)

    assert s == approx(-0.1220036263270905693e-7, 1e-18)

def test_s00a():
    s = S00a(2400000.5, 52541.0)

    assert s == approx(-0.1340684448919163584e-7, 1e-18)

def test_s00b():
    s = S00b(2400000.5, 52541.0)

    assert s == approx(-0.1340695782951026584e-7, 1e-18)

def test_s06():
    x = 0.5791308486706011000e-3
    y = 0.4020579816732961219e-4

    s = S00(2400000.5, 53736.0, x, y)

    assert s == approx(-0.1220032213076463117e-7, 1e-18)

def test_s06a():
    s = S06a(2400000.5, 52541.0)

    assert s == approx(-0.1340680437291812383e-7, 1e-18)

def test_sp00():
    assert Sp00(2400000.5, 52541.0) == approx(-0.6216698469981019309e-11, 1e-12)

def test_xy06():
    x, y = Xy06(2400000.5, 53736.0)

    assert x == approx(0.5791308486706010975e-3, 1e-15)
    assert y == approx(0.4020579816732958141e-4, 1e-16)

def test_xys00a():
    x, y, s = Xys00a(2400000.5, 53736.0)

    assert x == approx( 0.5791308472168152904e-3, 1e-14)
    assert y == approx( 0.4020595661591500259e-4, 1e-15)
    assert s == approx(-0.1220040848471549623e-7, 1e-18)

def test_xys00b():
    x, y, s = Xys00b(2400000.5, 53736.0)

    assert x == approx( 0.5791301929950208873e-3, 1e-14)
    assert y == approx( 0.4020553681373720832e-4, 1e-15)
    assert s == approx(-0.1220027377285083189e-7, 1e-18)

def test_xys06a():
    x, y, s = Xys06a(2400000.5, 53736.0)

    assert x == approx( 0.5791308482835292617e-3, 1e-14)
    assert y == approx( 0.4020580099454020310e-4, 1e-15)
    assert s == approx(-0.1220032294164579896e-7, 1e-18)

def test_ee00():

    epsa =  0.4090789763356509900
    dpsi = -0.9630909107115582393e-5

    ee = Ee00(2400000.5, 53736.0, epsa, dpsi)

    assert ee == approx(-0.8834193235367965479e-5, 1e-18)

def test_ee00a():
    ee = Ee00a(2400000.5, 53736.0)

    assert ee == approx(-0.8834192459222588227e-5, 1e-18)


def test_ee00b():
    ee = Ee00b(2400000.5, 53736.0)

    assert ee == approx(-0.8835700060003032831e-5, 1e-18)

def test_ee06a():
    ee = Ee06a(2400000.5, 53736.0)

    assert ee == approx(-0.8834195072043790156e-5, 1e-15)

def test_eect00():
    eect = Eect00(2400000.5, 53736.0)

    assert eect == approx(0.2046085004885125264e-8, 1e-20)

def test_eqeq94():
    eqeq = Eqeq94(2400000.5, 41234.0)

    assert eqeq == approx(0.5357758254609256894e-4, 1e-17)

def test_era00():
    era00 = Era00(2400000.5, 54388.0)

    assert era00 == approx(0.4022837240028158102, 1e-12)

def test_gmst00():
    theta = Gmst00(2400000.5, 53736.0, 2400000.5, 53736.0)

    assert theta == approx(1.754174972210740592, 1e-12)

def test_gmst06():
    theta = Gmst06(2400000.5, 53736.0, 2400000.5, 53736.0)    

    assert theta == approx(1.754174971870091203, 1e-12)

def test_gmst82():
    theta = Gmst82(2400000.5, 53736.0)

    assert theta == approx(1.754174981860675096, 1e-12)

def test_gst00a():
    theta = Gst00a(2400000.5, 53736.0, 2400000.5, 53736.0)

    assert theta == approx(1.754166138018281369, 1e-12)

def test_gst00b():
    theta = Gst00b(2400000.5, 53736.0)

    assert theta == approx(1.754166136510680589, 1e-12)

def test_gst06():
    rnpb = _np.array([[0.9999989440476103608, -0.1332881761240011518e-2, -0.5790767434730085097e-3],
                      [0.1332858254308954453e-2, 0.9999991109044505944, -0.4097782710401555759e-4],
                      [0.5791308472168153320e-3, 0.4020595661593994396e-4, 0.9999998314954572365]])

    theta = Gst06(2400000.5, 53736.0, 2400000.5, 53736.0, rnpb)

    assert theta == approx(1.754166138018167568, 1e-12)

def test_gst06a():
    theta = Gst06a(2400000.5, 53736.0, 2400000.5, 53736.0)

    assert theta == approx(1.754166137675019159, 1e-12)

def test_gst94():
    theta = Gst94(2400000.5, 53736.0)

    assert theta == approx(1.754166136020645203, 1e-12)

def test_d2dtf():
    iy, im, id, ihmsf = D2dtf('UTC', 5, 2400000.5, 49533.99999)

    assert iy == 1994
    assert im == 6
    assert id == 30
    assert ihmsf[0] == 23
    assert ihmsf[1] == 59
    assert ihmsf[2] == 60
    assert ihmsf[3] == 13599

def test_dat():
    deltat = Dat(2003, 6, 1, 0.0)
    assert deltat == 32.0

    deltat = Dat(2008, 1, 7, 0.0)
    assert deltat == 33.0

    deltat = Dat(2017, 9, 1, 0.0)
    assert deltat == 37.0

def test_dtdb():
    dtdb = Dtdb(2448939.5, 0.123, 0.76543, 5.0123, 5525.242, 3190.0)

    assert dtdb == approx(-0.1280368005936998991e-2, 1e-15)

def test_dtf2d():
    u1, u2 = Dtf2d("UTC", 1994, 6, 30, 23, 59, 60.13599)

    assert u1+u2 == approx(2449534.49999, 1e-6)

def test_taitt():
    t1, t2 = Taitt(2453750.5, 0.892482639)

    assert t1 == approx(2453750.5, 1e-6)
    assert t2 == approx(0.892855139, 1e-12)

def test_taiut1():
    u1, u2 = Taiut1(2453750.5, 0.892482639, -32.6659)

    assert u1 == approx(2453750.5, 1e-6)
    assert u2 == approx(0.8921045614537037037, 1e-12)

def test_taiutc():
    u1, u2 = Taiutc(2453750.5, 0.892482639)

    assert u1 == approx(2453750.5, 1e-6)
    assert u2 == approx(0.8921006945555555556, 1e-12)

def test_tcbtdb():
    b1, b2 = Tcbtdb(2453750.5, 0.893019599)

    assert b1 == approx(2453750.5, 1e-6)
    assert b2 == approx(0.8928551362746343397, 1e-12)

def test_tcgtt():
    t1, t2 = Tcgtt(2453750.5, 0.892862531)

    assert t1 == approx(2453750.5, 1e-6)
    assert t2 == approx(0.8928551387488816828, 1e-12)

def test_tdbtcb():
    b1, b2 = Tdbtcb(2453750.5, 0.892855137)

    assert b1 == approx(2453750.5, 1e-6)
    assert b2 == approx(0.8930195997253656716, 1e-12)

def test_tdbtt():
    t1, t2 = Tdbtt(2453750.5, 0.892855137, -0.000201)

    assert t1 == approx(2453750.5, 1e-6)
    assert t2 == approx(0.8928551393263888889, 1e-12)

def test_tttai():
    a1, a2 = Tttai(2453750.5, 0.892482639)

    assert a1 == approx(2453750.5, 1e-6)
    assert a2 == approx(0.892110139, 1e-12)

def test_tttcg():
    g1, g2 = Tttcg(2453750.5, 0.892482639)

    assert g1 == approx(2453750.5, 1e-6)
    assert g2 == approx(0.8924900312508587113, 1e-12)

def test_tttdb():
    b1, b2 = Tttdb(2453750.5, 0.892855139,  -0.000201)

    assert b1 == approx(2453750.5, 1e-6)
    assert b2 == approx(0.8928551366736111111, 1e-12)

def test_ttut1():
    u1, u2 = Ttut1(2453750.5, 0.892855139, 64.8499)

    assert u1 == approx(2453750.5, 1e-6)
    assert u2 == approx(0.8921045614537037037, 1e-12)

def test_ut1tai():
    a1, a2 = Ut1tai(2453750.5, 0.892104561, -32.6659)

    assert a1 == approx(2453750.5, 1e-6)
    assert a2 == approx(0.8924826385462962963, 1e-12)

def test_ut1utc():
    u1, u2 = Ut1utc(2453750.5, 0.892104561, 0.3341)

    assert u1 == approx(2453750.5, 1e-6)
    assert u2 == approx(0.8921006941018518519, 1e-12)

def test_utctai():
    u1, u2 = Utctai(2453750.5, 0.892100694)

    assert u1 == approx(2453750.5, 1e-6)
    assert u2 == approx(0.8924826384444444444, 1e-12)

def test_utcut1():
    u1, u2 = Utcut1(2453750.5, 0.892100694, 0.3341)

    assert u1 == approx(2453750.5, 1e-6)
    assert u2 == approx(0.8921045608981481481, 1e-12)

def test_rx():
    phi = 0.3456789

    r = _np.array([[2.0, 3.0, 2.0],
                  [3.0, 2.0, 3.0],
                  [3.0, 4.0, 5.0]])

    r = Rx(phi, r)

    assert r[0, 0] == approx(2.0, 0.0)
    assert r[0, 1] == approx(3.0, 0.0)
    assert r[0, 2] == approx(2.0, 0.0)

    assert r[1, 0] == approx(3.839043388235612460, 1e-12)
    assert r[1, 1] == approx(3.237033249594111899, 1e-12)
    assert r[1, 2] == approx(4.516714379005982719, 1e-12)

    assert r[2, 0] == approx(1.806030415924501684, 1e-12)
    assert r[2, 1] == approx(3.085711545336372503, 1e-12)
    assert r[2, 2] == approx(3.687721683977873065, 1e-12)

def test_ry():
    theta = 0.3456789

    r = _np.array([[2.0, 3.0, 2.0],
                   [3.0, 2.0, 3.0],
                   [3.0, 4.0, 5.0]])

    r = Ry(theta, r)

    assert r[0, 0] == approx(0.8651847818978159930, 1e-12)
    assert r[0, 1] == approx(1.467194920539316554, 1e-12)
    assert r[0, 2] == approx(0.1875137911274457342, 1e-12)

    assert r[1, 0] == approx(3, 1e-12)
    assert r[1, 1] == approx(2, 1e-12)
    assert r[1, 2] == approx(3, 1e-12)

    assert r[2, 0] == approx(3.500207892850427330, 1e-12)
    assert r[2, 1] == approx(4.779889022262298150, 1e-12)
    assert r[2, 2] == approx(5.381899160903798712, 1e-12)

def test_rz():
    psi = 0.3456789

    r = _np.array([[2.0, 3.0, 2.0],
                   [3.0, 2.0, 3.0],
                   [3.0, 4.0, 5.0]])

    r = Rz(psi, r)

    assert r[0, 0] == approx(2.898197754208926769, 1e-12)
    assert r[0, 1] == approx(3.500207892850427330, 1e-12)
    assert r[0, 2] == approx(2.898197754208926769, 1e-12)

    assert r[1, 0] == approx(2.144865911309686813, 1e-12)
    assert r[1, 1] == approx(0.865184781897815993, 1e-12)
    assert r[1, 2] == approx(2.144865911309686813, 1e-12)

    assert r[2, 0] == approx(3.0, 1e-12)
    assert r[2, 1] == approx(4.0, 1e-12)
    assert r[2, 2] == approx(5.0, 1e-12)