"""
Test all geodesic functions here.

This file compares results of geodesic functions with the
Black Hole Perturbation Toolkit written in Mathematica.
"""
import pytest
import numpy as np
from mpmath import mpf, mp, almosteq
from functions import radial_roots, polar_roots

digits = 100  # accuracy requested
eps = 10 ** (-digits)  # number of digits which may be different
mp.dps = digits  # set precision


def test_polar_sc_polar_circular():
    aa = 0
    slr = 6
    ecc = 0
    x = 0
    zp_ch, zm_ch = polar_roots(aa, slr, ecc, x, digits)
    zm = mpf("1")
    zp = mpf("0")
    assert almosteq(zm_ch, zm, eps)
    assert almosteq(zp_ch, zp, eps)


def test_polar_sc_inclined_circular():
    aa = 0
    slr = 6
    ecc = 0
    x = 0.5
    zp_ch, zm_ch = polar_roots(aa, slr, ecc, x, digits)
    zm = mpf(
        "0.86602540378443864676372317075293618347140262690519031402790348972596\
65084544000185405730933786242878378130707077033515149849725475"
    )
    zp = mpf(
        "3.46410161513775458705489268301174473388561050762076125611161395890386\
603381760007416229237351449715135125228283081340605993989019"
    )
    assert almosteq(zm_ch, zm, eps)
    assert almosteq(zp_ch, zp, eps)


def test_polar_sc_circular_equatorial():
    aa = 0
    slr = 6
    ecc = 0
    x = 1
    zp_ch, zm_ch = polar_roots(aa, slr, ecc, x, digits)
    zm = mpf("0")
    zp = mpf(
        "3.46410161513775458705489268301174473388561050762076125611161395890386\
603381760007416229237351449715135125228283081340605993989019"
    )
    assert almosteq(zm_ch, zm, eps)
    assert almosteq(zp_ch, zp, eps)


def test_polar_kerr_polar_circular():
    aa = 0.9
    slr = 6
    ecc = 0
    x = 0
    zp_ch, zm_ch = polar_roots(aa, slr, ecc, x, digits)
    zm = mpf("1")
    zp = mpf("0")
    assert almosteq(zm_ch, zm, eps)
    assert almosteq(zp_ch, zp, eps)


def test_polar_kerr_inclined_circular():
    aa = 0.9
    slr = 6
    ecc = 0
    x = 0.5
    zp_ch, zm_ch = polar_roots(aa, slr, ecc, x, digits)
    zm = mpf(
        "0.86602540378443864676372317075293618347140262690519031402790348972596\
65084544000185405730933786242878378130707077033515149849725475"
    )
    zp = mpf(
        "3.05543106615531533387676777007068231105836775333459432935070506527093\
51113541222589809215722610115925073534751033696119963632191"
    )
    assert almosteq(zm_ch, zm, eps)
    assert almosteq(zp_ch, zp, eps)


def test_polar_kerr_circular_equatorial():
    aa = 0.9
    slr = 6
    ecc = 0
    x = 1
    zp_ch, zm_ch = polar_roots(aa, slr, ecc, x, digits)
    zm = mpf("0")
    zp = mpf(
        "2.81576412442877462061915538563996125884423339096478120511885059807601\
627376072301848960531274422940925849602470358318845423602414"
    )
    assert almosteq(zm_ch, zm, eps)
    assert almosteq(zp_ch, zp, eps)


def test_polar_sc_polar_low_ecc():
    aa = 0
    slr = 6
    ecc = 0.1
    x = 0
    zp_ch, zm_ch = polar_roots(aa, slr, ecc, x, digits)
    zm = mpf("1")
    zp = mpf("0")
    assert almosteq(zm_ch, zm, eps)
    assert almosteq(zp_ch, zp, eps)


def test_polar_sc_inclined_low_ecc():
    aa = 0
    slr = 6
    ecc = 0.1
    x = 0.5
    zp_ch, zm_ch = polar_roots(aa, slr, ecc, x, digits)
    zm = mpf(
        "0.86602540378443864676372317075293618347140262690519031402790348972596\
65084544000185405730933786242878378130707077033515149849725475"
    )
    zp = mpf(
        "3.46988959179744133351400773640924026761245476148250978303554387894263\
507063490349890758982948765400825359194934322782659963821907"
    )
    assert almosteq(zm_ch, zm, eps)
    assert almosteq(zp_ch, zp, eps)


def test_polar_sc_equatorial_low_ecc():
    aa = 0
    slr = 6
    ecc = 0.1
    x = 1
    zp_ch, zm_ch = polar_roots(aa, slr, ecc, x, digits)
    zm = mpf("0")
    zp = mpf(
        "3.46988959179744133351400773640924026761245476148250978303554387894263\
507063490349890758982948765400825359194934322782659963821907"
    )
    assert almosteq(zm_ch, zm, eps)
    assert almosteq(zp_ch, zp, eps)


def test_polar_sc_equatorial_medium_ecc():
    aa = 0
    slr = 6
    ecc = 0.5
    x = 1
    zp_ch, zm_ch = polar_roots(aa, slr, ecc, x, digits)
    zm = mpf("0")
    zp = mpf(
        "3.61813613493316347176174480364074910973864204973384028770038052303616\
506466441146913692007177631707371674097219948449810066411613"
    )
    assert almosteq(zm_ch, zm, eps)
    assert almosteq(zp_ch, zp, eps)


def test_polar_sc_polar_high_ecc():
    aa = 0
    slr = 6
    ecc = 0.9999
    x = 0
    zp_ch, zm_ch = polar_roots(aa, slr, ecc, x, digits)
    zm = mpf("1")
    zp = mpf("0")
    assert almosteq(zm_ch, zm, eps)
    assert almosteq(zp_ch, zp, eps)


def test_polar_sc_high_ecc():
    aa = 0
    slr = 6
    ecc = 0.9999
    x = 0.5
    zp_ch, zm_ch = polar_roots(aa, slr, ecc, x, digits)
    zm = mpf(
        "0.86602540378443864676372317075293618347140262690519031402790348972596\
65084544000185405730933786242878378130707077033515149849725475"
    )
    zp = mpf(
        "4.24242858159851701578554748323491974415361125635078536064052378598976\
876818690433867661492273770095319175758399527157325757304661"
    )
    assert almosteq(zm_ch, zm, eps)
    assert almosteq(zp_ch, zp, eps)


def test_polar_sc_equatorial_high_ecc():
    aa = 0
    slr = 6
    ecc = 0.9999
    x = 1
    zp_ch, zm_ch = polar_roots(aa, slr, ecc, x, digits)
    zm = mpf("0")
    zp = mpf(
        "4.24242858159851701578554748323491974415361125635078536064052378598976\
876818690433867661492273770095319175758399527157325757304661"
    )
    assert almosteq(zm_ch, zm, eps)
    assert almosteq(zp_ch, zp, eps)


def test_polar_kerr_polar_high_ecc():
    aa = 0.9
    slr = 6
    ecc = 0.9999
    x = 0
    zp_ch, zm_ch = polar_roots(aa, slr, ecc, x, digits)
    zm = mpf("1")
    zp = mpf("0")
    assert almosteq(zm_ch, zm, eps)
    assert almosteq(zp_ch, zp, eps)


def test_polar_kerr_high_ecc():
    aa = 0.9
    slr = 6
    ecc = 0.9999
    x = 0.5
    zp_ch, zm_ch = polar_roots(aa, slr, ecc, x, digits)
    zm = mpf(
        "0.86602540378443864676372317075293618347140262690519031402790348972596\
65084544000185405730933786242878378130707077033515149849725475"
    )
    zp = mpf(
        "3.35955032142178191259340651151541677850998771597977762887568936095924\
34113586733182977383422102891310005387963394950"
    )
    assert almosteq(zm_ch, zm, eps)
    assert almosteq(zp_ch, zp, eps)


def test_radial_kerr_even_higher_ecc():
    aa = 0.9
    slr = 6
    ecc = 0.999999
    x = 0.5
    zp_ch, zm_ch = polar_roots(aa, slr, ecc, x, digits)
    zm = mpf(
        "0.86602540378443864676372317075293618347140262690519031402790348972596\
65084544000185405730933786242878378130707077033515149849725475"
    )
    zp = mpf(
        "3.35962035833558760981616667672183996231036459406120123792348715934372\
8405772076242436329420723358625696913331403315926505937159642756336278\
65960040873"
    )
    assert almosteq(zm_ch, zm, eps)
    assert almosteq(zp_ch, zp, eps)
