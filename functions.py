from mpmath import mp, mpf, cos
try:
    from constants.constants import calc_constants
    from geo_roots import radial_roots, polar_roots
    from frequencies import mino_freqs, find_omega, mino_freqs, boyer_freqs
    from coordinates.coords import calc_coords
except:
    from geodesic.constants.constants import calc_constants
    from geodesic.geo_roots import radial_roots, polar_roots
    from geodesic.frequencies import mino_freqs, find_omega, mino_freqs, boyer_freqs
    from geodesic.coordinates.coords import calc_coords
    from geodesic.coordinates.coords_gen import calc_gen_coords_mino


def calc_consts(aa, slr, ecc, x, digits):
    """
    Compute adiabatic constants.

    Parameters:
        aa (float): SMBH spin
        slr (float): semi-latus rectum
        ecc (float): eccentricity
        x (float): cos of the inclination
        digits (int): number of digits of accuracy requested

    Returns:
        En (mpf): energy
        Lz (mpf): angular momentum
        Q (mpf): Carter constant
    """
    try:
        mp.dps = digits
        prec = mp.prec
        mp.prec += 30
        aa = mpf(str(aa))
        slr = mpf(str(slr))
        ecc = mpf(str(ecc))
        x = mpf(str(x))
        En, Lz, Q = calc_constants(aa, slr, ecc, x)
        return En, Lz, Q
    finally:
        mp.prec = prec


def calc_radial_roots(aa, slr, ecc, x, digits):
    """
    Compute radial roots.

    Parameters:
        aa (float): SMBH spin
        slr (float): semi-latus rectum
        ecc (float): eccentricity
        x (float): cos of the inclination
        digits (int): number of digits of accuracy requested

    Returns:
        r1 (mpf): apastron
        r2 (mpf): periastron
        r3 (mpf): radial root 3
        r4 (mpf): radial root 4
    """
    try:
        mp.dps = digits
        prec = mp.prec
        mp.prec += 30
        aa = mpf(str(aa))
        slr = mpf(str(slr))
        ecc = mpf(str(ecc))
        x = mpf(str(x))
        En, Lz, Q = calc_constants(aa, slr, ecc, x)
        r1, r2, r3, r4 = radial_roots(En, Q, aa, slr, ecc)
        return r1, r2, r3, r4
    finally:
        mp.prec = prec


def calc_polar_roots(aa, slr, ecc, x, digits):
    """
    Compute polar roots.

    Parameters:
        aa (float): SMBH spin
        slr (float): semi-latus rectum
        ecc (float): eccentricity
        x (float): cos of the inclination
        digits (int): number of digits of accuracy requested

    Returns:
        zp (mpf): polar root
        zm (mpf): polar root
    """
    try:
        mp.dps = digits
        prec = mp.prec
        mp.prec += 30
        aa = mpf(str(aa))
        slr = mpf(str(slr))
        ecc = mpf(str(ecc))
        x = mpf(str(x))
        En, Lz, Q = calc_constants(aa, slr, ecc, x)
        zp, zm = polar_roots(En, Lz, aa, slr, x)
        return zp, zm
    finally:
        mp.prec = prec


def calc_mino_freqs(aa, slr, ecc, x, digits, M=1):
    """
    Compute Mino frequencies.

    Parameters:
        aa (float): SMBH spin
        slr (float): semi-latus rectum
        ecc (float): eccentricity
        x (float): cos of the inclination
        digits (int): number of digits of accuracy requested

    Returns:
        ups_r (mpf): radial Mino frequency
        ups_theta (mpf): polar Mino frequency
        ups_phi (mpf): azimuthal Mino frequency
        gamma (mpf): temporal Mino frequency
    """
    try:
        mp.dps = digits
        prec = mp.prec
        mp.prec += 30
        aa = mpf(str(aa))
        slr = mpf(str(slr))
        ecc = mpf(str(ecc))
        x = mpf(str(x))
        En, Lz, Q = calc_constants(aa, slr, ecc, x)
        r1, r2, r3, r4 = radial_roots(En, Q, aa, slr, ecc, M)
        ups_r, ups_theta, ups_phi, gamma = mino_freqs(
            r1, r2, r3, r4, En, Lz, Q, aa, slr, ecc, x
        )
        return ups_r, ups_theta, ups_phi, gamma
    finally:
        mp.prec = prec


def calc_boyer_freqs(aa, slr, ecc, x, digits):
    """
    Compute Boyer-Lindquist frequencies.

    Parameters:
        aa (float): SMBH spin
        slr (float): semi-latus rectum
        ecc (float): eccentricity
        x (float): cos of the inclination
        digits (int): number of digits of accuracy requested

    Returns:
        Omega_r (mpf): radial Boyer-Lindquist frequency
        Omega_theta (mpf): polar Boyer-Lindquist frequency
        Omega_phi (mpf): azimuthal Boyer-Lindquist frequency
    """
    mp.dps = digits
    aa = mpf(str(aa))
    slr = mpf(str(slr))
    ecc = mpf(str(ecc))
    x = mpf(str(x))
    En, Lz, Q = mp_calc_constants(aa, slr, ecc, x)
    r1, r2, r3, r4 = mp_radial_roots(En, Q, aa, slr, ecc, M)
    ups_r, ups_theta, ups_phi, gamma = mp_mino_freqs(
        r1, r2, r3, r4, En, Lz, Q, aa, slr, ecc, x
    )
    omega_r, omega_theta, omega_phi = mp_boyer_freqs(
        ups_r, ups_theta, ups_phi, gamma, aa, slr, ecc, x, M
    )
    return omega_r, omega_theta, omega_phi


def find_omega(en, em, kay, aa, slr, ecc, x, digits):
    """
    Compute gravitational wave frequency omega.

    Parameters:
        en (int): radial mode
        em (int): azimuthal mode
        kay (int): polar mode
        aa (float): SMBH spin
        slr (float): semi-latus rectum
        ecc (float): eccentricity
        x (float): cos of the inclination
        digits (int): number of digits of accuracy requested

    Returns:
        omega (mpf): gravitational wave frequency
    """
    mp.dps = digits
    aa = mpf(str(aa))
    slr = mpf(str(slr))
    ecc = mpf(str(ecc))
    x = mpf(str(x))
    En, Lz, Q = mp_calc_constants(aa, slr, ecc, x)
    r1, r2, r3, r4 = mp_radial_roots(En, Q, aa, slr, ecc, M)
    ups_r, ups_theta, ups_phi, gamma = mp_mino_freqs(
        r1, r2, r3, r4, En, Lz, Q, aa, slr, ecc, x
    )
    omega_r, omega_theta, omega_phi = mp_boyer_freqs(
        ups_r, ups_theta, ups_phi, gamma, aa, slr, ecc, x, M
    )
    omega = en * omega_r + em * omega_phi + kay * omega_theta
    return omega


def coordinates(psi, aa, slr, ecc, x, digits):
    """
    Compute coordinates of the orbit given two angles psi and chi.

    Parameters:
        psi (float): radial angle
        chi (float): polar angle
        aa (float): SMBH spin
        slr (float): semi-latus rectum
        ecc (float): eccentricity
        x (float): cos of the inclination
        digits (int): number of digits of accuracy requested

    Returns:
        t (mpf): time coordinate
        r (mpf): radial coordinate
        theta (mpf): theta coordinate
        phi (mpf): phi coordinate
    """
    if x ** 2 == 1:
        chi = 0  # equatorial orbits don't need chi
    try:
        mp.dps = digits
        prec = mp.prec
        mp.prec += 30
        psi = mpf(str(psi))
        # chi = mpf(str(chi))
        aa = mpf(str(aa))
        slr = mpf(str(slr))
        ecc = mpf(str(ecc))
        x = mpf(str(x))
        En, Lz, Q = calc_constants(aa, slr, ecc, x)
        r1, r2, r3, r4 = radial_roots(En, Q, aa, slr, ecc, M=1)
        zp, zm = polar_roots(En, Lz, aa, slr, x)
        ups_r, ups_theta, ups_phi, gamma = mino_freqs(
            r1, r2, r3, r4, En, Lz, Q, aa, slr, ecc, x
        )
        print(psi)
        # print(chi)
        t, r, theta, phi = calc_coords(
            psi,
            # chi,
            ups_r,
            ups_theta,
            ups_phi,
            gamma,
            r1,
            r2,
            r3,
            r4,
            zp,
            zm,
            En,
            Lz,
            Q,
            aa,
            slr,
            ecc,
            x,
        )
        return t, r, theta, phi
    finally:
        mp.prec = prec


def mino_coords(mino_t, aa, slr, ecc, x, digits):
    """
    Compute coordinates of the orbit given two angles psi and chi.

    Parameters:
        mino_t (float): mino time coordinate lambda
        aa (float): SMBH spin
        slr (float): semi-latus rectum
        ecc (float): eccentricity
        x (float): cos of the inclination
        digits (int): number of digits of accuracy requested

    Returns:
        t (mpf): time coordinate
        r (mpf): radial coordinate
        theta (mpf): theta coordinate
        phi (mpf): phi coordinate
    """
    try:
        mp.dps = digits
        prec = mp.prec
        mp.prec += 30
        mino_t = mpf(str(mino_t))
        aa = mpf(str(aa))
        slr = mpf(str(slr))
        ecc = mpf(str(ecc))
        x = mpf(str(x))
        En, Lz, Q = calc_constants(aa, slr, ecc, x)
        r1, r2, r3, r4 = radial_roots(En, Q, aa, slr, ecc, M=1)
        zp, zm = polar_roots(En, Lz, aa, slr, x)
        ups_r, ups_theta, ups_phi, gamma = mino_freqs(
            r1, r2, r3, r4, En, Lz, Q, aa, slr, ecc, x
        )
        t, r, theta, phi = calc_gen_coords_mino(
            mino_t,
            ups_r,
            ups_theta,
            ups_phi,
            gamma,
            r1,
            r2,
            r3,
            r4,
            zp,
            zm,
            En,
            Lz,
            Q,
            aa,
        )
        return t, r, theta, phi
    finally:
        mp.prec = prec