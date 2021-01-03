from mpmath import sqrt, mp
from sys import exit

# ------------------------------------------------------------------------------
#  Kerr spherical orbits (ecc = 0)
# ------------------------------------------------------------------------------


def spherical_energy(aa, slr, x):
    """
    Compute energy for the Kerr spherical case (ecc = 0).

    Parameters:
        aa (mpf): spin parameter (0, 1)
        slr (mpf): semi-latus rectum [6, inf)
        x (mpf): inclination value given by cos(theta_inc) (0, 1]
                   negative x -> retrograde
                   positive x -> prograde

    Returns:
        En (mpf): energy
    """
    try:
        prec = mp.prec
        mp.prec += 10

        slrm2 = (-2 + slr) ** 2
        slr2 = slr * slr
        slr3 = slr * slr2
        slr4 = slr * slr3
        slr5 = slr3 * slr2
        slr7 = slr5 * slr2
        x2 = x * x
        x3 = x2 * x
        aa2 = aa * aa
        aa3 = aa2 * aa
        aa4 = aa2 * aa2
        aa5 = aa4 * aa
        aa6 = aa4 * aa2
        x2m1 = x2 - 1
        x2m12 = x2m1 * x2m1

        mp.prec += 10
        return sqrt(
            (
                (-3 + slr) * slrm2 * slr5
                - 2 * aa5 * x * x2m1 * sqrt(slr3 + aa2 * slr * x2m1)
                + aa4 * slr2 * x2m1 * (4 - 5 * slr * (-1 + x2) + 3 * slr2 * x2m1)
                - aa6 * x2m12 * (x2 + slr2 * x2m1 - slr * (1 + 2 * x2))
                + aa2
                * slr3
                * (
                    4
                    - 4 * x2
                    + slr * (12 - 7 * x2)
                    - 3 * slr3 * (-1 + x2)
                    + slr2 * (-13 + 10 * x2)
                )
                + aa
                * (
                    -2 * slr ** 4.5 * x * sqrt(slr2 + aa2 * x2m1)
                    + 4 * slr3 * x * sqrt(slr3 + aa2 * slr * x2m1)
                )
                + 2
                * aa3
                * (
                    2 * slr * x * x2m1 * sqrt(slr3 + aa2 * slr * x2m1)
                    - x3 * sqrt(slr7 + aa2 * slr5 * x2m1)
                )
            )
            / (
                (slr2 - aa2 * x2m1)
                * (
                    (-3 + slr) ** 2 * slr4
                    - 2 * aa2 * slr2 * (3 + 2 * slr - 3 * x2 + slr2 * x2m1)
                    + aa4 * x2m1 * (-1 + x2 + slr2 * x2m1 - 2 * slr * (1 + x2))
                )
            )
        )
    finally:
        mp.prec = prec


def spherical_ang_momentum(En, aa, slr, x):
    """
    Compute angular momentum for the Kerr spherical case (ecc = 0).

    Parameters:
        En (mpf): energy
        aa (mpf): spin parameter (0, 1)
        slr (mpf): semi-latus rectum [6, inf)
        x (mpf): inclination value given by cos(theta_inc) (0, 1]
                   negative x -> retrograde
                   positive x -> prograde

    Returns:
        Lz (mpf): angular momentum
    """
    try:
        prec = mp.prec
        mp.prec += 10

        aa2 = aa * aa
        slr2 = slr * slr
        x2 = x * x
        mp.prec += 10
        g = 2 * aa * slr
        d = (aa2 + (-2 + slr) * slr) * (slr2 - aa2 * (-1 + x2))
        h = ((-2 + slr) * slr - aa2 * (-1 + x2)) / x2
        f = slr ** 4 + aa2 * (slr * (2 + slr) - (aa2 + (-2 + slr) * slr) * (-1 + x2))

        return (-(En * g) + sqrt((-(d * h) + En ** 2 * (g ** 2 + f * h)) / x2) * x) / h
    finally:
        mp.prec = prec


def gen_carter_const(En, Lz, aa, slr, ecc, x):
    """
    Compute Carter constant for generic orbit case.

    Parameters:
        En (mpf): energy
        Lz (mpf): angular momentum
        aa (mpf): spin parameter (0, 1)
        slr (mpf): semi-latus rectum [6, inf)
        ecc (mpf): eccentricity [0, 1)
        x (mpf): inclination value given by cos(theta_inc)
                   x < 0 -> retrograde
                   x > 0 -> prograde

    Returns:
        Q (mpf): Carter constant
    """
    try:
        prec = mp.prec
        mp.prec += 10
        zm = sqrt(1 - x * x)
        mp.prec += 10
        zm2 = zm * zm
        mp.prec += 10
        return zm2 * (aa * aa * (1 - En * En) + Lz * Lz / (1 - zm2))
    finally:
        mp.prec = prec


def calc_sph_constants(aa, slr, x):
    """
    Call spherical orbit (ecc = 0) functions.

    Parameters:
        aa (mpf): spin parameter (0, 1)
        slr (mpf): semi-latus rectum [6, inf)
        x (mpf): inclination value given by cos(theta_inc)
                   x < 0 -> retrograde
                   x > 0 -> prograde

    Returns:
        En (mpf): energy
        Lz (mpf): angular momentum
        Q (mpf): Carter constant
    """
    ecc = 0

    En = spherical_energy(aa, slr, x)
    Lz = spherical_ang_momentum(En, aa, slr, x)
    Q = gen_carter_const(En, Lz, aa, slr, ecc, x)

    return En, Lz, Q
