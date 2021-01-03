from mpmath import sqrt, mp
from sys import exit

# ------------------------------------------------------------------------------
#  SC orbits (a = 0)
# ------------------------------------------------------------------------------


def calc_sc_energy(slr, ecc, x):
    try:
        prec = mp.prec
        mp.prec += 10

        ecc2 = ecc * ecc
        slr2 = slr * slr
        x2 = x * x

        mp.prec += 10
        En = sqrt((-4 * ecc2 + (-2 + slr) ** 2) / (slr * (-3 - ecc2 + slr)))
        return En
    finally:
        mp.prec = prec


def calc_sc_carter(slr, ecc, x):
    try:
        prec = mp.prec
        mp.prec += 10

        ecc2 = ecc * ecc
        slr2 = slr * slr
        x2 = x * x

        mp.prec += 10
        Lz = (slr * x) / sqrt(-3 - ecc2 + slr)
        return Lz
    finally:
        mp.prec = prec


def calc_sc_ang_momentum(slr, ecc, x):
    try:
        prec = mp.prec
        mp.prec += 10

        ecc2 = ecc * ecc
        slr2 = slr * slr
        x2 = x * x

        mp.prec += 10
        Q = (slr2 * (-1 + x2)) / (3 + ecc2 - slr)
        return Q
    finally:
        mp.prec = prec


def calc_sc_constants(slr, ecc, x):
    """
    Energy, angular momentum, and carter constant calculation.

    Schwarzschild case (spin parameter a = 0)

    Parameters:
        slr (mpf): semi-latus rectum [6, inf)
        ecc (mpf): eccentricity [0, 1)
        x (mpf): inclination value given by cos(theta_inc) (0, 1]
                   negative x -> retrograde
                   positive x -> prograde

    Returns:
        En (mpf): energy
        Lz (mpf): angular momentum
        Q (mpf): Carter constant
    """
    try:
        prec = mp.prec
        mp.prec += 10

        ecc2 = ecc * ecc
        slr2 = slr * slr
        x2 = x * x

        mp.prec += 10
        En = sqrt((-4 * ecc2 + (-2 + slr) ** 2) / (slr * (-3 - ecc2 + slr)))
        Lz = (slr * x) / sqrt(-3 - ecc2 + slr)
        Q = (slr2 * (-1 + x2)) / (3 + ecc2 - slr)

        return En, Lz, Q
    finally:
        mp.prec = prec
