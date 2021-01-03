from mpmath import sqrt, mp
from sys import exit

# ------------------------------------------------------------------------------
#  Kerr equatorial orbits (x = +/- 1) -> sign is pro/retrograde motion
# ------------------------------------------------------------------------------


def eq_energy(aa, slr, ecc, x):
    """
    Energy for the Kerr equatorial case (inclination value x = 1).

    Parameters:
        aa (mpf): spin parameter (0, 1)
        slr (mpf): semi-latus rectum [~2, inf)
        ecc (mpf): eccentricity [0, 1)
        x (mpf): inclination value given by cos(theta_inc)
                   x = -1 -> retrograde
                   x = +1 -> prograde

    Returns:
        En (mpf): energy
    """
    try:
        prec = mp.prec
        mp.prec += 10

        ecc2 = ecc * ecc
        eta = ecc2 - 1
        eta2 = eta * eta
        aa2 = aa * aa
        aa4 = aa2 * aa2
        aa6 = aa4 * aa2
        mp.prec += 10
        return sqrt(
            1
            - (
                -eta
                * (
                    1
                    + (
                        eta
                        * (
                            aa2 * (1 + 3 * ecc2 + slr)
                            + slr
                            * (
                                -3
                                - ecc2
                                + slr
                                - 2
                                * sqrt(
                                    (
                                        aa6 * eta2
                                        + aa2 * (-4 * ecc2 + (-2 + slr) ** 2) * slr ** 2
                                        + 2 * aa4 * slr * (-2 + slr + ecc2 * (2 + slr))
                                    )
                                    / (slr ** 3 * x ** 2)
                                )
                                * x
                            )
                        )
                    )
                    / (-4 * aa2 * eta2 + (3 + ecc2 - slr) ** 2 * slr)
                )
            )
            / slr
        )
    finally:
        mp.prec = prec


def eq_ang_momentum(aa, slr, ecc, x):
    """
    Angular momentum for the Kerr equatorial case (inclination value x = 1).

    Parameters:
        aa (mpf): spin parameter (0, 1)
        slr (mpf): semi-latus rectum [6, inf)
        ecc (mpf): eccentricity [0, 1)
        x (mpf): inclination value given by cos(theta_inc)
                   x = -1 -> retrograde
                   x = +1 -> prograde

    Returns:
        Lz (mpf): angular momentum
    """
    try:
        prec = mp.prec
        mp.prec += 10

        ecc2 = ecc * ecc
        eta = ecc2 - 1
        eta2 = eta * eta
        aa2 = aa * aa
        aa4 = aa2 * aa2
        aa6 = aa4 * aa2
        slr2 = slr * slr
        slr3 = slr2 * slr
        x2 = x * x
        mp.prec += 10
        num_root = slr * (
            -3
            - ecc2
            + slr
            - 2
            * sqrt(
                (
                    aa6 * eta2
                    + aa2 * (-4 * ecc2 + (-2 + slr) ** 2) * slr2
                    + 2 * aa4 * slr * (-2 + slr + ecc2 * (2 + slr))
                )
                / (slr3 * x2)
            )
            * x
        )
        denom = -4 * aa2 * eta2 + (3 + ecc2 - slr) ** 2 * slr
        return slr * x * sqrt(
            (aa2 * (1 + 3 * ecc2 + slr) + num_root) / (denom * x2)
        ) + aa * sqrt(
            1
            - (-eta * (1 + (eta * (aa2 * (1 + 3 * ecc2 + slr) + num_root)) / denom))
            / slr
        )
    finally:
        mp.prec = prec


def calc_eq_constants(aa, slr, ecc, x):
    """
    Calculate equatorial constants in one function.

    Note that for the equatorial case, the Carter constant, Q = 0.

    Parameters:
        aa (mpf): spin parameter (0, 1)
        slr (mpf): semi-latus rectum [6, inf)
        ecc (mpf): eccentricity [0, 1)
        x (mpf): inclination value given by cos(theta_inc)
                   x = -1 -> retrograde
                   x = +1 -> prograde

    Returns:
        En (mpf): energy
        Lz (mpf): angular momentum
        Q (mpf): Carter constant
    """
    En = eq_energy(aa, slr, ecc, x)
    Lz = eq_ang_momentum(aa, slr, ecc, x)
    Q = 0  # equatorial case
    return En, Lz, Q
