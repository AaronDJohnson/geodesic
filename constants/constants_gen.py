from mpmath import sqrt, mp
from sys import exit

# ------------------------------------------------------------------------------
#  Generic orbit constant calculation
# ------------------------------------------------------------------------------


def calc_delta(r, aa):
    """
    Calculate ubiquitous function on Kerr spacetimes.

    Parameters:
        r (mpf): radius
        aa (mpf): spin parameter (0, 1)

    Returns:
        delta (mpf)
    """
    try:
        prec = mp.prec
        mp.prec += 20
        # return r * r - 2 * r + aa * aa
        return r * (r - 2) + aa * aa
    finally:
        mp.prec = prec


def calc_f(r, zm, aa):
    """
    DESCRIPTION TBD (FUNCTION USED IN GENERAL CONSTANTS).

    Parameters:
        r (mpf): radius
        zm (mpf): 1 - x * x
        aa (mpf): spin parameter (0, 1)

    Returns:
        f (mpf)
    """
    try:
        prec = mp.prec
        mp.prec += 10
        r2 = r * r
        mp.prec += 10
        r4 = r2 * r2
        aa2 = aa * aa
        zm2 = zm * zm

        delta = calc_delta(r, aa)
        mp.prec += 10
        return 2 * aa2 * r + aa2 * r2 + r4 + aa2 * zm2 * delta
    finally:
        mp.prec = prec


def calc_g(r, aa):
    """
    DESCRIPTION TBD (FUNCTION USED IN GENERAL CONSTANTS).

    Parameters:
        r (mpf): radius
        aa (mpf): spin parameter (0, 1)

    Returns:
        g (mpf)
    """
    try:
        prec = mp.prec
        mp.prec += 10
        return 2 * aa * r
    finally:
        mp.prec = prec


def calc_h(r, zm, aa):
    """
    DESCRIPTION TBD (FUNCTION USED IN GENERAL CONSTANTS).

    Parameters:
        r (mpf): radius
        zm (mpf): 1 - x * x
        aa (mpf): spin parameter (0, 1)

    Returns:
        h (mpf)
    """
    try:
        prec = mp.prec
        mp.prec += 10
        zm2 = zm * zm
        delta = calc_delta(r, aa)
        mp.prec += 10
        return (-2 + r) * r + (zm2 * delta) / (1 - zm2)
    finally:
        mp.prec = prec


def calc_d(r, zm, aa):
    """
    DESCRIPTION TBD (FUNCTION USED IN GENERAL CONSTANTS).

    Parameters:
        r (mpf): radius
        zm (mpf): 1 - x * x
        aa (mpf): spin parameter (0, 1)

    Returns:
        d (mpf)
    """
    try:
        prec = mp.prec
        mp.prec += 10
        r2 = r * r
        aa2 = aa * aa
        zm2 = zm * zm

        delta = calc_delta(r, aa)
        mp.prec += 10
        return (r2 + aa2 * zm2) * delta
    finally:
        mp.prec = prec


def gen_energy(zm, aa, slr, ecc, x):
    """
    Compute energy for generic orbit case.

    Parameters:
        zm (mpf): 1 - x * x
        aa (mpf): spin parameter (0, 1)
        slr (mpf): semi-latus rectum [~2, inf)
        ecc (mpf): eccentricity [0, 1)
        x (mpf): inclination value given by cos(theta_inc)
                   x < 0 -> retrograde
                   x > 0 -> prograde

    Returns:
        En (mpf): energy
    """
    try:
        prec = mp.prec
        mp.prec += 20
        r1 = slr / (1 - ecc)
        r2 = slr / (1 + ecc)

        dr1 = calc_d(r1, zm, aa)
        dr2 = calc_d(r2, zm, aa)
        gr1 = calc_g(r1, aa)
        gr2 = calc_g(r2, aa)
        hr1 = calc_h(r1, zm, aa)
        hr2 = calc_h(r2, zm, aa)
        fr1 = calc_f(r1, zm, aa)
        fr2 = calc_f(r2, zm, aa)

        mp.prec += 10
        kappa = dr1 * hr2 - hr1 * dr2
        epsilon = dr1 * gr2 - gr1 * dr2
        rho = fr1 * hr2 - hr1 * fr2
        eta = fr1 * gr2 - gr1 * fr2
        sigma = gr1 * hr2 - hr1 * gr2

        mp.prec += 10
        kappa2 = kappa * kappa
        epsilon2 = epsilon * epsilon
        rho2 = rho * rho
        x2 = x * x

        mp.prec += 10
        En = sqrt(
            (
                kappa * rho
                + 2 * epsilon * sigma
                - 2
                * sqrt(
                    (
                        sigma
                        * (-(eta * kappa2) + epsilon * kappa * rho + epsilon2 * sigma)
                    )
                    / x2
                )
                * x
            )
            / (rho2 + 4 * eta * sigma)
        )

        return En
    finally:
        mp.prec = prec


def gen_ang_momentum(En, aa, slr, ecc, x):
    """
    Compute angular momentum for generic orbit case.

    Parameters:
        En (mpf): energy
        aa (mpf): spin parameter (0, 1)
        slr (mpf): semi-latus rectum [6, inf)
        ecc (mpf): eccentricity [0, 1)
        x (mpf): inclination value given by cos(theta_inc)
                   x < 0 -> retrograde
                   x > 0 -> prograde

    Returns:
        z_ang_mom (mpf): angular momentum
    """
    try:
        prec = mp.prec
        mp.prec += 10
        r1 = slr / (1 - ecc)
        zm = sqrt(1 - x * x)

        # // delta_r1 = calc_delta(r1, aa)
        mp.prec += 10
        fr1 = calc_f(r1, zm, aa)
        gr1 = calc_g(r1, aa)
        hr1 = calc_h(r1, zm, aa)
        dr1 = calc_d(r1, zm, aa)

        mp.prec += 10
        En2 = En * En
        x2 = x * x
        mp.prec += 10
        Lz = (
            -(En * gr1) + x * sqrt((-(dr1 * hr1) + En2 * (gr1 * gr1 + fr1 * hr1)) / x2)
        ) / hr1
        return Lz
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


def calc_gen_constants(aa, slr, ecc, x):
    """
    Call generic orbit constant calculating functions in one function.

    Parameters:
        aa (mpf): spin parameter (0, 1)
        slr (mpf): semi-latus rectum [6, inf)
        ecc (mpf): eccentricity [0, 1)
        x (mpf): inclination value given by cos(theta_inc)
                   x < 0 -> retrograde
                   x > 0 -> prograde

    Returns:
        En (mpf): energy
        Lz (mpf): angular momentum
        Q (mpf): Carter constant
    """
    try:
        prec = mp.prec
        mp.prec += 30
        zm = sqrt(1 - x * x)

        En = gen_energy(zm, aa, slr, ecc, x)
        Lz = gen_ang_momentum(En, aa, slr, ecc, x)
        Q = gen_carter_const(En, Lz, aa, slr, ecc, x)
        return En, Lz, Q
    finally:
        mp.prec = prec
