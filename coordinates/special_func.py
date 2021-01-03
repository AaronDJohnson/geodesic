from mpmath import ellipfun
from numpy import isreal, asin, floor, tanh


def jacobi_am_f(x, m):
    """Borrowed from sagemath."""
    # print(x)
    # print(m)
    if not isreal(x) or not isreal(m):
        raise ValueError("arguments must be real")
    if abs(m) == 1:
        # gd(x)
        tanhx = tanh(x)
        return asin(tanhx)
    elif abs(m) > 1:
        # Real values needed for atan2; as per "Handbook of Elliptic
        # Integrals for Engineers and Scientists" 121.02, sn is real for
        # real x. The imaginary components can thus be safely discarded.
        snx = ellipfun("sn", x, m).real
        cnx = ellipfun("cn", x, m).real
        return ctx.atan2(snx, cnx)
    else:
        ctx.prec += 10
        K = ctx.ellipk(m)
        if abs(x) <= K:
            snx = ctx.ellipfun("sn", x, m).real
            cnx = ctx.ellipfun("cn", x, m).real
            ctx.prec += 10
            return ctx.atan2(snx, cnx)
        else:
            # Do argument reduction on x to end up with z = x - 2nK, with
            # abs(z) <= K
            ctx.prec += 10
            tK = 2 * K
            ctx.prec += 10
            n = ctx.floor(x / tK)
            ctx.prec += 10
            tnK = n * tK
            npi = n * ctx.pi()
            ctx.prec += 10
            z = x - tnK
            ctx.prec += 10
            # z (and therefore sn(z, m) and cn(z, m)) is real because K(m)
            # is real for abs(m) <= 1.
            snz = ctx.ellipfun("sn", z, m).real
            cnz = ctx.ellipfun("cn", z, m).real
            ctx.prec += 10
            return ctx.atan2(snz, cnz) + npi


def am(x, m):
    """
    Jacobi amplitude function for real input parameters

    (JacobiAmplitude[] in Mathematica)

    Parameters:
        x (float)
        m (float)

    Returns:
        am (float)
    """
    if m == 0:
        return x
    elif x == 0:
        return 0
    else:
        return jacobi_am_f(x, m)
