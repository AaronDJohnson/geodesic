# geodesic
 Computes geodesics of a stellar mass black hole around a massive black hole

* This work uses several of the routines from the Black Hole Perturbation Toolkit.
* Everything included is computed to machine precision using numpy or scipy where available.
* This still uses mpmath for elliptic pi because scipy doesn't implement this yet.

## Available functions:
 * coordinates of geodesics
 * adiabatic constants (energy, angular momentum, Carter constant)
 * Boyer Lindquist frequencies
 * Mino frequencies
