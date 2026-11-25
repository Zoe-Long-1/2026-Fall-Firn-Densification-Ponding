"""Firn densification with the Herron and Langway (1980) model.

For each layer (a Lagrangian parcel of firn moving down with burial):

    d rho / dt = c0 * b       * (rho_ice - rho)   for rho <  550 kg m^-3
    d rho / dt = c1 * b**0.5  * (rho_ice - rho)   for rho >= 550 kg m^-3

    c0 = 11  * exp(-10160 / (R T))      c1 = 575 * exp(-21400 / (R T))

with b the accumulation rate [m w.e. per year], T in kelvin, and the rate in
units of per year.  Over one time step the rate coefficient is held fixed
(operator splitting), so the ODE can be integrated exactly:

    rho_new = rho_ice - (rho_ice - rho) * exp(-C * dt)

which never overshoots ice density, whatever the step size.
"""

import numpy as np

from .params import R_GAS, RHO_ICE, RHO_WATER

RHO_CRITICAL = 550.0   # transition between the two densification stages [kg m^-3]


def rate_coefficient(T_celsius, rho, accumulation):
    """Herron-Langway rate coefficient C [1/yr], so that d rho/dt = C (rho_ice - rho)."""
    T_kelvin = np.asarray(T_celsius, dtype=float) + 273.15
    c0 = 11.0 * np.exp(-10160.0 / (R_GAS * T_kelvin))
    c1 = 575.0 * np.exp(-21400.0 / (R_GAS * T_kelvin))
    b = max(accumulation, 0.0)
    return np.where(np.asarray(rho) < RHO_CRITICAL, c0 * b, c1 * np.sqrt(b))


def densify(rho, T_celsius, accumulation, dt_years):
    """Return new layer densities after a step of ``dt_years``."""
    C = rate_coefficient(T_celsius, rho, accumulation)
    rho_new = RHO_ICE - (RHO_ICE - rho) * np.exp(-C * dt_years)
    return np.minimum(rho_new, RHO_ICE)


def steady_state_density(z, T_celsius, accumulation, rho_surface):
    """Analytical steady-state Herron-Langway density profile rho(z) [kg m^-3].

    With constant temperature and accumulation, firn moves down at
    v = b * rho_w / rho, and the steady profile satisfies
    v * d rho/dz = C (rho_ice - rho).  For each stage this gives

        d/dz ln( rho / (rho_ice - rho) ) = rho_ice * c0 / rho_w            (stage 1)
        d/dz ln( rho / (rho_ice - rho) ) = rho_ice * c1 / (rho_w sqrt(b))  (stage 2)

    i.e. a logistic curve in each stage, joined at 550 kg m^-3.
    """
    z = np.asarray(z, dtype=float)
    T_kelvin = T_celsius + 273.15
    c0 = 11.0 * np.exp(-10160.0 / (R_GAS * T_kelvin))
    c1 = 575.0 * np.exp(-21400.0 / (R_GAS * T_kelvin))
    slope1 = RHO_ICE * c0 / RHO_WATER
    slope2 = RHO_ICE * c1 / (RHO_WATER * np.sqrt(accumulation))

    def logit(r):
        return np.log(r / (RHO_ICE - r))

    def inverse_logit(y):
        return RHO_ICE * np.exp(y) / (1.0 + np.exp(y))

    z_critical = (logit(RHO_CRITICAL) - logit(rho_surface)) / slope1
    stage1 = inverse_logit(logit(rho_surface) + slope1 * z)
    stage2 = inverse_logit(logit(RHO_CRITICAL) + slope2 * (z - z_critical))
    return np.where(z < z_critical, stage1, stage2)
