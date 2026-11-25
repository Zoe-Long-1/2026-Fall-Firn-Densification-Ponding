"""Surface climate forcing: air temperature and surface melt."""

import math

from .params import DAYS_PER_YEAR


def air_temperature(t_days, p):
    """Air temperature [deg C] at time ``t_days`` since the start of the run.

    T_air(t) = T_mean + (warming/10) * t_years - T_amp * cos(omega * (t - coldest_day))

    with omega = 2*pi / (365 days).  The minimum falls on ``coldest_day``
    (mid-January by default) and the maximum half a year later.
    """
    omega = 2.0 * math.pi / DAYS_PER_YEAR
    t_years = t_days / DAYS_PER_YEAR
    return (p["T_mean"]
            + p["warming_per_decade"] / 10.0 * t_years
            - p["T_amp"] * math.cos(omega * (t_days - p["coldest_day"])))


def surface_melt(T_air, p, dt_days):
    """Surface melt during one step [kg m^-2] from the degree-day model.

    M = f_m * max(T_air, 0) * dt.  The degree-day factor is in mm w.e. per
    deg C per day, and 1 mm of water equivalent is 1 kg m^-2.
    """
    if not p["enable_melt"]:
        return 0.0
    return p["degree_day_factor"] * max(T_air, 0.0) * dt_days


def surface_freezing(T_air, p, dt_days):
    """Freezing of a surface pond exposed to cold air during one step [kg m^-2].

    Uses the same degree-day factor as melting, applied to negative
    degree-days.  This is a deliberately simple closure (see README).
    """
    if not p["enable_melt"]:
        return 0.0
    return p["degree_day_factor"] * max(-T_air, 0.0) * dt_days
