"""Physical constants and default model parameters.

All user-facing parameters live in one plain dictionary so the same values can
come from the web interface (as JSON), a test script, or the verification
notebook.  ``make_params`` fills in defaults and checks for typos.
"""

# ---------------------------------------------------------------------------
# Physical constants (SI units)
# ---------------------------------------------------------------------------
RHO_ICE = 917.0          # density of glacier ice                [kg m^-3]
RHO_WATER = 1000.0       # density of liquid water               [kg m^-3]
C_ICE = 2009.0           # specific heat capacity of ice          [J kg^-1 K^-1]
L_FUSION = 3.34e5        # latent heat of fusion of ice           [J kg^-1]
R_GAS = 8.314            # universal gas constant                 [J mol^-1 K^-1]
T_MELT = 0.0             # melting temperature                    [deg C]

SECONDS_PER_DAY = 86400.0
DAYS_PER_YEAR = 365.0
SECONDS_PER_YEAR = SECONDS_PER_DAY * DAYS_PER_YEAR

# ---------------------------------------------------------------------------
# Default parameters (units given in the comments)
# ---------------------------------------------------------------------------
DEFAULTS = {
    # Climate forcing: T_air(t) = T_mean + trend*t - T_amp*cos(2*pi*(t - coldest_day)/365)
    "T_mean": -16.0,               # mean annual air temperature          [deg C]
    "T_amp": 15.0,                 # seasonal half-range of air temperature [deg C]
    "warming_per_decade": 1.0,     # linear warming trend                 [deg C per decade]
    "coldest_day": 15.0,           # day of year of the temperature minimum [day]
    "accumulation": 0.3,           # snow accumulation rate               [m water eq. per year]
    "rho_snow": 350.0,             # density of fresh surface snow        [kg m^-3]

    # Meltwater
    "degree_day_factor": 5.0,      # melt per positive degree-day         [mm w.e. per deg C per day]
    "runoff_rate": 0.05,           # fraction of the surface pond draining per day [1/day]
    "irreducible_saturation": 0.05,  # fraction of pore space held by capillarity [-]
    "rho_impermeable": 830.0,      # density above which firn blocks water [kg m^-3]

    # Numerical grid and time stepping
    "depth": 30.0,                 # column depth H                       [m]
    "dz": 0.25,                    # target layer thickness               [m]
    "dt_days": 1.0,                # time step                            [day]
    "years": 40.0,                 # length of the run                    [year]
    "theta": 0.5,                  # 0.5 = Crank-Nicolson, 1.0 = backward Euler [-]

    # Switches (used to reduce the model to simple verification problems)
    "enable_melt": True,
    "enable_densification": True,
    "enable_accumulation": True,

    # Initial condition for density: "steady" = Herron-Langway steady-state
    # profile for T_mean and accumulation, or a number for a uniform density.
    "initial_density": "steady",
    "surface_density": 350.0,      # surface density for the steady-state profile [kg m^-3]

    # Output
    "snapshot_every_days": 30.0,   # interval between stored depth profiles [day]
}


def make_params(overrides=None, **kwargs):
    """Return a full parameter dictionary: defaults updated by ``overrides``.

    Unknown keys raise a ``KeyError`` so that a misspelled parameter cannot be
    silently ignored.
    """
    p = dict(DEFAULTS)
    updates = dict(overrides or {})
    updates.update(kwargs)
    for key, value in updates.items():
        if key not in DEFAULTS:
            raise KeyError(f"Unknown parameter '{key}'")
        p[key] = value

    # Coerce numeric types (JSON from the browser may give ints or strings).
    for key, default in DEFAULTS.items():
        if isinstance(default, bool):
            p[key] = bool(p[key])
        elif isinstance(default, float):
            p[key] = float(p[key])
    if p["initial_density"] != "steady":
        p["initial_density"] = float(p["initial_density"])

    if not 0.5 <= p["theta"] <= 1.0:
        raise ValueError("theta must be between 0.5 (Crank-Nicolson) and 1 (backward Euler)")
    if p["dz"] <= 0 or p["dt_days"] <= 0 or p["depth"] <= 2 * p["dz"]:
        raise ValueError("dz and dt_days must be positive and depth must exceed 2*dz")
    return p
