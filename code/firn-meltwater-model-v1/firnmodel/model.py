"""Time loop that couples heat conduction, densification and meltwater.

One time step (operator splitting: each process uses the latest values of
the others, which is accurate because density changes far more slowly than
temperature):

    1. accumulation   new snow is added at the surface
    2. surface BC     T_s = min(T_air, 0), or 0 deg C while a pond exists
    3. heat           theta-method (Crank-Nicolson) solve of the heat equation
    4. pond freezing  heat conducted from a pond into cold firn freezes pond
                      water onto the surface; cold air freezes it from above
    5. surface melt   degree-day melt removes ice from the top
    6. densification  Herron-Langway, exact exponential update per layer
    7. meltwater      refreezing, bucket percolation, perched water, ponding,
                      runoff of the pond
    8. regridding     split/merge the top layer, trim the base

Water mass balance (checked every run):

    melt = liquid stored + pond + refrozen + runoff + lateral + deep drainage
"""

import math
import time

import numpy as np

from . import densification, forcing, heat, hydrology
from .column import initial_column
from .params import (DAYS_PER_YEAR, L_FUSION, RHO_ICE, RHO_WATER,
                     SECONDS_PER_DAY, make_params)

ICE_DEPTH_WINDOW = 10.0   # depth range for the "ice in the upper firn" diagnostic [m]


def _diagnostics(col, p):
    """Column-integrated quantities used in the time series."""
    h = col.h
    edges = np.concatenate(([0.0], np.cumsum(h)))
    centers = 0.5 * (edges[:-1] + edges[1:])
    porosity = np.clip(1.0 - col.rho / RHO_ICE, 0.0, 1.0)
    is_ice = (col.rho >= p["rho_impermeable"]) & (centers < ICE_DEPTH_WINDOW)
    W_irr, _ = hydrology.storage_limits(col.m, col.rho, p["irreducible_saturation"])
    return {
        "firn_air_content": float(np.sum(porosity * h)),           # [m]
        "ice_upper": float(np.sum(h[is_ice])),                     # [m]
        "liquid_water": float(col.W.sum()),                        # [kg m^-2]
        "perched_water": float(np.sum(np.maximum(col.W - W_irr, 0.0))),  # [kg m^-2]
    }


def _snapshot(col, t_days):
    h = col.h
    return {
        "t_years": t_days / DAYS_PER_YEAR,
        "z_edges": col.depth_edges(),
        "T": col.T.copy(),
        "rho": col.rho.copy(),
        "water_percent": 100.0 * (col.W / RHO_WATER) / h,   # volumetric water content [%]
    }


def run(params=None, T_init=None, progress=None):
    """Run the firn model.

    Parameters
    ----------
    params : dict of parameters overriding ``params.DEFAULTS``
    T_init : optional initial temperature, an array (one value per layer) or a
             function of depth z [m]; default is uniform T_mean
    progress : optional callable progress(fraction) for the web interface

    Returns
    -------
    dict of numpy arrays and floats (see README for the full list).
    """
    wall_start = time.time()
    p = make_params(params)
    dt_days = p["dt_days"]
    dt = dt_days * SECONDS_PER_DAY
    dt_years = dt_days / DAYS_PER_YEAR
    n_steps = max(1, int(round(p["years"] * DAYS_PER_YEAR / dt_days)))
    dz = p["dz"]

    col = initial_column(p, T_init)
    initial = _snapshot(col, 0.0)

    pond = 0.0                                   # surface pond [kg m^-2]
    Ts_old = min(forcing.air_temperature(0.0, p), 0.0)
    totals = dict(melt=0.0, refrozen=0.0, runoff=0.0, lateral=0.0, deep=0.0,
                  accumulation=0.0)

    names = ["t_years", "T_air", "T_surface", "melt", "pond_depth", "liquid_water",
             "perched_water", "firn_air_content", "ice_upper",
             "cum_melt", "cum_refrozen", "cum_runoff"]
    series = {name: np.zeros(n_steps + 1) for name in names}

    def record(k, T_air, Ts, melt_step):
        d = _diagnostics(col, p)
        series["t_years"][k] = k * dt_days / DAYS_PER_YEAR
        series["T_air"][k] = T_air
        series["T_surface"][k] = Ts
        series["melt"][k] = melt_step / dt_days            # [mm w.e. per day]
        series["pond_depth"][k] = pond / RHO_WATER          # [m]
        series["liquid_water"][k] = d["liquid_water"]
        series["perched_water"][k] = d["perched_water"]
        series["firn_air_content"][k] = d["firn_air_content"]
        series["ice_upper"][k] = d["ice_upper"]
        series["cum_melt"][k] = totals["melt"]
        series["cum_refrozen"][k] = totals["refrozen"]
        series["cum_runoff"][k] = totals["runoff"] + totals["lateral"]

    record(0, forcing.air_temperature(0.0, p), Ts_old, 0.0)
    snapshots = [initial]
    next_snapshot = p["snapshot_every_days"]
    report_every = max(1, n_steps // 100)

    for step in range(1, n_steps + 1):
        t_new = step * dt_days
        T_air_new = forcing.air_temperature(t_new, p)
        T_air_mid = forcing.air_temperature(t_new - 0.5 * dt_days, p)

        # 1. Accumulation of new snow at the surface temperature
        if p["enable_accumulation"] and p["accumulation"] > 0.0:
            snow = p["accumulation"] * RHO_WATER * dt_years
            col.add_to_top(snow, p["rho_snow"], min(T_air_mid, 0.0), dz)
            totals["accumulation"] += snow

        # 2. Surface boundary condition
        Ts_new = 0.0 if pond > 0.0 else min(T_air_new, 0.0)

        # 3. Heat conduction (the PDE)
        col.T, q_surface = heat.heat_step(col.T, col.m, col.h, col.rho,
                                          Ts_old, Ts_new, dt, p["theta"])

        # 4. Pond freezing: from below by conduction into cold firn,
        #    and from above when the air is below freezing
        if pond > 0.0:
            frozen = min(pond, max(q_surface, 0.0) * dt / L_FUSION)
            frozen += min(pond - frozen, forcing.surface_freezing(T_air_mid, p, dt_days))
            pond -= frozen
            totals["refrozen"] += frozen
            col.add_to_top(frozen, RHO_ICE, 0.0, dz)

        # 5. Surface melt (degree-day) and melt of any layer warmer than 0 deg C
        melt_step = forcing.surface_melt(T_air_mid, p, dt_days)
        water_in = col.remove_from_top(melt_step) if melt_step > 0.0 else 0.0
        internal_melt = hydrology.melt_warm_layers(col)
        totals["melt"] += melt_step + internal_melt

        # 6. Densification
        if p["enable_densification"]:
            col.rho = densification.densify(col.rho, col.T, p["accumulation"], dt_years)
            # Liquid water cannot exceed the shrinking pore space; excess moves on.
            _, W_max = hydrology.storage_limits(col.m, col.rho, 0.0)
            squeezed = np.maximum(col.W - W_max, 0.0)
            col.W = col.W - squeezed
            water_in += float(squeezed.sum())

        # 7. Meltwater: refreeze stored water, then route new water down
        totals["refrozen"] += hydrology.refreeze(col)
        routed = hydrology.percolate(col, water_in + pond, p)
        totals["refrozen"] += routed["refrozen"]
        totals["lateral"] += routed["lateral"]
        totals["deep"] += routed["deep"]
        pond = routed["pond"]
        runoff = pond * (1.0 - math.exp(-p["runoff_rate"] * dt_days))
        pond -= runoff
        totals["runoff"] += runoff

        # 8. Regridding
        col.regrid_top(dz)
        totals["deep"] += col.regrid_bottom(p["depth"], dz)

        Ts_old = Ts_new
        record(step, T_air_new, Ts_new, melt_step + internal_melt)
        if t_new >= next_snapshot - 1e-9 or step == n_steps:
            snapshots.append(_snapshot(col, t_new))
            next_snapshot += p["snapshot_every_days"]
        if progress is not None and (step % report_every == 0 or step == n_steps):
            progress(step / n_steps)

    # ---------------------------------------------------------------- outputs
    stored = float(col.W.sum())
    water_out = (stored + pond + totals["refrozen"] + totals["runoff"]
                 + totals["lateral"] + totals["deep"])
    balance_error = abs(totals["melt"] - water_out)
    relative_error = balance_error / totals["melt"] if totals["melt"] > 0 else 0.0

    pond_series = series["pond_depth"]
    ponded = np.nonzero(pond_series > 1e-3)[0]
    first_pond_year = float(series["t_years"][ponded[0]]) if len(ponded) else None

    annual = _annual_totals(series, dt_days)

    melt_total = totals["melt"]
    summary = {
        "total_melt": melt_total,                                  # [kg m^-2]
        "total_refrozen": totals["refrozen"],
        "total_runoff": totals["runoff"] + totals["lateral"],
        "deep_drainage": totals["deep"],
        "liquid_water_final": stored,
        "pond_final": pond / RHO_WATER,                             # [m]
        "pond_max": float(pond_series.max()),                       # [m]
        "first_pond_year": first_pond_year,
        "refrozen_fraction": totals["refrozen"] / melt_total if melt_total > 0 else 0.0,
        "runoff_fraction": (totals["runoff"] + totals["lateral"]) / melt_total
                           if melt_total > 0 else 0.0,
        "ice_upper_final": float(series["ice_upper"][-1]),          # [m]
        "firn_air_content_initial": float(series["firn_air_content"][0]),
        "firn_air_content_final": float(series["firn_air_content"][-1]),
        "mass_balance_error": balance_error,                        # [kg m^-2]
        "mass_balance_relative_error": relative_error,
        "n_steps": n_steps,
        "n_layers_final": col.n,
        "runtime_seconds": time.time() - wall_start,
    }

    return {
        "params": p,
        "series": series,
        "annual": annual,
        "snapshots": snapshots,
        "initial": initial,
        "final": {
            "z_edges": col.depth_edges(),
            "z_centers": col.depth_centers(),
            "T": col.T.copy(),
            "rho": col.rho.copy(),
            "W": col.W.copy(),
            "h": col.h,
        },
        "summary": summary,
    }


def _annual_totals(series, dt_days):
    """Melt, refreezing and runoff in each model year [kg m^-2 = mm w.e.]."""
    t = series["t_years"]
    n_years = int(math.floor(t[-1] + 1e-9))
    years, melt, refrozen, runoff = [], [], [], []
    for y in range(n_years):
        i0 = int(np.searchsorted(t, y - 1e-9))
        i1 = int(np.searchsorted(t, y + 1 - 1e-9))
        i1 = min(i1, len(t) - 1)
        years.append(y + 1)
        melt.append(series["cum_melt"][i1] - series["cum_melt"][i0])
        refrozen.append(series["cum_refrozen"][i1] - series["cum_refrozen"][i0])
        runoff.append(series["cum_runoff"][i1] - series["cum_runoff"][i0])
    return {"year": np.array(years), "melt": np.array(melt),
            "refrozen": np.array(refrozen), "runoff": np.array(runoff)}


def to_jsonable(obj, decimals=4):
    """Convert results to plain Python lists and floats for JSON (web interface)."""
    if isinstance(obj, dict):
        return {k: to_jsonable(v, decimals) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [to_jsonable(v, decimals) for v in obj]
    if isinstance(obj, np.ndarray):
        return np.round(obj.astype(float), decimals).tolist()
    if isinstance(obj, (np.floating, float)):
        return float(obj) if math.isfinite(obj) else None
    if isinstance(obj, np.integer):
        return int(obj)
    return obj
