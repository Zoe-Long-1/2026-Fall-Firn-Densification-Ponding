"""Automated checks of the firnmodel package (run with: pytest)."""

import math
import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

import firnmodel  # noqa: E402
from firnmodel.densification import steady_state_density  # noqa: E402
from firnmodel.heat import conductivity, solve_tridiagonal  # noqa: E402
from firnmodel.params import C_ICE, DAYS_PER_YEAR, SECONDS_PER_DAY  # noqa: E402

RHO, T_MEAN, AMP, T0 = 500.0, -15.0, 10.0, 15.0


def exact_wave(z, t_days):
    kappa = float(conductivity(RHO)) / (RHO * C_ICE)
    omega = 2 * math.pi / (DAYS_PER_YEAR * SECONDS_PER_DAY)
    d = math.sqrt(2 * kappa / omega)
    return T_MEAN - AMP * np.exp(-z / d) * np.cos(omega * (t_days - T0) * SECONDS_PER_DAY - z / d)


def wave_error(dz, dt_days):
    r = firnmodel.run(dict(T_mean=T_MEAN, T_amp=AMP, warming_per_decade=0.0, coldest_day=T0,
                           enable_melt=False, enable_densification=False, enable_accumulation=False,
                           initial_density=RHO, depth=30.0, years=1.0, dz=dz, dt_days=dt_days),
                      T_init=lambda z: exact_wave(z, 0.0))
    f = r["final"]
    t_final = r["series"]["t_years"][-1] * DAYS_PER_YEAR
    return math.sqrt(np.mean((f["T"] - exact_wave(f["z_centers"], t_final)) ** 2))


def test_thomas_algorithm_matches_dense_solve():
    rng = np.random.default_rng(0)
    n = 12
    lower, upper = rng.random(n), rng.random(n)
    diag = 3.0 + rng.random(n)
    rhs = rng.random(n)
    A = np.diag(diag) + np.diag(lower[1:], -1) + np.diag(upper[:-1], 1)
    np.testing.assert_allclose(solve_tridiagonal(lower, diag, upper, rhs), np.linalg.solve(A, rhs))


def test_heat_solver_matches_seasonal_wave():
    assert wave_error(0.25, 1.0) < 5e-3


def test_heat_solver_second_order_in_space():
    e1, e2 = wave_error(0.4, 0.1), wave_error(0.2, 0.1)
    assert 1.8 < math.log2(e1 / e2) < 2.2


def test_densification_reaches_analytical_steady_state():
    r = firnmodel.run(dict(T_mean=-20.0, T_amp=0.0, warming_per_decade=0.0, enable_melt=False,
                           initial_density=350.0, dz=0.25, dt_days=10.0, years=150.0))
    f = r["final"]
    err = f["rho"] - steady_state_density(f["z_centers"], -20.0, 0.3, 350.0)
    assert math.sqrt(np.mean(err ** 2)) < 1.0


@pytest.mark.parametrize("overrides", [
    {},
    dict(T_mean=-14.0, warming_per_decade=0.0, years=30.0),
    dict(T_mean=-10.0, T_amp=14.0, warming_per_decade=0.0, years=30.0),
])
def test_water_mass_is_conserved(overrides):
    s = firnmodel.run(overrides)["summary"]
    assert s["total_melt"] > 0
    assert s["mass_balance_relative_error"] < 1e-10


def test_no_melt_means_no_water():
    s = firnmodel.run(dict(T_mean=-25.0, warming_per_decade=0.0, years=5.0))["summary"]
    assert s["total_melt"] == 0.0
    assert s["liquid_water_final"] == 0.0 and s["pond_max"] == 0.0


def test_temperatures_never_exceed_melting_point():
    r = firnmodel.run(dict(T_mean=-10.0, T_amp=14.0, warming_per_decade=0.0, years=10.0))
    assert max(np.max(s["T"]) for s in r["snapshots"]) <= 1e-9


def test_unknown_parameter_is_rejected():
    with pytest.raises(KeyError):
        firnmodel.make_params({"T_men": -10})
