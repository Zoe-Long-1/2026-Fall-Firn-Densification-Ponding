"""Heat conduction in the firn column: the PDE solved by the course method.

Governing equation (1D, depth z positive downward):

    rho * c * dT/dt = d/dz ( k(rho) * dT/dz ) + L * R

The latent-heat source L*R from refreezing is handled separately in
``hydrology.refreeze`` (operator splitting).  This module solves the
conduction part with the theta method on a finite-volume grid:

    theta = 0.5  -> Crank-Nicolson   (second order in time)
    theta = 1.0  -> backward Euler   (first order in time)

Each layer i is a control volume with ice mass m_i [kg m^-2], thickness h_i [m]
and temperature T_i at its centre.  Discretising the PDE over one layer:

    C_i (T_i^{n+1} - T_i^n) / dt = theta * F_i^{n+1} + (1 - theta) * F_i^n

    C_i = m_i * c                                   (heat capacity, J m^-2 K^-1)
    F_i = G_{i-1/2} (T_{i-1} - T_i) - G_{i+1/2} (T_i - T_{i+1})   (net heat flux, W m^-2)

G_{i+1/2} is the thermal conductance between neighbouring layer centres,
the two half-layers in series:

    G_{i+1/2} = 1 / ( h_i / (2 k_i) + h_{i+1} / (2 k_{i+1}) )

Boundary conditions:
    top    (z = 0): Dirichlet, T = T_s(t).  The surface sits half a layer
                    above the first centre, so G_s = 2 k_0 / h_0.
    bottom (z = H): zero heat flux, G_{N-1/2} = 0.

The resulting linear system is tridiagonal and is solved with the Thomas
algorithm below.
"""

import numpy as np

from .params import C_ICE


def conductivity(rho):
    """Effective thermal conductivity of firn [W m^-1 K^-1].

    Calonne et al. (2011): k = 2.5e-6 rho^2 - 1.23e-4 rho + 0.024, rho in kg m^-3.
    Gives about 0.26 W/m/K for fresh snow (350 kg/m3) and 2.0 for ice.
    """
    rho = np.asarray(rho, dtype=float)
    return 2.5e-6 * rho**2 - 1.23e-4 * rho + 0.024


def solve_tridiagonal(lower, diag, upper, rhs):
    """Solve a tridiagonal system with the Thomas algorithm.

    lower[i] multiplies x[i-1] (lower[0] unused), diag[i] multiplies x[i],
    upper[i] multiplies x[i+1] (upper[-1] unused).  The matrix here is
    diagonally dominant, so no pivoting is needed.
    """
    a = list(lower)
    b = list(diag)
    c = list(upper)
    d = list(rhs)
    n = len(b)
    c_prime = [0.0] * n
    d_prime = [0.0] * n

    # Forward elimination
    c_prime[0] = c[0] / b[0]
    d_prime[0] = d[0] / b[0]
    for i in range(1, n):
        denom = b[i] - a[i] * c_prime[i - 1]
        c_prime[i] = c[i] / denom
        d_prime[i] = (d[i] - a[i] * d_prime[i - 1]) / denom

    # Back substitution
    x = [0.0] * n
    x[n - 1] = d_prime[n - 1]
    for i in range(n - 2, -1, -1):
        x[i] = d_prime[i] - c_prime[i] * x[i + 1]
    return np.array(x)


def conductances(h, k):
    """Surface conductance G_s and interface conductances G_{i+1/2} [W m^-2 K^-1].

    Returns (G_s, G) where G has length N; G[i] couples layer i to i+1 and
    G[N-1] = 0 enforces the zero-flux bottom boundary.
    """
    G = np.zeros_like(h)
    G[:-1] = 1.0 / (0.5 * h[:-1] / k[:-1] + 0.5 * h[1:] / k[1:])
    G_s = 2.0 * k[0] / h[0]
    return G_s, G


def net_flux(T, Ts, G_s, G):
    """Net conductive heat flux into each layer, F_i [W m^-2]."""
    F = np.zeros_like(T)
    flux_down = G[:-1] * (T[:-1] - T[1:])   # flux from layer i to i+1
    F[:-1] -= flux_down
    F[1:] += flux_down
    F[0] += G_s * (Ts - T[0])               # flux from the surface into layer 0
    return F


def heat_step(T, m, h, rho, Ts_old, Ts_new, dt, theta=0.5):
    """Advance temperature one time step with the theta method.

    Parameters
    ----------
    T : array, layer temperatures at time n [deg C]
    m : array, ice mass per layer [kg m^-2]
    h : array, layer thickness [m]
    rho : array, layer density [kg m^-3] (sets conductivity)
    Ts_old, Ts_new : surface temperature at time n and n+1 [deg C]
    dt : time step [s]
    theta : 0.5 for Crank-Nicolson, 1.0 for backward Euler

    Returns
    -------
    T_new : array, temperatures at time n+1 [deg C]
    q_surface : mean heat flux from the surface into the firn over the step [W m^-2]
    """
    k = conductivity(rho)
    G_s, G = conductances(h, k)
    C = m * C_ICE                       # heat capacity of each layer [J m^-2 K^-1]

    G_above = np.concatenate(([G_s], G[:-1]))   # conductance to the layer above
    G_below = G                                  # conductance to the layer below

    # Implicit part (time n+1): tridiagonal matrix
    diag = C / dt + theta * (G_above + G_below)
    lower = -theta * G_above          # coefficient of T_{i-1}
    upper = -theta * G_below          # coefficient of T_{i+1}
    lower[0] = 0.0                    # surface handled through the right-hand side

    # Explicit part (time n) plus the known surface temperature at n+1
    rhs = C / dt * T + (1.0 - theta) * net_flux(T, Ts_old, G_s, G)
    rhs[0] += theta * G_s * Ts_new

    T_new = solve_tridiagonal(lower, diag, upper, rhs)

    q_surface = G_s * (theta * (Ts_new - T_new[0]) + (1.0 - theta) * (Ts_old - T[0]))
    return T_new, q_surface
