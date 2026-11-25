"""The firn column: a stack of Lagrangian layers.

Each layer i stores
    m[i]    ice mass                         [kg m^-2]
    rho[i]  density of the ice matrix        [kg m^-3]  (excludes liquid water)
    T[i]    temperature                      [deg C]
    W[i]    liquid water                     [kg m^-2]
and its thickness follows from h = m / rho.

Layers move down with the firn: new snow is added at the top, and layers are
removed at the bottom once the column is deeper than H.  Because the grid
moves with the material, the heat equation needs no advection term and
densification is an ODE for each layer.
"""

import numpy as np

from .densification import steady_state_density


class Column:
    def __init__(self, m, rho, T, W):
        self.m = np.asarray(m, dtype=float).copy()
        self.rho = np.asarray(rho, dtype=float).copy()
        self.T = np.asarray(T, dtype=float).copy()
        self.W = np.asarray(W, dtype=float).copy()

    # ------------------------------------------------------------------ geometry
    @property
    def n(self):
        return len(self.m)

    @property
    def h(self):
        """Layer thicknesses [m]."""
        return self.m / self.rho

    def depth_edges(self):
        """Depths of layer boundaries, from 0 (surface) to the column base [m]."""
        return np.concatenate(([0.0], np.cumsum(self.h)))

    def depth_centers(self):
        edges = self.depth_edges()
        return 0.5 * (edges[:-1] + edges[1:])

    # --------------------------------------------------------- top of the column
    def add_to_top(self, mass, rho_new, T_new, dz_target):
        """Add material (snow or refrozen pond ice) at the surface.

        Thick additions become a new layer; thin ones are merged into the top
        layer, conserving mass and heat content.
        """
        if mass <= 0.0:
            return
        h_new = mass / rho_new
        if h_new >= 0.5 * dz_target:
            self.m = np.insert(self.m, 0, mass)
            self.rho = np.insert(self.rho, 0, rho_new)
            self.T = np.insert(self.T, 0, T_new)
            self.W = np.insert(self.W, 0, 0.0)
            return
        h0 = self.m[0] / self.rho[0]
        m_total = self.m[0] + mass
        self.T[0] = (self.m[0] * self.T[0] + mass * T_new) / m_total
        self.rho[0] = m_total / (h0 + h_new)
        self.m[0] = m_total

    def remove_from_top(self, mass):
        """Melt ``mass`` of ice from the surface downward.

        Returns the liquid water released [kg m^-2]: the melted ice plus any
        liquid water held in layers that melt away completely.
        """
        released = 0.0
        while mass > 0.0 and self.n > 2:
            if self.m[0] > mass:
                self.m[0] -= mass          # density unchanged, layer gets thinner
                released += mass
                mass = 0.0
            else:
                released += self.m[0] + self.W[0]
                mass -= self.m[0]
                self._delete(0)
        return released

    def regrid_top(self, dz_target):
        """Split a top layer thicker than 2*dz, or merge one thinner than 0.3*dz."""
        h0 = self.m[0] / self.rho[0]
        if h0 > 2.0 * dz_target:
            half_m, half_W = 0.5 * self.m[0], 0.5 * self.W[0]
            self.m[0], self.W[0] = half_m, half_W
            self.m = np.insert(self.m, 0, half_m)
            self.rho = np.insert(self.rho, 0, self.rho[0])
            self.T = np.insert(self.T, 0, self.T[0])
            self.W = np.insert(self.W, 0, half_W)
        elif h0 < 0.3 * dz_target and self.n > 2:
            h1 = self.m[1] / self.rho[1]
            m_total = self.m[0] + self.m[1]
            self.T[1] = (self.m[0] * self.T[0] + self.m[1] * self.T[1]) / m_total
            self.rho[1] = m_total / (h0 + h1)
            self.m[1] = m_total
            self.W[1] += self.W[0]
            self._delete(0)

    # ------------------------------------------------------ bottom of the column
    def regrid_bottom(self, depth, dz_target):
        """Keep the column close to ``depth``.

        Layers that move below the base are removed; if melt has thinned the
        column, a layer copying the bottom one is appended.  Returns the
        liquid water carried out of the base [kg m^-2].
        """
        lost_water = 0.0
        h = self.h
        total = h.sum()
        while self.n > 3 and total - h[-1] >= depth:
            total -= h[-1]
            lost_water += self.W[-1]
            self._delete(-1)
            h = self.h
        while total < depth - 0.5 * dz_target:
            rho_b = self.rho[-1]
            self.m = np.append(self.m, rho_b * dz_target)
            self.rho = np.append(self.rho, rho_b)
            self.T = np.append(self.T, self.T[-1])
            self.W = np.append(self.W, 0.0)
            total += dz_target
        return lost_water

    def _delete(self, i):
        self.m = np.delete(self.m, i)
        self.rho = np.delete(self.rho, i)
        self.T = np.delete(self.T, i)
        self.W = np.delete(self.W, i)


def initial_column(p, T_init=None):
    """Build the initial column: uniform layers of thickness close to dz.

    Density is the Herron-Langway steady-state profile (or uniform if
    ``initial_density`` is a number).  Temperature is uniform at T_mean
    unless ``T_init`` is given as an array or as a function of depth.
    """
    n = max(3, int(round(p["depth"] / p["dz"])))
    h = np.full(n, p["depth"] / n)
    z_centers = np.cumsum(h) - 0.5 * h

    if p["initial_density"] == "steady" and p["accumulation"] > 0:
        rho = steady_state_density(z_centers, p["T_mean"], p["accumulation"],
                                   p["surface_density"])
    elif p["initial_density"] == "steady":
        rho = np.full(n, p["surface_density"])
    else:
        rho = np.full(n, float(p["initial_density"]))

    if T_init is None:
        T = np.full(n, p["T_mean"])
    elif callable(T_init):
        T = np.asarray(T_init(z_centers), dtype=float)
    else:
        T = np.asarray(T_init, dtype=float)
        if T.shape != (n,):
            raise ValueError(f"T_init must have {n} values, one per layer")

    return Column(m=rho * h, rho=rho, T=T, W=np.zeros(n))
