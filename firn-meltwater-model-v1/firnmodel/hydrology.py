"""Meltwater in the firn: refreezing, bucket percolation, perched water, ponding.

These are local rules, not a PDE.  Each time step, water from the surface is
passed down the column layer by layer ("bucket scheme"):

1. In each permeable layer the water first refreezes, limited by
   - the water available,
   - the cold content: freezing releases latent heat, which warms the layer;
     refreezing stops when the layer reaches 0 deg C,
   - the pore space: the layer cannot become denser than ice.
2. The layer then keeps up to its irreducible water content (water held by
   capillary forces) and passes the rest down.
3. A layer with density >= rho_impermeable (an ice layer or slab) blocks
   water.  Water arriving on top of it fills the layers above from the bottom
   up to saturation: this is perched water.
4. If the perched water table reaches the surface, the excess becomes a
   surface pond.

Water reaching the base of the column leaves as deep drainage.
"""

import numpy as np

from .params import C_ICE, L_FUSION, RHO_ICE, RHO_WATER


def storage_limits(m, rho, irreducible_saturation):
    """Irreducible and maximum (saturated) water per layer [kg m^-2].

    porosity  phi   = 1 - rho / rho_ice
    W_max           = phi * rho_w * h          (all pores filled)
    W_irr           = S_irr * W_max            (held against gravity)
    """
    h = m / rho
    porosity = np.clip(1.0 - rho / RHO_ICE, 0.0, 1.0)
    W_max = porosity * RHO_WATER * h
    return irreducible_saturation * W_max, W_max


def refreeze(col):
    """Refreeze liquid water in every layer as far as the limits allow.

    m_rf = min( W,  m c (0 - T) / L,  rho_ice h - m )

    Returns the total mass refrozen [kg m^-2].  The latent heat L*m_rf warms
    the layer; this is the source term L*R of the heat equation, applied as a
    separate step after the conduction solve (operator splitting).
    """
    h = col.h
    cold_limit = np.maximum(col.m * C_ICE * (0.0 - col.T), 0.0) / L_FUSION
    pore_limit = np.maximum(RHO_ICE * h - col.m, 0.0)
    m_rf = np.minimum(col.W, np.minimum(cold_limit, pore_limit))
    if not np.any(m_rf > 0.0):
        return 0.0
    m_new = col.m + m_rf
    col.T = (col.m * C_ICE * col.T + L_FUSION * m_rf) / (m_new * C_ICE)
    col.m = m_new
    col.rho = m_new / h            # thickness fixed, so density rises
    col.W = col.W - m_rf
    return float(m_rf.sum())


def melt_warm_layers(col):
    """Convert heat above 0 deg C into melt inside the column.

    The surface temperature never exceeds 0 deg C, but the Crank-Nicolson
    solution can overshoot very slightly.  Instead of discarding that energy,
    it melts ice: m_melt = m c T / L.  Returns the mass melted [kg m^-2].
    """
    warm = col.T > 0.0
    if not np.any(warm):
        return 0.0
    h = col.h
    m_melt = np.where(warm, np.minimum(col.m * C_ICE * col.T / L_FUSION, 0.5 * col.m), 0.0)
    col.m = col.m - m_melt
    col.rho = col.m / h
    col.W = col.W + m_melt
    col.T = np.where(warm, 0.0, col.T)
    return float(m_melt.sum())


def _refreeze_layer(m, T, W, h):
    """Scalar version of ``refreeze`` for one layer, used inside the percolation loop."""
    m_rf = min(W, max(m * C_ICE * (0.0 - T), 0.0) / L_FUSION, max(RHO_ICE * h - m, 0.0))
    if m_rf <= 0.0:
        return m, T, W, 0.0
    m_new = m + m_rf
    T_new = (m * C_ICE * T + L_FUSION * m_rf) / (m_new * C_ICE)
    return m_new, T_new, W - m_rf, m_rf


def percolate(col, water_in, p):
    """Route surface water through the column for one time step.

    Parameters
    ----------
    col : Column (modified in place)
    water_in : surface melt plus existing pond water [kg m^-2]
    p : parameter dictionary

    Returns
    -------
    dict with
        pond     water left at the surface                   [kg m^-2]
        lateral  water that could not be stored below a deeper ice layer and
                 is assumed to flow sideways                   [kg m^-2]
        deep     water draining out of the column base       [kg m^-2]
        refrozen water refrozen during percolation           [kg m^-2]
    """
    rho_imp = p["rho_impermeable"]
    s_irr = p["irreducible_saturation"]
    out = {"pond": 0.0, "lateral": 0.0, "deep": 0.0, "refrozen": 0.0}

    W_irr, _ = storage_limits(col.m, col.rho, s_irr)
    if water_in <= 0.0 and not np.any(col.W > W_irr + 1e-9):
        return out                 # dry, or only immobile water: nothing moves

    m = col.m.tolist()
    rho = col.rho.tolist()
    T = col.T.tolist()
    W = col.W.tolist()
    h = [mi / ri for mi, ri in zip(m, rho)]

    carry = water_in            # water moving down [kg m^-2]
    segment = []                # permeable layers since the last barrier
    first_segment = True        # the segment that starts at the surface

    for i in range(len(m)):
        if rho[i] >= rho_imp:
            # Impermeable layer: perched water fills the segment from the bottom up.
            for j in reversed(segment):
                porosity = max(1.0 - rho[j] / RHO_ICE, 0.0)
                room = porosity * RHO_WATER * h[j] - W[j]
                if room > 0.0:
                    add = min(room, carry)
                    W[j] += add
                    carry -= add
                if carry <= 0.0:
                    break
            if carry > 0.0:
                if first_segment:
                    out["pond"] += carry       # perched water reached the surface
                else:
                    out["lateral"] += carry
            carry = 0.0
            segment = []
            first_segment = False
            continue

        # Permeable layer: add water, refreeze, keep the irreducible part.
        W[i] += carry
        m[i], T[i], W[i], m_rf = _refreeze_layer(m[i], T[i], W[i], h[i])
        out["refrozen"] += m_rf
        rho[i] = m[i] / h[i]
        porosity = max(1.0 - rho[i] / RHO_ICE, 0.0)
        w_irr = s_irr * porosity * RHO_WATER * h[i]
        carry = max(W[i] - w_irr, 0.0)
        W[i] -= carry
        segment.append(i)

    out["deep"] += carry       # water leaving the base of the column

    col.m = np.array(m)
    col.rho = np.array(rho)
    col.T = np.array(T)
    col.W = np.array(W)
    return out
