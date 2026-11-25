"""Working example: run the default scenario without the web interface.

    python examples/run_example.py

Prints a summary (compare with examples/expected_output.txt) and saves
examples/example_output.png.  The default scenario is a firn column at
-16 deg C mean annual temperature that warms by 1 deg C per decade for
40 years: meltwater first refreezes in the firn, then an ice slab forms
and blocks it, and finally water ponds at the surface.
"""

import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

import firnmodel  # noqa: E402


def main():
    results = firnmodel.run()          # all default parameters
    s = results["summary"]
    p = results["params"]

    print("Default scenario:", f"T_mean = {p['T_mean']} degC,",
          f"warming = {p['warming_per_decade']} degC/decade,", f"{p['years']:.0f} years")
    print(f"Total surface melt        {s['total_melt'] / 1000:8.2f} m w.e.")
    print(f"Refrozen in the firn      {100 * s['refrozen_fraction']:8.1f} % of melt")
    print(f"Runoff                    {100 * s['runoff_fraction']:8.1f} % of melt")
    print(f"First surface pond        year {int(s['first_pond_year']) + 1}"
          if s["first_pond_year"] is not None else "First surface pond        none")
    print(f"Deepest pond              {s['pond_max']:8.2f} m")
    print(f"Ice in upper 10 m (end)   {s['ice_upper_final']:8.2f} m")
    print(f"Firn air content          {s['firn_air_content_initial']:8.2f} -> "
          f"{s['firn_air_content_final']:.2f} m")
    print(f"Water mass balance error  {s['mass_balance_relative_error']:8.1e} (relative)")

    # Figure: final profiles and the pond time series
    series = results["series"]
    final = results["final"]
    init = results["initial"]
    init_z = 0.5 * (init["z_edges"][1:] + init["z_edges"][:-1])

    fig, axes = plt.subplots(1, 3, figsize=(13, 4.5))
    axes[0].plot(init["rho"], init_z, "--", color="0.6", label="initial")
    axes[0].plot(final["rho"], final["z_centers"], "k-", label="after 40 years")
    axes[0].axvline(p["rho_impermeable"], color="tab:blue", ls=":", label="impermeable density")
    axes[0].set(xlabel="Density (kg/m³)", ylabel="Depth (m)", ylim=(p["depth"], 0))
    axes[0].legend(fontsize=8)

    axes[1].plot(init["T"], init_z, "--", color="0.6", label="initial")
    axes[1].plot(final["T"], final["z_centers"], "r-", label="after 40 years (midwinter)")
    axes[1].set(xlabel="Temperature (°C)", ylabel="Depth (m)", ylim=(p["depth"], 0))
    axes[1].legend(fontsize=8)

    axes[2].fill_between(series["t_years"], series["pond_depth"], color="tab:blue", alpha=0.4)
    axes[2].set(xlabel="Time (years)", ylabel="Surface pond depth (m)")
    ax2 = axes[2].twinx()
    ax2.plot(series["t_years"], series["ice_upper"], color="k", lw=1)
    ax2.set_ylabel("Ice in upper 10 m (m)")
    axes[2].set_title("Ponds (blue) appear after the ice slab (black) forms", fontsize=10)

    plt.tight_layout()
    out = os.path.join(HERE, "example_output.png")
    plt.savefig(out, dpi=120)
    print("Saved", out)


if __name__ == "__main__":
    main()
