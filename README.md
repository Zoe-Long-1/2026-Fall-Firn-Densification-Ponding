# Firn meltwater model

**When does firn stop absorbing meltwater and start ponding it?**

A browser-based 1D model of a firn column on an ice sheet. It solves the heat equation with
Crank–Nicolson finite volumes, coupled to firn densification, meltwater percolation and refreezing,
ice-slab formation, perched water and surface ponds. Python runs in the browser with Pyodide in a
Web Worker; there is no server.

- **Live app:** `https://<your-username>.github.io/<repository-name>/` *(replace after deploying)*
- **Verification notebook:** [`verification/verification.ipynb`](verification/verification.ipynb)
- **AI use statement:** [`AI_STATEMENT.md`](AI_STATEMENT.md)

## The geophysical problem

Much of the Greenland ice sheet is covered by firn, porous snow tens of metres thick. Firn acts as a
sponge: summer meltwater percolates into it and refreezes, so it never reaches the ocean. Each
refreezing event leaves ice behind, though, and after warm summers these ice layers can merge into
**ice slabs** metres thick that water cannot pass through. Meltwater then perches on the slab, fills
the firn above it, and finally ponds and runs off at the surface (Machguth et al., 2016;
MacFerrin et al., 2019). This model shows that transition in one column: from full retention, to
perched water, to slab growth, to surface ponds, as the climate warms.

## The model

### Governing PDE (parabolic)

$$\rho c \frac{\partial T}{\partial t} = \frac{\partial}{\partial z}\left(k(\rho)\frac{\partial T}{\partial z}\right) + L R$$

| Symbol | Meaning | Units |
|---|---|---|
| $T(z,t)$ | firn temperature | °C |
| $z$ | depth below the surface, positive down | m |
| $\rho$ | firn density (ice matrix, excluding liquid water) | kg m⁻³ |
| $c$ | specific heat of ice, 2009 | J kg⁻¹ K⁻¹ |
| $k(\rho)$ | thermal conductivity, $2.5\times10^{-6}\rho^2 - 1.23\times10^{-4}\rho + 0.024$ (Calonne et al., 2011) | W m⁻¹ K⁻¹ |
| $L$ | latent heat of fusion, $3.34\times10^5$ | J kg⁻¹ |
| $R$ | refreezing rate | kg m⁻³ s⁻¹ |

**Boundary conditions:** surface temperature $T_s = \min(T_{air}, 0\,°\text{C})$, or $0\,°$C while a pond
exists (Dirichlet); zero heat flux at the base $z = H$ (Neumann).
**Initial condition:** uniform $T = T_{mean}$; density from the Herron–Langway steady state.

**Forcing:** $T_{air}(t) = T_{mean} + \gamma t - A\cos\big(2\pi (t - t_0)/365\,\text{d}\big)$, with warming
trend $\gamma$ and the coldest day $t_0$ in mid-January.

### Numerical method (`firnmodel/heat.py`)

Finite volumes on layers that move down with the firn (a Lagrangian grid, so no advection term is
needed). For layer $i$ with ice mass $m_i$ and thickness $h_i$, the theta method gives

$$m_i c\,\frac{T_i^{n+1} - T_i^n}{\Delta t} = \theta F_i^{n+1} + (1-\theta)F_i^n, \qquad
F_i = G_{i-\frac12}(T_{i-1} - T_i) - G_{i+\frac12}(T_i - T_{i+1}),$$

with interface conductances $G_{i+\frac12} = \big(h_i/2k_i + h_{i+1}/2k_{i+1}\big)^{-1}$, surface
conductance $G_s = 2k_0/h_0$, and $G = 0$ at the base. $\theta = 0.5$ is Crank–Nicolson (default);
$\theta = 1$ is backward Euler. The tridiagonal system is solved with the Thomas algorithm.

### Coupled processes (operator splitting, each time step)

| Process | Module | Rule |
|---|---|---|
| Accumulation | `column.py` | snow at $\rho_{snow}$ added at the surface; thick top layers split, thin ones merge |
| Densification | `densification.py` | Herron–Langway (1980): $d\rho/dt = c_0 b(\rho_i - \rho)$ below 550 kg m⁻³, $c_1 b^{1/2}(\rho_i-\rho)$ above; integrated exactly over each step |
| Surface melt | `forcing.py` | degree-day: $M = f_m \max(T_{air}, 0)$ |
| Refreezing | `hydrology.py` | $m_{rf} = \min\big(W,\ m c(0 - T)/L,\ \rho_i h - m\big)$: limited by water, cold content and pore space; releases latent heat (the $LR$ term) |
| Percolation | `hydrology.py` | bucket scheme: each layer keeps its irreducible water $S_{irr}\,\phi\,\rho_w h$ and passes the rest down |
| Perched water | `hydrology.py` | layers with $\rho \ge \rho_{imp}$ block water, which saturates the firn above them from the bottom up |
| Ponds | `model.py` | water that reaches the surface ponds; drains at rate $r$; freezes from below by conduction into cold firn and from above at negative degree-days |

Density changes over years while temperature changes over days, so updating them in sequence within
each step is accurate: halving the time step from 1 to 0.5 day changes the default scenario's melt,
refrozen fraction and first pond year by less than 0.2%, and its final ice thickness by about 6%.

## Repository layout

```
firnmodel/            standalone simulation package (import it from any script)
  params.py           constants and default parameters
  forcing.py          air temperature and degree-day melt
  heat.py             heat-equation solver: Crank–Nicolson + Thomas algorithm
  densification.py    Herron–Langway densification and its analytical steady state
  column.py           Lagrangian layers, initial conditions, regridding
  hydrology.py        refreezing, percolation, perched water, ponding
  model.py            time loop coupling everything; run(params) -> results
web/                  user interface (no physics)
  index.html, style.css, main.js
  worker.js           Web Worker: loads Pyodide and runs firnmodel
verification/
  verification.ipynb  analytical verification, executed with outputs
tests/                pytest checks (also run by the deploy workflow)
examples/             working example script, expected output and figure
.github/workflows/deploy.yml   tests, then deploys to GitHub Pages
```

## Running it

### In the browser
Open the GitHub Pages link. The first load downloads Pyodide and NumPy (about 15 MB) and then runs the
default scenario; later runs take about 5–10 s for 40 model years. Choose a scenario or move the
sliders, press **Run simulation**, and use the **Time** slider to step through the depth profiles.

### Standalone Python
```bash
pip install -r requirements.txt
python examples/run_example.py          # working example (see below)
pytest -q tests                         # automated checks, about 15 s
jupyter notebook verification/verification.ipynb
```

```python
import firnmodel
results = firnmodel.run({"T_mean": -12.0, "warming_per_decade": 0.5, "years": 30})
print(results["summary"])
```

### Local copy of the website
```bash
mkdir -p _site/firnmodel && cp -r web/. _site/ && cp firnmodel/*.py _site/firnmodel/
python -m http.server --directory _site 8000    # then open http://localhost:8000
```
(Opening `index.html` directly from disk does not work, because the worker has to fetch the Python files.)

### Deploying to GitHub Pages
1. Push this repository to GitHub (branch `main`).
2. In **Settings → Pages**, set **Source** to **GitHub Actions**.
3. The workflow runs the tests, assembles `web/` and `firnmodel/` into a static site, and deploys it.
   The site URL appears in the workflow run and under Settings → Pages.

## Dependencies and versions

| Where | Package | Version |
|---|---|---|
| Browser | Pyodide (from cdn.jsdelivr.net) | 0.29.5, pinned in `web/worker.js` |
| Browser | NumPy (loaded by Pyodide) | 2.2.5, the version shipped with Pyodide 0.29.5 |
| Browser | Plotly.js (from cdn.plot.ly) | 2.35.2, pinned in `web/index.html` |
| Local | Python, NumPy, Matplotlib, Jupyter, pytest | see `requirements.txt`; tested with Python 3.13, NumPy 2.5.3 |

The simulation needs only NumPy, so the same code runs in Pyodide and in standard Python.

## Working example and expected results

`python examples/run_example.py` runs the default scenario (−16 °C mean annual temperature, 15 °C
seasonal amplitude, 0.3 m w.e./yr accumulation, warming 1 °C per decade, 40 years) and should print the
contents of [`examples/expected_output.txt`](examples/expected_output.txt):

```
Total surface melt            8.91 m w.e.
Refrozen in the firn          63.4 % of melt
Runoff                        36.6 % of melt
First surface pond        year 27
Deepest pond                  0.20 m
Ice in upper 10 m (end)       3.73 m
Firn air content             11.24 -> 6.30 m
```

The web app shows the same numbers for its default scenario. **Interpretation:** for the first ~20
years all meltwater refreezes in the firn, warming it through latent heat and filling its pores. Around
year 23 refrozen layers reach the impermeable density and merge into an ice slab 1–4 m below the
surface. From then on, water perches on the slab; from year 27 the perched water reaches the surface
every summer and ponds form, and a growing share of the melt runs off instead of being stored.

### App scenarios

| Scenario | What it shows |
|---|---|
| Warming: slab, then ponds (default) | the full transition described above |
| Dry snow zone | no melt; seasonal temperature wave and densification only |
| Percolation zone | all melt refreezes; firn air content falls but no slab forms in 30 years |
| Runoff zone with surface lakes | slab forms within a decade; most melt runs off; ponds every summer |

## Verification

The notebook [`verification/verification.ipynb`](verification/verification.ipynb) imports the same
`firnmodel` code as the app and shows:

1. **Heat equation vs. analytical seasonal wave** $T = T_{mean} - Ae^{-z/d}\cos(\omega(t-t_0) - z/d)$,
   with melt, densification and accumulation switched off: RMS error ≈ 3×10⁻³ °C at the app's
   resolution, **second-order convergence in Δz and Δt** for Crank–Nicolson, first order in Δt for
   backward Euler.
2. **Densification vs. the analytical Herron–Langway steady state**, built up from a uniform initial
   column: RMS error below 1 kg m⁻³.
3. **Water mass conservation** in the full model: relative error ~10⁻¹⁵.

## Limitations

- **1D:** a pond has a depth but no area; lateral flow is represented only by the runoff rate.
- **Bucket percolation:** water moves through permeable firn within one time step; the timing of
  wetting fronts and preferential (piping) flow are not represented. Richards-equation flow is not used.
- **Ice layers form only where refreezing raises a layer above the impermeable density**; thin ice
  lenses below the grid resolution are averaged into their layer.
- **Degree-day melt** stands in for a full surface energy balance, and there is no daily temperature cycle.
- **Pond freezing** uses the same degree-day factor; absorbed sunlight in ponds is ignored.
- **Herron–Langway** was calibrated for dry firn at mean annual temperatures; using it with seasonal
  temperatures and wet firn is a common but approximate choice.
- The initial temperature is uniform, so the first few years include spin-up.

## References

- Calonne, N., Flin, F., Morin, S., Lesaffre, B., Rolland du Roscoat, S., and Geindreau, C. (2011).
  Numerical and experimental investigations of the effective thermal conductivity of snow.
  *Geophysical Research Letters*, 38, L23501.
- Carslaw, H. S., and Jaeger, J. C. (1959). *Conduction of Heat in Solids*, 2nd ed. Oxford University Press.
- Coléou, C., and Lesaffre, B. (1998). Irreducible water saturation in snow: experimental results in a
  cold laboratory. *Annals of Glaciology*, 26, 64–68.
- Cuffey, K. M., and Paterson, W. S. B. (2010). *The Physics of Glaciers*, 4th ed. Academic Press.
- Herron, M. M., and Langway, C. C. (1980). Firn densification: an empirical model.
  *Journal of Glaciology*, 25(93), 373–385.
- MacFerrin, M., et al. (2019). Rapid expansion of Greenland's low-permeability ice slabs.
  *Nature*, 573, 403–407.
- Machguth, H., et al. (2016). Greenland meltwater storage in firn limited by near-surface ice formation.
  *Nature Climate Change*, 6, 390–393.
- Pyodide: https://pyodide.org · Plotly.js: https://plotly.com/javascript/

Check constants and citations against the original papers before relying on them.
