# AI use statement

> **Edit this file before submitting** so it describes what you actually did. The text below records
> how the code was produced; the section on checking is a starting point for your own account.

## Tools used

- **Claude (Anthropic), Opus 5.5 model, through claude.ai**, used as an AI coding agent.

## How the AI supported development

- **Project planning:** discussing feasibility, choosing the governing PDE (1D heat equation, parabolic),
  the numerical method (Crank–Nicolson finite volumes) and the verification test (seasonal temperature
  wave), and scoping the meltwater physics to a bucket scheme.
- **Code generation:** the AI wrote the first version of the `firnmodel` Python package, the web
  interface (`web/`), the Pyodide Web Worker integration, the GitHub Actions workflow, the tests and
  the verification notebook, following the course requirement to keep the simulation separate from the
  interface.
- **Testing and debugging:** the AI ran the model and tests, checked convergence rates against
  theory, and tested the web app in a headless browser with Pyodide, fixing the bugs this found.
- **Documentation:** README, code comments and this statement.

## What I checked myself

*(Replace with your own account, for example:)*

- I worked through the discretization in `firnmodel/heat.py` (control-volume balance, interface
  conductances, surface boundary at half a layer, Thomas algorithm) and confirmed it matches the
  equations in the README.
- I re-ran the verification notebook and confirmed the observed convergence orders.
- I checked the Herron–Langway and Calonne et al. constants against the original papers.
- I ran the app scenarios and checked that the results are physically sensible.

## External code and references

No external simulation code was used; the solver is written from the equations. Libraries: NumPy,
Matplotlib (notebook), Pyodide (Python in the browser), Plotly.js (browser plots). Scientific references
are listed in the README.
