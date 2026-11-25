// Web Worker: loads Pyodide and the firnmodel Python package, then runs
// simulations on request so the page stays responsive.
//
// Messages in:  {type: "run", params: {...}}
// Messages out: {type: "status", text}, {type: "ready", pyodideVersion},
//               {type: "progress", fraction}, {type: "result", json},
//               {type: "error", message}

const PYODIDE_VERSION = "0.29.5";
const DEFAULT_PYODIDE_URL = `https://cdn.jsdelivr.net/pyodide/v${PYODIDE_VERSION}/full/`;
// Optional override (for offline testing): worker.js?pyodide=<base url>
const PYODIDE_URL = new URL(self.location.href).searchParams.get("pyodide") || DEFAULT_PYODIDE_URL;

// The Python simulation files, copied next to the site by the deploy workflow.
const PY_MODULES = [
  "__init__.py", "params.py", "forcing.py", "heat.py",
  "densification.py", "column.py", "hydrology.py", "model.py",
];

importScripts(PYODIDE_URL + "pyodide.js");

async function setup() {
  self.postMessage({ type: "status", text: "Loading Python (Pyodide)…" });
  const pyodide = await loadPyodide({ indexURL: PYODIDE_URL });

  self.postMessage({ type: "status", text: "Loading NumPy…" });
  await pyodide.loadPackage("numpy");

  self.postMessage({ type: "status", text: "Loading the firn model…" });
  pyodide.FS.mkdirTree("/home/pyodide/firnmodel");
  for (const name of PY_MODULES) {
    const response = await fetch(`firnmodel/${name}`, { cache: "no-cache" });
    if (!response.ok) throw new Error(`Could not load firnmodel/${name} (${response.status})`);
    pyodide.FS.writeFile(`/home/pyodide/firnmodel/${name}`, await response.text());
  }

  // Thin glue: JSON in, JSON out.  All physics lives in the firnmodel package.
  pyodide.runPython(`
import sys, json
sys.path.insert(0, "/home/pyodide")
import firnmodel

def run_json(params_json, progress_callback):
    params = json.loads(params_json)
    results = firnmodel.run(params, progress=progress_callback)
    return json.dumps(firnmodel.to_jsonable(results))
`);
  self.postMessage({ type: "ready", pyodideVersion: pyodide.version });
  return pyodide;
}

const pyodideReady = setup().catch((err) => {
  self.postMessage({ type: "error", message: String(err) });
  throw err;
});

self.onmessage = async (event) => {
  const msg = event.data;
  if (msg.type !== "run") return;
  try {
    const pyodide = await pyodideReady;
    const runJson = pyodide.globals.get("run_json");
    const progress = (fraction) => self.postMessage({ type: "progress", fraction });
    const json = runJson(JSON.stringify(msg.params), progress);
    runJson.destroy();
    self.postMessage({ type: "result", json });
  } catch (err) {
    self.postMessage({ type: "error", message: String(err) });
  }
};
