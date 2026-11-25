// User interface: builds the controls, talks to the Web Worker, draws plots.
// No physics here; the simulation is the Python package in firnmodel/.

"use strict";

// ---------------------------------------------------------------- controls
const CONTROL_GROUPS = [
  {
    title: "Climate",
    controls: [
      { key: "T_mean", label: "Mean annual air temperature", unit: "°C", min: -30, max: -2, step: 0.5, decimals: 1 },
      { key: "T_amp", label: "Seasonal amplitude", unit: "°C", min: 0, max: 25, step: 0.5, decimals: 1,
        hint: "Summer peak = mean + amplitude" },
      { key: "warming_per_decade", label: "Warming trend", unit: "°C/decade", min: 0, max: 2, step: 0.1, decimals: 1 },
      { key: "accumulation", label: "Snow accumulation", unit: "m w.e./yr", min: 0.05, max: 1.0, step: 0.05, decimals: 2 },
    ],
  },
  {
    title: "Meltwater",
    controls: [
      { key: "degree_day_factor", label: "Degree-day melt factor", unit: "mm w.e./°C/day", min: 1, max: 10, step: 0.5, decimals: 1 },
      { key: "runoff_rate", label: "Pond runoff rate", unit: "1/day", min: 0, max: 0.5, step: 0.01, decimals: 2,
        hint: "Fraction of the pond draining per day; stands in for surface slope" },
      { key: "irreducible_saturation", label: "Irreducible water saturation", unit: "–", min: 0.01, max: 0.1, step: 0.01, decimals: 2 },
      { key: "rho_impermeable", label: "Impermeable density", unit: "kg/m³", min: 780, max: 900, step: 5, decimals: 0 },
    ],
  },
  {
    title: "Run and numerics",
    controls: [
      { key: "years", label: "Run length", unit: "years", min: 1, max: 60, step: 1, decimals: 0 },
      { key: "dz", label: "Layer thickness", unit: "m", min: 0.1, max: 1.0, step: 0.05, decimals: 2 },
      { key: "dt_days", label: "Time step", unit: "days", min: 0.25, max: 5, step: 0.25, decimals: 2 },
      { key: "theta", label: "Time stepping", type: "select",
        options: [["0.5", "Crank–Nicolson (θ = 0.5)"], ["1", "Backward Euler (θ = 1)"]] },
    ],
  },
];

const BASE = {
  T_mean: -16, T_amp: 15, warming_per_decade: 1.0, accumulation: 0.3,
  degree_day_factor: 5, runoff_rate: 0.05, irreducible_saturation: 0.05, rho_impermeable: 830,
  years: 40, dz: 0.25, dt_days: 1, theta: 0.5,
};

const PRESETS = [
  { name: "Warming: slab, then ponds", values: {} },
  { name: "Dry snow zone (no melt)", values: { T_mean: -25, warming_per_decade: 0, years: 30 } },
  { name: "Percolation zone (all refreezes)", values: { T_mean: -14, warming_per_decade: 0, years: 30 } },
  { name: "Runoff zone with surface lakes", values: { T_mean: -10, T_amp: 14, warming_per_decade: 0, years: 30 } },
];

const state = { params: { ...BASE }, results: null, busy: false, ready: false };

function buildControls() {
  const container = document.getElementById("control-groups");
  for (const group of CONTROL_GROUPS) {
    const section = document.createElement("section");
    const h = document.createElement("h2");
    h.textContent = group.title;
    section.appendChild(h);
    for (const c of group.controls) {
      const wrap = document.createElement("div");
      wrap.className = "control";
      const label = document.createElement("label");
      label.htmlFor = `ctl-${c.key}`;
      const name = document.createElement("span");
      name.textContent = c.unit ? `${c.label} (${c.unit})` : c.label;
      label.appendChild(name);
      let input;
      if (c.type === "select") {
        input = document.createElement("select");
        for (const [value, text] of c.options) {
          const opt = document.createElement("option");
          opt.value = value;
          opt.textContent = text;
          input.appendChild(opt);
        }
        input.addEventListener("change", () => { state.params[c.key] = Number(input.value); });
      } else {
        const out = document.createElement("output");
        out.id = `out-${c.key}`;
        label.appendChild(out);
        input = document.createElement("input");
        input.type = "range";
        input.min = c.min; input.max = c.max; input.step = c.step;
        input.addEventListener("input", () => {
          state.params[c.key] = Number(input.value);
          out.textContent = Number(input.value).toFixed(c.decimals);
        });
      }
      input.id = `ctl-${c.key}`;
      wrap.appendChild(label);
      wrap.appendChild(input);
      if (c.hint) {
        const hint = document.createElement("p");
        hint.className = "hint";
        hint.textContent = c.hint;
        wrap.appendChild(hint);
      }
      section.appendChild(wrap);
    }
    container.appendChild(section);
  }

  const presetBox = document.getElementById("presets");
  for (const preset of PRESETS) {
    const b = document.createElement("button");
    b.type = "button";
    b.textContent = preset.name;
    b.addEventListener("click", () => { setParams({ ...BASE, ...preset.values }); runSimulation(); });
    presetBox.appendChild(b);
  }
}

function setParams(values) {
  state.params = { ...values };
  for (const group of CONTROL_GROUPS) {
    for (const c of group.controls) {
      const input = document.getElementById(`ctl-${c.key}`);
      input.value = String(values[c.key]);
      const out = document.getElementById(`out-${c.key}`);
      if (out) out.textContent = Number(values[c.key]).toFixed(c.decimals);
    }
  }
}

// ------------------------------------------------------------------ worker
const pyodideOverride = new URLSearchParams(window.location.search).get("pyodide");
const worker = new Worker(pyodideOverride ? `worker.js?pyodide=${encodeURIComponent(pyodideOverride)}` : "worker.js");
const statusEl = document.getElementById("status");
const progressBar = document.getElementById("progress-bar");
const runButton = document.getElementById("run-button");

function setStatus(text, isError = false) {
  statusEl.textContent = text;
  statusEl.classList.toggle("error", isError);
}

worker.onmessage = (event) => {
  const msg = event.data;
  if (msg.type === "status") {
    setStatus(msg.text);
  } else if (msg.type === "ready") {
    state.ready = true;
    setStatus(`Python ready (Pyodide ${msg.pyodideVersion}).`);
    runSimulation();
  } else if (msg.type === "progress") {
    progressBar.style.width = `${(100 * msg.fraction).toFixed(1)}%`;
    setStatus(`Running… ${(100 * msg.fraction).toFixed(0)}%`);
  } else if (msg.type === "result") {
    state.busy = false;
    runButton.classList.remove("busy");
    runButton.textContent = "Run simulation";
    state.results = JSON.parse(msg.json);
    const s = state.results.summary;
    setStatus(`Done: ${s.n_steps} time steps in ${s.runtime_seconds.toFixed(1)} s.`);
    renderAll();
  } else if (msg.type === "error") {
    state.busy = false;
    runButton.classList.remove("busy");
    runButton.textContent = "Run simulation";
    setStatus(`Error: ${msg.message}`, true);
  }
};

function runSimulation() {
  if (state.busy) return;
  if (!state.ready) { setStatus("Python is still loading; the run will start when it's ready."); return; }
  state.busy = true;
  runButton.classList.add("busy");
  runButton.textContent = "Running…";
  progressBar.style.width = "0%";
  worker.postMessage({ type: "run", params: { ...state.params } });
}

runButton.addEventListener("click", runSimulation);

// ------------------------------------------------------------------- plots
const PLOT_CONFIG = { responsive: true, displaylogo: false };
const FONT = { family: "system-ui, -apple-system, Segoe UI, Roboto, sans-serif", size: 12 };
const COLORS = { temp: "#d62728", rho: "#444444", water: "#1f77b4", ice: "#7f7f7f", pond: "#1f77b4",
                 refrozen: "#9467bd", runoff: "#ff7f0e", melt: "#2ca02c", air: "#d62728" };
const MONTHS = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"];

function stepProfile(zEdges, values) {
  // Layers are uniform inside: draw each as a vertical segment.
  const x = [], y = [];
  for (let i = 0; i < values.length; i++) {
    x.push(values[i], values[i]);
    y.push(zEdges[i], zEdges[i + 1]);
  }
  return { x, y };
}

function timeLabel(tYears) {
  const year = Math.floor(tYears + 1e-9);
  const day = (tYears - year) * 365;
  const month = MONTHS[Math.min(11, Math.floor(day / (365 / 12)))];
  return `Year ${year + 1}, ${month} (t = ${tYears.toFixed(2)} yr)`;
}

function renderSummary() {
  const s = state.results.summary;
  const pct = (x) => `${(100 * x).toFixed(0)}%`;
  const tiles = [
    ["Total surface melt", `${(s.total_melt / 1000).toFixed(2)} m w.e.`, "over the whole run"],
    ["Refrozen in the firn", pct(s.refrozen_fraction), "of melt, incl. frozen ponds"],
    ["Runoff", pct(s.runoff_fraction), "of melt, lost from the column"],
    ["Deepest pond", `${s.pond_max.toFixed(2)} m`, s.first_pond_year === null ? "no pond formed"
      : `first pond in year ${Math.floor(s.first_pond_year) + 1}`],
    ["Ice in upper 10 m", `${s.ice_upper_final.toFixed(2)} m`, "at the end; ≥ impermeable density"],
    ["Firn air content", `${s.firn_air_content_initial.toFixed(1)} → ${s.firn_air_content_final.toFixed(1)} m`,
      "pore space available for water"],
    ["Water mass balance", s.total_melt > 0 ? `${s.mass_balance_relative_error.toExponential(1)}` : "no melt",
      "relative error; should be ~1e-15"],
  ];
  const box = document.getElementById("summary");
  box.innerHTML = "";
  for (const [label, value, note] of tiles) {
    const div = document.createElement("div");
    div.className = "tile";
    for (const [cls, text] of [["label", label], ["value", value], ["note", note]]) {
      const el = document.createElement("div");
      el.className = cls;
      el.textContent = text;
      div.appendChild(el);
    }
    box.appendChild(div);
  }
}

function profileRanges() {
  const snaps = state.results.snapshots;
  let zMax = 0, tMin = 0, wMax = 1;
  for (const s of snaps) {
    zMax = Math.max(zMax, s.z_edges[s.z_edges.length - 1]);
    tMin = Math.min(tMin, ...s.T);
    wMax = Math.max(wMax, ...s.water_percent);
  }
  return { zMax, tMin: Math.floor(tMin - 1), wMax: Math.ceil(wMax * 1.1) };
}

function renderProfiles(index) {
  const r = state.results;
  const snap = r.snapshots[index];
  const init = r.initial;
  const ranges = state.ranges;
  const yaxis = { title: "Depth (m)", range: [ranges.zMax, 0], zeroline: false };
  const margin = { l: 55, r: 10, t: 30, b: 45 };
  const panels = [
    { id: "plot-T", key: "T", color: COLORS.temp, title: "Temperature", xtitle: "Temperature (°C)",
      xrange: [ranges.tMin, 1] },
    { id: "plot-rho", key: "rho", color: COLORS.rho, title: "Density", xtitle: "Density (kg/m³)",
      xrange: [300, 940] },
    { id: "plot-water", key: "water_percent", color: COLORS.water, title: "Liquid water",
      xtitle: "Volumetric water content (%)", xrange: [0, ranges.wMax] },
  ];
  for (const panel of panels) {
    const now = stepProfile(snap.z_edges, snap[panel.key]);
    const start = stepProfile(init.z_edges, init[panel.key]);
    const traces = [
      { ...start, mode: "lines", name: "Initial", line: { color: "#aaaaaa", dash: "dash", width: 1.5 },
        hovertemplate: "%{x:.2f} at %{y:.2f} m<extra>initial</extra>" },
      { ...now, mode: "lines", name: "Selected time", line: { color: panel.color, width: 2 },
        hovertemplate: "%{x:.2f} at %{y:.2f} m<extra></extra>" },
    ];
    const shapes = panel.key === "rho" ? [{ type: "line", x0: r.params.rho_impermeable, x1: r.params.rho_impermeable,
      y0: 0, y1: ranges.zMax, line: { color: COLORS.water, dash: "dot", width: 1.5 } }] : [];
    Plotly.react(panel.id, traces, {
      title: { text: panel.title, font: { size: 14 } }, font: FONT, margin, showlegend: false,
      xaxis: { title: panel.xtitle, range: panel.xrange }, yaxis, shapes,
    }, { ...PLOT_CONFIG, displayModeBar: false });
  }
  document.getElementById("time-label").textContent = timeLabel(snap.t_years);
  const marker = (t) => [{ type: "line", xref: "x", yref: "paper", x0: t, x1: t, y0: 0, y1: 1,
                            line: { color: "#888", dash: "dot", width: 1 } }];
  for (const id of ["plot-water-series", "plot-firn-series", "plot-forcing"]) {
    Plotly.relayout(id, { shapes: marker(snap.t_years) });
  }
}

function renderSeries() {
  const r = state.results;
  const s = r.series;
  const layoutBase = (ytitle, y2title) => {
    const layout = {
      font: FONT, margin: { l: 60, r: 60, t: 10, b: 45 }, hovermode: "x unified",
      legend: { orientation: "h", y: 1.12 },
      xaxis: { title: "Time (years)" },
      yaxis: { title: ytitle, rangemode: "tozero" },
    };
    if (y2title) layout.yaxis2 = { title: y2title, overlaying: "y", side: "right", rangemode: "tozero" };
    return layout;
  };

  Plotly.react("plot-water-series", [
    { x: s.t_years, y: s.pond_depth, name: "Surface pond depth (m)", line: { color: COLORS.pond, width: 1.5 },
      fill: "tozeroy", fillcolor: "rgba(31,119,180,0.2)" },
    { x: s.t_years, y: s.perched_water, name: "Perched water in firn (mm w.e.)", yaxis: "y2",
      line: { color: COLORS.refrozen, width: 1.5 } },
  ], layoutBase("Pond depth (m)", "Perched water (mm w.e.)"), PLOT_CONFIG);

  Plotly.react("plot-firn-series", [
    { x: s.t_years, y: s.firn_air_content, name: "Firn air content (m)", line: { color: COLORS.rho, width: 1.5 } },
    { x: s.t_years, y: s.ice_upper, name: "Ice in upper 10 m (m)", yaxis: "y2",
      line: { color: COLORS.water, width: 1.5 } },
  ], layoutBase("Firn air content (m)", "Ice thickness (m)"), PLOT_CONFIG);

  const a = r.annual;
  Plotly.react("plot-annual", [
    { x: a.year, y: a.melt, type: "bar", name: "Melt", marker: { color: COLORS.melt } },
    { x: a.year, y: a.refrozen, type: "bar", name: "Refrozen", marker: { color: COLORS.refrozen } },
    { x: a.year, y: a.runoff, type: "bar", name: "Runoff", marker: { color: COLORS.runoff } },
  ], { font: FONT, margin: { l: 60, r: 20, t: 10, b: 45 }, barmode: "group",
       legend: { orientation: "h", y: 1.12 },
       xaxis: { title: "Model year" }, yaxis: { title: "Water (mm w.e. per year)" } }, PLOT_CONFIG);

  Plotly.react("plot-forcing", [
    { x: s.t_years, y: s.T_air, name: "Air temperature (°C)", line: { color: COLORS.air, width: 1 } },
    { x: s.t_years, y: s.T_surface, name: "Surface temperature (°C)", line: { color: "#333", width: 1 } },
  ], { ...layoutBase("Temperature (°C)"), yaxis: { title: "Temperature (°C)" } }, PLOT_CONFIG);
}

function renderAll() {
  const r = state.results;
  state.ranges = profileRanges();
  renderSummary();
  renderSeries();
  const slider = document.getElementById("time-slider");
  slider.max = String(r.snapshots.length - 1);
  // Start at the warmest snapshot of the final year, when water is most visible.
  slider.value = String(defaultSnapshot(r));
  renderProfiles(Number(slider.value));
}

function defaultSnapshot(r) {
  const snaps = r.snapshots;
  const lastYear = Math.floor(snaps[snaps.length - 1].t_years - 1e-9);
  let best = snaps.length - 1, bestWater = -1;
  snaps.forEach((s, i) => {
    if (Math.floor(s.t_years - 1e-9) !== lastYear) return;
    const water = s.water_percent.reduce((a, b) => a + b, 0);
    if (water > bestWater) { bestWater = water; best = i; }
  });
  return best;
}

let pending = false;
document.getElementById("time-slider").addEventListener("input", (e) => {
  if (!state.results || pending) return;
  pending = true;
  requestAnimationFrame(() => {
    pending = false;
    renderProfiles(Number(document.getElementById("time-slider").value));
  });
});

// Link to the GitHub repository when served from GitHub Pages.
(function setRepoLink() {
  const host = window.location.hostname;
  if (host.endsWith("github.io")) {
    const user = host.split(".")[0];
    const repo = window.location.pathname.split("/").filter(Boolean)[0] || "";
    document.getElementById("repo-link").href = `https://github.com/${user}/${repo}`;
  }
})();

buildControls();
setParams(BASE);
