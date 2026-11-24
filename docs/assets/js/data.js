// Loads the single results file that every number on the site comes from.
// docs/assets/results.json is written by `python scripts/run_all.py`.
// If it does not exist yet (or the page is opened with ?mock) we fall back to
// docs/assets/mock_results.json, and the page shows a red MOCK banner.
let cache = null;

export function loadResults() {
  if (!cache) cache = fetchResults();
  return cache;
}

async function fetchResults() {
  const wantMock = new URLSearchParams(location.search).has("mock");
  if (!wantMock) {
    try {
      const r = await fetch("assets/results.json", { cache: "no-cache" });
      if (r.ok) return await r.json();
    } catch (e) {
      /* fall through to the mock */
    }
  }
  const r = await fetch("assets/mock_results.json", { cache: "no-cache" });
  if (!r.ok) throw new Error(`could not load results (HTTP ${r.status})`);
  const d = await r.json();
  d.meta = Object.assign({}, d.meta, { mock: true });
  return d;
}

// ---------- small helpers shared by viewer.js and report.js ----------
export const fmt = (v, d = 1) =>
  v === undefined || v === null || Number.isNaN(+v) ? "—" : Number(v).toFixed(d);

// Parameter values may be plain numbers or {value, unit}.
export function param(data, name) {
  const p = data.params && data.params[name];
  if (p === undefined || p === null) return undefined;
  return typeof p === "object" ? p.value : p;
}

export function get(obj, path) {
  return path.split(".").reduce((o, k) => (o === undefined || o === null ? undefined : o[k]), obj);
}

// Constraint rows of one design, plus the thickness bounds checked from d.
export function constraintRows(data, design) {
  const rows = (design.constraints || []).map((c) => ({
    name: c.label || prettyName(c.name),
    key: c.name,
    value: c.value_C ?? c.value,
    trueMax: c.true_max_C,
    limit: c.limit_C ?? c.limit,
    unit: c.value_C !== undefined ? "°C" : c.unit || "",
    sense: c.sense || "<=",
    pass: c.pass,
  }));
  const tvMin = param(data, "tv_min"), tcMin = param(data, "tc_min");
  const dd = design.d || {};
  if (tvMin !== undefined && dd.vane_thickness_mm)
    rows.push(boundRow("Minimum vane thickness", Math.min(...dd.vane_thickness_mm), tvMin * 1e3));
  if (tcMin !== undefined && dd.cheek_thickness_mm)
    rows.push(boundRow("Minimum cheek thickness", Math.min(...dd.cheek_thickness_mm), tcMin * 1e3));
  return rows;
}

function boundRow(name, v, lim) {
  return { name, key: name, value: v, limit: lim, unit: "mm", sense: ">=", pass: v >= lim - 1e-9, bound: true };
}

export function prettyName(n) {
  const map = {
    peak_face_T: "Peak friction-face temperature",
    peak_hub_T: "Peak hub / bell-mount temperature",
  };
  return map[n] || String(n).replace(/_/g, " ");
}

export const DESIGN_LABEL = {
  initial: "Initial",
  optimized: "Optimized",
  optimized_linear: "Optimized without nonlinearity",
};
