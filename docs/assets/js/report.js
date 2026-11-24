// Fills every number, table and plot in the report from results.json, so no
// number on the page is typed by hand. Markup used in index.md:
//   <span class="rj" data-k="designs.optimized.mass_kg" data-d="2"></span>   value at a JSON path
//   <span class="rj" data-calc="mass_saving_pct" data-d="1"></span>          derived value (CALC below)
//   <div class="rj-table" data-table="results"></div>                        table (TABLES below)
//   <div class="rj-plot" data-plot="taylor"></div>                           SVG plot (PLOTS below)
// Anything missing from results.json shows as "pending" instead of breaking the page.
import { loadResults, fmt, get, param, constraintRows, prettyName, DESIGN_LABEL } from "./data.js";

const PENDING = '<span class="rj-pending" title="not in results.json yet">pending</span>';

loadResults().then(fill).catch((e) => {
  console.error(e);
  document.querySelectorAll(".rj, .rj-table, .rj-plot").forEach((n) => (n.innerHTML = PENDING));
});

const num = (v, d) => (v === undefined || v === null || Number.isNaN(+v) ? PENDING : fmt(v, d));
const sci = (v, d = 1) => (v === undefined || v === null ? PENDING : Number(v) === 0 ? "0" : Number(v).toExponential(d));
const passMark = (p) => (p === undefined ? "" : p ? '<span class="rv-pass">✓ pass</span>' : '<span class="rv-fail">✗ fail</span>');

const CALC = {
  mass_saving_pct: (r) => pct(r.designs?.initial?.mass_kg, r.designs?.optimized?.mass_kg),
  mass_saving_kg: (r) => diff(r.designs?.initial?.mass_kg, r.designs?.optimized?.mass_kg),
  face_peak_initial: (r) => con(r, "initial", "peak_face_T")?.true_max_C ?? con(r, "initial", "peak_face_T")?.value_C,
  face_peak_optimized: (r) => con(r, "optimized", "peak_face_T")?.true_max_C ?? con(r, "optimized", "peak_face_T")?.value_C,
  hub_peak_optimized: (r) => con(r, "optimized", "peak_hub_T")?.true_max_C ?? con(r, "optimized", "peak_hub_T")?.value_C,
  face_limit: (r) => con(r, "initial", "peak_face_T")?.limit_C,
  hub_limit: (r) => con(r, "initial", "peak_hub_T")?.limit_C,
  adjoint_speedup: (r) => {
    const c = r.sensitivities?.cost;
    return c && c.fd_s && c.adjoint_s ? c.fd_s / c.adjoint_s : undefined;
  },
  n_nodes: (r) => r.designs?.initial?.mesh?.nodes?.length,
  n_elems: (r) => r.designs?.initial?.mesh?.hexes?.length,
  n_frames: (r) => r.designs?.initial?.temperature?.frames?.length,
  linear_mass_delta_pct: (r) => pct(r.designs?.optimized?.mass_kg, r.designs?.optimized_linear?.mass_kg),
};
function pct(a, b) { return a && b !== undefined ? (100 * (a - b)) / a : undefined; }
function diff(a, b) { return a !== undefined && b !== undefined ? a - b : undefined; }
function con(r, k, name) { return (r.designs?.[k]?.constraints || []).find((c) => c.name === name); }

// Accept the alternative key names / shapes the optimizer thread writes.
function normalize(r) {
  const S = r.sensitivities || {};
  if (!S.taylor && S.taylor_plot) S.taylor = S.taylor_plot;
  if (S.cost_view) S.cost = S.cost_view; // compact view: {n_vars, forward_s, adjoint_s, fd_s, speedup}
  for (const st of r.optimization?.starts || []) {
    if (Array.isArray(st.history)) {
      const h = {};
      for (const row of st.history) for (const [k, v] of Object.entries(row)) (h[k] = h[k] || []).push(v);
      st.history = h;
    }
    if (st.history && !st.history.iter && st.history.mass_kg) st.history.iter = st.history.mass_kg.map((_, i) => i);
  }
  return r;
}
const startLabel = (name) => (name.includes(":") ? "start " + name.split(":")[0] : name);

// Any object or list as tables: lists of rows become tables, scalars a key/value table.
function anyTable(v) {
  if (v === undefined || v === null) return null;
  if (Array.isArray(v)) return genericTable(v.filter((x) => x && typeof x === "object")) || cell(v);
  if (typeof v !== "object") return cell(v);
  const scalars = [], parts = [];
  for (const [k, x] of Object.entries(v)) {
    if (Array.isArray(x) && x.length && typeof x[0] === "object") parts.push(`<h4>${k.replace(/_/g, " ")}</h4>` + genericTable(x));
    else if (x && typeof x === "object" && !Array.isArray(x)) parts.push(`<h4>${k.replace(/_/g, " ")}</h4>` + anyTable(x));
    else if (!(Array.isArray(x) && x.length > 6)) scalars.push([k.replace(/_/g, " "), typeof x === "string" ? x : cell(x, k)]);
  }
  return (scalars.length ? table(["Quantity", "Value"], scalars) : "") + parts.join("");
}

function fill(r) {
  normalize(r);
  if (r.meta?.mock) {
    const b = document.getElementById("mock-banner");
    if (b) { b.hidden = false; b.innerHTML = "<b>MOCK DATA.</b> Every number on this page is a placeholder from <code>scripts/make_mock_results.py</code>. Run <code>python scripts/run_all.py</code> to replace them."; }
  }
  document.querySelectorAll(".rj").forEach((n) => {
    const d = n.dataset.d !== undefined ? +n.dataset.d : 2;
    let v;
    if (n.dataset.param) v = param(r, n.dataset.param);
    else if (n.dataset.calc) v = CALC[n.dataset.calc]?.(r);
    else if (n.dataset.k) v = get(r, n.dataset.k);
    if (v !== undefined && n.dataset.scale) v = v * +n.dataset.scale;
    if (v !== undefined && n.dataset.offset) v = v + +n.dataset.offset;
    if (typeof v === "string") n.textContent = v;
    else n.innerHTML = n.dataset.fmt === "sci" ? sci(v, d) : num(v, d);
  });
  document.querySelectorAll(".rj-table").forEach((n) => {
    const f = TABLES[n.dataset.table];
    try { n.innerHTML = (f && f(r)) || PENDING; } catch (e) { console.error(e); n.innerHTML = PENDING; }
  });
  document.querySelectorAll(".rj-plot").forEach((n) => {
    const f = PLOTS[n.dataset.plot];
    try { n.innerHTML = (f && f(r)) || PENDING; } catch (e) { console.error(e); n.innerHTML = PENDING; }
  });
}

// ---------------------------------------------------------------- tables
const table = (head, rows, cls = "") =>
  `<div class="rj-scroll"><table class="${cls}"><thead><tr>${head.map((h) => `<th>${h}</th>`).join("")}</tr></thead>` +
  `<tbody>${rows.map((row) => `<tr>${row.map((c) => `<td>${c}</td>`).join("")}</tr>`).join("")}</tbody></table></div>`;

function genericTable(list) {
  if (!Array.isArray(list) || !list.length) return null;
  const keys = [...new Set(list.flatMap((o) => Object.keys(o)))];
  return table(keys.map((k) => k.replace(/_/g, " ")), list.map((o) => keys.map((k) => cell(o[k], k))));
}
function cell(v, k = "") {
  if (v === undefined || v === null) return "";
  if (typeof v === "boolean") return k === "pass" || k === "converged" ? passMark(v) : String(v);
  if (typeof v === "number") return Math.abs(v) !== 0 && (Math.abs(v) < 1e-2 || Math.abs(v) >= 1e5) ? sci(v, 2) : fmt(v, Math.abs(v) >= 100 ? 1 : 3);
  if (Array.isArray(v)) return v.length > 6 ? `[${v.length} values]` : v.map((x) => cell(x)).join(", ");
  if (typeof v === "object") return JSON.stringify(v);
  return String(v);
}

function designTable(r, want) {
  const keys = want.filter((k) => r.designs?.[k]);
    if (!keys.length) return null;
  const base = constraintRows(r, r.designs[keys[0]]);
  const head = ["Quantity", ...keys.map((k) => DESIGN_LABEL[k]), "Limit"];
  const rows = [["Rotor mass (objective)", ...keys.map((k) => `${num(r.designs[k].mass_kg, 2)} kg`), "minimize"]];
  base.forEach((c0, i) => {
    rows.push([c0.name, ...keys.map((k) => {
      const c = constraintRows(r, r.designs[k])[i];
      const tm = c.trueMax !== undefined ? ` <small>(peak ${fmt(c.trueMax, 1)})</small>` : "";
      return `${num(c.value, c.unit === "mm" ? 2 : 1)} ${c.unit}${tm} ${passMark(c.pass)}`;
    }), `${c0.sense === ">=" ? "≥" : "≤"} ${fmt(c0.limit, c0.unit === "mm" ? 1 : 0)} ${c0.unit}`]);
  });
  return table(head, rows, "rj-results");
}

const TABLES = {
  results(r) { return designTable(r, ["initial", "optimized", "optimized_linear"]); },
  nonlinear(r) { return r.designs?.optimized_linear ? designTable(r, ["optimized", "optimized_linear"]) : null; },
  design_vars(r) {
    const i = r.designs?.initial?.d, o = r.designs?.optimized?.d;
    if (!i) return null;
    const rows = i.station_radius_mm.map((rad, k) => [k + 1, fmt(rad, 1), fmt(i.vane_thickness_mm[k], 2), o ? fmt(o.vane_thickness_mm[k], 2) : PENDING,
      fmt(i.cheek_thickness_mm[k], 2), o ? fmt(o.cheek_thickness_mm[k], 2) : PENDING]);
    return table(["Station", "Radius (mm)", "Vane t<sub>v</sub> initial (mm)", "Vane t<sub>v</sub> optimized (mm)", "Cheek t<sub>c</sub> initial (mm)", "Cheek t<sub>c</sub> optimized (mm)"], rows);
  },
  verification(r) {
    const v = r.verification;
    if (!Array.isArray(v) || !v.length) return null;
    return table(["ID", "Test", "Ours", "Exact", "Rel. error", "Status", "Note"],
      v.map((t) => [t.id, t.name ?? t.test, cell(t.ours ?? t.computed), cell(t.exact), sci(t.rel_error ?? t.error, 1), passMark(t.pass), t.note || ""]));
  },
  gradient_check(r) {
    const g = r.sensitivities?.gradient_check;
    if (!Array.isArray(g) || !g.length) return null;
    if (!("adjoint" in g[0])) return genericTable(g);
    const hasModel = g.some((t) => t.model !== undefined);
    return table(["Function", "Variable", ...(hasModel ? ["Model"] : []), "Adjoint", "Finite difference", "Rel. error"],
      g.map((t) => [t.function, t.variable, ...(hasModel ? [t.model ?? ""] : []), Number(t.adjoint).toPrecision(8), Number(t.fd).toPrecision(8), sci(t.rel_error, 1)]));
  },
  cost(r) {
    const c = r.sensitivities?.cost;
    if (!c) return null;
    return table(["Design variables", "One forward solve (s)", "Adjoint gradients (s)", "Finite-difference gradients (s)", "Speed-up"],
      [[c.n_vars, fmt(c.forward_s, 2), fmt(c.adjoint_s, 2), fmt(c.fd_s, 1), c.fd_s && c.adjoint_s ? `${fmt(c.fd_s / c.adjoint_s, 1)}×` : ""]]);
  },
  starts(r) {
    const s = r.optimization?.starts;
    if (!Array.isArray(s) || !s.length) return null;
    return table(["Start", "Iterations", "Final mass (kg)", "Converged"],
      s.map((t) => [t.name, t.iterations ?? (t.history?.iter ? t.history.iter.length - 1 : ""), fmt(t.final_mass_kg, 3),
        t.converged === undefined ? "" : t.converged ? '<span class="rv-pass">✓ yes</span>' : `✗ no${t.message ? ` (${t.message})` : ""}`]));
  },
  multipliers(r) {
    const a = r.optimization?.active;
    if (!Array.isArray(a) || !a.length) return null;
    const rows = a.map((t) => [prettyName(t.constraint), `${cell(t.multiplier)} ${t.unit || ""}`, t.meaning || ""]);
    const b = r.optimization?.active_bounds;
    if (b && (b.lower?.length || b.upper?.length)) {
      const lst = (v) => v.join(", ");
      if (b.lower?.length) rows.push(["Lower thickness bounds", `${b.lower.length} active`, `at the minimum: ${lst(b.lower)}`]);
      if (b.upper?.length) rows.push(["Upper thickness bounds", `${b.upper.length} active`, `at the maximum: ${lst(b.upper)}`]);
    }
    return table(["Active constraint", "Multiplier", "What it means"], rows);
  },
  n2(r) { return anyTable(r.optimization?.n2_comparison); },
  robust(r) { return anyTable(r.optimization?.robust); },
  params(r) {
    if (!r.params) return null;
    const rows = Object.entries(r.params).map(([k, v]) => [`<code>${k}</code>`, typeof v === "object" && v ? cell(v.value) : cell(v), typeof v === "object" && v ? v.unit || "" : ""]);
    return table(["Parameter", "Value", "Unit"], rows);
  },
  meta(r) {
    const m = r.meta || {};
    return table(["Generated", "Git commit", "Command"], [[m.generated || "", `<code>${m.git_commit || ""}</code>`, `<code>${m.command || ""}</code>`]]);
  },
};

// ---------------------------------------------------------------- plots
const COLORS = ["#c2410c", "#2563eb", "#6b7280", "#15803d", "#7c3aed"];

function plot({ series, xlog = false, ylog = false, xlabel = "", ylabel = "", hlines = [], width = 640, height = 300 }) {
  series = series.filter((s) => s && s.x && s.y && s.x.length);
  if (!series.length) return null;
  const tx = (v) => (xlog ? Math.log10(v) : v), ty = (v) => (ylog ? Math.log10(Math.max(v, 1e-300)) : v);
  let x0 = Infinity, x1 = -Infinity, y0 = Infinity, y1 = -Infinity;
  for (const s of series) s.x.forEach((x, i) => {
    const X = tx(x), Y = ty(s.y[i]);
    if (!Number.isFinite(X) || !Number.isFinite(Y)) return;
    x0 = Math.min(x0, X); x1 = Math.max(x1, X); y0 = Math.min(y0, Y); y1 = Math.max(y1, Y);
  });
  for (const h of hlines) { y0 = Math.min(y0, ty(h.y)); y1 = Math.max(y1, ty(h.y)); }
  if (x1 === x0) x1 = x0 + 1;
  if (y1 === y0) y1 = y0 + 1;
  const pad = (y1 - y0) * 0.06; y0 -= pad; y1 += pad;
  const L = 62, R = 14, T = 14, B = 44, W = width, H = height, pw = W - L - R, ph = H - T - B;
  const X = (v) => L + ((tx(v) - x0) / (x1 - x0)) * pw, Y = (v) => T + ph * (1 - (ty(v) - y0) / (y1 - y0));
  const ticks = (a, b, log) => {
    if (log) { const out = []; for (let e = Math.ceil(a); e <= Math.floor(b); e++) out.push(10 ** e); const st = Math.ceil(out.length / 7); return out.filter((_, i) => i % st === 0); }
    const span = b - a, mag = 10 ** Math.floor(Math.log10(span / 5)), stp = [1, 2, 5, 10].map((m) => m * mag).find((s) => span / s <= 6);
    const out = []; for (let v = Math.ceil(a / stp) * stp; v <= b + 1e-12; v += stp) out.push(+v.toPrecision(10)); return out;
  };
  const lab = (v, log) => (log ? `10<tspan dy="-6" font-size="9">${Math.round(Math.log10(v))}</tspan>` : String(v));
  let svg = `<svg viewBox="0 0 ${W} ${H}" class="rj-svg" role="img" aria-label="${ylabel} versus ${xlabel}">`;
  for (const v of ticks(y0, y1, ylog)) svg += `<line x1="${L}" x2="${L + pw}" y1="${Y(v)}" y2="${Y(v)}" class="grid"/><text x="${L - 6}" y="${Y(v) + 4}" text-anchor="end">${lab(v, ylog)}</text>`;
  for (const v of ticks(x0, x1, xlog)) svg += `<line y1="${T}" y2="${T + ph}" x1="${X(v)}" x2="${X(v)}" class="grid"/><text x="${X(v)}" y="${T + ph + 16}" text-anchor="middle">${lab(v, xlog)}</text>`;
  svg += `<rect x="${L}" y="${T}" width="${pw}" height="${ph}" class="frame"/>`;
  for (const h of hlines) svg += `<line x1="${L}" x2="${L + pw}" y1="${Y(h.y)}" y2="${Y(h.y)}" class="ref"/><text x="${L + pw - 4}" y="${Y(h.y) - 4}" text-anchor="end">${h.label}</text>`;
  series.forEach((s, k) => {
    const c = s.color || COLORS[k % COLORS.length];
    const pts = s.x.map((x, i) => [X(x), Y(s.y[i])]).filter(([a, b]) => Number.isFinite(a) && Number.isFinite(b));
    svg += `<polyline points="${pts.map((p) => p.join(",")).join(" ")}" fill="none" stroke="${c}" stroke-width="2" ${s.dash ? 'stroke-dasharray="5 4"' : ""}/>`;
    if (s.markers !== false) for (const [a, b] of pts) svg += `<circle cx="${a}" cy="${b}" r="2.6" fill="${c}"/>`;
  });
  let lx = L + 8, ly = T + 12;
  series.forEach((s, k) => {
    if (!s.label) return;
    const c = s.color || COLORS[k % COLORS.length];
    const wlab = 30 + 6.2 * s.label.length;
    if (lx + wlab > L + pw && lx > L + 8) { lx = L + 8; ly += 15; } // wrap the legend
    svg += `<rect x="${lx - 3}" y="${ly - 9}" width="${wlab - 4}" height="14" class="legbg"/>`;
    svg += `<line x1="${lx}" x2="${lx + 16}" y1="${ly}" y2="${ly}" stroke="${c}" stroke-width="2" ${s.dash ? 'stroke-dasharray="4 3"' : ""}/><text x="${lx + 20}" y="${ly + 4}">${s.label}</text>`;
    lx += wlab;
  });
  svg += `<text x="${L + pw / 2}" y="${H - 6}" text-anchor="middle" class="axis">${xlabel}</text>`;
  svg += `<text transform="translate(14 ${T + ph / 2}) rotate(-90)" text-anchor="middle" class="axis">${ylabel}</text></svg>`;
  return svg;
}

// reference line of a given slope through the first point of a series (log-log)
function slopeRef(x, y, p, label) {
  if (!x || !x.length) return null;
  return { x: [x[0], x[x.length - 1]], y: [y[0], y[0] * (x[x.length - 1] / x[0]) ** p], label, dash: true, markers: false, color: "#9ca3af" };
}

const PLOTS = {
  taylor(r) {
    const t = r.sensitivities?.taylor;
    if (!t || !t.eps) return null;
    const keys = Object.keys(t).filter((k) => k !== "eps" && Array.isArray(t[k]));
    const names = { err_fd: "|f(x+εp) − f(x)|  (should fall as ε)", err_adjoint: "|f(x+εp) − f(x) − ε ∇fᵀp|  (should fall as ε²)" };
    return plot({ xlog: true, ylog: true, xlabel: "step size ε", ylabel: "Taylor remainder",
      series: [...keys.map((k) => ({ x: t.eps, y: t[k], label: names[k] || prettyName(k) })), slopeRef(t.eps, t[keys[keys.length - 1]], 2, "slope 2")].filter(Boolean) });
  },
  history(r) {
    const s = r.optimization?.starts;
    if (!Array.isArray(s) || !s.length) return null;
    return plot({ xlabel: "SLSQP iteration", ylabel: "rotor mass (kg)",
      series: s.map((st) => ({ x: st.history.iter, y: st.history.mass_kg, label: startLabel(st.name) })) });
  },
  history_g(r) {
    const s = r.optimization?.starts;
    if (!Array.isArray(s) || !s.length || !s[0].history.g_face_C) return null;
    const lim = con(r, "initial", "peak_face_T")?.limit_C;
    return plot({ xlabel: "SLSQP iteration", ylabel: "KS face temperature (°C)", hlines: lim ? [{ y: lim, label: `limit ${lim} °C` }] : [],
      series: s.map((st) => ({ x: st.history.iter, y: st.history.g_face_C, label: startLabel(st.name) })) });
  },
  mms(r) {
    const c = r.convergence;
    if (!c) return null;
    const NAMES = { mms_space: "Manufactured solution, space", mms_time: "Manufactured solution, time",
      kirchhoff: "Nonlinear k(T) slab (Kirchhoff), space", radiation_time: "Radiation-cooled lump, time" };
    const out = [];
    let k = 0;
    for (const [name, v] of Object.entries(c)) {
      if (!v || typeof v !== "object") continue;
      const xk = ["h", "h_rel", "dt"].find((q) => Array.isArray(v[q]));
      const yk = ["err_L2", "err"].find((q) => Array.isArray(v[q]));
      if (!xk || !yk) continue;
      const p = xk === "dt" ? 1 : 2;
      const xl = xk === "dt" ? "time step Δt (s)" : xk === "h_rel" ? "relative element size h" : "element size h (m)";
      out.push(`<figure>${plot({ xlog: true, ylog: true, xlabel: xl, ylabel: "error", width: 330, height: 260,
        series: [{ x: v[xk], y: v[yk], label: v.slope !== undefined ? `measured slope ${fmt(v.slope, 2)}` : "error", color: COLORS[k++ % COLORS.length] },
                 slopeRef(v[xk], v[yk], p, `slope ${p}`)] })}<figcaption>${NAMES[name] || name}</figcaption></figure>`);
    }
    return out.length ? `<div class="rj-pair">${out.join("")}</div>` : null;
  },
  newton(r) {
    const n = r.convergence?.newton?.residuals;
    if (!n) return null;
    return plot({ ylog: true, xlabel: "Newton iteration", ylabel: "residual norm", width: 420, height: 240,
      series: [{ x: n.map((_, i) => i), y: n, label: "one braking time step" }] });
  },
  thickness(r) {
    const keys = ["initial", "optimized", "optimized_linear"].filter((k) => r.designs?.[k]?.d);
    if (!keys.length) return null;
    const s = [];
    keys.forEach((k, i) => {
      const d = r.designs[k].d;
      s.push({ x: d.station_radius_mm, y: d.cheek_thickness_mm, label: `${DESIGN_LABEL[k]} cheek`, color: COLORS[i] });
      s.push({ x: d.station_radius_mm, y: d.vane_thickness_mm, label: `${DESIGN_LABEL[k]} vane`, color: COLORS[i], dash: true });
    });
    return plot({ xlabel: "radius (mm)", ylabel: "thickness (mm)", series: s, height: 320 });
  },
  peaks(r) {
    const keys = ["initial", "optimized", "optimized_linear"].filter((k) => r.designs?.[k]?.history);
    if (!keys.length) return null;
    const lim = con(r, "initial", "peak_face_T")?.limit_C;
    return plot({ xlabel: "time (s)", ylabel: "peak friction-face T (°C)", hlines: lim ? [{ y: lim, label: `limit ${lim} °C` }] : [],
      series: keys.map((k, i) => ({ x: r.designs[k].history.times, y: r.designs[k].history.peak_face_C, label: DESIGN_LABEL[k], color: COLORS[i], markers: false })) });
  },
};
