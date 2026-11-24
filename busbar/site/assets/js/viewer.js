// 3D initial-vs-optimized busbar viewer. Everything comes from results.json.
// Handout viewer requirements, and where each is met:
//  1 side by side, linked cameras ............ one shared camera, two OrbitControls
//  2 true 3D geometry ......................... half model mirrored about y = 0; optional, labelled thickness exaggeration
//  3 single-hue temperature scale + colour bar; faces above their limit are hatched (not colour alone)
//  4 field on the geometry with the element's shape functions (bilinear on each hex face)
//  5 animation: heating during the current pulse and the recovery afterwards (play / scrub)
//  6 hover or tap: element, surface, position, temperature now and peak, current density, limit margin
//  7 metrics panel per design: mass and every constraint with pass / fail
//  8 loads results.json, responsive to phone width, renders in about a second
import * as THREE from "three";
import { OrbitControls } from "three/addons/controls/OrbitControls.js";
import { loadResults, fmt, constraintRows, DESIGN_LABEL } from "./data.js";

const ROOT = document.getElementById("rotor-viewer");
const SUB = 3;
const HEX_FACES = [[0, 3, 2, 1], [4, 5, 6, 7], [0, 1, 5, 4], [2, 3, 7, 6], [0, 4, 7, 3], [1, 2, 6, 5]];
const SURF = ["contact_in", "contact_out", "bottom", "top", "side", "ends", "other"];
const SURF_LABEL = {
  contact_in: "inlet joint face (bolted to the module terminal)",
  contact_out: "outlet joint face (bolted to the module terminal)",
  bottom: "underside", top: "top face", side: "edge", ends: "end face", other: "surface",
};
const RAMP = ["#fff5eb", "#fee6ce", "#fdd0a2", "#fdae6b", "#fd8d3c", "#f16913", "#d94801", "#a63603", "#7f2704"];
const RAMP_LIN = RAMP.map((h) => new THREE.Color(h));
const el = (tag, cls, html) => { const e = document.createElement(tag); if (cls) e.className = cls; if (html !== undefined) e.innerHTML = html; return e; };

if (ROOT) {
  loadResults()
    .then((data) => { ROOT.innerHTML = ""; build(data); })
    .catch((err) => { console.error(err); ROOT.innerHTML = `<p class="rv-error">The 3D viewer could not load its data: ${err.message}</p>`; });
}

function rampColor(u, out) {
  u = Math.min(1, Math.max(0, u)) * (RAMP_LIN.length - 1);
  const i = Math.min(RAMP_LIN.length - 2, Math.floor(u));
  const f = u - i, a = RAMP_LIN[i], b = RAMP_LIN[i + 1];
  out[0] = a.r + (b.r - a.r) * f; out[1] = a.g + (b.g - a.g) * f; out[2] = a.b + (b.b - a.b) * f;
}

function limitFor(design) {
  let bar = Infinity, joint = Infinity;
  for (const c of design.constraints || []) {
    if (c.name === "peak_bar_T") bar = c.limit_C;
    if (c.name === "peak_joint_T") joint = c.limit_C;
  }
  return SURF.map((s) => (s.startsWith("contact") ? Math.min(joint, bar) : bar));
}

function build(data) {
  const D = data.designs;
  const rightKeys = ["optimized", "optimized_linear", "optimized_robust"].filter((k) => D[k]);
  const leftKey = D.initial ? "initial" : rightKeys.shift();
  let rightKey = rightKeys[0];

  let tMin = Infinity, tMax = -Infinity;
  for (const k of Object.keys(D)) for (const fr of D[k].temperature.frames) for (const v of fr) { if (v < tMin) tMin = v; if (v > tMax) tMax = v; }
  const step = 5;
  tMin = Math.floor(tMin / step) * step;
  const lims = Object.values(D).flatMap((d) => limitFor(d)).filter(Number.isFinite);
  const cMax = Math.ceil(Math.max(tMax, Math.min(...lims) + step) / step) * step;

  const bar = el("div", "rv-toolbar");
  const playBtn = el("button", "rv-btn", "▶ Play"); playBtn.type = "button";
  playBtn.setAttribute("aria-label", "Play or pause the current-pulse animation");
  const times = D[leftKey].temperature.times;
  const tStart = times[0], tEnd = times[times.length - 1];
  const slider = el("input", "rv-slider");
  Object.assign(slider, { type: "range", min: tStart, max: tEnd, step: 0.1, value: tStart });
  slider.setAttribute("aria-label", "Time");
  const tLabel = el("span", "rv-tlabel");
  const exSel = el("select", "rv-select");
  exSel.setAttribute("aria-label", "Thickness exaggeration");
  for (const [v, t] of [["1", "True thickness"], ["4", "Thickness × 4 (exaggerated)"]])
    exSel.appendChild(Object.assign(document.createElement("option"), { value: v, textContent: t }));
  bar.append(playBtn, slider, tLabel, exSel);
  let cmpSel = null;
  if (rightKeys.length > 1) {
    cmpSel = el("select", "rv-select");
    cmpSel.setAttribute("aria-label", "Design shown on the right");
    for (const k of rightKeys) cmpSel.appendChild(Object.assign(document.createElement("option"), { value: k, textContent: "Right: " + DESIGN_LABEL[k] }));
    bar.appendChild(cmpSel);
  }
  ROOT.appendChild(bar);
  const panels = el("div", "rv-panels"); ROOT.appendChild(panels);
  const legend = el("div", "rv-legend"); ROOT.appendChild(legend);
  drawColorbar(legend, tMin, cMax, D[leftKey]);
  const chartWrap = el("div", "rv-chart"); const chart = el("canvas"); chartWrap.appendChild(chart); ROOT.appendChild(chartWrap);
  const note = el("p", "rv-note"); ROOT.appendChild(note);
  const setNote = (ex) => (note.innerHTML =
    `${ex > 1 ? `<b>Thickness exaggerated ${ex}×</b> (lengths and widths true). ` : "True geometry, no exaggeration. "}` +
    `The solver models half of the bar (mirror plane through its centre line); the viewer mirrors it. The bolted joints ` +
    `are the 25 mm pads on the underside at each end (rotate to see them). Drag to rotate (both views move together), scroll or pinch to zoom, ` +
    `hover or tap a face to inspect it. Hatched faces are above their temperature limit.`);
  setNote(1);
  const tip = el("div", "rv-tip"); tip.hidden = true; document.body.appendChild(tip);

  const camera = new THREE.PerspectiveCamera(35, 1, 1, 5000);
  camera.up.set(0, 0, 1);
  const target = new THREE.Vector3();
  let dirty = true, ex = 1;
  const left = makeView(leftKey);
  const right = rightKey ? makeView(rightKey) : null;
  const views = [left, right].filter(Boolean);

  function makeView(key) {
    const panel = el("div", "rv-panel"), head = el("div", "rv-head"), wrap = el("div", "rv-canvas");
    const renderer = new THREE.WebGLRenderer({ antialias: true, alpha: true });
    renderer.setPixelRatio(Math.min(window.devicePixelRatio, 2));
    wrap.appendChild(renderer.domElement);
    const metrics = el("div", "rv-metrics");
    panel.append(head, wrap, metrics); panels.appendChild(panel);
    const scene = new THREE.Scene();
    scene.add(new THREE.HemisphereLight(0xffffff, 0xaaaaaa, 2.2));
    const sun = new THREE.DirectionalLight(0xffffff, 0.8); sun.position.set(1, -1, 2); scene.add(sun);
    const hl = new THREE.LineLoop(new THREE.BufferGeometry().setFromPoints([0, 1, 2, 3].map(() => new THREE.Vector3())),
      new THREE.LineBasicMaterial({ color: 0x1f6feb, depthTest: false }));
    hl.renderOrder = 2; hl.visible = false; scene.add(hl);
    const controls = new OrbitControls(camera, renderer.domElement);
    controls.target = target; controls.enablePan = false; controls.minDistance = 30; controls.maxDistance = 2000;
    controls.addEventListener("change", () => (dirty = true));
    const v = { panel, head, wrap, renderer, scene, hl, controls, metrics };
    setDesign(v, key);
    return v;
  }

  function setDesign(v, key) {
    if (v.mesh) { v.scene.remove(v.mesh); v.mesh.geometry.dispose(); }
    const d = D[key];
    v.key = key; v.d = d;
    v.geo = buildGeometry(d);
    v.limits = limitFor(d);
    v.mesh = new THREE.Mesh(v.geo.geometry, hatchMaterial());
    v.mesh.position.x = -v.geo.xMid;
    v.mesh.scale.z = ex;
    v.scene.add(v.mesh);
    v.head.innerHTML = `<strong>${DESIGN_LABEL[key] || key}</strong> <span>${fmt(1000 * d.mass_kg, 1)} g</span>`;
    fillMetrics(v.metrics, data, d, key);
  }

  function setView() {
    const Lb = left.geo.xMid * 2;
    target.set(0, 0, 0);
    camera.position.set(0.35 * Lb, -0.9 * Lb, 0.55 * Lb);
    views.forEach((v) => v.controls.update());
  }

  let tNow = tStart, playing = false, lastMs = 0;
  const PLAY_SECONDS = 15;
  function frameAt(ts, t) {
    let i = 0;
    while (i < ts.length - 2 && ts[i + 1] < t) i++;
    const f = Math.min(1, Math.max(0, (t - ts[i]) / (ts[i + 1] - ts[i] || 1)));
    return [i, f];
  }
  const H = D[leftKey].history || {};
  const hT = H.times || times, cur = H.current_A || null;
  const Ibase = cur ? Math.min(...cur) : 0;
  function setTime(t) {
    tNow = Math.min(tEnd, Math.max(tStart, t));
    slider.value = tNow;
    let txt = `t = ${fmt(tNow, 0)} s`;
    if (cur) { const [i] = frameAt(hT, tNow); const I = cur[Math.min(i + 1, cur.length - 1)]; txt += ` · ${fmt(I, 0)} A · ${I > Ibase ? "boost pulse" : "continuous rating"}`; }
    tLabel.textContent = txt;
    for (const v of views) colorize(v);
    drawChart();
    dirty = true;
  }
  function colorize(v) {
    const fr = v.d.temperature.frames, ts = v.d.temperature.times;
    const [i, f] = frameAt(ts, tNow);
    const A = fr[i], B = fr[i + 1] || fr[i];
    const { corner, weight, color, over, vRegion, Tv } = v.geo;
    const c = [0, 0, 0], span = cMax - tMin;
    for (let k = 0, n = Tv.length; k < n; k++) {
      let T = 0;
      for (let j = 0; j < 4; j++) { const id = corner[4 * k + j]; T += weight[4 * k + j] * (A[id] + (B[id] - A[id]) * f); }
      Tv[k] = T;
      rampColor((T - tMin) / span, c);
      color[3 * k] = c[0]; color[3 * k + 1] = c[1]; color[3 * k + 2] = c[2];
      over[k] = T > v.limits[vRegion[k]] ? 1 : 0;
    }
    v.geo.geometry.attributes.color.needsUpdate = true;
    v.geo.geometry.attributes.over.needsUpdate = true;
  }

  const SERIES = { initial: "#6b7280", optimized: "#c2410c", optimized_linear: "#2563eb", optimized_robust: "#15803d" };
  const jointLimit = left.limits[0], barLimit = left.limits[3];
  function drawChart() {
    const W = chartWrap.clientWidth, Hh = 170, dpr = Math.min(window.devicePixelRatio, 2);
    if (!W) return;
    chart.width = W * dpr; chart.height = Hh * dpr; chart.style.height = Hh + "px";
    const g = chart.getContext("2d");
    g.setTransform(dpr, 0, 0, dpr, 0, 0);
    const css = getComputedStyle(ROOT);
    const ink = css.getPropertyValue("--rv-ink").trim() || "#222", muted = css.getPropertyValue("--rv-muted").trim() || "#888";
    const L = 44, R = 44, T = 22, B = 26, pw = W - L - R, ph = Hh - T - B;
    const X = (t) => L + ((t - tStart) / (tEnd - tStart)) * pw;
    const Y = (v) => T + ph * (1 - (v - tMin) / (cMax - tMin));
    const imax = cur ? Math.max(...cur) * 1.1 : 1, Yi = (v) => T + ph * (1 - v / imax);
    g.clearRect(0, 0, W, Hh); g.font = "11px system-ui, sans-serif";
    if (cur) {
      g.fillStyle = "rgba(127,127,127,0.14)";
      let a = null;
      for (let i = 0; i < hT.length; i++) {
        const on = cur[i] > Ibase;
        if (on && a === null) a = i > 0 ? hT[i - 1] : hT[i];
        if ((!on || i === hT.length - 1) && a !== null) { g.fillRect(X(a), T, X(hT[on ? i : i - 1]) - X(a), ph); a = null; }
      }
    }
    g.strokeStyle = muted; g.lineWidth = 1; g.strokeRect(L, T, pw, ph);
    g.fillStyle = muted; g.textAlign = "right";
    for (const v of [tMin, (tMin + cMax) / 2, cMax]) g.fillText(fmt(v, 0), L - 4, Y(v) + 4);
    g.textAlign = "center";
    for (const t of [tStart, (tStart + tEnd) / 2, tEnd]) g.fillText(fmt(t, 0), X(t), T + ph + 13);
    g.fillText("time (s)", L + pw / 2, Hh - 3);
    g.save(); g.translate(11, T + ph / 2); g.rotate(-Math.PI / 2); g.fillText("temperature (°C)", 0, 0); g.restore();
    if (cur) {
      g.textAlign = "left";
      for (const v of [0, Math.max(...cur)]) g.fillText(String(Math.round(v)), L + pw + 4, Yi(v) + 4);
      g.save(); g.translate(W - 5, T + ph / 2); g.rotate(Math.PI / 2); g.textAlign = "center"; g.fillText("current (A)", 0, 0); g.restore();
      g.strokeStyle = muted; g.setLineDash([2, 3]);
      g.beginPath(); hT.forEach((t, i) => (i ? g.lineTo(X(t), Yi(cur[i])) : g.moveTo(X(t), Yi(cur[i])))); g.stroke(); g.setLineDash([]);
    }
    for (const [lim, lab] of [[jointLimit, "joint limit"], [barLimit, "insulation limit"]]) {
      if (!Number.isFinite(lim) || lim > cMax) continue;
      g.strokeStyle = ink; g.setLineDash([5, 3]); g.beginPath(); g.moveTo(L, Y(lim)); g.lineTo(L + pw, Y(lim)); g.stroke(); g.setLineDash([]);
      g.fillStyle = ink; g.textAlign = "right"; g.fillText(`${lab} ${fmt(lim, 0)} °C`, L + pw - 4, Y(lim) - 4);
    }
    let lx = L + 4; g.textAlign = "left";
    for (const k of [leftKey, rightKey].filter(Boolean)) {
      const h = D[k].history; if (!h) continue;
      g.strokeStyle = SERIES[k] || ink; g.lineWidth = 2;
      g.beginPath(); h.times.forEach((t, i) => (i ? g.lineTo(X(t), Y(h.peak_bar_C[i])) : g.moveTo(X(t), Y(h.peak_bar_C[i])))); g.stroke();
      g.setLineDash([4, 3]);
      g.beginPath(); h.times.forEach((t, i) => (i ? g.lineTo(X(t), Y(h.peak_joint_C[i])) : g.moveTo(X(t), Y(h.peak_joint_C[i])))); g.stroke();
      g.setLineDash([]);
      g.fillStyle = SERIES[k] || ink; g.fillRect(lx, 8, 14, 3);
      g.fillStyle = ink; g.fillText(DESIGN_LABEL[k] || k, lx + 18, 13);
      lx += 30 + g.measureText(DESIGN_LABEL[k] || k).width;
    }
    if (lx < W - 190) { g.fillStyle = muted; g.fillText("— bar peak  ╌ joint  ┈ current", lx, 13); }
    g.strokeStyle = "#1f6feb"; g.lineWidth = 1.5;
    g.beginPath(); g.moveTo(X(tNow), T); g.lineTo(X(tNow), T + ph); g.stroke();
    g.lineWidth = 1;
    chart._x2t = (x) => tStart + ((x - L) / pw) * (tEnd - tStart);
  }
  const scrub = (e) => { const r = chart.getBoundingClientRect(); stop(); setTime(chart._x2t(e.clientX - r.left)); };
  chart.addEventListener("pointerdown", (e) => { scrub(e); chart.setPointerCapture(e.pointerId); });
  chart.addEventListener("pointermove", (e) => { if (e.buttons) scrub(e); });

  const ray = new THREE.Raycaster(), ndc = new THREE.Vector2(), tri = new THREE.Triangle(), bary = new THREE.Vector3();
  for (const v of views) {
    const dom = v.renderer.domElement;
    dom.addEventListener("pointermove", (e) => inspect(v, e));
    dom.addEventListener("pointerdown", (e) => inspect(v, e));
    dom.addEventListener("pointerleave", () => { tip.hidden = true; views.forEach((w) => (w.hl.visible = false)); dirty = true; });
  }
  function inspect(v, e) {
    const r = v.renderer.domElement.getBoundingClientRect();
    ndc.set(((e.clientX - r.left) / r.width) * 2 - 1, -((e.clientY - r.top) / r.height) * 2 + 1);
    camera.aspect = r.width / r.height; camera.updateProjectionMatrix();
    ray.setFromCamera(ndc, camera);
    const hit = ray.intersectObject(v.mesh, false)[0];
    views.forEach((w) => (w.hl.visible = false)); dirty = true;
    if (!hit) { tip.hidden = true; return; }
    const g = v.geo, d = v.d;
    const qg = Math.floor(hit.faceIndex / (2 * SUB * SUB));
    const q = qg % g.nq, copy = Math.floor(qg / g.nq);
    const reg = g.quadRegion[q], elem = g.quadElem[q];
    const P = g.geometry.attributes.position;
    const ia = hit.face.a, ib = hit.face.b, ic = hit.face.c;
    const local = v.mesh.worldToLocal(hit.point.clone());
    tri.set(new THREE.Vector3().fromBufferAttribute(P, ia), new THREE.Vector3().fromBufferAttribute(P, ib), new THREE.Vector3().fromBufferAttribute(P, ic));
    tri.getBarycoord(local, bary);
    const T = bary.x * g.Tv[ia] + bary.y * g.Tv[ib] + bary.z * g.Tv[ic];
    const nodes = g.quads.slice(4 * q, 4 * q + 4);
    const peak = Math.max(...nodes.map((n) => g.nodePeak[n]));
    const lim = v.limits[reg];
    const jd = d.current_density ? d.current_density.elem[elem] : undefined;
    const status = T > lim ? `<b class="rv-fail">✗ OVER LIMIT by ${fmt(T - lim, 1)} °C</b> (limit ${fmt(lim, 0)} °C)` : `✓ ${fmt(lim - T, 1)} °C below the ${fmt(lim, 0)} °C limit`;
    tip.innerHTML = `<b>${DESIGN_LABEL[v.key] || v.key}</b> · ${SURF_LABEL[SURF[reg]]}<br>` +
      `element ${elem} (${copy ? "mirror half" : "model half"}) · x = ${fmt(local.x, 1)} mm<br>` +
      `T now <b>${fmt(T, 1)} °C</b> · peak over the history ${fmt(peak, 1)} °C<br>` +
      (jd !== undefined ? `current density at the hottest moment ${fmt(jd, 2)} A/mm²<br>` : "") + status;
    tip.hidden = false;
    const tw = tip.offsetWidth, th = tip.offsetHeight;
    let x = e.clientX + 14, y = e.clientY + 14;
    if (x + tw > window.innerWidth - 8) x = e.clientX - tw - 14;
    if (y + th > window.innerHeight - 8) y = e.clientY - th - 14;
    tip.style.left = Math.max(8, x) + "px"; tip.style.top = Math.max(8, y) + "px";
    const sy = copy ? -1 : 1;
    v.hl.geometry.setFromPoints(nodes.map((n) => v.mesh.localToWorld(new THREE.Vector3(g.X[3 * n], sy * g.X[3 * n + 1], g.X[3 * n + 2]))));
    v.hl.visible = true;
  }

  function stop() { playing = false; playBtn.textContent = "▶ Play"; }
  playBtn.addEventListener("click", () => {
    playing = !playing;
    if (playing && tNow >= tEnd - 1e-6) setTime(tStart);
    playBtn.textContent = playing ? "❚❚ Pause" : "▶ Play";
    lastMs = performance.now();
  });
  slider.addEventListener("input", () => { stop(); setTime(+slider.value); });
  exSel.addEventListener("change", () => { ex = +exSel.value; views.forEach((v) => (v.mesh.scale.z = ex)); setNote(ex); dirty = true; });
  if (cmpSel) cmpSel.addEventListener("change", () => { rightKey = cmpSel.value; setDesign(right, rightKey); setTime(tNow); });

  function resize() { for (const v of views) v.renderer.setSize(v.wrap.clientWidth, v.wrap.clientHeight, false); drawChart(); dirty = true; }
  new ResizeObserver(resize).observe(ROOT);
  function loop(ms) {
    requestAnimationFrame(loop);
    if (playing) {
      const nt = tNow + (((ms - lastMs) / 1000) * (tEnd - tStart)) / PLAY_SECONDS;
      lastMs = ms;
      if (nt >= tEnd) { setTime(tEnd); playing = false; playBtn.textContent = "▶ Replay"; } else setTime(nt);
    }
    if (!dirty) return;
    dirty = false;
    for (const v of views) {
      const w = v.wrap.clientWidth, h = v.wrap.clientHeight;
      if (!w || !h) continue;
      camera.aspect = w / h; camera.updateProjectionMatrix();
      v.renderer.render(v.scene, camera);
    }
  }
  setView(); resize();
  const hh = D[rightKey || leftKey].history;
  if (hh && hh.peak_bar_C) setTime(hh.times[hh.peak_bar_C.indexOf(Math.max(...hh.peak_bar_C))]); else setTime(tStart);
  window.__rotorViewer = { ready: true };
  requestAnimationFrame(loop);
}

// Boundary faces of the half-bar hex mesh (sym_y hidden), mirrored about y = 0.
// Each face is drawn as SUB x SUB sub-quads from the face's bilinear shape functions.
function buildGeometry(d) {
  const M = d.mesh, nn = M.nodes.length;
  const X = new Float32Array(3 * nn);
  M.nodes.forEach((p, i) => { X[3 * i] = p[0] * 1e3; X[3 * i + 1] = p[1] * 1e3; X[3 * i + 2] = p[2] * 1e3; });
  let xMax = 0;
  for (let i = 0; i < nn; i++) xMax = Math.max(xMax, X[3 * i]);
  const quads = [], quadRegion = [], quadElem = [];
  const key = (q) => [...q].sort((a, b) => a - b).join(",");
  const elemOf = new Map();
  M.hexes.forEach((h, e) => { for (const f of HEX_FACES) elemOf.set(key(f.map((i) => h[i])), e); });
  for (const [name, list] of Object.entries(M.faces)) {
    if (name.startsWith("sym")) continue;
    const ri = SURF.indexOf(name);
    for (const q of list) { quads.push(...q); quadRegion.push(ri < 0 ? 6 : ri); quadElem.push(elemOf.get(key(q)) ?? -1); }
  }
  const nq = quadRegion.length, copies = [1, -1];
  const per = (SUB + 1) * (SUB + 1), nv = 2 * nq * per;
  const pos = new Float32Array(3 * nv), color = new Float32Array(3 * nv), over = new Float32Array(nv);
  const corner = new Int32Array(4 * nv), weight = new Float32Array(4 * nv), vRegion = new Uint8Array(nv);
  const index = new Uint32Array(2 * nq * SUB * SUB * 6);
  let k = 0, ii = 0;
  copies.forEach((sy) => {
    for (let q = 0; q < nq; q++) {
      const n = quads.slice(4 * q, 4 * q + 4), base = k;
      for (let b = 0; b <= SUB; b++) for (let a = 0; a <= SUB; a++) {
        const s = a / SUB, t = b / SUB, w = [(1 - s) * (1 - t), s * (1 - t), s * t, (1 - s) * t];
        for (let c = 0; c < 3; c++) pos[3 * k + c] = (c === 1 ? sy : 1) * (w[0] * X[3 * n[0] + c] + w[1] * X[3 * n[1] + c] + w[2] * X[3 * n[2] + c] + w[3] * X[3 * n[3] + c]);
        for (let j = 0; j < 4; j++) { corner[4 * k + j] = n[j]; weight[4 * k + j] = w[j]; }
        vRegion[k] = quadRegion[q]; k++;
      }
      for (let b = 0; b < SUB; b++) for (let a = 0; a < SUB; a++) {
        const i0 = base + b * (SUB + 1) + a, i1 = i0 + 1, i3 = i0 + SUB + 1, i2 = i3 + 1;
        if (sy < 0) index.set([i0, i2, i1, i0, i3, i2], ii); else index.set([i0, i1, i2, i0, i2, i3], ii);
        ii += 6;
      }
    }
  });
  const geometry = new THREE.BufferGeometry();
  geometry.setAttribute("position", new THREE.BufferAttribute(pos, 3));
  geometry.setAttribute("color", new THREE.BufferAttribute(color, 3).setUsage(THREE.DynamicDrawUsage));
  geometry.setAttribute("over", new THREE.BufferAttribute(over, 1).setUsage(THREE.DynamicDrawUsage));
  geometry.setIndex(new THREE.BufferAttribute(index, 1));
  geometry.computeVertexNormals(); geometry.computeBoundingSphere();
  const nodePeak = new Float32Array(nn).fill(-Infinity);
  for (const fr of d.temperature.frames) for (let i = 0; i < nn; i++) if (fr[i] > nodePeak[i]) nodePeak[i] = fr[i];
  return { geometry, corner, weight, color, over, vRegion, Tv: new Float32Array(nv), nq, quads, quadRegion, quadElem, X, xMid: xMax / 2, nodePeak };
}

function hatchMaterial() {
  const m = new THREE.MeshLambertMaterial({ vertexColors: true, side: THREE.DoubleSide });
  m.onBeforeCompile = (sh) => {
    sh.vertexShader = sh.vertexShader.replace("#include <common>", "#include <common>\nattribute float over;\nvarying float vOver;")
      .replace("#include <begin_vertex>", "#include <begin_vertex>\nvOver = over;");
    sh.fragmentShader = sh.fragmentShader.replace("#include <common>", "#include <common>\nvarying float vOver;")
      .replace("#include <color_fragment>", "#include <color_fragment>\nif (vOver > 0.5 && mod(gl_FragCoord.x + gl_FragCoord.y, 9.0) < 3.0) diffuseColor.rgb = vec3(0.02);");
  };
  return m;
}

function fillMetrics(box, data, d, key) {
  const init = data.designs.initial;
  const dm = init && key !== "initial" ? ` <small>(${fmt((100 * (d.mass_kg - init.mass_kg)) / init.mass_kg, 1)} %)</small>` : "";
  const rows = [`<tr><th scope="row">Objective: busbar mass</th><td>${fmt(1000 * d.mass_kg, 1)} g${dm}</td><td>minimize</td><td></td></tr>`];
  for (const c of constraintRows(data, d)) {
    const dig = c.unit === "mm" ? 2 : 1;
    rows.push(`<tr><th scope="row">${c.name}</th><td>${fmt(c.value, dig)} ${c.unit}</td><td>${c.sense === ">=" ? "≥" : "≤"} ${fmt(c.limit, dig)}</td>` +
      `<td>${c.pass ? '<span class="rv-pass">✓ pass</span>' : '<span class="rv-fail">✗ FAIL</span>'}</td></tr>`);
  }
  box.innerHTML = `<table><thead><tr><th>Quantity</th><th>Value</th><th>Limit</th><th>Status</th></tr></thead><tbody>${rows.join("")}</tbody></table>` +
    `<p class="rv-small">Temperatures are true peaks over the whole history; the optimizer used their KS aggregates (a few kelvin above).</p>`;
}

function drawColorbar(box, tMin, tMax, design) {
  const stops = RAMP.map((c, i) => `${c} ${(100 * i) / (RAMP.length - 1)}%`).join(",");
  const L = limitFor(design);
  const marks = [[L[0], "joint limit"], [L[3], "insulation limit"]].filter(([v]) => Number.isFinite(v) && v >= tMin && v <= tMax);
  const ticks = [0, 0.25, 0.5, 0.75, 1].map((u) => `<span style="left:${100 * u}%">${Math.round(tMin + u * (tMax - tMin))}</span>`).join("");
  const lims = marks.map(([v, n]) => { const u = (100 * (v - tMin)) / (tMax - tMin); return `<i style="left:${u}%"><b class="${u > 60 ? "rv-lab-l" : ""}">${n} ${fmt(v, 0)}</b></i>`; }).join("");
  box.innerHTML = `<div class="rv-cb-title">Temperature (°C), same scale for every design</div>
    <div class="rv-cb" style="background:linear-gradient(90deg,${stops})">${lims}</div>
    <div class="rv-cb-ticks">${ticks}</div>
    <div class="rv-cb-key"><span class="rv-hatch"></span> hatched = above its limit (joint faces: joint limit; everywhere else: insulation limit)</div>`;
}
