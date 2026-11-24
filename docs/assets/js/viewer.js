// 3D initial-vs-optimized brake rotor viewer. Everything comes from results.json.
// Handout viewer requirements, and where each is met:
//  1 side by side, linked cameras ............ one shared camera, two OrbitControls
//  2 true 3D geometry ......................... quarter sector mirrored + rotated, no exaggeration
//  3 single-hue temperature scale + colour bar, over-limit faces hatched (not colour alone)
//  4 field on the geometry with the element's shape functions (bilinear on each hex face)
//  5 animation: heating and cooling through the braking stops (play / scrub)
//  6 hover or tap: element id, surface, radius, temperature now and peak, limit margin
//  7 metrics panel per design: mass and every constraint with pass / fail
//  8 loads results.json, responsive to phone width, renders in about a second
import * as THREE from "three";
import { OrbitControls } from "three/addons/controls/OrbitControls.js";
import { loadResults, fmt, constraintRows, DESIGN_LABEL } from "./data.js";

const ROOT = document.getElementById("rotor-viewer");
const SUB = 3; // each FE face drawn as SUB x SUB sub-quads from its bilinear shape functions
const HEX_FACES = [[0, 3, 2, 1], [4, 5, 6, 7], [0, 1, 5, 4], [2, 3, 7, 6], [0, 4, 7, 3], [1, 2, 6, 5]];
const SURF = ["friction", "face_free", "channel", "rim", "hub", "other"];
const SURF_LABEL = {
  friction: "friction face (pad rubbing band)",
  face_free: "cheek face inside the pad band",
  channel: "vane / cooling channel wall",
  rim: "outer rim",
  hub: "hub / bell-mount face",
  other: "surface",
};
// ColorBrewer "Oranges": single hue, light = cool, dark = hot.
const RAMP = ["#fff5eb", "#fee6ce", "#fdd0a2", "#fdae6b", "#fd8d3c", "#f16913", "#d94801", "#a63603", "#7f2704"];
const RAMP_LIN = RAMP.map((h) => new THREE.Color(h));

const el = (tag, cls, html) => {
  const e = document.createElement(tag);
  if (cls) e.className = cls;
  if (html !== undefined) e.innerHTML = html;
  return e;
};

if (ROOT) {
  ROOT.innerHTML = '<p class="rv-loading">Loading 3D comparison…</p>';
  loadResults()
    .then((data) => { ROOT.innerHTML = ""; build(data); })
    .catch((err) => {
      console.error(err);
      ROOT.innerHTML = `<p class="rv-error">The 3D viewer could not load its data: ${err.message}</p>`;
    });
}

function rampColor(u, out) {
  u = Math.min(1, Math.max(0, u)) * (RAMP_LIN.length - 1);
  const i = Math.min(RAMP_LIN.length - 2, Math.floor(u));
  const f = u - i, a = RAMP_LIN[i], b = RAMP_LIN[i + 1];
  out[0] = a.r + (b.r - a.r) * f; out[1] = a.g + (b.g - a.g) * f; out[2] = a.b + (b.b - a.b) * f;
}

function limitFor(design) {
  const lim = {};
  for (const c of design.constraints || []) {
    const L = c.limit_C ?? c.limit;
    if (c.name === "peak_face_T") lim.friction = L;
    if (c.name === "peak_hub_T") lim.hub = L;
  }
  return SURF.map((s) => (lim[s] === undefined ? Infinity : lim[s]));
}

function build(data) {
  const D = data.designs;
  const rightKeys = ["optimized", "optimized_linear"].filter((k) => D[k]);
  const leftKey = D.initial ? "initial" : rightKeys.shift();
  let rightKey = rightKeys[0];

  // ---------- colour range shared by every design ----------
  let tMin = Infinity, tMax = -Infinity;
  for (const k of Object.keys(D)) for (const fr of D[k].temperature.frames) for (const v of fr) {
    if (v < tMin) tMin = v; if (v > tMax) tMax = v;
  }
  const step = tMax - tMin > 300 ? 50 : 10;
  tMin = Math.floor(tMin / step) * step;
  const cMax = Math.ceil(tMax / step) * step;

  // ---------- DOM ----------
  if (data.meta && data.meta.mock)
    ROOT.appendChild(el("div", "rv-mock", "MOCK DATA: placeholders until <code>python scripts/run_all.py</code> writes results.json."));
  const bar = el("div", "rv-toolbar");
  const playBtn = el("button", "rv-btn", "▶ Play");
  playBtn.type = "button";
  playBtn.setAttribute("aria-label", "Play or pause the braking animation");
  const times = D[leftKey].temperature.times;
  const tStart = times[0], tEnd = times[times.length - 1];
  const slider = el("input", "rv-slider");
  Object.assign(slider, { type: "range", min: tStart, max: tEnd, step: 0.01, value: tStart });
  slider.setAttribute("aria-label", "Time in the braking sequence");
  const tLabel = el("span", "rv-tlabel");
  const viewSel = el("select", "rv-select");
  viewSel.setAttribute("aria-label", "How much of the rotor to show");
  for (const [v, t] of [["cut", "Cutaway (see the vanes)"], ["full", "Full rotor"], ["sector", "One vane pitch"]])
    viewSel.appendChild(Object.assign(document.createElement("option"), { value: v, textContent: t }));
  bar.append(playBtn, slider, tLabel, viewSel);
  let cmpSel = null;
  if (rightKeys.length > 1) {
    cmpSel = el("select", "rv-select");
    cmpSel.setAttribute("aria-label", "Design shown on the right");
    for (const k of rightKeys) cmpSel.appendChild(Object.assign(document.createElement("option"), { value: k, textContent: "Right: " + DESIGN_LABEL[k] }));
    bar.appendChild(cmpSel);
  }
  ROOT.appendChild(bar);
  const panels = el("div", "rv-panels");
  ROOT.appendChild(panels);
  const legend = el("div", "rv-legend");
  ROOT.appendChild(legend);
  drawColorbar(legend, tMin, cMax, D[leftKey]);
  const chartWrap = el("div", "rv-chart");
  const chart = el("canvas");
  chartWrap.appendChild(chart);
  ROOT.appendChild(chartWrap);
  const nV = D[leftKey].mesh.n_vanes;
  ROOT.appendChild(el("p", "rv-note",
    `True geometry, no exaggeration. The solver models one quarter of a vane pitch; the viewer mirrors it about the vane ` +
    `mid-plane and the disc mid-plane, then repeats it ${nV} times around the axis. Drag to rotate (both discs move together), ` +
    `scroll or pinch to zoom, hover or tap a face to inspect it. Hatched faces are above their temperature limit. ` +
    `Lighting darkens some faces a little; the tooltip gives the exact value.`));
  const tip = el("div", "rv-tip");
  tip.hidden = true;
  document.body.appendChild(tip);

  // ---------- shared camera = linked views ----------
  const camera = new THREE.PerspectiveCamera(35, 1, 1, 5000);
  camera.up.set(0, 0, 1);
  const target = new THREE.Vector3();
  let dirty = true;

  const left = makeView(leftKey);
  const right = rightKey ? makeView(rightKey) : null;
  const views = [left, right].filter(Boolean);

  function makeView(key) {
    const panel = el("div", "rv-panel");
    const head = el("div", "rv-head");
    const wrap = el("div", "rv-canvas");
    const renderer = new THREE.WebGLRenderer({ antialias: true, alpha: true });
    renderer.setPixelRatio(Math.min(window.devicePixelRatio, 2));
    wrap.appendChild(renderer.domElement);
    const metrics = el("div", "rv-metrics");
    panel.append(head, wrap, metrics);
    panels.appendChild(panel);
    const scene = new THREE.Scene();
    scene.add(new THREE.HemisphereLight(0xffffff, 0xaaaaaa, 2.2));
    const sun = new THREE.DirectionalLight(0xffffff, 0.8);
    sun.position.set(1, -1, 2);
    scene.add(sun);
    const hl = new THREE.LineLoop(new THREE.BufferGeometry().setFromPoints([0, 1, 2, 3].map(() => new THREE.Vector3())),
      new THREE.LineBasicMaterial({ color: 0x1f6feb, depthTest: false }));
    hl.renderOrder = 2; hl.visible = false;
    scene.add(hl);
    const controls = new OrbitControls(camera, renderer.domElement);
    controls.target = target; // the same Vector3 object for both views
    controls.enablePan = false;
    controls.minDistance = 60;
    controls.maxDistance = 1500;
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
    v.mesh = new THREE.InstancedMesh(v.geo.geometry, hatchMaterial(), d.mesh.n_vanes);
    v.scene.add(v.mesh);
    v.head.innerHTML = `<strong>${DESIGN_LABEL[key] || key}</strong> <span>${fmt(d.mass_kg, 2)} kg</span>`;
    fillMetrics(v.metrics, data, d, key);
  }

  // ---------- view modes (which vane pitches are drawn) ----------
  let mode = "cut";
  function layout() {
    const m = new THREE.Matrix4();
    for (const v of views) {
      const n = v.d.mesh.n_vanes, dTh = (2 * Math.PI) / n;
      let ids = [...Array(n).keys()];
      if (mode === "sector") ids = [0];
      if (mode === "cut") { const gap = Math.round(n / 4); ids = ids.filter((i) => i < n - gap - 1 || i === n - 1); }
      ids.forEach((s, slot) => v.mesh.setMatrixAt(slot, m.makeRotationZ(s * dTh - Math.PI / 4 + Math.PI / n)));
      v.mesh.count = ids.length;
      v.mesh.instanceMatrix.needsUpdate = true;
      v.mesh.computeBoundingSphere();
      v.slotSector = ids;
      v.hl.visible = false;
    }
    dirty = true;
  }
  function setMode(m) {
    mode = m;
    layout();
    const rO = left.geo.rMax;
    if (mode === "sector") {
      const a = -Math.PI / 4 + Math.PI / left.d.mesh.n_vanes;
      target.set(0.8 * rO * Math.cos(a), 0.8 * rO * Math.sin(a), 0);
      camera.position.set(target.x + 0.9 * rO, target.y - 0.9 * rO, 0.6 * rO);
    } else {
      target.set(0, 0, 0);
      camera.position.set(1.6 * rO, -2.0 * rO, 1.45 * rO);
    }
    views.forEach((v) => v.controls.update());
  }

  // ---------- time and colouring ----------
  let tNow = tStart, playing = false, lastMs = 0;
  const PLAY_SECONDS = 15;
  function frameAt(ts, t) {
    let i = 0;
    while (i < ts.length - 2 && ts[i + 1] < t) i++;
    const f = Math.min(1, Math.max(0, (t - ts[i]) / (ts[i + 1] - ts[i] || 1)));
    return [i, f];
  }
  const H = D[leftKey].history || {};
  const hT = H.times || times, spd = H.speed_kmh || null;
  const isBraking = (i) => (H.braking ? !!H.braking[i + 1] : spd[i + 1] < spd[i]);
  function setTime(t) {
    tNow = Math.min(tEnd, Math.max(tStart, t));
    slider.value = tNow;
    let txt = `t = ${fmt(tNow, 1)} s`;
    if (spd) {
      const [i, f] = frameAt(hT, tNow);
      const s = spd[i] + (spd[i + 1] - spd[i]) * f;
      txt += ` · ${fmt(s, 0)} km/h · ${isBraking(i) ? "braking" : "accelerating / cooling"}`;
    }
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
      for (let j = 0; j < 4; j++) {
        const id = corner[4 * k + j];
        T += weight[4 * k + j] * (A[id] + (B[id] - A[id]) * f);
      }
      Tv[k] = T;
      rampColor((T - tMin) / span, c);
      color[3 * k] = c[0]; color[3 * k + 1] = c[1]; color[3 * k + 2] = c[2];
      over[k] = T > v.limits[vRegion[k]] ? 1 : 0;
    }
    v.geo.geometry.attributes.color.needsUpdate = true;
    v.geo.geometry.attributes.over.needsUpdate = true;
  }

  // ---------- chart: peak face temperature and speed through the sequence ----------
  const SERIES = { initial: "#6b7280", optimized: "#c2410c", optimized_linear: "#2563eb" };
  const faceLimit = left.limits[0];
  function drawChart() {
    const W = chartWrap.clientWidth, Hh = 160, dpr = Math.min(window.devicePixelRatio, 2);
    if (!W) return;
    chart.width = W * dpr; chart.height = Hh * dpr; chart.style.height = Hh + "px";
    const g = chart.getContext("2d");
    g.setTransform(dpr, 0, 0, dpr, 0, 0);
    const css = getComputedStyle(ROOT);
    const ink = css.getPropertyValue("--rv-ink").trim() || "#222";
    const muted = css.getPropertyValue("--rv-muted").trim() || "#888";
    const L = 44, R = 40, T = 22, B = 26, pw = W - L - R, ph = Hh - T - B;
    const X = (t) => L + ((t - tStart) / (tEnd - tStart)) * pw;
    const Y = (v) => T + ph * (1 - (v - tMin) / (cMax - tMin));
    const vmax = spd ? Math.max(...spd) * 1.05 : 1;
    const Yv = (v) => T + ph * (1 - v / vmax);
    g.clearRect(0, 0, W, Hh);
    g.font = "11px system-ui, sans-serif";
    if (spd) {
      g.fillStyle = "rgba(127,127,127,0.14)";
      for (let i = 0; i < hT.length - 1; i++) if (isBraking(i)) g.fillRect(X(hT[i]), T, X(hT[i + 1]) - X(hT[i]) + 0.5, ph);
    }
    g.strokeStyle = muted; g.lineWidth = 1; g.strokeRect(L, T, pw, ph);
    g.fillStyle = muted; g.textAlign = "right";
    for (const v of [tMin, (tMin + cMax) / 2, cMax]) g.fillText(fmt(v, 0), L - 4, Y(v) + 4);
    g.textAlign = "center";
    for (const t of [tStart, (tStart + tEnd) / 2, tEnd]) g.fillText(fmt(t, 0), X(t), T + ph + 13);
    g.fillText("time (s)", L + pw / 2, Hh - 3);
    g.save(); g.translate(11, T + ph / 2); g.rotate(-Math.PI / 2); g.fillText("peak face T (°C)", 0, 0); g.restore();
    if (spd) {
      g.textAlign = "left";
      for (const v of [0, Math.round(vmax / 1.05)]) g.fillText(String(v), L + pw + 4, Yv(v) + 4);
      g.save(); g.translate(W - 5, T + ph / 2); g.rotate(Math.PI / 2); g.textAlign = "center"; g.fillText("speed (km/h)", 0, 0); g.restore();
      g.strokeStyle = muted; g.setLineDash([2, 3]);
      g.beginPath(); hT.forEach((t, i) => (i ? g.lineTo(X(t), Yv(spd[i])) : g.moveTo(X(t), Yv(spd[i])))); g.stroke();
      g.setLineDash([]);
    }
    if (Number.isFinite(faceLimit)) {
      g.strokeStyle = ink; g.setLineDash([5, 3]);
      g.beginPath(); g.moveTo(L, Y(faceLimit)); g.lineTo(L + pw, Y(faceLimit)); g.stroke(); g.setLineDash([]);
      g.fillStyle = ink; g.textAlign = "right"; g.fillText(`face limit ${fmt(faceLimit, 0)} °C`, L + pw - 4, Y(faceLimit) - 4);
    }
    const shown = [leftKey, rightKey].filter(Boolean);
    let lx = L + 4;
    g.textAlign = "left";
    for (const k of shown) {
      const h = D[k].history;
      if (!h || !h.peak_face_C) continue;
      g.strokeStyle = SERIES[k] || ink; g.lineWidth = 2;
      g.beginPath(); h.times.forEach((t, i) => (i ? g.lineTo(X(t), Y(h.peak_face_C[i])) : g.moveTo(X(t), Y(h.peak_face_C[i])))); g.stroke();
      g.fillStyle = SERIES[k] || ink; g.fillRect(lx, 8, 14, 3);
      g.fillStyle = ink; g.fillText(DESIGN_LABEL[k] || k, lx + 18, 13);
      lx += 30 + g.measureText(DESIGN_LABEL[k] || k).width;
    }
    if (spd && lx < W - 120) { g.fillStyle = muted; g.fillText("┈ speed   ▒ braking", lx, 13); }
    g.lineWidth = 1;
    g.strokeStyle = "#1f6feb"; g.lineWidth = 1.5;
    g.beginPath(); g.moveTo(X(tNow), T); g.lineTo(X(tNow), T + ph); g.stroke();
    chart._x2t = (x) => tStart + ((x - L) / pw) * (tEnd - tStart);
  }
  const scrub = (e) => {
    const r = chart.getBoundingClientRect();
    stop();
    setTime(chart._x2t(e.clientX - r.left));
  };
  chart.addEventListener("pointerdown", (e) => { scrub(e); chart.setPointerCapture(e.pointerId); });
  chart.addEventListener("pointermove", (e) => { if (e.buttons) scrub(e); });

  // ---------- hover / tap inspection ----------
  const ray = new THREE.Raycaster(), ndc = new THREE.Vector2(), tri = new THREE.Triangle(), bary = new THREE.Vector3();
  const inst = new THREE.Matrix4(), instInv = new THREE.Matrix4();
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
    views.forEach((w) => (w.hl.visible = false));
    dirty = true;
    if (!hit) { tip.hidden = true; return; }
    const g = v.geo, d = v.d;
    const qg = Math.floor(hit.faceIndex / (2 * SUB * SUB));
    const q = qg % g.nq, copy = Math.floor(qg / g.nq);
    const reg = g.quadRegion[q], elem = g.quadElem[q];
    v.mesh.getMatrixAt(hit.instanceId, inst);
    instInv.copy(inst).invert();
    const local = hit.point.clone().applyMatrix4(instInv);
    const P = g.geometry.attributes.position;
    const ia = hit.face.a, ib = hit.face.b, ic = hit.face.c;
    tri.set(new THREE.Vector3().fromBufferAttribute(P, ia), new THREE.Vector3().fromBufferAttribute(P, ib), new THREE.Vector3().fromBufferAttribute(P, ic));
    tri.getBarycoord(local, bary);
    const T = bary.x * g.Tv[ia] + bary.y * g.Tv[ib] + bary.z * g.Tv[ic];
    const nodes = g.quads.slice(4 * q, 4 * q + 4);
    const peak = Math.max(...nodes.map((n) => g.nodePeak[n]));
    const lim = v.limits[reg];
    const rad = Math.hypot(local.x, local.y);
    let status = "no temperature limit on this surface";
    if (Number.isFinite(lim)) status = T > lim
      ? `<b class="rv-fail">✗ OVER LIMIT by ${fmt(T - lim, 0)} °C</b> (limit ${fmt(lim, 0)} °C)`
      : `✓ ${fmt(lim - T, 0)} °C below the ${fmt(lim, 0)} °C limit`;
    tip.innerHTML = `<b>${DESIGN_LABEL[v.key] || v.key}</b> · ${SURF_LABEL[SURF[reg]]}<br>` +
      `element ${elem} of the model sector (vane pitch ${v.slotSector[hit.instanceId] + 1}/${d.mesh.n_vanes}, mirror copy ${copy + 1}/4)<br>` +
      `radius ${fmt(rad, 1)} mm · T now <b>${fmt(T, 1)} °C</b> · peak over the sequence ${fmt(peak, 1)} °C<br>${status}`;
    tip.hidden = false;
    const tw = tip.offsetWidth, th = tip.offsetHeight;
    let x = e.clientX + 14, y = e.clientY + 14;
    if (x + tw > window.innerWidth - 8) x = e.clientX - tw - 14;
    if (y + th > window.innerHeight - 8) y = e.clientY - th - 14;
    tip.style.left = Math.max(8, x) + "px";
    tip.style.top = Math.max(8, y) + "px";
    const sgn = g.copies[copy];
    v.hl.geometry.setFromPoints(nodes.map((n) => new THREE.Vector3(g.X[3 * n], sgn[0] * g.X[3 * n + 1], sgn[1] * g.X[3 * n + 2]).applyMatrix4(inst)));
    v.hl.visible = true;
  }

  // ---------- controls ----------
  function stop() { playing = false; playBtn.textContent = "▶ Play"; }
  playBtn.addEventListener("click", () => {
    playing = !playing;
    if (playing && tNow >= tEnd - 1e-6) setTime(tStart);
    playBtn.textContent = playing ? "❚❚ Pause" : "▶ Play";
    lastMs = performance.now();
  });
  slider.addEventListener("input", () => { stop(); setTime(+slider.value); });
  viewSel.addEventListener("change", () => setMode(viewSel.value));
  if (cmpSel) cmpSel.addEventListener("change", () => {
    rightKey = cmpSel.value;
    setDesign(right, rightKey);
    layout(); setTime(tNow);
  });

  function resize() {
    for (const v of views) v.renderer.setSize(v.wrap.clientWidth, v.wrap.clientHeight, false);
    drawChart();
    dirty = true;
  }
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

  setMode("cut");
  resize();
  // open on the hottest moment of the design on the right, so the first view says something
  const hk = rightKey || leftKey;
  const hh = D[hk].history;
  if (hh && hh.peak_face_C) setTime(hh.times[hh.peak_face_C.indexOf(Math.max(...hh.peak_face_C))]);
  else setTime(tStart);
  window.__rotorViewer = { ready: true }; // for automated checks
  requestAnimationFrame(loop);
}

// Boundary faces of the quarter-sector hex mesh, mirrored into a full vane pitch.
// Each face is drawn as SUB x SUB sub-quads whose positions and temperatures come
// from the face's bilinear shape functions (the trilinear hex restricted to the face).
function buildGeometry(d) {
  const M = d.mesh;
  const nn = M.nodes.length;
  const X = new Float32Array(3 * nn);
  M.nodes.forEach((p, i) => { X[3 * i] = p[0] * 1e3; X[3 * i + 1] = p[1] * 1e3; X[3 * i + 2] = p[2] * 1e3; });
  const key = (q) => [...q].sort((a, b) => a - b).join(",");
  // face set membership (preferred) ...
  const setOf = new Map();
  const FS = M.faces || M.face_sets; // named boundary-face sets from the solver
  if (FS) for (const [name, list] of Object.entries(FS)) {
    if (name.startsWith("sym")) continue; // symmetry / periodic planes lie inside the full disc
    for (const q of list) setOf.set(key(q), name);
  }
  // ... boundary faces by counting
  const count = new Map(), first = new Map();
  M.hexes.forEach((h, e) => {
    for (const f of HEX_FACES) {
      const q = f.map((i) => h[i]), k = key(q);
      count.set(k, (count.get(k) || 0) + 1);
      if (!first.has(k)) first.set(k, [q, e]);
    }
  });
  let rMax = 0;
  for (let i = 0; i < nn; i++) rMax = Math.max(rMax, Math.hypot(X[3 * i], X[3 * i + 1]));
  const halfPitch = Math.PI / M.n_vanes;
  const quads = [], quadRegion = [], quadElem = [];
  for (const [k, c] of count) {
    if (c !== 1) continue;
    const [q, e] = first.get(k);
    let name = setOf.get(k);
    if (!name) {
      if (FS) continue; // not in any drawn set => symmetry / periodic plane, hidden inside the disc
      // no face sets: drop faces lying on z=0, theta=0 or theta=pi/n
      const zs = q.map((n) => Math.abs(X[3 * n + 2])), ths = q.map((n) => Math.atan2(X[3 * n + 1], X[3 * n]));
      if (Math.max(...zs) < 1e-6 || Math.max(...ths.map(Math.abs)) < 1e-6 || Math.max(...ths.map((t) => Math.abs(t - halfPitch))) < 1e-6) continue;
      name = "other";
    }
    const ri = SURF.indexOf(name);
    quads.push(...q); quadRegion.push(ri < 0 ? 5 : ri); quadElem.push(e);
  }
  const nq = quadRegion.length;
  const copies = [[1, 1], [-1, 1], [1, -1], [-1, -1]]; // (y sign, z sign): theta and z mirrors
  const per = (SUB + 1) * (SUB + 1), nv = 4 * nq * per;
  const pos = new Float32Array(3 * nv), color = new Float32Array(3 * nv), over = new Float32Array(nv);
  const corner = new Int32Array(4 * nv), weight = new Float32Array(4 * nv), vRegion = new Uint8Array(nv);
  const index = new Uint32Array(4 * nq * SUB * SUB * 6);
  let k = 0, ii = 0;
  copies.forEach(([sy, sz]) => {
    const flip = sy * sz < 0; // one mirror reverses the winding
    for (let q = 0; q < nq; q++) {
      const n = quads.slice(4 * q, 4 * q + 4), base = k;
      for (let b = 0; b <= SUB; b++) for (let a = 0; a <= SUB; a++) {
        const s = a / SUB, t = b / SUB;
        const w = [(1 - s) * (1 - t), s * (1 - t), s * t, (1 - s) * t];
        for (let c = 0; c < 3; c++) {
          const sg = c === 1 ? sy : c === 2 ? sz : 1;
          pos[3 * k + c] = sg * (w[0] * X[3 * n[0] + c] + w[1] * X[3 * n[1] + c] + w[2] * X[3 * n[2] + c] + w[3] * X[3 * n[3] + c]);
        }
        for (let j = 0; j < 4; j++) { corner[4 * k + j] = n[j]; weight[4 * k + j] = w[j]; }
        vRegion[k] = quadRegion[q];
        k++;
      }
      for (let b = 0; b < SUB; b++) for (let a = 0; a < SUB; a++) {
        const i0 = base + b * (SUB + 1) + a, i1 = i0 + 1, i3 = i0 + SUB + 1, i2 = i3 + 1;
        if (flip) index.set([i0, i2, i1, i0, i3, i2], ii); else index.set([i0, i1, i2, i0, i2, i3], ii);
        ii += 6;
      }
    }
  });
  const geometry = new THREE.BufferGeometry();
  geometry.setAttribute("position", new THREE.BufferAttribute(pos, 3));
  geometry.setAttribute("color", new THREE.BufferAttribute(color, 3).setUsage(THREE.DynamicDrawUsage));
  geometry.setAttribute("over", new THREE.BufferAttribute(over, 1).setUsage(THREE.DynamicDrawUsage));
  geometry.setIndex(new THREE.BufferAttribute(index, 1));
  geometry.computeVertexNormals();
  geometry.computeBoundingSphere();
  const nodePeak = new Float32Array(nn).fill(-Infinity);
  for (const fr of d.temperature.frames) for (let i = 0; i < nn; i++) if (fr[i] > nodePeak[i]) nodePeak[i] = fr[i];
  return { geometry, corner, weight, color, over, vRegion, Tv: new Float32Array(nv), nq, quads, quadRegion, quadElem, copies, X, rMax, nodePeak };
}

// Lambert shading with vertex colours, plus screen-space hatching where the
// interpolated "over" flag is set: over-limit areas never rely on colour alone.
function hatchMaterial() {
  const m = new THREE.MeshLambertMaterial({ vertexColors: true, side: THREE.DoubleSide });
  m.onBeforeCompile = (sh) => {
    sh.vertexShader = sh.vertexShader
      .replace("#include <common>", "#include <common>\nattribute float over;\nvarying float vOver;")
      .replace("#include <begin_vertex>", "#include <begin_vertex>\nvOver = over;");
    sh.fragmentShader = sh.fragmentShader
      .replace("#include <common>", "#include <common>\nvarying float vOver;")
      .replace("#include <color_fragment>",
        "#include <color_fragment>\nif (vOver > 0.5 && mod(gl_FragCoord.x + gl_FragCoord.y, 9.0) < 3.0) diffuseColor.rgb = vec3(0.02);");
  };
  return m;
}

function fillMetrics(box, data, d, key) {
  const init = data.designs.initial;
  const dm = init && key !== "initial" ? ` <small>(${fmt((100 * (d.mass_kg - init.mass_kg)) / init.mass_kg, 1)} %)</small>` : "";
  const rows = [`<tr><th scope="row">Objective: rotor mass</th><td>${fmt(d.mass_kg, 2)} kg${dm}</td><td>minimize</td><td></td></tr>`];
  for (const c of constraintRows(data, d)) {
    const val = c.trueMax !== undefined ? `${fmt(c.value, 1)} <small>(true peak ${fmt(c.trueMax, 1)})</small>` : fmt(c.value, c.unit === "mm" ? 2 : 1);
    rows.push(`<tr><th scope="row">${c.name}</th><td>${val} ${c.unit}</td><td>${c.sense === ">=" ? "≥" : "≤"} ${fmt(c.limit, c.unit === "mm" ? 1 : 0)}</td>` +
      `<td>${c.pass ? '<span class="rv-pass">✓ pass</span>' : '<span class="rv-fail">✗ FAIL</span>'}</td></tr>`);
  }
  box.innerHTML = `<table><thead><tr><th>Quantity</th><th>Value</th><th>Limit</th><th>Status</th></tr></thead><tbody>${rows.join("")}</tbody></table>` +
    `<p class="rv-small">Temperature constraints use the KS aggregate over all surface nodes and time steps (always at or above the true peak).</p>`;
}

function drawColorbar(box, tMin, tMax, design) {
  const stops = RAMP.map((c, i) => `${c} ${(100 * i) / (RAMP.length - 1)}%`).join(",");
  const L = limitFor(design);
  const marks = [[L[0], "face limit"], [L[4], "hub limit"]].filter(([v]) => Number.isFinite(v) && v >= tMin && v <= tMax);
  const ticks = [0, 0.25, 0.5, 0.75, 1].map((u) => `<span style="left:${100 * u}%">${Math.round(tMin + u * (tMax - tMin))}</span>`).join("");
  const lims = marks.map(([v, n]) => {
    const u = (100 * (v - tMin)) / (tMax - tMin);
    return `<i style="left:${u}%"><b class="${u > 60 ? "rv-lab-l" : ""}">${n} ${fmt(v, 0)}</b></i>`;
  }).join("");
  box.innerHTML = `<div class="rv-cb-title">Temperature (°C), same scale for every design</div>
    <div class="rv-cb" style="background:linear-gradient(90deg,${stops})">${lims}</div>
    <div class="rv-cb-ticks">${ticks}</div>
    <div class="rv-cb-key"><span class="rv-hatch"></span> hatched = above its limit (rubbing face: face limit, hub face: hub limit)</div>`;
}
