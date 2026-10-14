// Front end of the alignment GUI. Talks to the JSON API in app.py; all
// chemistry happens on the server, this file only shows and edits state.

const token = new URLSearchParams(location.search).get("token") || "";
const $ = (id) => document.getElementById(id);

const COLORS = { ref: "#8d96a0", refCore: "#6c7681", core: "#e8871e", rest: "#2aa6b8", overlay: "#7aa7b0" };

let state = null;          // last /api/state
let selected = null;       // name of the analog in the list
let view = null;           // /api/analogs/{name} for the selected analog
let refView = null;        // /api/reference
let overlay = null;        // /api/overlay, fetched when the overlay is switched on
let pollTimer = null;
let needZoom = true;
let viewer = null;
let noWebGL = false;
let colorReference = readSetting("colorReference") === "1";

// ---------------------------------------------------------------- helpers

// A remembered preference; the page works without storage (private windows, blocked cookies).
function readSetting(key) {
  try { return localStorage.getItem(`msldprep.${key}`); } catch (_) { return null; }
}
function saveSetting(key, value) {
  try { localStorage.setItem(`msldprep.${key}`, value); } catch (_) { /* not available */ }
}

async function api(path, { method = "GET", body } = {}) {
  const headers = { "X-Token": token };
  if (body !== undefined) headers["Content-Type"] = "application/json";
  const res = await fetch(path, { method, headers, body: body === undefined ? undefined : JSON.stringify(body) });
  if (!res.ok) {
    let detail = res.statusText;
    try {
      detail = (await res.json()).detail ?? detail;
    } catch (_) { /* not json */ }
    if (typeof detail !== "string") detail = detail.map((d) => d.msg).join("; ");
    throw new Error(detail);
  }
  return res.json();
}

function say(message, isError = false) {
  const banner = $("banner");
  banner.textContent = message;
  banner.classList.toggle("error", isError);
}

async function guarded(action) {
  try {
    await action();
  } catch (err) {
    say(err.message, true);
  }
}

const el = (tag, className, text) => {
  const node = document.createElement(tag);
  if (className) node.className = className;
  if (text !== undefined) node.textContent = text;
  return node;
};

// ------------------------------------------------------------------ render

function renderState() {
  const ref = state.reference;
  $("reference-info").replaceChildren();
  if (ref) {
    const strong = el("strong", null, ref.name);
    const sign = (n) => (n > 0 ? `+${n}` : String(n));
    let note = ` · ${ref.atoms} atoms (with hydrogens) · net charge ${sign(ref.net_charge)}`;
    if (ref.loaded_net_charge !== ref.net_charge) note += ` (file: ${sign(ref.loaded_net_charge)})`;
    if (ref.charged_atoms) note += ` · ${ref.charged_atoms} charged atom${ref.charged_atoms > 1 ? "s" : ""}`;
    $("reference-info").append(strong, note);
  }

  const analogs = state.analogs;
  $("count").textContent = analogs.length ? `(${analogs.length})` : "";
  const list = $("list");
  list.replaceChildren();
  for (const a of analogs) {
    const li = el("li");
    li.tabIndex = 0;
    li.setAttribute("role", "option");
    li.setAttribute("aria-selected", String(a.name === selected));
    li.dataset.name = a.name;
    li.append(el("span", "name", a.name));
    const badges = el("span");
    const status = a.decision !== "undecided" ? a.decision : a.status;
    badges.append(el("span", `badge ${status}`, status));
    li.append(badges);
    const meta = a.status === "aligned"
      ? `RMSD ${a.rmsd.toFixed(3)} Å · ${a.core_atoms} core atoms · seed ${a.seed}`
      : a.status === "failed" ? a.error : "not aligned yet";
    li.append(el("span", "meta", meta));
    li.addEventListener("click", () => guarded(() => select(a.name)));
    li.addEventListener("keydown", (e) => {
      if (e.key === "Enter" || e.key === " ") { e.preventDefault(); guarded(() => select(a.name)); }
    });
    list.append(li);
  }

  const job = state.job;
  const progress = $("progress");
  progress.hidden = !job.running && job.done === 0;
  progress.max = Math.max(job.total, 1);
  progress.value = job.done;
  $("align").disabled = job.running || !ref || analogs.length === 0;
  $("align").textContent = job.running ? `Aligning ${job.done}/${job.total}…` : "Align all";
  if (!noWebGL) $("empty").hidden = Boolean(ref);
  updateDownload();
}

function syncParams() {
  const p = state.params;
  $("p-complete").checked = p.complete_rings_only;
  $("p-ringring").checked = p.ring_matches_ring_only;
  $("p-anybond").checked = p.ignore_bond_orders;
  $("p-timeout").value = p.timeout;
  $("p-seed").value = p.seed;
  $("neutralize").checked = state.neutralize_reference;
  $("t-color").checked = colorReference;
}

function renderDetail() {
  const box = $("detail");
  box.replaceChildren();
  $("actions").hidden = !view || view.status === "pending";
  if (!view) {
    box.append(el("p", "hint", "Select an analog from the list."));
    $("detail-title").textContent = "Details";
    return;
  }
  $("detail-title").textContent = view.name;
  const dl = el("dl");
  const add = (term, value) => { dl.append(el("dt", null, term), el("dd", null, value)); };
  add("Status", view.decision !== "undecided" ? `${view.status}, ${view.decision}` : view.status);
  add("SMILES", view.smiles);
  if (view.status === "aligned") {
    add("RMSD of core", `${view.rmsd.toFixed(3)} Å`);
    add("Core atoms", String(view.core_atoms.length));
    add("Seed", String(view.seed));
    $("r-seed").value = view.seed + 1;
  }
  box.append(dl);
  if (view.error) box.append(el("p", "error", view.error));
  if (view.svg) {
    const svg = el("div", "svg");
    svg.innerHTML = view.svg; // drawn by RDKit on the server from the molecule only
    box.append(svg);
    box.append(el("p", "hint", "Matched core in blue."));
  }
  $("accept").disabled = $("reject").disabled = view.status !== "aligned";
}

function updateDownload() {
  const any = state && state.analogs.some((a) => a.status === "aligned");
  $("download").disabled = !any;
  const mol2 = $("d-format").value === "mol2";
  $("d-hint").textContent = (mol2
    ? "mol2 as in the notebook: Sybyl atom types, new atom names, no charges. "
    : "sdf keeps each analog's RMSD and seed. ")
    + "Nothing is stored on the server.";
}

async function download() {
  const query = new URLSearchParams({
    format: $("d-format").value,
    which: $("d-which").value,
    reference: $("d-reference").checked,
    per_file: $("d-perfile").checked,
  });
  const res = await fetch(`/api/export?${query}`, { headers: { "X-Token": token } });
  if (!res.ok) {
    let detail = res.statusText;
    try { detail = (await res.json()).detail ?? detail; } catch (_) { /* not json */ }
    throw new Error(typeof detail === "string" ? detail : detail.map((d) => d.msg).join("; "));
  }
  const name = /filename="([^"]+)"/.exec(res.headers.get("Content-Disposition") || "")?.[1] || "aligned";
  const url = URL.createObjectURL(await res.blob());
  const link = el("a");
  link.href = url;
  link.download = name;
  document.body.append(link);
  link.click();
  link.remove();
  URL.revokeObjectURL(url);
  say(`Downloaded ${name}.`);
}

// ---------------------------------------------------------------- 3D view

function ensureViewer() {
  if (!viewer && !noWebGL) {
    try {
      viewer = $3Dmol.createViewer("viewer", { backgroundAlpha: 0 });
    } catch (err) {
      const empty = $("empty");
      empty.textContent = "The 3D view needs WebGL, which this browser does not offer. The lists, drawings and downloads still work.";
      empty.hidden = false;
      viewer = null;
      noWebGL = true;
    }
  }
  return viewer;
}

function jmolHex(symbol) {
  const value = $3Dmol.elementColors.Jmol[symbol];
  return value === undefined ? "#888888" : `#${value.toString(16).padStart(6, "0")}`;
}

function renderElementLegend() {
  const box = $("elements");
  box.replaceChildren();
  box.classList.toggle("on", colorReference && Boolean(refView));
  if (!colorReference || !refView) return;
  const shown = new Set(refView.elements.filter((e) => e !== "H" || $("t-hydrogens").checked));
  for (const symbol of refView.elements) {
    if (!shown.has(symbol)) continue;
    const chip = el("span", "el", symbol);
    chip.style.setProperty("--c", jmolHex(symbol));
    box.append(chip);
  }
}

// The reference's 2D drawing in panel 1; with an aligned analog selected, the atoms it
// matched are marked.
function renderReference2D() {
  const box = $("ref2d");
  box.hidden = !refView;
  if (!refView) return;
  const matched = view && view.status === "aligned";
  $("ref2d-svg").innerHTML = matched ? view.reference_svg : refView.svg; // drawn by RDKit on the server
  $("ref2d-hint").textContent = matched ? `Matched core of ${view.name} in blue.` : "";
}

function draw() {
  renderReference2D();
  renderElementLegend();
  if (!refView) return;
  const v = ensureViewer();
  if (!v) return;
  v.removeAllModels();
  const ref = v.addModel(refView.molblock, "sdf");
  // Gray, or Jmol element colors (carbon stays gray) when the toggle is on.
  const refStick = (radius, gray) =>
    ({ stick: colorReference ? { radius, colorscheme: "Jmol" } : { radius, color: gray } });
  ref.setStyle({}, refStick(0.09, COLORS.ref));

  if (overlay && $("t-overlay").checked) {
    for (const item of overlay) {
      if (view && item.name === view.name) continue;
      v.addModel(item.molblock, "sdf").setStyle({}, { stick: { radius: 0.07, color: COLORS.overlay, opacity: 0.5 } });
    }
  }

  if (view && view.status === "aligned") {
    const refAtoms = ref.selectedAtoms({});
    for (const i of view.reference_core_atoms) {
      ref.setStyle({ serial: refAtoms[i].serial }, refStick(0.13, COLORS.refCore));
    }
    const mol = v.addModel(view.molblock, "sdf");
    mol.setStyle({}, { stick: { radius: 0.15, color: COLORS.rest } });
    const atoms = mol.selectedAtoms({});
    for (const i of view.core_atoms) {
      mol.setStyle({ serial: atoms[i].serial }, { stick: { radius: 0.2, color: COLORS.core } });
    }
  }

  if (!$("t-hydrogens").checked) {
    for (const m of v.getModelList()) m.setStyle({ elem: "H" }, {});
  }
  if (needZoom) { v.zoomTo(); needZoom = false; }
  v.render();
}

async function ensureOverlay() {
  overlay = $("t-overlay").checked ? await api("/api/overlay") : null;
}

// ----------------------------------------------------------------- actions

async function refreshState() {
  state = await api("/api/state");
  renderState();
}

async function select(name) {
  selected = name;
  view = await api(`/api/analogs/${encodeURIComponent(name)}`);
  renderState();
  renderDetail();
  draw();
}

async function loadReference(file) {
  const text = await file.text();
  state = await api("/api/reference", { method: "PUT", body: { name: file.name.replace(/\.[^.]+$/, ""), text } });
  refView = await api("/api/reference");
  view = null; selected = null; overlay = null; needZoom = true;
  say(`Reference ${state.reference.name} loaded.`);
  renderState(); renderDetail(); draw();
}

async function loadAnalogs(text) {
  state = await api("/api/analogs", { method: "PUT", body: { text } });
  view = null; selected = null; overlay = null;
  const failed = state.analogs.filter((a) => a.status === "failed").length;
  say(`${state.analogs.length} analogs read` + (failed ? `, ${failed} with a SMILES that could not be parsed.` : "."), failed > 0);
  renderState(); renderDetail(); draw();
}

async function align() {
  await api("/api/params", {
    method: "PUT",
    body: {
      complete_rings_only: $("p-complete").checked,
      ring_matches_ring_only: $("p-ringring").checked,
      ignore_bond_orders: $("p-anybond").checked,
      timeout: Number($("p-timeout").value) || 0,
      seed: Number($("p-seed").value) || 0,
    },
  });
  state = await api("/api/align", { method: "POST" });
  view = null; selected = null; overlay = null;
  say("Aligning…");
  renderState(); renderDetail(); draw();
  clearInterval(pollTimer);
  pollTimer = setInterval(() => guarded(poll), 400);
}

async function poll() {
  await refreshState();
  if (state.job.running) return;
  clearInterval(pollTimer);
  const ok = state.analogs.filter((a) => a.status === "aligned").length;
  const bad = state.analogs.length - ok;
  say(`${ok} aligned` + (bad ? `, ${bad} failed (see the list).` : "."), bad > 0);
  const first = state.analogs.find((a) => a.status === "aligned") || state.analogs[0];
  if (first) await select(first.name);
}

async function decide(decision) {
  if (!view || view.status !== "aligned") return;
  const next = view.decision === decision ? "undecided" : decision;
  state = await api(`/api/analogs/${encodeURIComponent(view.name)}/decision`, { method: "POST", body: { decision: next } });
  view.decision = next;
  overlay = null;
  await ensureOverlay();
  renderState(); renderDetail(); draw();
}

async function realign() {
  if (!view) return;
  say(`Aligning ${view.name} again…`);
  view = await api(`/api/analogs/${encodeURIComponent(view.name)}/realign`, { method: "POST", body: { seed: Number($("r-seed").value) || 0 } });
  await refreshState();
  await ensureOverlay();
  renderDetail(); draw();
  say(view.status === "aligned" ? `${view.name} aligned again with seed ${view.seed}.` : `${view.name}: ${view.error}`, view.status !== "aligned");
}

function move(step) {
  if (!state || !state.analogs.length) return;
  const names = state.analogs.map((a) => a.name);
  const i = Math.min(Math.max(names.indexOf(selected) + step, 0), names.length - 1);
  guarded(() => select(names[i]));
}

// ------------------------------------------------------------------ wiring

$("reference-file").addEventListener("change", (e) => {
  const file = e.target.files[0];
  if (file) guarded(() => loadReference(file));
  e.target.value = "";
});
$("analogs-file").addEventListener("change", async (e) => {
  const file = e.target.files[0];
  if (file) $("analogs-text").value = await file.text();
  e.target.value = "";
});
$("analogs-load").addEventListener("click", () => guarded(() => loadAnalogs($("analogs-text").value)));
$("align").addEventListener("click", () => guarded(align));
$("accept").addEventListener("click", () => guarded(() => decide("accepted")));
$("reject").addEventListener("click", () => guarded(() => decide("rejected")));
$("realign").addEventListener("click", () => guarded(realign));
$("d-format").addEventListener("change", () => state && updateDownload());
$("download").addEventListener("click", () => guarded(download));
$("t-hydrogens").addEventListener("change", draw);
$("ref2d").open = readSetting("ref2dOpen") === "1";
$("ref2d").addEventListener("toggle", (e) => saveSetting("ref2dOpen", e.target.open ? "1" : "0"));
$("t-color").addEventListener("change", (e) => {
  colorReference = e.target.checked;
  saveSetting("colorReference", colorReference ? "1" : "0");
  draw();
});
$("neutralize").addEventListener("change", (e) => guarded(async () => {
  const wanted = e.target.checked;
  try {
    state = await api("/api/reference/neutralize", { method: "PUT", body: { neutralize: wanted } });
  } catch (err) {
    e.target.checked = !wanted;
    throw err;
  }
  view = null; selected = null; overlay = null;
  if (state.reference) {
    refView = await api("/api/reference");
    const ref = state.reference;
    say(wanted
      ? `Charges removed from the reference: net charge ${ref.loaded_net_charge} to ${ref.net_charge}. Align again.`
      : "Reference charges kept as in the file. Align again.");
  }
  renderState(); renderDetail(); draw();
}));
$("t-overlay").addEventListener("change", () => guarded(async () => { await ensureOverlay(); draw(); }));
$("t-reset").addEventListener("click", () => { needZoom = true; draw(); });
$("reset").addEventListener("click", () => guarded(async () => {
  if (!confirm("Forget the reference, the analogs and all results?")) return;
  state = await api("/api/reset", { method: "POST" });
  refView = null; view = null; selected = null; overlay = null; needZoom = true;
  if (viewer) { viewer.removeAllModels(); viewer.render(); }
  $("analogs-text").value = "";
  syncParams(); renderState(); renderDetail();
  say("Started over.");
}));
window.addEventListener("resize", () => { if (viewer) { viewer.resize(); viewer.render(); } });
document.addEventListener("keydown", (e) => {
  if (e.target.closest("input, textarea, select")) return;
  if (e.key === "ArrowDown") { e.preventDefault(); move(1); }
  else if (e.key === "ArrowUp") { e.preventDefault(); move(-1); }
  else if (e.key === "a") guarded(() => decide("accepted"));
  else if (e.key === "r") guarded(() => decide("rejected"));
});

guarded(async () => {
  await refreshState();
  syncParams();
  if (state.reference) refView = await api("/api/reference");
  const first = state.analogs.find((a) => a.status === "aligned");
  if (first) await select(first.name); else { renderDetail(); draw(); }
  if (state.job.running) pollTimer = setInterval(() => guarded(poll), 400);
});
