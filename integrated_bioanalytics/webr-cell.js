/*
 * webr-cell.js — drop-in runnable R cells + self-checks for the Integrated Bioanalytics book.
 *
 * Turns every <div class="webr-cell"> on a page into a live R cell (edit + Run + output), and every
 * <div class="selfcheck"> into a checked multiple-choice question. All cells on a page share ONE WebR
 * instance and ONE R session (boot once, lazily on the first Run) — so a variable made in one cell is
 * available in the next, exactly like a notebook. R runs entirely in the browser; nothing is installed.
 *
 * Page config (packages/datasets/setup) comes from a global set before this script loads:
 *   <script>window.WEBR_CELL_CONFIG = { packages:[...], datasets:[{name,url}], setup:"..." };</script>
 *
 * Authoring (raw HTML in an .Rmd, or any static page):
 *   <div class="webr-cell"><textarea class="webr-code">head(algae_data)</textarea></div>
 *   <div class="selfcheck" data-answer="1" data-explain="geom_col draws bars.">
 *     <p class="scq">Which geom draws bars?</p>
 *     <label><input type="radio"> geom_point()</label>
 *     <label><input type="radio"> geom_col()</label>
 *   </div>
 */
import { WebR } from "https://webr.r-wasm.org/latest/webr.mjs";

let webR = null, ready = false, booting = null;
const statusEls = new Set();                       // every cell's status line — the boot message shows on all
const setBootStatus = m => statusEls.forEach(s => { s.textContent = m; });

async function boot() {
  if (ready) return;
  if (booting) return booting;
  booting = (async () => {
    const cfg = window.WEBR_CELL_CONFIG || {};
    setBootStatus("Starting R in your browser… (first time ~20–40s)");
    webR = new WebR({ interactive: false });
    await webR.init();
    const pkgs = cfg.packages || [];
    if (pkgs.length) { setBootStatus("Installing R packages: " + pkgs.join(", ") + " …"); await webR.installPackages(pkgs, { quiet: true }); }
    for (const ds of (cfg.datasets || [])) {
      setBootStatus("Loading data: " + ds.name + " …");
      const resp = await fetch(ds.url);
      if (!resp.ok) throw new Error("could not fetch " + ds.url);
      const bytes = new Uint8Array(await resp.arrayBuffer());
      const p = "/home/web_user/" + ds.name + ".csv";
      await webR.FS.writeFile(p, bytes);
      await webR.evalRVoid(`${ds.name} <- readr::read_csv("${p}", show_col_types = FALSE)`);
    }
    if (cfg.setup) { setBootStatus("Preparing session…"); await webR.evalRVoid(cfg.setup); }
    ready = true; setBootStatus("");
  })();
  return booting;
}

async function runCode(code, outEl) {
  outEl.innerHTML = "";
  if (!ready) await boot();
  if (!ready) return;
  const shelter = await new webR.Shelter();
  try {
    const result = await shelter.captureR(code, { withAutoprint: true, captureStreams: true, captureGraphics: { width: 720, height: 460 } });
    const text = result.output.filter(o => o.type === "stdout" || o.type === "stderr").map(o => o.data).join("\n");
    if (text.trim().length) { const pre = document.createElement("pre"); pre.className = "webr-out"; pre.textContent = text; outEl.appendChild(pre); }
    for (const img of (result.images || [])) {
      const c = document.createElement("canvas"); c.width = img.width; c.height = img.height; c.className = "webr-plot";
      c.getContext("2d").drawImage(img, 0, 0); outEl.appendChild(c);
    }
    if (!text.trim().length && !(result.images || []).length) { const pre = document.createElement("pre"); pre.className = "webr-out muted"; pre.textContent = "(no output)"; outEl.appendChild(pre); }
  } catch (err) {
    const pre = document.createElement("pre"); pre.className = "webr-out err"; pre.textContent = "Error: " + (err && err.message ? err.message : err); outEl.appendChild(pre);
  } finally { shelter.purge(); }
}

function initCells() {
  document.querySelectorAll(".webr-cell").forEach(cell => {
    if (cell.dataset.wired) return; cell.dataset.wired = "1";
    const ta = cell.querySelector("textarea");
    let bar = cell.querySelector(".webr-bar"); if (!bar) { bar = document.createElement("div"); bar.className = "webr-bar"; cell.appendChild(bar); }
    let btn = bar.querySelector("button.webr-run"); if (!btn) { btn = document.createElement("button"); btn.className = "webr-run"; btn.textContent = "▶ Run"; bar.appendChild(btn); }
    let stat = bar.querySelector(".webr-status"); if (!stat) { stat = document.createElement("span"); stat.className = "webr-status"; bar.appendChild(stat); }
    let out = cell.querySelector(".webr-output"); if (!out) { out = document.createElement("div"); out.className = "webr-output"; cell.appendChild(out); }
    statusEls.add(stat);
    const run = async () => {
      btn.disabled = true; const t = btn.textContent; btn.textContent = "running…";
      try { await runCode(ta.value, out); } finally { btn.disabled = false; btn.textContent = t; if (ready) stat.textContent = ""; }
    };
    btn.onclick = run;
    ta.addEventListener("keydown", e => { if ((e.metaKey || e.ctrlKey) && e.key === "Enter") { e.preventDefault(); run(); } });
  });
}

function initSelfChecks() {
  document.querySelectorAll(".selfcheck").forEach((sc, i) => {
    if (sc.dataset.wired) return; sc.dataset.wired = "1";
    const inputs = [...sc.querySelectorAll('input[type="radio"]')];
    inputs.forEach(r => { r.name = sc.dataset.name || ("selfcheck-" + i); });
    let btn = sc.querySelector("button.sccheck"); if (!btn) { btn = document.createElement("button"); btn.className = "sccheck"; btn.textContent = "Check"; sc.appendChild(btn); }
    let fb = sc.querySelector(".scfb"); if (!fb) { fb = document.createElement("div"); fb.className = "scfb"; sc.appendChild(fb); }
    btn.onclick = () => {
      const pick = inputs.findIndex(r => r.checked);
      if (pick < 0) { fb.className = "scfb"; fb.textContent = "Pick an answer first."; return; }
      const correct = pick === Number(sc.dataset.answer);
      fb.className = "scfb " + (correct ? "ok" : "no");
      fb.textContent = (correct ? "✓ Correct. " : "✗ Not quite. ") + (sc.dataset.explain || "");
    };
  });
}

function init() { initCells(); initSelfChecks(); }
if (document.readyState !== "loading") init(); else document.addEventListener("DOMContentLoaded", init);
window.WebRCell = { init, run: runCode };   // re-callable if a page injects cells later
