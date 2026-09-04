/*
 * webr-cell.js — drop-in runnable R cells + self-checks for the Integrated Bioanalytics book.
 *
 * Turns every <div class="webr-cell"> on a page into a live R cell (edit + Run + output), and every
 * <div class="selfcheck"> into a checked multiple-choice question. All cells on a page share ONE WebR
 * instance and ONE R session (boot once, lazily on the first Run) — so a variable made in one cell is
 * available in the next, exactly like a notebook. R runs entirely in the browser; nothing is installed.
 *
 * The WebR half is NOT implemented here. It is `WebRConsole`, shared with the escape rooms and the
 * sandbox (see the import below). This file used to carry a copy-pasted second implementation of boot +
 * run; the two drifted (the escape-room copy handled a failed boot, this one did not), so they were
 * consolidated on 2026-09-04. Keep this file to CELL UI only — anything about running R belongs in
 * `escape_rooms/shared/webr-console.js` so all three surfaces get it at once.
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
 *
 * Students get `view(x)` for free in every cell — the RStudio-style table viewer, defined by the shared
 * console. It works on derived data too: view(algae_data %>% filter(...)).
 *
 * CROSS-REPO IMPORT: the shared console lives in the escape_rooms repo and is fetched by its public
 * URL, exactly as sandbox.html already does. Its `?v=` token must be bumped by hand in lockstep with
 * escape_rooms/shared/ — nothing validates it from this side. See this dir's AGENTS.md.
 */
import { WebRConsole } from "/escape_rooms/shared/webr-console.js?v=88";

let rconsole = null;
function session() {
  if (!rconsole) {
    const cfg = window.WEBR_CELL_CONFIG || {};
    rconsole = new WebRConsole(
      { packages: cfg.packages, datasets: cfg.datasets, setup: cfg.setup },
      {}                                   // no single output/status: each cell supplies its own
    );
    // The shared console's "R is ready" chatter belongs in a standalone console, not under every cell
    // in a chapter — swap it for a blank line once the boot succeeds.
    const origSet = rconsole.setStatus.bind(rconsole);
    rconsole.setStatus = m => origSet(rconsole.ready ? "" : m);
  }
  return rconsole;
}

function initCells() {
  document.querySelectorAll(".webr-cell").forEach(cell => {
    if (cell.dataset.wired) return; cell.dataset.wired = "1";
    const ta = cell.querySelector("textarea");
    let bar = cell.querySelector(".webr-bar"); if (!bar) { bar = document.createElement("div"); bar.className = "webr-bar"; cell.appendChild(bar); }
    let btn = bar.querySelector("button.webr-run"); if (!btn) { btn = document.createElement("button"); btn.className = "webr-run"; btn.textContent = "▶ Run"; bar.appendChild(btn); }
    let stat = bar.querySelector(".webr-status"); if (!stat) { stat = document.createElement("span"); stat.className = "webr-status"; bar.appendChild(stat); }
    let out = cell.querySelector(".webr-output"); if (!out) { out = document.createElement("div"); out.className = "webr-output"; cell.appendChild(out); }
    // One session, many cells: register this cell's status line so the boot message shows wherever
    // the student actually clicked.
    session().addStatusEl(stat);
    const run = async () => {
      btn.disabled = true; const t = btn.textContent; btn.textContent = "running…";
      try { await session().runFrom(ta, out); } finally { btn.disabled = false; btn.textContent = t; if (session().ready) stat.textContent = ""; }
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
// re-callable if a page injects cells later; `session` exposed for console debugging in class
window.WebRCell = { init, session };
