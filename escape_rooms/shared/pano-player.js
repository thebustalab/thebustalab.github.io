/*
 * pano-player.js — the shared pseudo-360 room player, chapter-agnostic.
 *
 * A chapter shell is tiny: link the Pannellum CSS + ../shared/pano-player.css,
 * load the Pannellum script + the chapter data (a plain script that sets
 * window.CHAPTER), then `<script type="module" src="../shared/pano-player.js">`.
 * This module injects its own DOM and drives the chapter.
 *
 * A chapter is a chain of pseudo-360 rooms. Each room has a `panorama` (+ optional
 * `panoramaOpen`), a `wrap` {haov,vaov,hfov,vOffset,pitch}, and `hotspots` with
 * boxes [x0,y0,x1,y1] (fractions of the flat scene) + per-type content:
 *   puzzle → { starterCode, question:{prompt,options,correct,maxAttempts,feedback} }
 *   clue   → { body }
 *   door   → (box only; swaps to panoramaOpen on solve and advances)
 * Nav is ‹ › arrows (yaw only). WebR boots once and mounts into the puzzle modal.
 * The submission codec is not wired here yet (add shared/codec.js when needed).
 *
 * window.CHAPTER: { title, subtitle, story, enterLabel?, done?:{title,body},
 *   packages, datasets:[{name,url}], setup, rooms:[…] }
 */
import { WebRConsole } from "./webr-console.js";

const CH = window.CHAPTER;
if (!CH) {
  document.body.innerHTML = '<p style="color:#ff9b9b;font:15px system-ui;padding:24px">' +
    "pano-player: window.CHAPTER is not defined — load the chapter data script before this module.</p>";
  throw new Error("pano-player: window.CHAPTER is not defined");
}
document.title = CH.title || "Escape room";

// ---- inject the player DOM ----
const root = document.createElement("div");
root.innerHTML = `
  <section id="screen1" class="screen active">
    <div class="intro">
      <h1 id="s1title"></h1>
      <div class="sub" id="s1sub"></div>
      <p id="s1story"></p>
      <div class="x500">
        <input id="x500" placeholder="your x500 (e.g. bust0037)" autocomplete="off" />
        <button id="enter"></button>
      </div>
    </div>
  </section>
  <section id="screen2" class="screen">
    <div id="pano"></div>
    <div id="hud"><span id="hudroom"></span> · <b>use ‹ › to look around</b> · click the glowing markers</div>
    <button class="arrow l" id="prev">‹</button>
    <button class="arrow r" id="next">›</button>
    <div id="toast"></div>
    <div id="modal">
      <div class="mbox">
        <button class="ghost mback" id="mback">← Back to the room</button>
        <div class="mtitle" id="mtitle"></div>
        <div id="mbody"></div>
      </div>
    </div>
    <div id="done">
      <div class="card">
        <h2 id="doneTitle"></h2>
        <p id="doneBody"></p>
        <button class="ghost" id="replay">Play again</button>
      </div>
    </div>
  </section>
  <div id="console-holder">
    <div id="console-block" class="console">
      <div id="webr-status">R console</div>
      <textarea id="code-input" spellcheck="false"></textarea>
      <div><button id="run-btn" disabled>▶ Run</button></div>
      <div id="webr-output"></div>
    </div>
  </div>`;
document.body.appendChild(root);

const $ = s => document.querySelector(s);
let viewer = null, rconsole = null, roomIdx = 0, room = null, solved = false;

// screen 1
$("#s1title").textContent = CH.title || "";
$("#s1sub").textContent = CH.subtitle || "";
$("#s1story").textContent = CH.story || "";
$("#enter").textContent = CH.enterLabel || "Begin →";
$("#doneTitle").textContent = (CH.done && CH.done.title) || ((CH.title || "Chapter") + " — complete");
$("#doneBody").textContent = (CH.done && CH.done.body) || "Nice work — you've finished this chapter.";
$("#replay").onclick = () => location.reload();
$("#x500").addEventListener("keydown", e => { if (e.key === "Enter") $("#enter").click(); });
$("#enter").onclick = () => {
  const id = $("#x500").value.trim();
  if (!/\S/.test(id)) { $("#x500").focus(); return; }
  window.__x500 = id;
  $("#screen1").classList.remove("active");
  $("#screen2").classList.add("active");
  bootConsole();
  startRoom(0);
};

// box [x0,y0,x1,y1] fractions -> yaw/pitch across the room's wrap coverage
function boxToYP(box, c) {
  const cx = (box[0] + box[2]) / 2, cy = (box[1] + box[3]) / 2;
  return { yaw: cx * c.haov - c.haov / 2, pitch: c.vaov / 2 - cy * c.vaov + (c.vOffset || 0) };
}

function startRoom(i) {
  roomIdx = i; room = CH.rooms[i]; solved = false;
  $("#hudroom").textContent = room.title || "";
  buildViewer(room.panorama);
}

function buildViewer(img) {
  if (viewer) { try { viewer.destroy(); } catch (e) {} }
  const c = room.wrap;
  const p = c.pitch || 0, f = c.hfov || 110;
  // Build hotspots up front and pass them in the config (the reliable path on a
  // static, non-draggable viewer — addHotSpot-after-load leaves them unpositioned).
  const hotSpots = (room.hotspots || []).map(h => {
    const { yaw, pitch } = boxToYP(h.box, c);
    let cssClass = "hsmark " + h.type;
    if (h.type === "door") cssClass += solved ? " open" : " locked";
    if (h.type === "puzzle" && solved) cssClass += " done";
    return { id: h.id, yaw, pitch, cssClass, clickHandlerFunc: onHotspot, clickHandlerArgs: h };
  });
  // No up/down / no zoom is enforced by disabling drag + zoom and only moving yaw
  // via the arrows — NOT by min===max pitch/hfov locks (those break hotspot projection).
  viewer = pannellum.viewer("pano", {
    type: "equirectangular", panorama: img,
    haov: c.haov, vaov: c.vaov, vOffset: c.vOffset || 0,
    hfov: f, pitch: p, yaw: 0,
    autoLoad: true, showControls: false, autoRotate: 0,
    draggable: false, mouseZoom: false, doubleClickZoom: false,
    keyboardZoom: false, disableKeyboardCtrl: true,
    hotSpots: hotSpots,
    backgroundColor: [0.02, 0.05, 0.09],
  });
}
const TURN = 45;
$("#prev").onclick = () => viewer && viewer.setYaw(viewer.getYaw() - TURN, 600);
$("#next").onclick = () => viewer && viewer.setYaw(viewer.getYaw() + TURN, 600);

function onHotspot(evt, h) {   // Pannellum calls clickHandlerFunc(event, clickHandlerArgs)
  try {
    if (h.type === "clue") return openClue(h);
    if (h.type === "puzzle") return openPuzzle(h);
    if (h.type === "door") return solved ? goThrough() : toast("The door won't budge — solve the puzzle first.");
  } catch (e) { console.error("hotspot handler error", e); }  // Pannellum swallows handler throws
}

// ---- modal helpers ----
function openModal(title, node) {
  $("#mtitle").textContent = title;
  const body = $("#mbody"); body.innerHTML = ""; body.appendChild(node);
  $("#modal").classList.add("open");
}
$("#mback").onclick = closeModal;
function closeModal() { unmountConsole(); $("#modal").classList.remove("open"); }
function unmountConsole() { $("#console-holder").appendChild($("#console-block")); }

function openClue(h) {
  const d = document.createElement("div");
  d.innerHTML = `<p>${h.body || ""}</p>`;
  openModal(h.label || "Clue", d);
}

function openPuzzle(h) {
  // Set code/output while the console block is still in the document — querying
  // #code-input after moving it into a detached div would return null.
  $("#code-input").value = h.starterCode || "";
  $("#webr-output").innerHTML = "";
  const cb = $("#console-block");
  const grid = document.createElement("div"); grid.className = "qa";
  const left = document.createElement("div");   // console pane
  const right = document.createElement("div");   // question pane
  left.appendChild(cb);                          // move the live console in
  right.appendChild(buildQuestion(h.question, () => { closeModal(); solveRoom(); }));
  grid.appendChild(left); grid.appendChild(right);
  openModal(h.label || "Puzzle", grid);
}

// multiple-choice card — gate is the *product* of running the analysis
function buildQuestion(q, onSolved) {
  const maxA = q.maxAttempts || 4; let attempts = 0, sel = -1;
  const card = document.createElement("div"); card.className = "qcard";
  const grp = "q" + Math.round(performance.now());
  card.innerHTML =
    `<div class="qprompt">${q.prompt}</div><div class="qopts"></div>
     <div class="qfeedback"></div><button class="qsubmit" disabled>Submit answer</button>`;
  const opts = card.querySelector(".qopts"), fb = card.querySelector(".qfeedback"), sub = card.querySelector(".qsubmit");
  q.options.forEach((o, i) => {
    const l = document.createElement("label"); l.className = "qopt";
    l.innerHTML = `<input type="radio" name="${grp}"><span>${o}</span>`;
    l.querySelector("input").addEventListener("change", () => { sel = i; sub.disabled = false; });
    opts.appendChild(l);
  });
  sub.addEventListener("click", () => {
    if (sel < 0) return;
    attempts++;
    if (sel === q.correct) {
      fb.className = "qfeedback ok"; fb.innerHTML = q.feedback.correct;
      sub.disabled = true; opts.style.pointerEvents = "none";
      setTimeout(onSolved, 900);
    } else if (attempts >= maxA) {
      fb.className = "qfeedback out"; fb.innerHTML = q.feedback.reveal || "Out of attempts.";
      sub.disabled = true; opts.style.pointerEvents = "none";
    } else {
      fb.className = "qfeedback no";
      const hints = q.feedback.wrong || [];
      fb.innerHTML = (hints[Math.min(attempts - 1, hints.length - 1)] || "Not quite.") +
        ` <span class="attempts">(attempt ${attempts} of ${maxA})</span>`;
    }
  });
  return card;
}

// solving swaps to the open-door panorama; hotspots re-place with the door live
function solveRoom() {
  solved = true;
  buildViewer(room.panoramaOpen || room.panorama);
  toast("The door is open. Look for the way through.");
}

function goThrough() {
  if (roomIdx + 1 < CH.rooms.length) startRoom(roomIdx + 1);
  else $("#done").classList.add("open");
}

let toastT = null;
function toast(msg) {
  const t = $("#toast"); t.textContent = msg; t.classList.add("show");
  clearTimeout(toastT); toastT = setTimeout(() => t.classList.remove("show"), 2600);
}

// ---- WebR (boot once) ----
function bootConsole() {
  if (rconsole) return;
  rconsole = new WebRConsole(
    { packages: CH.packages, datasets: CH.datasets, setup: CH.setup },
    { status: $("#webr-status"), output: $("#webr-output") }
  );
  const runBtn = $("#run-btn");
  rconsole.init().then(() => { runBtn.disabled = false; })
    .catch(e => { $("#webr-status").textContent = "R failed to start: " + (e.message || e); });
  runBtn.addEventListener("click", () => rconsole.run($("#code-input").value));
}
