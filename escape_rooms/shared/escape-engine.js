/*
 * escape-engine.js — drives a two-screen escape-room scenario.
 *
 * Reads window.SCENARIO (see alaska/scenario.js), wires the story screen to
 * the workroom screen, runs the multiple-choice steps one at a time with
 * feedback and attempt-tracking, and produces a Canvas submission code at the
 * end via window.EscapeCodec.
 */
import { WebRConsole } from "./webr-console.js";

const SECRET = "chem5725-noatak-2026"; // change per course; obfuscation only
const VERSION = 1;
const PASTE_WARN_CHARS = 40;

const S = window.SCENARIO;
const state = { studentId: "", steps: [], current: 0 };

document.addEventListener("DOMContentLoaded", () => {
  buildScreen1();
});

/* ---------- Screen 1: the story ---------- */
function buildScreen1() {
  document.getElementById("scenario-title").textContent = S.title;
  const img = document.getElementById("scene-image");
  img.src = S.screen1.image;
  img.onerror = () => { img.style.display = "none"; };
  document.getElementById("scene-heading").textContent = S.screen1.title;
  document.getElementById("scene-story").innerHTML = S.screen1.story;
  const enterBtn = document.getElementById("enter-btn");
  enterBtn.textContent = S.screen1.enterLabel || "Enter →";

  const idInput = document.getElementById("student-id");
  enterBtn.addEventListener("click", () => {
    const id = idInput.value.trim();
    if (!id) {
      document.getElementById("id-error").textContent =
        "Please enter your x500 first — it ties your code to you.";
      return;
    }
    state.studentId = id;
    startMusic(); // triggered by this click, so autoplay is allowed
    enterWorkroom();
  });
}

/* ---------- Ambience ---------- */
function startMusic() {
  if (!S.music) return;
  const audio = document.getElementById("ambience");
  audio.src = S.music;
  audio.volume = typeof S.musicVolume === "number" ? S.musicVolume : 0.35;
  audio.loop = true;
  audio.play().catch(() => { /* some browsers still block; the toggle recovers it */ });

  const toggle = document.getElementById("music-toggle");
  toggle.classList.add("show");
  toggle.addEventListener("click", () => {
    if (audio.paused) { audio.play(); toggle.textContent = "♪ on"; }
    else { audio.pause(); toggle.textContent = "♪ off"; }
  });
}

/* ---------- Screen 2: the workroom ---------- */
async function enterWorkroom() {
  document.getElementById("screen1").classList.remove("active");
  document.getElementById("screen2").classList.add("active");
  window.scrollTo(0, 0);

  if (S.keepBackground && S.screen1.image) {
    const bg = document.getElementById("screen2-bg");
    bg.style.backgroundImage = "url('" + S.screen1.image + "')";
    bg.classList.add("show");
  }

  document.getElementById("briefing").innerHTML = S.briefing;

  const input = document.getElementById("code-input");
  input.value = S.starterCode || "";
  wirePasteWarning(input);

  const consoleUI = {
    status: document.getElementById("webr-status"),
    output: document.getElementById("webr-output"),
  };
  const rconsole = new WebRConsole(
    { packages: S.packages, datasets: S.datasets, setup: S.setup },
    consoleUI
  );

  const runBtn = document.getElementById("run-btn");
  runBtn.disabled = true;
  rconsole.init()
    .then(() => { runBtn.disabled = false; })
    .catch((e) => { consoleUI.status.textContent = "R failed to start: " + e.message; });

  runBtn.addEventListener("click", () => rconsole.run(input.value));

  renderStep();
}

function wirePasteWarning(input) {
  const banner = document.getElementById("paste-banner");
  input.addEventListener("paste", (e) => {
    const text = (e.clipboardData || window.clipboardData).getData("text");
    if (text && text.length > PASTE_WARN_CHARS) {
      banner.classList.add("show");
      clearTimeout(banner._t);
      banner._t = setTimeout(() => banner.classList.remove("show"), 8000);
    }
    // paste is allowed — we only nudge.
  });
  document.getElementById("paste-dismiss").addEventListener("click", () => {
    banner.classList.remove("show");
  });
}

/* ---------- Steps ---------- */
function renderStep() {
  const container = document.getElementById("steps");
  const idx = state.current;

  if (idx >= S.steps.length) {
    finish();
    return;
  }
  const step = S.steps[idx];
  const maxAttempts = step.maxAttempts || 4;
  let attempts = 0;

  const card = document.createElement("div");
  card.className = "step-card";
  card.innerHTML = `
    <div class="step-number">Question ${idx + 1} of ${S.steps.length}</div>
    <div class="step-prompt">${step.prompt}</div>
    <div class="step-options"></div>
    <div class="step-feedback"></div>
    <button class="submit-btn" disabled>Submit answer</button>
  `;
  container.appendChild(card);
  card.scrollIntoView({ behavior: "smooth", block: "start" });

  const optionsEl = card.querySelector(".step-options");
  const feedbackEl = card.querySelector(".step-feedback");
  const submitBtn = card.querySelector(".submit-btn");
  let selected = -1;

  step.options.forEach((opt, i) => {
    const label = document.createElement("label");
    label.className = "option";
    label.innerHTML = `<input type="radio" name="step${idx}" value="${i}"> <span>${opt}</span>`;
    label.querySelector("input").addEventListener("change", () => {
      selected = i;
      submitBtn.disabled = false;
    });
    optionsEl.appendChild(label);
  });

  function lockAndAdvance(answerIndex) {
    optionsEl.querySelectorAll("input").forEach((el) => (el.disabled = true));
    submitBtn.disabled = true;
    submitBtn.style.display = "none";
    state.steps.push({ answer: answerIndex, attempts: attempts });
    state.current += 1;
    setTimeout(renderStep, 700);
  }

  submitBtn.addEventListener("click", () => {
    if (selected < 0) return;
    attempts += 1;
    if (selected === step.correct) {
      feedbackEl.className = "step-feedback ok";
      feedbackEl.innerHTML = step.feedback.correct;
      lockAndAdvance(selected);
    } else if (attempts >= maxAttempts) {
      feedbackEl.className = "step-feedback out";
      const reveal = step.feedback.reveal || "Moving on.";
      feedbackEl.innerHTML = reveal;
      lockAndAdvance(selected); // records their last (wrong) pick + attempts
    } else {
      feedbackEl.className = "step-feedback no";
      const hints = step.feedback.wrong || [];
      const hint = hints[Math.min(attempts - 1, hints.length - 1)] ||
        "Not quite — take another look at the data and try again.";
      feedbackEl.innerHTML =
        hint + ` <span class="attempts-left">(attempt ${attempts} of ${maxAttempts})</span>`;
    }
  });
}

/* ---------- Finish: build the code ---------- */
function finish() {
  const code = window.EscapeCodec.encode({
    version: VERSION,
    scenarioId: S.id,
    steps: state.steps,
    studentId: state.studentId,
    secret: SECRET,
  });

  const panel = document.getElementById("finish");
  panel.classList.add("show");
  document.getElementById("final-code").textContent = code;
  document.getElementById("finish-body").innerHTML = S.finishMessage || "";
  panel.scrollIntoView({ behavior: "smooth" });

  document.getElementById("copy-code").addEventListener("click", () => {
    navigator.clipboard.writeText(code).then(() => {
      document.getElementById("copy-code").textContent = "Copied ✓";
    });
  });
}
