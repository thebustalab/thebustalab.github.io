#!/usr/bin/env node
/*
 * Regression tests for buildLeaderboard() in apps_script.gs.js.
 *
 * Failure mode these pin (2026-08-19): the in-tool board and the newsletter's
 * per-term prize board disagreed. This board scored EVERY newsletter form
 * submission at 1 point and counted ALL TIME; the newsletter weighted by type
 * (data update 10, abstract 5) and scoped to the current term. A student with
 * two data updates saw 2 here and 20 in the newsletter, which reads as being
 * short-changed. The two implementations must now agree row for row — see
 * AGENTS.md → Leaderboard, and newsletter/AGENTS.md → Leaderboard Scoring.
 *
 * ALSO pinned here — two HIGH security findings fixed 2026-09-08:
 *
 *  1. STORED XSS in the leaderboard (index.html). The board rendered labeller
 *     names into `tbody.innerHTML` via a template literal. Names are not
 *     trusted input: the ?k= link is shared lab-wide and anyone can POST the
 *     Apps Script endpoint directly, bypassing the roster dropdown entirely —
 *     so a `labeller` of `<img src=x onerror=…>` was stored in the Sheet and
 *     then executed in every labeller's browser on every board refresh. Fixed
 *     by rendering each cell with textContent (renderLeaderboardRows), with a
 *     single escapeHtml() for the one remaining HTML-string interpolation.
 *
 *  2. FORMULA INJECTION in apps_script.gs.js. Free-text submission fields were
 *     written straight into cells, so a value starting with = + - @ tab or CR
 *     was stored as a LIVE spreadsheet formula (=IMPORTRANGE/=HYPERLINK can
 *     exfiltrate the sheet). Fixed by sheetSafe(), applied to every cell at the
 *     appendRows() choke point so no handler can forget it.
 *
 * The XSS cases below load the two functions out of index.html by name and run
 * them against a minimal DOM stub whose innerHTML setter THROWS — so a future
 * edit that goes back to string-built markup fails here rather than shipping.
 *
 * Run: node test_leaderboard.js
 */
const fs = require("fs");
const path = require("path");
const assert = require("assert");

// ─── Fixture sheets, keyed by tab name (the stub ignores spreadsheet ids) ───
// Dates are derived from "now" so the tests don't rot when the term rolls over.
const now = new Date();
const m = now.getMonth() + 1;
const termStartMonth = m <= 4 ? 1 : m <= 8 ? 5 : 9;
const inTerm = new Date(now.getFullYear(), termStartMonth - 1, 2, 12, 0, 0).toISOString();
const beforeTerm = new Date(now.getFullYear(), termStartMonth - 2, 15, 12, 0, 0).toISOString();

const SHEETS = {
  labels: [
    ["timestamp", "labeller", "abstract_id", "abstract_title"],
    // 7 triples on one abstract — capped at MAX_ROWS_PER_ABSTRACT (5) × 1 pt.
    ...Array.from({ length: 7 }, () => [inTerm, "Ada", "abs1", "t"]),
    [inTerm, "Ada", "abs2", "t"],                 // +1
    [beforeTerm, "Ada", "abs3", "t"],             // last term — excluded
    ["", "Ada", "abs4", "t"],                     // undated — excluded
  ],
  enzyme_labels: [
    ["timestamp", "labeller", "abstract_id", "abstract_title"],
    [inTerm, "Bo", "e1", "t"],                    // 3 pts
    [inTerm, "Bo", "e1", "t"],                    // same abstract, under cap: +3
    [beforeTerm, "Bo", "e2", "t"],                // excluded
  ],
  escape_feedback: [
    ["timestamp", "labeller", "scenario", "scenario_title"],
    [inTerm, "Cy", "data_vis/alaska", "Signal in the Cold"],   // 10
    [inTerm, "Cy", "data_vis/alaska", "Signal in the Cold"],   // re-submission: capped
    [inTerm, "Cy", "data_vis/hawaii", "Saltwater Intrusion"],  // +10
    [beforeTerm, "Cy", "data_vis2/hospital", "Vital Signs"],   // excluded
  ],
  "Form Responses 1": [
    ["Timestamp", "What sort of submission is this?", "Your name"],
    [inTerm, "Data Update Submission", "Di"],     // 10, not 1
    [inTerm, "Abstract Submission", "Di"],        // 5
    [inTerm, "Photo Submission", "Di"],           // 1
    [inTerm, "Reminder Submission", "Di"],        // 1
    [inTerm, "Photo Submission", ""],             // no name — ignored
    [beforeTerm, "Data Update Submission", "Di"], // excluded
  ],
};

const WRITES = [];   // every appendRows() write, captured by the stub below

global.SpreadsheetApp = {
  openById() {
    return {
      getSheetByName(name) {
        const v = SHEETS[name];
        if (!v) return null;
        const width = Math.max(...v.map((r) => r.length));
        const grid = v.map((r) => {
          const c = r.slice();
          while (c.length < width) c.push("");
          return c;
        });
        return {
          getLastRow: () => grid.length,
          getLastColumn: () => width,
          getRange: (r, c, nr, nc) => ({
            getValues: () => grid.slice(r - 1, r - 1 + nr).map((row) => row.slice(c - 1, c - 1 + nc)),
            // Records what a handler would write, for the sanitisation cases.
            setValues: (vals) => { WRITES.push({ tab: name, rows: vals }); },
          }),
        };
      },
    };
  },
};

global.ContentService = {
  MimeType: { JSON: "application/json" },
  createTextOutput: (text) => ({ getContent: () => text, setMimeType() { return this; } }),
};

// Load the Apps Script source into this scope (it has no module system).
const src = fs.readFileSync(path.join(__dirname, "apps_script.gs.js"), "utf8");
eval(src);

const board = buildLeaderboard();
const points = Object.fromEntries(board.map((r) => [r.name, r.count]));

// Pathogen: 5 (capped) + 1. Out-of-term and undated rows must not count.
assert.strictEqual(points.Ada, 6, "pathogen cap / term window");
// Enzyme: 2 rows on one abstract, under the cap, × ENZYME_POINTS_PER_ROW.
assert.strictEqual(points.Bo, 6, "enzyme weight");
// Escape: 10 per scenario, one per scenario per person, this term only.
assert.strictEqual(points.Cy, 20, "escape weight + per-scenario cap");
// Form: the regression itself — weighted by type, not 1 apiece, this term only.
assert.strictEqual(points.Di, 17, "form submissions weighted by type");

// Ordering: descending by count, ties alphabetical.
assert.deepStrictEqual(board.map((r) => r.name), ["Cy", "Di", "Ada", "Bo"], "sort order");

// The weights themselves, so a silent retune on one side trips a test here.
assert.strictEqual(submissionPoints("Data Update Submission"), 10);
assert.strictEqual(submissionPoints("Abstract Submission"), 5);
assert.strictEqual(submissionPoints("Photo Submission"), 1);
assert.strictEqual(submissionPoints(""), 1);
assert.strictEqual(ENZYME_POINTS_PER_ROW, 3);
assert.strictEqual(ESCAPE_POINTS_PER_ROW, 10);
assert.strictEqual(MAX_ROWS_PER_ABSTRACT, 5);

// Term split: Spring Jan–Apr, Summer May–Aug, Fall Sep–Dec, [start, end).
const summer = currentSeasonWindow(new Date(2026, 6, 15));
assert.strictEqual(summer.start.getTime(), new Date(2026, 4, 1).getTime());
assert.strictEqual(summer.end.getTime(), new Date(2026, 8, 1).getTime());
const fall = currentSeasonWindow(new Date(2026, 11, 31));
assert.strictEqual(fall.end.getTime(), new Date(2027, 0, 1).getTime());

// Timestamps arrive as ISO strings (labelling tabs) or Date cells (form sheet).
const win = currentSeasonWindow(new Date(2026, 6, 15));
assert.ok(inWindow("2026-07-01T00:00:00.000Z", win), "ISO string in window");
assert.ok(inWindow(new Date(2026, 6, 1), win), "Date cell in window");
assert.ok(!inWindow("2026-04-30T00:00:00.000Z", win), "previous term excluded");
assert.ok(!inWindow("", win), "undated row excluded");
assert.ok(!inWindow("not a date", win), "unparseable row excluded");

// ───────────────────────────────────────────────────────────────────────────
// Fix 1 (2026-09-08) — STORED XSS: labeller names must render inert.
// The two client-side renderers are lifted out of index.html by name and run
// against a DOM stub that THROWS if innerHTML is touched.
// ───────────────────────────────────────────────────────────────────────────
const indexSrc = fs.readFileSync(path.join(__dirname, "index.html"), "utf8");

function extractFunction(src, name) {
  const start = src.indexOf("function " + name + "(");
  assert.ok(start !== -1, name + "() not found in index.html");
  let depth = 0;
  for (let j = src.indexOf("{", start); j < src.length; j++) {
    if (src[j] === "{") depth++;
    else if (src[j] === "}" && --depth === 0) return src.slice(start, j + 1);
  }
  throw new Error("unbalanced braces reading " + name + "() from index.html");
}

function stubEl(tag) {
  return {
    tagName: tag,
    className: "",
    value: "",
    children: [],
    _text: "",
    get textContent() { return this._text; },
    set textContent(v) { this._text = String(v); this.children.length = 0; },
    appendChild(child) { this.children.push(child); return child; },
    set innerHTML(v) {
      throw new Error("innerHTML used to render untrusted values in <" + tag + ">");
    },
  };
}

const documentStub = { createElement: (t) => stubEl(t) };
const loadFn = (name) =>
  new Function("document", extractFunction(indexSrc, name) + "; return " + name + ";")(documentStub);

const renderLeaderboardRows = loadFn("renderLeaderboardRows");
const escapeHtml = loadFn("escapeHtml");
const optionEl = loadFn("optionEl");

// The renderer itself must contain no innerHTML at all.
assert.ok(
  !extractFunction(indexSrc, "renderLeaderboardRows").includes("innerHTML"),
  "renderLeaderboardRows must not use innerHTML"
);

// A labeller name straight off an unauthenticated POST, carrying a payload.
const XSS_NAME = '<img src=x onerror="alert(1)"><script>alert(2)</script>';
const tbody = stubEl("tbody");
renderLeaderboardRows(tbody, [[XSS_NAME, 5], ["Ada", 6]], "Ada");
assert.strictEqual(tbody.children.length, 2, "one row per scorer");
const nameCell = tbody.children[0].children[0];
assert.strictEqual(nameCell.textContent, XSS_NAME, "name stored verbatim as TEXT");
assert.strictEqual(nameCell.children.length, 0, "markup in a name is never parsed into nodes");
assert.strictEqual(tbody.children[0].children[1].textContent, "5", "count rendered as text");
assert.strictEqual(tbody.children[1].className, "me", "current labeller row still highlighted");

// The one remaining HTML-string interpolation (the article link) is escaped.
assert.strictEqual(escapeHtml("<script>alert(1)</script>"),
  "&lt;script&gt;alert(1)&lt;/script&gt;", "escapeHtml neutralises tags");
assert.strictEqual(escapeHtml('" onmouseover="alert(1)'),
  "&quot; onmouseover=&quot;alert(1)", "escapeHtml closes attribute breakout");
assert.ok(indexSrc.includes("href=\"${escapeHtml(url)}\""), "article link href is escaped");

// Product names come back from the enzyme sheet; they become option VALUES.
const opt = optionEl("<script>alert(1)</script>", "");
assert.strictEqual(opt.tagName, "option");
assert.strictEqual(opt.value, "<script>alert(1)</script>", "product set as a value, not markup");

// ───────────────────────────────────────────────────────────────────────────
// Fix 2 (2026-09-08) — FORMULA INJECTION: no submitted string may land in a
// cell as a live formula.
// ───────────────────────────────────────────────────────────────────────────
assert.strictEqual(sheetSafe("=SUM(A1:A9)"), "'=SUM(A1:A9)", "= is neutralised");
assert.strictEqual(sheetSafe("+1+1"), "'+1+1", "+ is neutralised");
assert.strictEqual(sheetSafe("-1+1"), "'-1+1", "- is neutralised");
assert.strictEqual(sheetSafe("@A1"), "'@A1", "@ is neutralised");
assert.strictEqual(sheetSafe("\t=1"), "'\t=1", "leading tab is neutralised");
assert.strictEqual(sheetSafe("\r=1"), "'\r=1", "leading CR is neutralised");
assert.strictEqual(sheetSafe("Ada Lovelace"), "Ada Lovelace", "ordinary text untouched");
assert.strictEqual(sheetSafe("a = b"), "a = b", "an = mid-string is not a formula");
assert.strictEqual(sheetSafe(true), true, "booleans pass through as booleans");
assert.strictEqual(sheetSafe(42), 42, "numbers pass through as numbers");
assert.strictEqual(sheetSafe(undefined), "", "undefined becomes an empty cell");

// End-to-end through a real handler: every string cell must be inert.
WRITES.length = 0;
handlePathogenPost({
  labeller: '=IMPORTRANGE("1AbC","A1")',
  abstract_id: "absX",
  abstract_title: "=1+1",
  abstract_text: "plain text",
  submission_id: "subX",
  notes: "-payload",
  client_timestamp: "2026-09-08T10:00:00.000Z",
  triples: [{ compound: '=HYPERLINK("https://evil.example","click")',
              pathogen: "@evil", direction: "inhibits" }],
});
assert.strictEqual(WRITES.length, 1, "one batched append");
const written = WRITES[0].rows[0];
for (const cell of written) {
  if (typeof cell === "string") {
    assert.ok(!/^[=+\-@\t\r]/.test(cell),
      "cell would be a live formula: " + JSON.stringify(cell));
  }
}
assert.strictEqual(written[1], '\'=IMPORTRANGE("1AbC","A1")', "labeller sanitised");
assert.strictEqual(written[3], "'=1+1", "title sanitised");
assert.strictEqual(written[5], '\'=HYPERLINK("https://evil.example","click")', "compound sanitised");
assert.strictEqual(written[6], "'@evil", "pathogen sanitised");
assert.strictEqual(written[8], false, "boolean flag column unchanged");
assert.strictEqual(written[12], "'-payload", "notes sanitised");

console.log(`ok — ${board.length} people scored; leaderboard + XSS + formula-injection assertions passed`);
