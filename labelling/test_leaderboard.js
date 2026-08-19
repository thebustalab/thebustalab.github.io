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
          }),
        };
      },
    };
  },
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

console.log(`ok — ${board.length} people scored; all leaderboard assertions passed`);
