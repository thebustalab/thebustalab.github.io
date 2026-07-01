/*
 * Apps Script backend for the abstract-labelling tool.
 *
 * ===========================================================================
 * REDEPLOY CHECKLIST — READ THIS FIRST
 * ---------------------------------------------------------------------------
 * This file lives in git with PLACEHOLDER constants (no secrets committed).
 * Pasting it over Code.gs in the cloud editor therefore WIPES the live values.
 * After every paste, re-set the CONFIG constants below before deploying:
 *
 *   SHEET_ID                = the labelling Sheet's ID (from its URL)
 *   SHEET_NAME              = "labels"   (the live pathogen-labels tab)
 *   ENZYME_SHEET_ID         = the Characterised Enzymes Sheet's ID (from its URL;
 *                             may be the same spreadsheet as SHEET_ID with a
 *                             different tab, or a separate spreadsheet)
 *   ENZYME_SHEET_NAME       = "enzyme_labels"  (the enzyme data tab)
 *   ACCESS_KEY              = the token in the labeller link's ?k=...  (the
 *                             real value lives in the link / password manager,
 *                             NOT here — keep it out of git)
 *   FORM_RESPONSES_SHEET_ID = the newsletter form-response Sheet's ID
 *
 * Then: Deploy → Manage deployments → Edit (pencil) → Version: New version →
 * Deploy. Re-using the same deployment keeps the /exec URL stable.
 *
 * Sanity check (incognito, signed out): open
 *   <EXEC_URL>?k=<ACCESS_KEY>&action=leaderboard
 * Expect a JSON array. "bad token" => ACCESS_KEY mismatch; a Google Drive
 * "unable to open the file" page => deployment access isn't set to "Anyone".
 * ===========================================================================
 *
 * Deployment (one-time):
 *   1. Create a Google Sheet with a tab matching SHEET_NAME below; first row
 *      is headers
 *        timestamp | labeller | abstract_id | abstract_title | abstract_text |
 *        compound | pathogen | direction | no_relationships_found |
 *        submission_id | elapsed_seconds | client_timestamp | notes
 *      Also create the enzyme tab (ENZYME_SHEET_NAME, in ENZYME_SHEET_ID's
 *      spreadsheet); first row is headers
 *        timestamp | labeller | abstract_id | abstract_title | abstract_text |
 *        doi | enzyme_name | product | sequence | not_characterization |
 *        submission_id | elapsed_seconds | client_timestamp | notes
 *   2. Extensions → Apps Script. Paste this file's contents into Code.gs.
 *   3. Set the CONFIG constants below (see REDEPLOY CHECKLIST).
 *   4. Deploy → New deployment → type: Web app → execute as Me, access "Anyone".
 *   5. Copy the /exec URL into the labelling page's CONFIG.APPS_SCRIPT_URL.
 *   6. Generate the high-entropy access key (e.g. `openssl rand -hex 16`) and put
 *      it into ACCESS_KEY below AND into the labeller link you send out as
 *      ?k=<same value>.
 *
 * When labelling is done: Deploy → Manage deployments → Archive the web app
 * to stop accepting submissions. The Sheet stays intact.
 */

// ───────────────────────────────────────────────────────────────────────────
// CONFIG
// ───────────────────────────────────────────────────────────────────────────
const SHEET_ID   = "REPLACE_WITH_SHEET_ID";   // Pathogen labels sheet ID, from its URL
const SHEET_NAME = "labels";          // live pathogen-labels tab name
const ACCESS_KEY = "REPLACE_WITH_HIGH_ENTROPY_TOKEN";  // same value the client sends as ?k=…

// Characterised-enzyme labels destination. Set in the cloud editor; keep the
// placeholder in git. May point at the same spreadsheet as SHEET_ID (just a
// different tab) or a separate spreadsheet — either works. Enzyme-task
// submissions route here; pathogen submissions route to SHEET_ID/SHEET_NAME.
const ENZYME_SHEET_ID   = "REPLACE_WITH_ENZYME_SHEET_ID";
const ENZYME_SHEET_NAME = "enzyme_labels";

// Newsletter form-response sheet — feeds the leaderboard's "newsletter
// submission" tally. The Apps Script owner must have read access to this
// sheet. Set in the cloud editor; keep this repo file with the placeholder.
const FORM_RESPONSES_SHEET_ID = "REPLACE_WITH_FORM_RESPONSES_SHEET_ID";
const FORM_RESPONSES_SHEET_NAME = "Form Responses 1";

// ───────────────────────────────────────────────────────────────────────────
// HTTP handlers
// ───────────────────────────────────────────────────────────────────────────

function doGet(e) {
  const params = e.parameter || {};
  if (params.k !== ACCESS_KEY) {
    return jsonResponse({ status: "error", error: "bad token" });
  }
  if (params.action === "state") {
    return jsonResponse(buildState());
  }
  if (params.action === "leaderboard") {
    return jsonResponse(buildLeaderboard());
  }
  if (params.action === "products") {
    return jsonResponse(buildProducts());
  }
  return jsonResponse({ status: "error", error: "unknown action" });
}

function doPost(e) {
  let payload;
  try {
    payload = JSON.parse(e.postData.contents);
  } catch (err) {
    return jsonResponse({ status: "error", error: "bad JSON" });
  }
  if (payload.k !== ACCESS_KEY) {
    return jsonResponse({ status: "error", error: "bad token" });
  }
  if (!payload.labeller || !payload.abstract_id || !payload.submission_id) {
    return jsonResponse({ status: "error", error: "missing required fields" });
  }

  // Route by task type. Untyped payloads are pathogen submissions (back-compat).
  if (payload.task_type === "enzyme") {
    return handleEnzymePost(payload);
  }
  return handlePathogenPost(payload);
}

// ─── Pathogen submissions → SHEET_ID / SHEET_NAME ───────────────────────────
function handlePathogenPost(payload) {
  const ts = new Date().toISOString();
  const rows = [];
  const title = payload.abstract_title || "";
  const text = payload.abstract_text || "";
  const notes = payload.notes || "";  // per-submission, denormalised into every row

  if (payload.no_relationships_found || !payload.triples || payload.triples.length === 0) {
    rows.push([
      ts, payload.labeller, payload.abstract_id, title, text,
      "",   // compound
      "",   // pathogen
      "",   // direction
      true, // no_relationships_found
      payload.submission_id,
      payload.elapsed_seconds || "",
      payload.client_timestamp || "",
      notes,
    ]);
  } else {
    for (const t of payload.triples) {
      rows.push([
        ts, payload.labeller, payload.abstract_id, title, text,
        t.compound || "",
        t.pathogen || "",
        t.direction || "",
        false,
        payload.submission_id,
        payload.elapsed_seconds || "",
        payload.client_timestamp || "",
        notes,
      ]);
    }
  }

  appendRows(SHEET_ID, SHEET_NAME, rows);
  return jsonResponse({ status: "ok", appended: rows.length });
}

// ─── Enzyme submissions → ENZYME_SHEET_ID / ENZYME_SHEET_NAME ───────────────
// One row per (submission, enzyme). A "not a functional characterisation"
// verdict writes a single row with empty enzyme fields and the flag set.
function handleEnzymePost(payload) {
  const ts = new Date().toISOString();
  const rows = [];
  const title = payload.abstract_title || "";
  const text = payload.abstract_text || "";
  const doi = payload.doi || "";
  const notes = payload.notes || "";
  const enzymes = payload.enzymes || [];

  if (payload.not_characterization || enzymes.length === 0) {
    rows.push([
      ts, payload.labeller, payload.abstract_id, title, text,
      doi,
      "",   // enzyme_name
      "",   // product
      "",   // sequence
      true, // not_characterization
      payload.submission_id,
      payload.elapsed_seconds || "",
      payload.client_timestamp || "",
      notes,
    ]);
  } else {
    for (const en of enzymes) {
      rows.push([
        ts, payload.labeller, payload.abstract_id, title, text,
        doi,
        en.enzyme_name || "",
        en.product || "",
        en.sequence || "",
        false,
        payload.submission_id,
        payload.elapsed_seconds || "",
        payload.client_timestamp || "",
        notes,
      ]);
    }
  }

  appendRows(ENZYME_SHEET_ID, ENZYME_SHEET_NAME, rows);
  return jsonResponse({ status: "ok", appended: rows.length });
}

// Batched append to minimise quota use.
function appendRows(sheetId, sheetName, rows) {
  if (!rows || rows.length === 0) return;
  const sheet = SpreadsheetApp.openById(sheetId).getSheetByName(sheetName);
  const startRow = sheet.getLastRow() + 1;
  sheet.getRange(startRow, 1, rows.length, rows[0].length).setValues(rows);
}

// ───────────────────────────────────────────────────────────────────────────
// State summary
// ───────────────────────────────────────────────────────────────────────────

function buildState() {
  // Returns { item_id: [distinct labeller names], ... } merged across BOTH the
  // pathogen and enzyme sheets. Item ids are unique across the whole pool, so a
  // single merged map drives eligibility for every task type. Both sheets share
  // the column layout timestamp(0) | labeller(1) | abstract_id(2).
  const map = {};
  collectLabellers(SHEET_ID, SHEET_NAME, map);
  if (ENZYME_SHEET_ID && String(ENZYME_SHEET_ID).indexOf("REPLACE_") !== 0) {
    collectLabellers(ENZYME_SHEET_ID, ENZYME_SHEET_NAME, map);
  }
  // Serialise sets to arrays for JSON
  const out = {};
  for (const k in map) out[k] = Array.from(map[k]);
  return out;
}

function collectLabellers(sheetId, sheetName, map) {
  try {
    const sheet = SpreadsheetApp.openById(sheetId).getSheetByName(sheetName);
    if (!sheet) return;
    const lastRow = sheet.getLastRow();
    if (lastRow < 2) return;
    const values = sheet.getRange(2, 1, lastRow - 1, 3).getValues();
    for (const row of values) {
      const labeller = String(row[1]).trim();
      const itemId = String(row[2]).trim();
      if (!labeller || !itemId) continue;
      if (!map[itemId]) map[itemId] = new Set();
      map[itemId].add(labeller);
    }
  } catch (err) {
    // A missing/unconfigured enzyme sheet must not break pathogen state.
  }
}

// ───────────────────────────────────────────────────────────────────────────
// Leaderboard: combined tally of labelling triples + newsletter submissions.
// Returns [{name, count}, ...] sorted descending. The client seeds zeros for
// the full roster, so missing names just don't appear here.
// ───────────────────────────────────────────────────────────────────────────

function buildLeaderboard() {
  const counts = {};

  // Labelling sheets — one point per row (one triple, or one enzyme, per row).
  tallyLabellerColumn(SHEET_ID, SHEET_NAME, counts);
  if (ENZYME_SHEET_ID && String(ENZYME_SHEET_ID).indexOf("REPLACE_") !== 0) {
    tallyLabellerColumn(ENZYME_SHEET_ID, ENZYME_SHEET_NAME, counts);
  }

  // Newsletter form-response sheet — one point per submission (any type).
  // Submitter name lives in the "Your name" column added in May 2026.
  try {
    if (FORM_RESPONSES_SHEET_ID && !FORM_RESPONSES_SHEET_ID.startsWith("REPLACE_")) {
      const formSheet = SpreadsheetApp.openById(FORM_RESPONSES_SHEET_ID).getSheetByName(FORM_RESPONSES_SHEET_NAME);
      if (formSheet) {
        const lastRow = formSheet.getLastRow();
        const lastCol = formSheet.getLastColumn();
        if (lastRow >= 2 && lastCol >= 1) {
          const headers = formSheet.getRange(1, 1, 1, lastCol).getValues()[0];
          let nameIdx = -1;
          for (let i = 0; i < headers.length; i++) {
            if (String(headers[i]).trim().toLowerCase().indexOf("your name") !== -1) {
              nameIdx = i;
              break;
            }
          }
          if (nameIdx >= 0) {
            const rows = formSheet.getRange(2, nameIdx + 1, lastRow - 1, 1).getValues();
            for (const row of rows) {
              const name = String(row[0]).trim();
              if (name) counts[name] = (counts[name] || 0) + 1;
            }
          }
        }
      }
    }
  } catch (err) {
    // Same — don't fail the leaderboard if the form sheet is unreachable.
  }

  const out = [];
  for (const name in counts) out.push({ name: name, count: counts[name] });
  out.sort(function (a, b) { return b.count - a.count || a.name.localeCompare(b.name); });
  return out;
}

// Add one point per data row's labeller (column 2) to `counts`.
function tallyLabellerColumn(sheetId, sheetName, counts) {
  try {
    const sheet = SpreadsheetApp.openById(sheetId).getSheetByName(sheetName);
    if (!sheet) return;
    const lastRow = sheet.getLastRow();
    if (lastRow < 2) return;
    const col = sheet.getRange(2, 2, lastRow - 1, 1).getValues();
    for (const row of col) {
      const name = String(row[0]).trim();
      if (name) counts[name] = (counts[name] || 0) + 1;
    }
  } catch (err) {
    // Don't fail the leaderboard if a sheet is missing/unconfigured.
  }
}

// ───────────────────────────────────────────────────────────────────────────
// Product typeahead — distinct enzyme products already recorded, for the
// client datalist so labeller spellings converge on existing entries.
// ───────────────────────────────────────────────────────────────────────────

function buildProducts() {
  const set = {};
  try {
    if (!ENZYME_SHEET_ID || String(ENZYME_SHEET_ID).indexOf("REPLACE_") === 0) return [];
    const sheet = SpreadsheetApp.openById(ENZYME_SHEET_ID).getSheetByName(ENZYME_SHEET_NAME);
    if (!sheet) return [];
    const lastRow = sheet.getLastRow();
    if (lastRow < 2) return [];
    // product is column 8 in the enzyme schema
    const col = sheet.getRange(2, 8, lastRow - 1, 1).getValues();
    for (const row of col) {
      const p = String(row[0]).trim();
      if (p) set[p] = true;
    }
  } catch (err) {
    return [];
  }
  return Object.keys(set).sort();
}


// ───────────────────────────────────────────────────────────────────────────
// Response helper — JSON with the right content type
// ───────────────────────────────────────────────────────────────────────────

function jsonResponse(obj) {
  return ContentService
    .createTextOutput(JSON.stringify(obj))
    .setMimeType(ContentService.MimeType.JSON);
}
