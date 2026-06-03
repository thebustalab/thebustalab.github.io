/*
 * Apps Script backend for the abstract-labelling tool.
 *
 * ===========================================================================
 * REDEPLOY CHECKLIST — READ THIS FIRST
 * ---------------------------------------------------------------------------
 * This file lives in git with PLACEHOLDER constants (no secrets committed).
 * Pasting it over Code.gs in the cloud editor therefore WIPES the live values.
 * After every paste, re-set the four CONFIG constants below before deploying:
 *
 *   SHEET_ID                = the labelling Sheet's ID (from its URL)
 *   SHEET_NAME              = "labels"   (the live data tab)
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
 *   1. Create a Google Sheet with one tab matching SHEET_NAME below; first row
 *      is headers
 *        timestamp | labeller | abstract_id | abstract_title | abstract_text |
 *        compound | pathogen | direction | no_relationships_found |
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
const SHEET_ID   = "REPLACE_WITH_SHEET_ID";   // Labels sheet ID, from its URL
const SHEET_NAME = "labels";          // live data tab name
const ACCESS_KEY = "REPLACE_WITH_HIGH_ENTROPY_TOKEN";  // same value the client sends as ?k=…

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

  const sheet = SpreadsheetApp.openById(SHEET_ID).getSheetByName(SHEET_NAME);
  const ts = new Date().toISOString();
  const rows = [];

  const title = payload.abstract_title || "";
  const text = payload.abstract_text || "";
  const notes = payload.notes || "";  // per-submission, denormalised into every row

  if (payload.no_relationships_found || !payload.triples || payload.triples.length === 0) {
    rows.push([
      ts,
      payload.labeller,
      payload.abstract_id,
      title,
      text,
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
        ts,
        payload.labeller,
        payload.abstract_id,
        title,
        text,
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

  // Append all rows in one batched setValues call to minimise quota use
  if (rows.length > 0) {
    const startRow = sheet.getLastRow() + 1;
    sheet.getRange(startRow, 1, rows.length, rows[0].length).setValues(rows);
  }

  return jsonResponse({ status: "ok", appended: rows.length });
}

// ───────────────────────────────────────────────────────────────────────────
// State summary
// ───────────────────────────────────────────────────────────────────────────

function buildState() {
  // Returns { abstract_id: [distinct labeller names], ... }
  const sheet = SpreadsheetApp.openById(SHEET_ID).getSheetByName(SHEET_NAME);
  const lastRow = sheet.getLastRow();
  if (lastRow < 2) return {};

  const values = sheet.getRange(2, 1, lastRow - 1, 3).getValues();
  // columns: timestamp(0), labeller(1), abstract_id(2)
  const map = {};
  for (const row of values) {
    const labeller = String(row[1]).trim();
    const abstractId = String(row[2]).trim();
    if (!labeller || !abstractId) continue;
    if (!map[abstractId]) map[abstractId] = new Set();
    map[abstractId].add(labeller);
  }
  // Serialise sets to arrays for JSON
  const out = {};
  for (const k in map) out[k] = Array.from(map[k]);
  return out;
}

// ───────────────────────────────────────────────────────────────────────────
// Leaderboard: combined tally of labelling triples + newsletter submissions.
// Returns [{name, count}, ...] sorted descending. The client seeds zeros for
// the full roster, so missing names just don't appear here.
// ───────────────────────────────────────────────────────────────────────────

function buildLeaderboard() {
  const counts = {};

  // Labels sheet — one point per row (already one triple per row).
  try {
    const labelsSheet = SpreadsheetApp.openById(SHEET_ID).getSheetByName(SHEET_NAME);
    const lastRow = labelsSheet.getLastRow();
    if (lastRow >= 2) {
      const labellerCol = labelsSheet.getRange(2, 2, lastRow - 1, 1).getValues();
      for (const row of labellerCol) {
        const name = String(row[0]).trim();
        if (name) counts[name] = (counts[name] || 0) + 1;
      }
    }
  } catch (err) {
    // Don't fail the whole leaderboard if labels sheet is missing.
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


// ───────────────────────────────────────────────────────────────────────────
// Response helper — JSON with the right content type
// ───────────────────────────────────────────────────────────────────────────

function jsonResponse(obj) {
  return ContentService
    .createTextOutput(JSON.stringify(obj))
    .setMimeType(ContentService.MimeType.JSON);
}
