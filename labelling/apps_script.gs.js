/*
 * Apps Script backend for the abstract-labelling tool.
 *
 * Deployment (one-time):
 *   1. Create a Google Sheet with one tab named "labels"; first row is headers
 *        timestamp | labeller | abstract_id | compound | pathogen | direction |
 *        no_relationships_found | submission_id | elapsed_seconds | client_timestamp
 *   2. Extensions → Apps Script. Paste this file's contents into Code.gs.
 *   3. Set the CONFIG constants below (SHEET_ID, ACCESS_KEY).
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
const SHEET_ID   = "REPLACE_WITH_SHEET_ID";   // Sheet ID from its URL
const SHEET_NAME = "labels";
const ACCESS_KEY = "REPLACE_WITH_HIGH_ENTROPY_TOKEN";  // same value the client sends as ?k=…

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

  if (payload.no_relationships_found || !payload.triples || payload.triples.length === 0) {
    rows.push([
      ts,
      payload.labeller,
      payload.abstract_id,
      "",   // compound
      "",   // pathogen
      "",   // direction
      true, // no_relationships_found
      payload.submission_id,
      payload.elapsed_seconds || "",
      payload.client_timestamp || "",
    ]);
  } else {
    for (const t of payload.triples) {
      rows.push([
        ts,
        payload.labeller,
        payload.abstract_id,
        t.compound || "",
        t.pathogen || "",
        t.direction || "",
        false,
        payload.submission_id,
        payload.elapsed_seconds || "",
        payload.client_timestamp || "",
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
// Response helper — JSON with the right content type
// ───────────────────────────────────────────────────────────────────────────

function jsonResponse(obj) {
  return ContentService
    .createTextOutput(JSON.stringify(obj))
    .setMimeType(ContentService.MimeType.JSON);
}
