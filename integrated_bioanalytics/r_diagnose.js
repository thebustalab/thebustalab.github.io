/*
 * r_diagnose.js — turn R's error output into something a student can act on.
 *
 * Editing this file? Read integrated_bioanalytics/_agent_reference/webr_surfaces.md first — the design decisions, the failure histories and the anti-"improvement" guards live there.
 *
 * WHY THIS EXISTS. WebR hands back only the FIRST LINE of an R error. R itself produced more: for a
 * mistyped column name, R said "Problem while computing aesthetics." AND "Caused by error: object
 * 'watertemp' not found" — and webR drops the second half, which is the half naming the actual mistake.
 * The clean fixes (a global calling handler, rlang::global_entrace) are both refused by webR, because its
 * capture machinery already has handlers on the stack. So instead of recovering R's text we reconstruct
 * the MEANING from the headline, which survives intact and is a reliable signpost.
 *
 * Two independent halves:
 *   structuralProblems(src)  — parse-level faults found by walking our own token stream: unclosed
 *                              brackets, unterminated strings, a line left dangling on %>% or +. This
 *                              beats R's parse error, which reports where the parser GAVE UP rather than
 *                              where the mistake is: leave a `(` open on line 1 of a 20-line script and R
 *                              points at line 20. We can point at line 1.
 *   explainError(text, ctx)  — a plain-English gloss of a runtime error, plus the next thing to TYPE.
 *
 * The house style for a hint: say what R meant, then name a command that moves the student forward. Never
 * rewrite their code for them — pointing at `colnames(d)` teaches a habit, pasting the right column name
 * does not.
 *
 * PURE — no DOM, no WebR. Tested in Node (`tests/r_diagnose.test.mjs`).
 */
import { tokenizeR } from "./r_highlight.js?v=1";

const OPENERS = { "(": ")", "[": "]", "{": "}" };
const CLOSERS = { ")": "(", "]": "[", "}": "{" };

// 1-based line number of a character offset.
export function lineAt(src, offset) {
  let line = 1;
  for (let i = 0; i < offset && i < src.length; i++) if (src[i] === "\n") line++;
  return line;
}

/*
 * Structural faults, most-useful-first. Returns [] for code that is structurally fine — which includes
 * plenty of code that will still fail at runtime; this half only knows about shape.
 *
 * Only the FIRST bracket fault is reported. One missing bracket cascades into a pile of downstream
 * complaints, and a student handed five problems does not know which one is real.
 */
export function structuralProblems(src) {
  const problems = [];
  const tokens = tokenizeR(src);

  // An unterminated string swallows the rest of the file, so everything after it is noise — report it
  // alone and stop.
  const openString = tokens.find(t => t.unterminated);
  if (openString) {
    const q = openString.quote === "`" ? "backtick" : `${openString.quote} quote`;
    problems.push({
      kind: "unterminated-string",
      line: lineAt(src, openString.start),
      message: openString.quote === "`"
        ? `The backtick name starting on line ${lineAt(src, openString.start)} is never closed.`
        : `The text starting on line ${lineAt(src, openString.start)} is never closed — it needs a matching ${q}.`,
    });
    return problems;
  }

  const stack = [];
  for (const t of tokens) {
    if (t.type !== "op") continue;           // brackets inside strings/comments are not brackets
    const v = t.value;
    if (OPENERS[v]) { stack.push(t); continue; }
    if (!CLOSERS[v]) continue;
    const open = stack.pop();
    if (!open) {
      problems.push({
        kind: "unexpected-close",
        line: lineAt(src, t.start),
        message: `There is a ${v} on line ${lineAt(src, t.start)} with no matching ${CLOSERS[v]} before it.`,
      });
      return problems;
    }
    if (OPENERS[open.value] !== v) {
      problems.push({
        kind: "mismatched-bracket",
        line: lineAt(src, open.start),
        message: `The ${open.value} opened on line ${lineAt(src, open.start)} is closed by a ${v} on line ${lineAt(src, t.start)}.`,
      });
      return problems;
    }
  }
  if (stack.length) {
    const open = stack[0];
    problems.push({
      kind: "unclosed-bracket",
      line: lineAt(src, open.start),
      message: `The ${open.value} opened on line ${lineAt(src, open.start)} is never closed — it needs a ${OPENERS[open.value]}.`,
    });
    return problems;
  }

  // Code that ends mid-expression. R's "unexpected end of input" is the same symptom as an unclosed
  // bracket, and students hit this constantly by leaving a trailing %>% or + on the last line.
  const meaningful = tokens.filter(t => t.type !== "ws" && t.type !== "co");
  const last = meaningful[meaningful.length - 1];
  if (last && last.type === "op" && !CLOSERS[last.value] && ![")", "]", "}"].includes(last.value)) {
    problems.push({
      kind: "dangling-operator",
      line: lineAt(src, last.start),
      message: `Your code ends with ${last.value} on line ${lineAt(src, last.start)}, so R is still waiting for what comes next.`,
    });
  }
  return problems;
}

// Edit distance, capped — only used to ask "did you mean this existing name?".
function distance(a, b) {
  if (Math.abs(a.length - b.length) > 3) return 99;
  const prev = Array.from({ length: b.length + 1 }, (_, i) => i);
  for (let i = 1; i <= a.length; i++) {
    let carry = prev[0];
    prev[0] = i;
    for (let j = 1; j <= b.length; j++) {
      const tmp = prev[j];
      prev[j] = Math.min(prev[j] + 1, prev[j - 1] + 1, carry + (a[i - 1] === b[j - 1] ? 0 : 1));
      carry = tmp;
    }
  }
  return prev[b.length];
}

/*
 * The closest known name, if one is close enough to be worth suggesting. Case-insensitive, because
 * `View`/`view` and `TRUE`/`True` are exactly the sort of slip this should catch.
 */
export function didYouMean(name, candidates) {
  let best = null, bestD = Infinity;
  for (const c of candidates || []) {
    const d = distance(name.toLowerCase(), c.toLowerCase());
    if (d < bestD) { bestD = d; best = c; }
  }
  // Allow 1 edit for short names, 2 for longer ones — enough for a typo, not enough to guess wildly.
  const limit = name.length <= 4 ? 1 : 2;
  return bestD <= limit ? best : null;
}

const PACKAGE_LIST = "dplyr, ggplot2, tidyr, readr and stringr";

/*
 * Explain a runtime error. `text` is what the console rendered; `ctx` is {source, objects, packages}.
 * Returns {message, try} or null when we have nothing useful to add — and saying nothing is the right
 * answer more often than it looks. A wrong hint is worse than no hint: it sends the student off to fix
 * something that was never broken.
 */
export function explainError(text, ctx) {
  const c = ctx || {};
  const objects = c.objects || [];
  const src = c.source || "";
  const hint = (message, tryThis) => ({ message, try: tryThis || null });

  // Parse errors: our own structural analysis is strictly better than R's position, so prefer it.
  if (/unexpected (end of input|symbol|string constant|'.*')|unexpected INCOMPLETE/i.test(text)) {
    const problems = structuralProblems(src);
    if (problems.length) return hint(problems[0].message);
    return hint("R could not read your code as written — the usual causes are a missing comma, a missing operator between two things, or a bracket in the wrong place.");
  }

  let m = text.match(/object '([^']+)' not found/);
  if (m) {
    const name = m[1];
    const near = didYouMean(name, objects);
    if (near) return hint(`R has nothing called ${name}. There is one called ${near} — is that the one you meant?`);
    return hint(
      `R has nothing called ${name}. Check the spelling, and check the Environment pane — if it is not listed there, it does not exist yet.`,
      `colnames(alaska_lake_data)  # if ${name} is meant to be a column, it only works inside the data`);
  }

  m = text.match(/could not find function "([^"]+)"/);
  if (m) {
    const near = didYouMean(m[1], ["filter", "select", "mutate", "summarise", "group_by", "arrange",
      "ggplot", "aes", "geom_point", "geom_bar", "geom_line", "mean", "median", "sum", "length",
      "nrow", "ncol", "colnames", "head", "view", "str_wrap", "pivot_longer", "pivot_wider"]);
    if (near) return hint(`R has no function called ${m[1]}. Did you mean ${near}()?`);
    return hint(`R has no function called ${m[1]}. Check the spelling, or it may live in a package this sandbox does not load (it has ${PACKAGE_LIST}).`);
  }

  m = text.match(/there is no package called ['‘]([^'’]+)['’]/);
  if (m) {
    return hint(`The package ${m[1]} is not available in the browser — packages cannot be installed here. This sandbox has ${PACKAGE_LIST}.`);
  }

  // The headline webR leaves us when a name in aes() is wrong; R's own "Caused by: object 'x' not
  // found" is the part webR threw away, so reconstruct the meaning instead of the text.
  if (/Problem while computing aesthetics|Problem while mapping/i.test(text)) {
    return hint(
      "One of the names inside aes() does not match a column in your data. R usually names the culprit here, but the browser drops that part of the message — so check the spelling of each name in aes() against the columns.",
      "colnames(your_data)");
  }

  if (/We detected a named input|named argument/i.test(text)) {
    return hint('You used = where R wants ==. Inside filter(), lake == "Lava_Lake" TESTS whether they are equal; lake = "Lava_Lake" tries to ASSIGN, which is not allowed there.');
  }

  if (/Discrete value supplied to a continuous scale/i.test(text)) {
    return hint("That column holds text or categories, but the scale you asked for expects numbers. Either drop the scale_*_continuous() line, or plot a numeric column on that axis.");
  }

  if (/Continuous value supplied to a discrete scale/i.test(text)) {
    return hint("That column holds numbers, but the scale you asked for expects categories.");
  }

  if (/non-numeric argument to binary operator/i.test(text)) {
    return hint("You are doing arithmetic on something that is not a number — often a text column, or a whole data frame where a single column was meant.");
  }

  if (/argument "([^"]+)" is missing, with no default/.test(text)) {
    const arg = text.match(/argument "([^"]+)" is missing/)[1];
    return hint(`The function needs its ${arg} argument and you have not given it one.`);
  }

  if (/\$ operator is invalid for atomic vectors/i.test(text)) {
    return hint("You used $ on something that is not a data frame or list — often a single column that has already been pulled out.");
  }

  if (/undefined columns selected|subscript out of bounds/i.test(text)) {
    return hint("You asked for a column or position that is not there.", "colnames(your_data)");
  }

  // dplyr wraps failures from inside a verb like this, and again webR drops the "Caused by" detail.
  m = text.match(/In argument: `([^`]+)`/);
  if (m) {
    return hint(`Something went wrong inside ${m[1]}. A common cause is a bare word where text was meant — R reads Lava_Lake as the name of an object, while "Lava_Lake" is the text.`);
  }

  return null;
}

/*
 * The silent failure: `ggplot(...)` then `geom_point()` on the next line WITHOUT a trailing +. R runs two
 * separate expressions, prints the internals of the second, and draws nothing. No error at all, which is
 * why it needs its own check — the student sees output and reasonably assumes it worked.
 */
export function looksLikeOrphanLayer(outputText) {
  if (!outputText) return null;
  const t = outputText.trim();
  if (!/^(mapping:|geom_|stat_|position_|<ggproto)/m.test(t)) return null;
  if (!/^(mapping:|geom_|stat_|position_|<ggproto)/.test(t.split("\n")[0])) return null;
  return {
    message: "That looks like a single ggplot layer printed on its own rather than a plot. Check that the line before it ends with a + — every layer has to be joined to the one above.",
  };
}
