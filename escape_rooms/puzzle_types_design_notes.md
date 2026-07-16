---
authority: intent
---

# Reusable console puzzle types — design spec

Status: design (2026-07-15). Three reusable puzzle "types" for the
Integrated Bioanalytics escape rooms, each tied to real coding skills from the
book and authored once, reused per chapter by swapping dataset + check. Grounded
in `teaching/CHEM5725/exercises.csv` (every existing exercise is a Type 1).

## The shared primitive — the console `check`

Today the engine grades a multiple-choice **index**: `makeQuestionCard` compares
`selected === q.correct`. The live WebR session (`WebRConsole.webR`) is never
inspected — the console is a scratchpad the student reads an answer *out of*,
then clicks a radio button.

Every type below leans on one new capability: **grade on the R session state.**
A node/step may carry a `check`:

```js
check: {
  requires: ["answer"],                    // vars that must exist first
  expr: "round(answer, 2) == 6.83",        // R expr → single logical, run in the student's session
  hint: "Assign your result to `answer`."  // shown if requires missing or expr errors
}
```

Engine support (one focused change, unlocks all three):

1. Expose the booted `webR` handle from the console to the question card.
2. Add a **Check my answer** button on a console-checked card. On click:
   - verify each `requires` var exists (`exists("answer")`); if not, show `hint`;
   - eval `expr` via `webR.evalRBoolean(expr)`; `TRUE` → solved, `FALSE` →
     wrong-attempt (same feedback ladder as MCQ).
3. **Codec:** a console-checked step encodes `answer = solved ? 1 : 0` plus
   `attempts` — fits the existing 5-bit/3-bit byte with no schema change. This
   also **retires the 32-option ceiling worry**: compound answers are checked by
   expression, not enumerated as options, so "class + count" tuples no longer
   need a 35-way option cross-product.

Everything else (`starterCode` per node, `maxAttempts`, two-dataset loading) is
already supported. MCQ grading stays available; a node picks `check` **or**
`correct`. Console-checked is the stronger, recommended default.

---

## Type 1 — Compute-the-Key (the backbone)

**Skill:** produce an analysis; its result *is* the answer. This is what all nine
existing quizzes already are (Q1/Q3). The craft is a dataset where the answer is
**un-eyeballable** — you must run the technique.

**Shape:** student writes a pipeline, assigns the result to a named variable,
hits Check; engine verifies the variable(s). Supports **compound** answers via a
conjunction in `expr` (the corpus's "class + count", "strain + significance").

**Worked example — wrangling (ch 6), `beer_components`** (real Q1):

```js
{
  key: "esters",
  technique: "group_by + summarise + arrange",
  prompt: "Rank the three most abundant Aliphatic_ester compounds in hops by " +
          "decreasing mean abundance. Assign their names, in order, to `top3`.",
  starterCode: "# hops esters → mean per compound → arrange desc → names → top3",
  check: {
    requires: ["top3"],
    expr: 'identical(top3[1:3], c("Methyl_6_methylheptanoate",' +
          '"Methylheptanoate","Methyl_2_methylheptanoate"))',
    hint: "top3 should be a character vector of compound names, highest mean first."
  }
}
```

**Reuse:** every chapter's Q1/Q3.
- PCA (ch 8): `which.max(abs(loadings[,"PC1"]))` name == biomarker.
- Comparing means (ch 10): a computed difference-in-means + significance flag.
- Data-vis (ch 3–5): the outlier only visible once faceted/log-transformed.

---

## Type 2 — Classify-the-Unknown

**Skill:** build a reference from labelled data, place a **mystery sample**
against it, return a verdict — often binary, often one-shot. Console-forced: the
verdict depends entirely on where the unknown lands. This is the dominant
structure in the back half of the corpus (guilty/innocent, poisonous/safe, which
patient) and *is* the escape-room narrative — there's always an unknown.

**Shape:** two datasets (reference + unknown); the student fits/clusters/embeds
the reference, projects the unknown, assigns a verdict variable. `oneShot: true`
→ `maxAttempts: 1`.

**Worked example — hierarchical clustering (ch 7), `wood_smoke` + `unknown_smoke`**
(real Q1):

```js
{
  key: "smoke",
  technique: "hierarchical clustering",
  oneShot: true,
  datasets: ["wood_smoke", "unknown_smoke"],
  prompt: "Cluster the reference smokes, place the jacket sample, and decide: " +
          "Red Oak or Paper Birch? Assign 'red_oak' or 'paper_birch' to `match`.",
  check: { requires: ["match"], expr: 'match == "red_oak"' }
}
```

**Reuse:** ch 7 dendrogram membership · ch 8 nearest in PCA space · ch 9 which
k-means/dbscan cluster · ch 11 model `predict()` class · embeddings nearest
neighbour by cosine. Same node shape; swap the technique that produces the
placement. Engine needs only the `check` primitive (two-dataset load already
works).

---

## Type 3 — Repair-the-Pipeline

**Skill:** read and debug real code — the one thing a conclusion-MCQ never tests
(a student can reach a conclusion with sloppy or AI-pasted code; they cannot fix
a broken pivot without understanding it). New to the corpus, deliberately.

**Shape:** `starterCode` pre-loads a nearly-right block that errors or returns a
telltale-wrong result. The student edits and reruns until the `check` on the
corrected output passes. The **error message itself teaches** — this is where the
live console earns its keep. (An MCQ "which edit fixes it" flavour exists but is
more guessable; console-checked is recommended.)

**Worked example — comparing means (ch 10), `algae_data`** (formula reversed):

```js
{
  key: "aov_fix",
  technique: "ANOVA",
  prompt: "This test should compare omega-3 abundance ACROSS harvesting regimes, " +
          "but the model is backwards and the p-value is nonsense. Fix it so " +
          "`fit` models abundance as a function of regime, then Check.",
  starterCode: "fit <- aov(harvest ~ omega_3_polyunsaturated_Fas, data = algae_data)\nsummary(fit)",
  check: {
    requires: ["fit"],
    expr: 'all.vars(formula(fit))[1] == "omega_3_polyunsaturated_Fas"',
    hint: "In aov(), the response goes on the LEFT of ~, the predictor on the right."
  }
}
```

**Reuse:** every technique has a canonical mistake — a ggplot missing its `aes()`,
a `pivot_longer` with the wrong `names_from`, a pipe with verbs out of order, a
`kmeans` with the wrong `centers`. Author the broken version + a check on the
corrected result.

---

## Build order (phased)

1. **Console-check primitive** + Check button + codec `solved` byte. Unlocks all
   three. Retrofit one existing `datavis1` MCQ to console-checked as the pilot.
2. **Type 1** — the backbone; convert the proven Q1/Q3 pattern.
3. **Type 2** — highest narrative payoff; the mystery-sample rooms.
4. **Type 3** — new pedagogy, smallest surface.

## Open decisions

- Keep MCQ as a fallback per node, or go console-check only once proven?
- Boss figure (Q2/Q4 image upload) stays hand-graded, or gets a lightweight
  `check` on the plot object's structure?
- Chapter 2 (data-vis II) room: which second dataset — reuse `solvents` (its real
  Q3) or pick one that hides a facet/scale reveal?
