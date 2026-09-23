---
authority: intent
---

# New chapter — networks

**Status: EXTRACTION + RENUMBERING EXECUTED 2026-08-27.** Decided 2026-08-11 (Lucas): networks becomes
a **standalone chapter** teaching **both** halves — network *visualisation* (derived similarity networks)
and network *data* (relational, observed edges). The `networks/subway` and `networks/beacons` escape-room
scenarios are its exercise sets.

**What was done 2026-08-27:**
- `chapters/7_networks.Rmd` **created**. Part 1 authored (the ch.5 network section lifted, plus the two
  improvements this plan asked for: the distance→similarity **inversion** warning and the
  **threshold-is-the-analysis** section with a cutoff sweep). Part 3 authored (layout is arbitrary; edge
  meaning determines hub meaning; when a network is a hairball hiding a table). **Part 2 is an authoring
  stub** — prose is in place, code chunks are not; see *Part 2 data — published as US airports* below.
- `chapters/5_datavis_3.Rmd` — `### network plots {-}` **removed** (65 lines) and the chapter intro
  paragraph rewritten so it no longer promises similarity networks.
- **Renumbered** old 7–18 → 8–19 (12 files), `index.Rmd` child chunks updated and the new
  `7_networks.Rmd` chunk inserted after wrangling in `# (PART) STATISTICAL METHODS`.
- `AGENTS.md` chapter table + tree entry updated; in-prose `ch.N` cross-references swept in both this repo
  and `../escape_rooms/` (25 targeted replacements — **not** a blanket regex, which corrupts strings like
  "bir**ch 7**.42" in `spa/notes.md`).
- **`_template.Rmd`** "chapter 8 has a worked live example" → chapter 9.

**Not done / still open:** Part 2 code (blocked, below); a header image for the chapter
(`images/networks.png` — every other chapter has one); a **full book render has NOT been run**, so the new
chapter is unverified end-to-end.

**Unrelated defect noticed while doing this:** `15_generative_language_models.Rmd` (old 14) has **no
`child=` chunk in `index.Rmd`** and therefore is not in the rendered book. Pre-existing, not caused by the
renumbering. Recorded in the `AGENTS.md` chapter table.

## Why a standalone chapter

Networks currently live as one section inside `chapters/5_datavis_3.Rmd`, a grab-bag of seven
unrelated plot types (3D scatter, similarity networks, marginal summaries, distributions,
Venn, ternary, maps) with no progression between them. Two problems follow:

1. **Only half the technique is taught.** Ch.5 teaches network *drawing*: compute a distance
   matrix, threshold it, plot nodes and edges. Every edge is derived. The chapter never
   touches a network whose edges are **observed facts** — and that is where the actual
   analysis lives (degree, components, betweenness, articulation points, flow).
2. **No ladder is possible.** A scenario is supposed to track its chapter's technique
   sequence, but ch.5 has no sequence — just seven unrelated plots.

The organising idea for the new chapter, and the thing worth teaching:

> **What a hub *means* depends entirely on what an edge means.**
> In a similarity network an edge means "these two are alike", so a hub is a *mixture* —
> something that belongs to nothing. In a transit or interaction network an edge means "these
> two are connected", so a hub is the most important place in the system. Same topology,
> opposite significance. Reading a network without knowing what its edges mean is how people
> get networks wrong.

## Placement — chapter 7, CONFIRMED (Lucas, 2026-08-11) — **DONE 2026-08-27**

Insert **between wrangling (6) and hierarchical clustering (7)**, renumbering current
7–18 → 8–19. Agreed; **executed 2026-08-27**.

Rationale, and it is a genuine dependency chain rather than a preference:
- The relational half **needs `group_by`/`summarise`** (aggregating journeys to stations), so
  networks must come *after* wrangling (6).
- The derived half **introduces distance matrices**, which hierarchical clustering then builds
  on directly — so networks → hclust reads as one continuous thread about similarity.

**Cost:** renumbering 12 chapter files, updating the `index.Rmd` child chunks (chapters are
merged by explicit child chunks in reading order, *not* auto-collected), and updating the
chapter table in `AGENTS.md`. Any in-prose "see chapter N" cross-references need a sweep.
Escape-room links are absolute URLs into the `escape_rooms` repo and are **not** affected.

Alternative if renumbering is unwanted: append as a high number and place it correctly via the
`index.Rmd` child order. Ordering is explicit there, so the filename number and the reading
position need not agree — but that leaves the filenames misleading, so it is not recommended.

## Chapter outline (draft)

**Part 1 — networks you build (derived / visualisation).** Mostly a lift-and-improve of the
existing ch.5 network section.
- what a network is: nodes, edges, an edge list
- from a data matrix to a distance matrix (`runMatrixAnalysis(analysis = "dist")`)
- melting a distance matrix to an edge list; the order-invariant pair key; dropping
  one-direction links
- turning distance into similarity: `edgeweight = 1/(1 + distance) * 100` — **the inversion
  students get backwards**, plotting the most *dissimilar* things as the most connected
- **the threshold is the analysis.** Too low and everything connects (a hairball that says
  nothing); too high and it falls to dust. Sweep it and show the structure changing. The
  clusters a reader sees are partly an artefact of where the cutoff was put.
- `buildNetwork()` + plotting; joining node attributes and colouring by them

**Part 2 — networks you're given (relational / analysis).** New material.
- edges as observed facts: an edge list you are handed, not one you compute
  (`buildNetwork()` already accepts a bare two-column edge list — weights optional, **no new
  tooling needed**)
- directed vs undirected; edge weights as counts; asymmetry (an origin-destination table is
  not symmetric, and the asymmetry is the finding)
- node-level summaries: degree, and **strength/volume** via `group_by`/`summarise` over both
  ends — the bridge back to ch.6
- **volume is not importance.** The busiest node can be structurally irrelevant; a quiet node
  can hold the whole graph together. Demonstrated on the subway data below.
- components and communities: who exchanges with whom
- critical nodes: articulation points / betweenness — which node's removal breaks the network

**Part 3 — reading networks honestly.** Short closing section.
- layout is arbitrary; two plots of the same graph can look nothing alike
- what an edge means determines what a hub means (the organising idea above)
- when a network is the right picture and when it is a hairball hiding a table

## Worked examples — both verified, ready to lift

**Part 1 (derived) — gene co-expression.** The escape-room scenario's shipping dataset,
`../escape_rooms/rooms/networks/subway/data/lichen_expression.csv` (24 genes × 200 samples).
Demonstrates the chapter's central claim on real biology: **abundance is not connectivity.**
The most abundant gene in the transcriptome (`XAN_3514`, a cell-wall structural gene, 5644
units) has degree 3 and is a bystander. The hub of the pigment module is `XAN_1192`, the
**lowest-expressed gene in the entire dataset** (72 units, rank 24 of 24), connected to every
other member of its module.
You cannot find regulators by ranking expression level — which is the reason co-expression
analysis exists.

Two constraints discovered while building it, kept as **design notes for us, not chapter
content** (Lucas, 2026-08-11) — full statement in the scenario's `notes.md`: correlation
modules become **cliques** unless module membership is graded, and a gene cannot correlate
above **1/√2 ≈ 0.707** with two orthogonal programmes at once, which makes connector genes and
therefore articulation points **impossible** on a correlation network.

**Part 2 (relational) — passenger movement.** Now the US airports dataset (`passenger_flows` +
`us_airports`); numbers and design in *Part 2 data — published as US airports* below. It replaces the
20-station subway prototype (`../escape_rooms/rooms/networks/subway/_scratch/superseded/passenger_flows.csv`,
local-only). Same lesson: the busiest node can be deleted with zero structural effect while a
near-quietest one holds the network together. Note this lesson is only available on **relational** data — see the
0.707 constraint above for why a correlation network cannot produce it.

A third verified dataset, `superseded/tunnel_dust.csv` (24 stations × 11 analytes of settled
tunnel dust), is available as an alternative Part 1 example if a chemistry framing is wanted
instead of genes.

## Exercises

`../escape_rooms/rooms/networks/subway/` — "the subway", puzzle phase complete. Its `notes.md`
carries the verified ladder and the escape design.

`../escape_rooms/rooms/networks/beacons/` — **"Line of Sight", puzzle + story phases complete
(2026-08-13).** This is the exercise set that closes the *"Part 2 needs its own exercises"* item.

Its ladder tracks Part 2's sequence 1:1 on a handed-over directed edge list of **garrison dispatch
traffic** (`data/dispatch_ledger.csv`, 20 nodes, 184 edges, seeded generator + full verification in
`_scratch/build_dispatch_ledger.py`):

| rung | chapter beat | answer |
|---|---|---|
| 1 | edges as observed facts | heaviest single link — a quarry to its forge, 3400 (40% clear of any other pair's best) |
| 2 | direction; *the asymmetry is the finding* | the frontier watch that sends 7.55× what it receives (500% clear) |
| 3 | node summaries over both ends — **the shortcut** | the depot, busiest by 131% |
| boss | volume is not importance; critical nodes | the gate fort: sole articulation point, **20th of 20 by volume** |

Removing the busiest garrison leaves the network completely unchanged; removing the quietest splits it
into 11 and 8. **Bearing on open item 4** ("how far to take Part 2"): the boss is deliberately
answerable with **components alone** — remove a node, recount components — so the scenario does **not**
depend on whether betweenness and articulation points get taught properly or only demonstrated.
Betweenness reproduces the same answer (5.6× the next node) for anyone who reaches for it.

The escape is data-free and sits outside the dataset entirely: a sightline network between beacon posts,
collected by eye and never tabulated. Verified in `_scratch/verify_chain.py`; diagram in `_scratch/draw_chain.py`.

It splits the chapter's two halves across its two objectives: the **graded rooms teach Part 1**
(derived co-expression networks — modules, threshold, hubs, abundance-is-not-connectivity), and
the **escape has the player inhabit Part 2** (a relational rail network whose edges they rode).
The two datasets are deliberately unrelated, so nothing about the genes implies anything about
the stations.

That means the scenario grades Part 1 only. **Part 2 needs its own exercises** — and that is now
answered: **`../escape_rooms/rooms/networks/beacons/`** (agreed 2026-08-12) is the partner
scenario, and it grades Part 2. A signal-beacon chain along a mountain range, where **an edge means
"I can see you"** and weather *deletes* edges without moving a node — the chapter's organising idea
made physical. It reskins the `passenger_flows` data above, and satisfies the corpus's
two-scenarios-per-chapter pairing convention. Design notes only; no ladder yet.

## WebR

No blocker. Checked against the wasm repo for R 4.4 and 4.5: `igraph`, `network`, `ggnetwork`,
`sna`, `statnet.common`, `ggrepel` all have builds. Unlike hierarchical clustering (blocked on
`ggtree`), nothing here needs swapping out. `buildNetwork()` uses igraph only for the
force-directed layout and falls back to `kamadakawai` without it.

## Migration out of ch.5 — DONE, and ch.5 was dissolved entirely (2026-08-27)

The plan said: *"Worth asking at that point whether what is left still justifies a chapter."* It was
asked, and the answer was no. `5_datavis_3.Rmd` is **retired** —
`z_archive/retired_chapters/5_datavis_3.Rmd.retired_2026-08-27`, with a README recording the mapping.
Its seven sections were **distributed, not deleted**:

| former ch.5 section | new home |
|---|---|
| `### network plots {-}` | `chapters/7_networks.Rmd` — Part 1 |
| `### marginal summaries {-}` | `chapters/11_comparing_means.Rmd` — new `## looking at the distribution first {-}` |
| `### representing distributions {-}` | same section as above |
| `### 3D scatter plots {-}` | `chapters/z_specialized_plots.Rmd` (appendix) |
| `### venn diagrams {-}` | appendix |
| `### ternary plots {-}` | appendix |
| `## map data {-}` | appendix |
| `## further reading {-}` | appendix (all three links were map/ternary links) |

**Why the distribution material went to comparing means, not the appendix.** Showing a distribution in
full before you test it is the habit that chapter is already trying to build — it makes students check
normality and homogeneity of variance a section later. The raincloud plot belongs next to the t-test.

**Maps were considered for a chapter of their own** and parked as an appendix section instead (Lucas,
2026-08-27). The argument for promoting them later still stands and is worth recording: projections
genuinely distort, choosing one is an analytical decision, and the existing shoreline code already
computes **haversine distances** on a sphere. That is real analysis, unlike the other appendix plots.
If it is ever promoted, the blocking dependency is `pfas_data_private.csv` (a Mac-only absolute path
that already breaks cold renders) — it would need a redacted public version or a substitute dataset.

**Consequences worth knowing:**
- **There is now no chapter 5.** The numbering runs 3, 4, 6, 7 … 19. Closing the gap would mean a second
  renumbering sweep immediately after the first, which would invalidate every `ch.N` reference just
  updated across both repos — so it was left open **pending Lucas's call**.
- `images/datavis3.png` is now **unreferenced**.
- The **ggtern leak hazard** and the **private-PFAS cold-render failure**, both previously attributed to
  `5_datavis_3.Rmd` in `AGENTS.md`, now live in `z_specialized_plots.Rmd`. Both notes were updated.
- Data vis III had **no escape-room scenario** and now never needs one — one fewer chapter to pair.

## Open items

1. ~~**Confirm placement + renumbering**~~ **DONE 2026-08-27** — chapter 7, old 7–18 → 8–19, executed.
2. ~~**Possible bug in the existing ch.5 network example.**~~ **RESOLVED — and the underlying hazard was
   removed at the root 2026-09-22.** The original worry came from two similarly-named functions returning
   **different types** for the same `analysis = "dist"`. There is now **one** function:
   `runMatrixAnalysis()`. The two return shapes survive as an explicit, named choice rather than a name
   collision — the default gives a base R `dist` object, and `output_format = "long"` gives the
   one-row-per-pair frame. Chapter 7 passes `output_format = "long"` at all three call sites. Both paths
   are now **executed**, in a real browser, by `escape_rooms/tests/matrix_analysis_smoke.mjs` — which also
   answers the old caveat that this had only ever been settled by reading the source.
3. **Ch.5's exercises block is wrong regardless** — `5_datavis_3.Rmd:387–397` is commented out
   and its text is about normality tests and t-tests, copy-pasted from comparing means.
4. ~~**How far to take Part 2.**~~ **SETTLED 2026-09-22 — by demonstration, not by formula.** Degree and
   strength are taught with `group_by()`/`summarize()`. Components and articulation points are shown by
   *deleting a node with `filter()` and re-plotting* — the `### node removal {-}` section — so the only
   new idea is the deletion itself, and no new function is introduced. The term **articulation point** is
   named at the moment it is demonstrated, and `igraph` is pointed at in further reading for anyone who
   wants the numbers. Betweenness is not taught at all; it earned no place once the picture made the
   point.
5. ~~**Does the chapter need a helper function** for degree/components?~~ **SETTLED 2026-09-22 — no
   helper, and the one-line shortcut once floated here is a TRAP.** `table(c(edges$from, edges$to))`
   double-counts whenever every route appears in both directions, which is true of the flight ledger AND
   of any long-format distance matrix, before and after thresholding — so the same line silently returns
   different answers depending on a property of the data the student was never asked to check. The
   chapter groups each column separately instead, which is longer and never asks anyone to guess. The
   trap itself is now written up in the chapter as a caution.

## Duplicate edges and scaling (2026-09-18, session "network")

- **Duplicate edges fixed in the toolkit, not the chapter.** The long distance matrix lists every pair in
  both directions with identical distances, so Part 1 was drawing every edge twice (darker alpha, doubled
  layout weight, 380 edges for 20 lakes instead of 190). Lucas chose an all-or-nothing rule: if every
  repeated pair agrees on its edge values, `buildNetwork()` collapses to one edge per pair; if any differ,
  every edge is kept; `directed =` overrides. Behaviour detail lives in `../phylochemistry/AGENTS.md`.
  Lucas's manual `table(duplicated(...))` check was removed from the chapter; the thresholds paragraph and
  the Part 2 intro now explain the rule.
- **Part 2 consequence:** the dispatch/journey edge lists are directed. When their chunks are written,
  pass `directed = TRUE`, or reciprocal pairs with equal (or no) values will be merged.
- **Scaling explained here, not in ch.8.** The `## distances and scaling {-}` section, after the unscaled
  `dist()`: the unscaled Alaska matrix correlates 0.99 with a chloride-only distance (Cl sd ~78 mg/L vs P
  sd ~0.001), then `scale()`, then when not to scale. `scale_variance` is set explicitly on the
  `runMatrixAnalysis(analysis = "dist")` calls rather than relied on as a default.
- The two stale fragments once listed here — the thresholds paragraph naming "`wood_1`, `wood_2`" and the
  unexplained hard-coded `filter(distance < 3.672035, ...)` — are both gone from the chapter.

## Part 2 data — published as US airports (2026-09-18)

The old blocker (the subway `passenger_flows.csv` lived only in a gitignored `_scratch/` folder) is
resolved by a **new** dataset rather than by moving the old one. At Lucas's request it uses **US
airports**, so students can plot the map beside the force layout and see that layout is not geography.

- Files: `../phylochemistry/sample_data/passenger_flows.csv` + `us_airports.csv`, registered in
  `modules/datasets.R` as `passenger_flows` and `us_airports`. Generator + assertions:
  `notes/build_passenger_flows.R` (phylochem conda R).
- Simulated counts, real codes and coordinates, simplified routes (Alaska reaches the lower 48 only via
  Ketchikan; a Fairbanks–Juneau link was added so Anchorage is not a second articulation point).
- Verified: ATL busiest (1.90M, 36% over LAX), removal changes nothing (1 component stays 1); KTN 27th
  of 29 by volume, sole articulation point, removal gives 2 components; KTN tops betweenness but only
  1.2x SEA; Hawaii airports top inflow/outflow (1.11–1.18). `buildNetwork()` keeps all 192 rows (directed).
- Layout lesson, checked by drawing it: with passenger weights the lower 48 collapses into a tight
  ball and Alaska/Hawaii swing to different sides on each run; without weights the layout is readable
  and HNL sits beside SFO/LAX, DLH beside DTW, KTN between PDX and JNU. Layout-to-map distance
  correlation is ~0.5, because routes partly follow geography.
- **Delivered 2026-09-22:** the Part 2 code chunks are written, and the old subway placeholder numbers
  (20 stations, 40,000, 301%, 19th of 20, 2→3 components, 8x betweenness) are replaced with the airport
  ones, each re-derived from the CSV rather than copied. The TODO comment that carried them is gone.

## Terminology, the function merge, and a local render path (2026-09-22, session "networks")

**The chapter's spine is now `similarity networks` vs `interaction networks`.** Earlier drafts paired
derived/relational, then derived/observed. Both mixed two independent axes: *provenance* (computed vs
handed to you) and *direction* (does the relation have an orientation). The pair that survived names what
an **edge means**, which is what the chapter's closing warning turns on, and it lets the co-expression vs
protein-protein-interaction confusion be stated directly — a confusion worth pre-empting because students
meet a co-expression network in the exercises. "Observed" was rejected late: the distances section already
opens by calling the lake measurements observations, so the word collides with itself.

**Direction is NOT the provenance axis** and now lives in `### directed and undirected {-}` under
interpretation. The settled account: a similarity network built from distances is *always* undirected
because distance is symmetric; an interaction network may be either; and a directed network that comes out
**balanced** is still directed — the balance is a result, destroyed by collapsing it. That last case is why
`buildNetwork()` cannot decide for the user: identical reciprocal values mean "listed twice" in a distance
matrix and "genuinely balanced" in a flight ledger, and the numbers look the same either way. Nice accident:
the airport data *is* a balanced digraph — every airport has equal in- and out-degree with unequal passenger
counts — so the example proves its own point.

**`runMatrixAnalyses` (plural) is retired; `runMatrixAnalysis` (singular) absorbed it.** They were two
independent implementations (~550 and ~680 lines) of the same pipeline, differing mainly in that the plural
was wide-input-only. Three capabilities were ported into the singular *first*: `tree_method`, `gower` via
`cluster::daisy`, and the long-format `dist` return. Two traps surfaced doing it, both recorded in
`../phylochemistry/AGENTS.md`: `tree_method` had to be added to the argument-name **whitelist** (which
tests every formal, not just what was passed, so omitting it fails *every* call), and `tree_method` defaults
to `linkage_dendrogram` **not** the plural's `neighbor_joining`, because escape-room trees are graded
against the linkage route. A named regression guard in the smoke test pins that default. Gower was a latent
bug: it was an accepted value the singular never implemented, so asking for it silently returned euclidean.

**Voice: expository prose is impersonal, book-wide.** Lucas asked for "you"/"your" out of the explanation.
Applied to ch.7 by hand and then to ch.3-4, 6, 8-19. Convention and the four exempt categories are canon in
`chapters/AGENTS.md` -> *Voice*. The thing worth remembering is the method: **recast, don't substitute** —
across the whole book the conversion needed an actor noun only seven times, so "the analyst" appearing often
would be the smell of a pass done badly. `2_installation.Rmd` and `index.Rmd` were deliberately left whole,
being procedure and course policy.

**This box can render the book.** The `phylochem` conda env already had the full stack; the only missing
piece was the harness, which `_agent_reference/build_harness.md` had recorded as lost with its paths
"unknown -- ask Lucas". They were all findable on the box. `_local_harness.R` + `_render_check.R` are the
result; detail in that reference file. Rendering earned its keep immediately by catching a bug no amount of
reading would have: **`S4Vectors::rename()` masks `dplyr::rename()`** in a full-scope book session, so a
bare `rename()` in a chapter dies with a message about arguments needing to be character. Invisible in WebR
(no Bioconductor) and in a plain tidyverse session. Canon in `AGENTS.md`.

**Figure captions and sizes** are now on every figure-producing chunk in ch.7, following the caption guide
in `index.Rmd`. Steps 5 and 6 of the chapter-enhancement checklist in `chapters/AGENTS.md` were added for
this, including the point that plot sizing is a **render-time** check that cannot be done by reading source.

**Still genuinely owed for ch.7:** nothing blocking. Chapters 8-10 have many figures still at the 4x3 global
default -- i.e. never explicitly sized -- which `_render_check.R` flags by design.
