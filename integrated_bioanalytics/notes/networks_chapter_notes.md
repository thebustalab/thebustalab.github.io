---
authority: intent
---

# New chapter — networks

**Status:** intent / planning. Decided 2026-08-11 (Lucas): networks becomes a **standalone
chapter** teaching **both** halves — network *visualisation* (derived similarity networks)
and network *data* (relational, observed edges). The `networks/subway` escape-room scenario
is its exercise set.

Nothing has been written or renumbered yet. This document is the plan.

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

## Placement — chapter 7, CONFIRMED (Lucas, 2026-08-11)

Insert **between wrangling (6) and hierarchical clustering (7)**, renumbering current
7–18 → 8–19. Agreed; the renumbering work is not yet done.

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

**Part 2 (relational) — passenger movement.**

Passenger-movement data, 20 stations, origin/destination/journeys —
`../escape_rooms/rooms/networks/subway/_scratch/superseded/passenger_flows.csv`. Built for an
earlier version of the scenario and set aside there, but it is the cleanest demonstration of
**volume versus criticality** we have:

- **Busiest station:** ~40,000 journeys, **301%** more than the runner-up.
- **Remove it:** the network is completely unchanged. 2 components before, 2 after.
- **The critical station** ranks **19th of 20** by volume — nearly the quietest in the system.
- **Remove it:** 2 components → 3. Sole articulation point, betweenness 8× the next station's.

The busiest station can be deleted with zero structural effect; the nineteenth-busiest holds
the whole thing together. Note this lesson is only available on **relational** data — see the
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

## Migration out of ch.5

When this chapter lands, `5_datavis_3.Rmd` loses its `### network plots {-}` section. The
remaining six plot types (3D scatter, marginal summaries, distributions, Venn, ternary, maps)
stay. Worth asking at that point whether what is left still justifies a chapter or should be
folded into ch.4.

## Open items

1. **Confirm placement + renumbering** (chapter 7, current 7–18 → 8–19). Lucas's call.
2. **Possible bug in the existing ch.5 network example.** `runMatrixAnalysis(analysis = "dist")`
   returns a **long data frame** (`sample_1`, `sample_2`, `distance`), but
   `5_datavis_3.Rmd:84` calls `as.data.frame(as.table(as.matrix(wood_dist)))` as though it
   were a dist object. Verify before lifting that code into the new chapter — it determines
   what students are handed. Unverified; do not assume either way.
3. **Ch.5's exercises block is wrong regardless** — `5_datavis_3.Rmd:387–397` is commented out
   and its text is about normality tests and t-tests, copy-pasted from comparing means.
4. **How far to take Part 2.** Degree and components are clearly in scope. Betweenness and
   articulation points are the payoff but are heavier; decide whether they are taught properly
   or introduced by demonstration only.
5. **Does the chapter need a helper function** for degree/components, or is base R plus
   `igraph` enough? (`table(c(edges$from, edges$to))` gives degree in one line.)
