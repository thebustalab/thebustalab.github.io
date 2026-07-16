---
authority: intent
---

# Hawai‘i Aquifers — room sequence design notes

Escape-room case for the **Data Visualization** chapter (first data-vis
assignment). Three practice rooms then a boss room. The boss is the existing
CHEM 5725 exercise (`teaching/CHEM5725/exercises.csv`, Data Visualization Q3 —
Kiana & Dr Kamaka, saltwater intrusion). The three lead-up rooms scaffold the
skills the boss needs — **filtering rows and building plots** — each adding
exactly one new move, all on the same dataset so it's familiar by the boss.

## The dataset (`hawaii_aquifers`)

`https://thebustalab.github.io/phylochemistry/sample_data/hawaii_aquifers.csv`

- **Long format**, one row per well per analyte. Columns: `aquifer_code`,
  `well_name`, `longitude`, `latitude`, `analyte`, `abundance`.
  (longitude/latitude are mostly `NA`.)
- 106 wells across **10 aquifer codes** (`aquifer_1` … `aquifer_10`).
- **9 analytes:** Ca, Cl, HCO3, K, Mg, Na, SiO2, SO4, `dissolved_solids`.
- `dissolved_solids` reaches the hundreds; every other analyte sits under ~100.
  Cl max = 99, Na max = 90.

## Why reverse-engineer from the boss

The boss requires a student to stack several moves:
`filter(analyte %in% c("Na","Cl"))` → `filter(abundance > 50)` →
`facet_grid(analyte ~ .)` → colour/label by `well_name` → spot the joint
outliers → map wells to a region. That's 4–5 ideas at once. So each practice
room introduces **one** of those moves, on the same data.

## The ladder (no filter → one filter → multi-value filter + facet → full combo)

### Room 1 — first plot, no filtering
- **Skill:** plot grammar only — pick data, map x/y, choose a geom.
- **Code shape:** `ggplot(hawaii_aquifers) + geom_point(aes(x = abundance, y = analyte))`
- **Question:** which analyte reaches by far the highest values?
- **Answer:** `dissolved_solids` (hundreds vs everything else < ~100).
  Un-eyeballable without plotting; needs no filter. An easy first win that
  familiarises them with the columns.

### Room 2 — one filter
- **Skill:** `filter()` with a single `==` condition, then read an outlier.
- **Code shape:** `filter(analyte == "Cl")` → plot Cl across wells.
- **Question:** which well has the most chloride?
- **Answer:** `LALAMILO_D` (aquifer_7, Cl = 99). The core boss move in miniature.

### Room 3 — two-value filter + facet
- **Skill:** `%in%` for multiple values, and `facet_grid` small multiples.
- **Code shape:** `filter(analyte %in% c("Na","Cl"))` + `facet_grid(analyte ~ .)`
- **Question:** which aquifer stands out as high in **both** sodium and chloride?
- **Answer:** aquifer_7 (Lalamilo wells top both; aquifer_2 / Kamaile is the
  runner-up). This is the boss minus the threshold, the colour-by-well, and the
  geography step.

### Boss — the CSV question as-is
- Add `abundance > 50`, colour/label by `well_name`, identify the specific
  wells, google them to a region.
- Reference code (from exercises.csv):
  ```r
  hawaii_aquifers %>% filter(analyte %in% c("Cl","Na")) %>%
    ggplot() + geom_point(aes(x = abundance, y = aquifer_code)) +
    facet_grid(analyte ~ .)
  hawaii_aquifers %>% filter(analyte %in% c("Cl","Na"), abundance > 50) %>%
    ggplot() + geom_point(aes(x = abundance, y = aquifer_code, color = well_name)) +
    facet_grid(analyte ~ .)
  ```

## Room look & feel (decided)

**Format: pseudo-360 panorama** (the `alaska_pano/` direction), not the flat
two-screen style. Each room is a gpt-image-2 scene you look around, with
clickable hotspots that pop the puzzle and a door/scene swap on solve.
**The flat rooms (`alaska/`, `datavis1/`, `demo_hub/`) are to be deleted** —
superseded by the pano approach. (Deletion not yet done; awaiting the go-ahead.)

**Narrative:** one continuous case following Kiana and Dr Kamaka (from the
exercises.csv story), same two characters throughout, a distinct backdrop per
room that tracks the fieldwork:

- **Room 1 — the Honolulu lab bench.** All the water samples spread out for a
  first look. Pairs with the overview plot (no filter).
- **Room 2 — a coastal wellhead site.** Out in the field, zeroing in on
  chloride. Pairs with the one-filter step.
- **Room 3 — back at the bench.** Sodium and chloride compared side by side.
  Pairs with the two-analyte + facet step.
- **Boss — the big-island map room.** Deciding which community to warn. Pairs
  with the full analysis + geography step.

Visual continuity across all four: same field-station world, dusk/night,
teal-and-amber palette, painterly cinematic, no people, no text. Each scene
carries natural hotspot candidates (a laptop showing the plot, the island map,
shelves of labelled sample bottles).

## Format & mechanic decisions

- **Puzzle mechanic — phased.** Phase 1: build all four rooms **multiple-choice**
  (what the engine does today) so the whole chain is playable end-to-end and the
  ladder is proven. Phase 2: upgrade to the **console-check** mechanic from
  `../../puzzle_types_design_notes.md` (student writes the pipeline, assigns to a
  named variable, hits Check, engine grades on the live R session) — clearly the
  better pedagogy for a filter-and-plot sequence, but needs the one engine change
  built first. Biggest wins first; validate the sequence before the polish.
- **Engine mode:** the chain is a `flow: "journey"` scenario (rooms worked in
  order, solving each unlocks the next), already supported by `escape-engine.js`.

## Scene prompts (authoring)

Draft gpt-image-2 prompts for the four scenes live as the **default column
prompts in `alaska_pano/harness_gpt.html`** (`ROOM_PROMPTS`, one per column,
tags `room1`…`room4`). Generate/tune them there; copy the winners back here or
into the scenario when the room is built.

**Each scene includes a CLOSED DOOR** — the swap portal that opens on solve — and
each has a matching **open-door modifier** (`DOOR_PROMPTS`, same index in the
harness), so the columns pair scene↔open all the way down. The closed doors:
Room 1 a heavy wooden door on a side wall; Room 2 a weathered plank supply door
in the shelter's solid wall; Room 3 a wooden door beside the island map; Boss a
tall door beneath the wall map. The door hotspot's masked `/api/dooropen` edit
repaints only that door's box, so opening one door leaves the rest of the scene
untouched.

## Authoring tooling — four-column harness (BUILT 2026-07-15)

The gpt-image-2 harness (`authoring/harness_server.py` + `alaska_pano/harness_gpt.html`,
`harness_ui` tmux on :8751) now generates **the whole room series in one pass**:
four side-by-side generate+pick columns, one per room, each with its own tag,
prompt and candidate grid, running gpt-image-2 concurrently (per-slot job state,
tag-namespaced filenames `gpt_<tag>_NNN.png`, atomic index reservation). A single
shared wrap + hotspot editor below acts on whichever base was last clicked. Full
detail in `../../AGENTS.md` → "Four-column generation".

**Still open on the pano side** (from `../../AGENTS.md`): wiring a hotspot click
to a real WebR MC puzzle (the viewer currently pops a placeholder modal) is the
next integration step before a pano room is actually playable.

## OPEN — the region answer key (affects Room 3 + boss)

The wells genuinely high in **both** Na and Cl in this data are the **Lalamilo**
wells (aquifer_7, Big Island South Kohala) and the **Kamaile** wells (aquifer_2,
Leeward O‘ahu) — neither obviously "Kona". Yet the CSV key and the existing
`datavis1/scenario.js` both point at **Kona**, and that scenario's own comment
says "adjust here if the key differs". So the well→region mapping is unresolved.
**Decide before building:** derive the correct region from the actual outlier
wells, or confirm a reading that lands on Kona. This also fixes Room 3's answer.

## TODO next time — wire the submission code + boss figure into the pano player

The pano rooms are now genuinely playable: the shared **`shared/pano-player.js`**
(+ `pano-player.css`) drives any chapter's `window.CHAPTER` — pseudo-360 rooms with a
WebR-editor puzzle + multiple-choice gate, door swap on solve, ‹ › nav. (This supersedes
the "wiring a hotspot click to a real WebR MC puzzle is the next step" note above — that
part is done.) Two pieces are still **not** wired, and they're what the real Canvas
assignment needs:

1. **Submission code.** `shared/codec.js` already mints the per-student code the R decoder
   grades, but `pano-player.js` doesn't call it yet. When a chapter's rooms are chained,
   feed each room's `{answer, attempts}` (the MC result the player already tracks) + a boss
   byte into `codec.js` `encode()` at the finish screen, and show the code to copy into
   Canvas. Keep the JS/R codec contract in sync (see top-level `AGENTS.md` → "codec contract";
   version-then-scenario-id arg order, base32 precision). Give this chapter a unique scenario
   `id` and add its decoder key in `decoder/decode_codes.R`.
2. **Boss figure deliverable.** The boss room ends in a **figure** the student builds in the
   WebR console (the plot + `buildCaption()`), then **downloads as a PNG watermarked with
   their x500 + the code** and uploads to Canvas (Lucas grades it by hand). `pano-player.js`
   already captures the x500 (`window.__x500`); add a "Download your figure" button that
   exports the WebR plot canvas to PNG with the watermark baked in. Verify `buildCaption()`
   actually runs in WebR first (bustalab function — may pull heavy deps).
