# Disciplinary Expertise Post — Empirical Support to Add on Revision

*Moved out of root `todo.md` on 2026-09-14 (to-do triage). This is verified reference material for
revising the published disciplinary-expertise post, and for the "exponential gap" framing in the
university-training piece (`notes_university_piece.md`). Confidence tiers follow
`.claude/skills/citation_verification/SKILL.md`.*

---

## The core empirical claim: expertise, not profession

**Anthropic, "How Claude Code is used in practice"** — CONFIRMED tier.
Official page: https://www.anthropic.com/research/claude-code-expertise

Privacy-preserving analysis of **~400,000 Claude Code sessions**, **Oct 2025 – Apr 2026**,
**~235,000 people**. Citable numbers, consistent across the Anthropic page and multiple secondary
write-ups:

- Software engineers reached verified success in **~34%** of code-producing sessions vs **~29%** for
  non-software professionals — **a ~5-point gap**. "Every major occupation succeeds at nearly the
  same rate as software engineers."
- Expertise-level gradient: **novices 15%** verified / 77% at-least-partial success;
  **intermediate 28%** / 91%; **expert 33%** / 92%.
- Output quality scales with expertise: expert users trigger **~12 Claude actions and ~3,200 words
  of output per prompt** vs novices' **~5 actions and ~600 words** — a 2.4× action / 5× output
  multiplier.

Headline for the post: **coding agents don't substitute for domain expertise — they reward it.**

Route in via Mollick, "The Twilight of the Chatbots" (*One Useful Thing*, 2026-06-30), which is where
this surfaced.

### ⚠ Nuance to preserve — do not conflate two studies

The "domain expertise, not profession, predicts success and quality" claim is the **Anthropic**
400K-session study above. Mollick's specific line about legal, HR and other non-tech functions having
**adopted** agents at nearly the same rate appears to lean on a **separate OpenAI + academic-economists
adoption study** — that is an *adoption-rate* finding, not the *success-rate* finding.

For the disciplinary-expertise argument, Anthropic is the right primary citation. If a draft also wants
the legal/HR adoption sentence, source it to the OpenAI study separately. Verify which claim you are
making before pinning it to a source.

---

## Acceleration numbers (verified 2026-07-14)

Three distinct sources. Keep them straight — **there are two different "14-hour" figures and they must
not be conflated.**

### METR task-completion time horizons — CONFIRMED
metr.org/time-horizons; "Time Horizon 1.1", 2026-01-29.

The task duration (in human-expert time) an agent can complete at 50% reliability has grown
**exponentially for ~6 years, doubling roughly every 7 months — accelerating to every ~4 months over
2024–25**. Latest frontier models sit at **~16–20 h on the 50% horizon (~3–4 h at 80%)**, with METR's
own caveat that **measurements above ~16 h are unreliable with the current task suite**.

Best single line for a post: *"the length of task an AI can do autonomously has been doubling every
four to seven months."*

### UK AISI test-time-compute / cyber horizon — CONFIRMED (**14-hour figure #1**)
aisi.gov.uk/blog, "More compute, more capability".

On AISI's cyber-CTF suite, horizons have **doubled every ~4.7 months since late 2024**. Crucially,
raising the per-task budget from **2.5M → 50M tokens lifts the frontier horizon from ~2 h to ~14 h**.
This is the "throw more compute, get more capability" point — **a compute-scaling result, NOT a single
model run.**

### Epoch/METR "MirrorCode" — **14-hour figure #2**; headline CONFIRMED, exact stats TENTATIVE

MirrorCode tests rebuilding a full program with no access to the original source. Per secondary
write-ups (AlphaSignal; the-decoder): **Claude Opus 4.7 leads at ~56% solve rate** and reimplemented
`gotree` (~16,000 lines of Go, 40+ commands) autonomously in **~14 h for ~$251** — a job a human would
need **~2–17 weeks** for. A separate MirrorCode task reportedly ran **19 days for ~$2,600**.

**These precise numbers come from secondary sources — confirm against Epoch's own MirrorCode/model page
before putting them in a post.**

### ⚠ Framing caution

The METR and AISI numbers are **capability-horizon measures under favourable conditions** (50%
reliability, large compute budgets) — not "AI reliably does 14-hour jobs unaided." State them as
**trajectory or ceiling, not typical performance**, or a reviewer will rightly push back.
