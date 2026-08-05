# Disciplinary Expertise in an AI-Assisted World

*Working outline — companion post to "When Should a Scientist Use AI?"*

---

## Central question

AI is increasingly capable of tasks that used to require disciplinary training: summarising literature, generating code, drafting analyses, even suggesting experimental designs. If AI has (or will soon have) its own form of domain knowledge, why does deep human expertise still matter — and in what specific ways?

## Thesis

Disciplinary expertise becomes *more* important with AI, not less. AI shifts the bottleneck from generation to judgement. The skills that matter most are precisely the ones that are hardest to delegate: knowing what question to ask, recognising what matters in messy data, validating whether an answer is right for the right reasons, and maintaining the capacity to do all of this over a career.

---

## Proposed structure

### 1. The question people actually ask

- "Won't AI just learn my field?" / "Isn't domain expertise becoming obsolete?"
- Acknowledge the genuine capability: frontier models can now produce competent-looking work in many scientific domains
- But competent-looking and correct-for-the-right-reasons are different things — and telling them apart is what expertise does

### 2. What expertise actually is (and what AI approximates)

- Expertise is not a body of facts — it is the ability to make judgements under uncertainty, recognise what is surprising or anomalous, and know which details matter
- AI has access to vast training data and can pattern-match effectively — this is a real form of knowledge, and it would be dishonest to dismiss it
- But AI's "expertise" is statistical: it can tell you what is typical, but it struggles with what is *important* — a distinction that requires understanding the field's open questions, its history of dead ends, and its current tensions
- The gap: AI can produce a plausible methods section but cannot tell you whether this is the right experiment to run. It can summarise a paper's claims but may miss the buried caveat that an experienced reader would flag immediately.

**[DRAFT — added 2026-08-04, reading-list session, from Goedecke "LLMs reward expertise"; review before use]**
- Terence Tao's public ChatGPT conversation about the counterexample to the Jacobian Conjecture is a concrete, citable illustration of this gap. Tao's prompts are short and to-the-point; because he signals expertise, the model shifts into dense "talking-to-mathematicians" mode rather than "explaining-to-amateurs" mode. He pushes back on weak answers without flatly contradicting ("this looks more complex than I was hoping for"), and — critically — he drives the direction of the conversation himself, rarely taking the model's suggested next step.
- The mechanism isn't a prompting trick you can copy; it's that Tao actually understands the mathematics well enough to pull the right thread out of a multi-paragraph response, propose reformulations, and recognise what "looks weird." That's the expertise doing the work, not the phrasing.

### 3. The five things expertise does that AI cannot (yet) do

Frame these as specific, concrete capabilities rather than abstract virtues:

**a. Asking the right question**
- The most consequential intellectual act in research is choosing what to investigate
- AI can generate hypotheses from patterns in literature, but selecting which hypothesis is worth pursuing requires understanding of what the field needs, what is technically feasible, and what would actually change how we think — this is judgement shaped by years of immersion
- [Could use an example from chemistry/plant biology here]
- **[DRAFT — added 2026-08-04, reading-list session]** The Tao/Jacobian Conjecture example from §2 also fits here: Tao "almost never takes the model's advice about where to go next" — he's the one asking "does X work here?" and "given Y and Z, why A?" A chem/plant-bio-specific analogue would strengthen this (a reader may not know who Tao is), but the Tao case works as a citable placeholder in the meantime.

**b. Recognising what matters in ambiguous data**
- Experienced scientists notice the anomaly, the unexpected peak, the result that doesn't fit — and they know when that anomaly is an artefact versus a discovery
- AI can flag statistical outliers but cannot distinguish between "interesting" and "broken instrument"
- The mass spec example from the decision framework notes is strong here

**c. Validation that goes beyond checking**
- From the decision framework: "the code running is not the same as the analysis being correct"
- Validation requires not just comparing output against a reference, but understanding whether the approach was appropriate in the first place
- This is the core of the "human in the loop must mean something" argument
- Reference: the Patterns paper on "systematic calibration" as irreducible intellectual labour

**d. Mentoring and knowledge transfer**
- Training the next generation of scientists requires a kind of knowledge that is embodied, contextual, and relational — how to read a gel, when to trust a measurement, how to present a negative result honestly
- This is irreducibly interpersonal. AI can answer questions about technique, but it cannot model scientific judgement in the way a mentor does through years of working alongside a trainee
- The time reclaimed by AI could go *towards* more mentoring — or it could erode the very experiences through which mentoring knowledge develops

**e. Ethical and professional judgement in context**
- Knowing when a result is ready to publish, when a reviewer's comment reveals a genuine flaw versus a misunderstanding, when a collaboration is productive versus extractive
- These are judgements about the social and professional fabric of science that require being embedded in it

### 4. The deskilling risk: what happens when expertise isn't practised

- The expertise erosion paradox (from Patterns paper): GenAI can perform the foundational activities through which expertise *develops* — reading, writing, designing, coding
- If trainees delegate these activities, they never build the judgement that makes validation possible
- Concrete evidence: the Lancet colonoscopy study (AI-assisted endoscopists performed worse without it); the ACM "deskilling paradox"; the DeepMind proposal for AI to assign humans tasks to prevent atrophy
- This is not hypothetical hand-wringing — it is a measurable phenomenon
- The implication for training: we need to be deliberate about which tasks students do themselves and which they delegate, and the decision should be based on which activities build expertise, not just which are fastest

### 5. But AI *does* have domain knowledge — so what changes?

- Take the question seriously rather than dismissing it
- AI's domain knowledge is genuine but partial: it lacks the ability to weigh, prioritise, and contextualise
- **[DRAFT — added 2026-08-04, reading-list session, from Goedecke "LLMs reward expertise"; review before use]** A useful rebuttal to "but the model already has the domain knowledge baked in": expertise isn't removed from the loop, it's relocated. Goedecke's example — OpenAI's models producing mathematical discoveries — still depends on a team of expert mathematicians checking and filtering the model's output before anything is trusted. The expertise didn't disappear; it moved from generation to validation. That's a clean, concrete instance of this post's own §3c argument (validation as irreducible labour), reinforcing rather than undercutting the thesis.
- What changes is the *distribution* of where human time goes. The mechanical parts of expertise (recall, formatting, first-pass organisation) are increasingly handled by AI. The judgement parts (what to pursue, what to trust, what to question) become the primary locus of human contribution
- This is a rebalancing, not a replacement. But it demands that training programmes and career development adapt — if we keep training scientists the same way, AI will handle the things they spent years learning to do and leave them unprepared for the things that now matter most
- The "centaur scientist" framing: the most effective researcher is one who brings deep domain judgement and uses AI to extend their reach, not one who uses AI as a substitute for the judgement they never developed

### 6. What this means for how we train scientists

- Practical implications for graduate education and mentoring
- The argument for *more* rigorous disciplinary training, not less — but with a shift in emphasis from rote execution toward judgement, experimental design, and critical evaluation
- The case for deliberate practice: some tasks should be done without AI specifically to build the expertise that makes AI use safe and effective later
- Connection back to the decision framework: the validateability axis only works if the person in the loop has genuine expertise. Investing in that expertise is not a luxury — it is a prerequisite for responsible AI use

### 7. Closing

- Not a summary — a forward-looking statement
- The question is not whether AI will continue to get better at domain-specific tasks (it will). The question is whether we continue to develop the human judgement that makes those capabilities useful rather than dangerous
- Expertise is not competing with AI. Expertise is what makes AI worth using.

---

## Interactive element (TBD)

Each post in the series has a companion interactive. Ideas for this one:
- A "validation challenge" — readers are shown AI-generated outputs (a methods section, a data interpretation, a literature summary) and asked to spot the subtle error that requires domain expertise to catch
- A decision tree: "What kind of expertise does this task require?" mapping different scientific tasks to the specific expertise capabilities from section 3
- Something else?

---

## Key references to draw on

- "Recalibrating academic expertise in the age of generative AI" (Cell/Patterns, 2025) — expertise erosion paradox, three meta-skills framework
- "AI cannot automate science" (The Conversation) — philosopher's argument for irreducibly human aspects of research
- "AI has supercharged scientists — but may have shrunk science" (Science/AAAS) — AI expanding individual impact while narrowing collective focus
- "The AI Deskilling Paradox" (Communications of the ACM)
- Deskilling in medicine (Lancet colonoscopy study; Springer mixed-method review)
- Google DeepMind (2026) — AI assigning tasks to humans to prevent skill atrophy
- "AI Inverts the Disciplinary Hierarchy" — humanities/judgement skills become more important when generation is cheap
- Stack Overflow "Domain expertise still wanted" (March 2026)
- Lubars & Tan (2019) — task delegability, trust as key factor
- "Delegation and Verification Under AI" (2026) — contract-first verification principle
- Segal & Rachitsky, "How tech workers are feeling in 2026: a workforce splitting in two" (Lenny's Newsletter, 7 Jul 2026) — second annual tech-worker sentiment survey (n≈5,332). Empirical, worker-voice evidence for two threads: **(§4 deskilling)** the comprehension-debt quotes — "my brain is rotting," "I just follow Claude, I don't fully understand what I merge," "I don't think hard enough anymore"; and **(§6 training / the on-ramp)** the disappearing-entry-level finding — NPS –39 for recommending one's own field, "the water's fine; don't come in," early-career ICs the most destabilised, and the leadership warning "don't let the bottom rung rot" (AI doing the entry-level work juniors learned on = starving next year's senior talent). Directly supports the student-facing argument that we must be deliberate about which tasks trainees do themselves. Full text: `to_read/z_done/workforce.txt`.
- Sean Goedecke, "LLMs reward expertise" (seangoedecke.com, 24 Jul 2026) — argues domain expertise, not generic prompting skill, is what lets people get more out of LLMs. Terence Tao's public ChatGPT conversation on the Jacobian Conjecture counterexample as the worked example (short expert-mode prompts, model shifts register, Tao drives direction rather than following the model). Supports **§2/§3a** (a concrete illustration of the expertise gap, candidate for the "[Could use an example...]" placeholder) and **§5** (the "expertise is relocated to validation, not removed" rebuttal to the "but AI already knows this" objection). Full text: `to_read/z_done/expertise.txt`.
- David Brooks, "The People Who Will Thrive in the AI Age" (*The Atlantic*, June 2026) — sorts users by "need for cognition" into three groups (low: erode fastest; medium: intend to resist, get pulled in anyway; high: actively use AI to increase rather than replace agency) and argues AI reveals volition, not intelligence, as the scarce differentiating trait. Gives concrete, citable empirical backing for **§4 (deskilling)** that the draft currently lacks: MIT Media Lab found brain connectivity drops up to 55% during ChatGPT use; a Wharton study found people accepted AI's deliberately-planted wrong answers 80% of the time; endoscopists' precancerous-lesion detection rate dropped from 28.4% to 22.4% after AI access was removed post-adoption. Strengthens the draft's existing "hollowing out the junior bundle" argument (`draft.md:83`) with hard numbers rather than a purely economic framing. Full text: `to_read/z_done/Atlantic.txt`.
- Decision framework post (cross-reference: validateability axis, human-in-the-loop, worked examples)

---

## Cross-references to other posts in the series

- Post 1 (decision framework): the validateability axis, "human in the loop has to mean something," worked examples showing expertise as the differentiator
- Post 3 (environmental impacts): brief nod — expertise also helps scientists make informed choices about *when* the computational cost of AI is justified
