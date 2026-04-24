# When to Use AI in the Natural Sciences — Framework Draft

A subchapter for the digital online textbook. Working notes and emerging framework.

---

## Ground rule

Anyone who uses AI to produce content is fully responsible for that content. This means there must be a human in the loop — a person cannot be responsible for something they are not meaningfully involved in. Everything below flows from this premise.

---

## Decision flow overview

The framework is structured as a decision flow, not a single evaluation step. Some nodes are universal gates (they apply regardless of context), others are conditional. Paths can diverge and reconverge.

```
┌─────────────────────────────────────┐
│  1. DATA CONFIDENTIALITY GATE       │  Universal — applies in all modes
│  Does this task involve sensitive,   │
│  confidential, or pre-publication   │
│  data?                              │
│                                     │
│  YES ──► Can you use a local model? │
│          YES ──► proceed ───────┐   │
│          NO ──► Can you          │   │
│            de-identify safely?   │   │
│            YES ──► proceed ──┐  │   │
│            NO ──► STOP       │  │   │
│                              │  │   │
│  NO ──► proceed ─────────────┤  │   │
│    (also sets tooling:       ▼  ▼   │
│     cloud / local / de-id)          │
└──────────────────────────┬──────────┘
                           │
                           ▼
┌─────────────────────────────────────┐
│  2. MODE: EXPLORATORY vs CONSIDERED │
│  Am I exploring/learning, or doing  │
│  professional work?                 │
│                                     │
│  EXPLORATORY ──► lighter evaluation │
│    (but responsibility still        │
│     applies — see node 5)           │
│                                     │
│  CONSIDERED ──► full evaluation     │
│    (nodes 3 + 4 below)             │
└──────────┬──────────┬───────────────┘
           │          │
     exploratory   considered
           │          │
           │          ▼
           │  ┌───────────────────────┐
           │  │  3. TASK              │
           │  │  DECOMPOSITION        │
           │  │  Break the task into  │
           │  │  subtasks. Evaluate   │
           │  │  each one separately. │
           │  └───────────┬───────────┘
           │              │
           │              ▼
           │  ┌───────────────────────┐
           │  │  4a. VALIDATEABILITY  │
           │  │  Can you efficiently  │
           │  │  verify the output?   │
           │  │                       │
           │  │  4b. AI CAPABILITY    │
           │  │  FIT                  │
           │  │  Is AI well-suited to │
           │  │  this task?           │
           │  │                       │
           │  │  (Go / reshape / no-  │
           │  │   go decision per     │
           │  │   subtask)            │
           │  └───────────┬───────────┘
           │              │
           ▼              ▼
┌─────────────────────────────────────┐
│  5. RESPONSIBILITY                  │  Universal — applies in all modes
│  Ethical implications, reproduci-   │
│  bility, transparency, perception,  │
│  professional consequences.         │
│  Everything passes through here.    │
└─────────────────────────────────────┘
```

**Note on node interactions:** Node 1 can constrain node 4b. If the confidentiality gate routes you to a local model, that model may be less capable than a cloud-hosted one, which affects AI capability fit downstream. This interaction should be considered when evaluating subtasks at node 4.

---

## Node 1: Data confidentiality gate

This is the first question, before anything else. It applies in both exploratory and considered modes — you don't upload confidential data to a third-party model regardless of your purpose.

- **Does the task involve sensitive, confidential, or pre-publication data?** This includes unpublished research data, student records (FERPA), health data (HIPAA), pre-publication findings, proprietary datasets, or anything the researcher wouldn't want ingested into a third-party training set.
- **If yes — can you use a local model?** Local models (running on your own hardware) keep data on your machine. Data never leaves. This is a viable path for many tasks, though local models may be less capable than cloud-hosted ones.
- **If not — can you de-identify sufficiently?** De-identification must be robust: even de-identified data can pose re-identification risks if the dataset contains enough variables that, combined with other publicly available data, could identify individuals.
- **If neither local models nor adequate de-identification are possible — stop.** AI is not appropriate for this task in its current form.
- **If no sensitive data is involved — proceed.**

The output of this node is not just "go/no-go" but also *which tool to use*: cloud-hosted AI, local model, or AI with de-identified data. This tooling choice carries forward through the rest of the flow.

---

## Node 2: Exploratory vs. considered mode

Two modes of AI use, which affect how rigorously the subsequent evaluation needs to be:

1. **Exploratory mode** — trying out various tasks, seeing how well AI handles them, learning what works and what doesn't. This is valuable for building intuition and is how most people start. A lot of genuine learning happens here.

2. **Considered mode** — deliberate, framework-guided use of AI in professional or research contexts. The shift from exploratory to considered use is a mark of maturing AI practice.

In exploratory mode, you might try things that are hard to validate or that AI isn't well suited for — because learning is the point. The decomposition and axis evaluation (nodes 3–4) can be lighter. But node 5 (responsibility) still applies fully: exploratory use doesn't exempt you from ethical obligations, transparency, or reproducibility concerns.

In considered mode, nodes 3 and 4 are where the serious evaluation happens.

---

## Node 3: Task decomposition

Before evaluating a task on the axes, break it into subtasks. This is a formal step, not an afterthought — multiple tasks only become good AI candidates once decomposed.

Most real-world tasks can be broken into subtasks, some of which are good AI candidates and some of which are not. Each subtask should be evaluated independently through node 4.

**Example: writing a recommendation letter**
- Researching the target institution's mission, programme details, and values — good AI subtask (text-based search, easy to verify against source)
- Drafting the letter itself — poor AI subtask (personal voice essential, generic output counterproductive)
- Formatting and proofreading the final letter — good AI subtask (mechanical, verifiable)

**Example: preparing a conference talk**
- Outlining the talk structure and drafting bullet points — good AI subtask (text-based, easy to review)
- Building the actual slides in PowerPoint/Keynote — poor AI subtask (requires GUI interaction)
- Writing speaker notes in the presenter's voice — poor AI subtask (personal voice, delivery style)

The decomposition step often reveals that the answer to "should I use AI for this?" is not yes or no, but "yes for these parts, no for those parts."

---

## Node 4: Evaluation axes (considered mode)

In considered mode, these two axes determine whether each subtask (from node 3) is a good candidate for AI assistance, or whether it needs to be reshaped first.

### 4a. Validateability

Can the human in the loop efficiently verify the AI's output?

- If validation is quick and reliable, AI offers genuine time and energy savings — the AI does the work, the human confirms it's correct, and the task is done.
- If validation is difficult or time-consuming, the benefit erodes. The time spent checking may exceed the time it would have taken to do the task directly.
- At the extreme, if validation requires essentially redoing the task, the AI contribution was pointless.

Factors that feed into this axis:
- **Domain expertise of the validator** — a PI reviewing an AI-generated summary in their own field can validate quickly and confidently. The same PI reviewing an AI summary of regulatory requirements in an unfamiliar jurisdiction cannot, even though both are "text summaries." This makes a strong argument for disciplinary-specific training: expertise is what makes validation meaningful rather than superficial.
- **Whether verifiable ground truth exists** — reformatting a citation to a specific style has a clear right answer, easy to check. Interpreting an ambiguous dataset has no single correct answer; "validation" becomes a judgement call, which is fundamentally harder.

### 4b. AI capability fit

Is the task within AI's current practical capabilities?

- AI excels at text-in, text-out tasks: drafting, summarising, reformatting, searching, coding.
- AI struggles with tasks that require navigating institutional web interfaces, managing credentials, working through complex multi-step browser workflows, or interacting with systems that have no API or text interface.
- This axis will shift over time as AI capabilities evolve, but at any given moment there is a clear boundary between what AI can practically access and what it cannot.

Factors that feed into this axis:
- **Whether personal voice is required** — a recommendation letter needs to sound like *you* and convey genuine personal knowledge of the candidate. AI-generated letters risk sounding generic, which can actively harm the candidate. This is a capability limitation: AI is poorly set up for authentic personal expression.
- **Tooling constraints from node 1** — if the data confidentiality gate routed you to a local model, that model may be less capable than a cloud-hosted alternative. A subtask that would be a good fit for a frontier cloud model might be a poor fit for a smaller local one. This interaction between nodes should be considered explicitly.

### Decision outcomes at node 4

- **Both axes favourable** — proceed to node 5 (responsibility check).
- **One or both axes unfavourable** — consider reshaping the subtask using strategies (see below) before proceeding, or decide that AI is not the right tool for this subtask.

---

## Node 5: Responsibility (universal)

This node applies to everything — exploratory and considered, cloud and local, regardless of how the task scored on validateability and capability fit. It is the floor that all AI use sits on.

### Ethical implications

- Ethics represents a category of cost that people tend to underweight because it is less tangible than a factual error. The cost of an ethical misstep — loss of trust from colleagues, violations of journal or institutional policy, damage to the broader credibility of the field — is not always immediate or obvious, but can be severe.
- Some ethical costs arise regardless of output quality. Using AI to draft peer review comments without disclosing it, or generating preliminary data descriptions that end up in a publication without proper attribution of AI's role, are integrity issues even if the content is factually correct.
- For readers evaluating responsibility: consider not just "what if the output is wrong?" but also "what are the ethical implications of using AI for this task at all, even if the output is right?"

### Reproducibility

- If AI is part of a research workflow that produces published results, can the process be reproduced? Prompts, model versions, and outputs are often ephemeral. This is particularly relevant to natural scientists who already operate within reproducibility norms.
- Irreproducible work has serious professional and scientific costs. Documenting AI use (prompts, model versions, outputs) is a practical step towards reproducibility.

### Human perception and the "AI-made" signal

- If a recipient perceives that a task was done by AI, they may react differently than if they believe it was done by a human. This is a real and potentially significant cost, distinct from whether the output is objectively correct.
- This argues for upfront transparency in AI use — but it also means that for some tasks, the cost of a mistake is genuinely hard to judge, because it depends on how others will perceive and interpret the output's provenance.
- The perception cost may be especially acute in contexts like recommendation letters, grant reviews, or any communication where the recipient expects authentic human engagement.

### Transparency

- Full responsibility means ethical responsibility, not just factual accuracy. Transparency about AI use is both an ethical obligation and a practical protection.

### Responsibility checklist for readers

When passing through this node, ask:
- Have I been transparent about AI's role in this work?
- Could my AI use be perceived as violating institutional or journal policies?
- If this involves published research, is the AI-assisted workflow documented and reproducible?
- Would the recipient of this output react differently if they knew AI was involved? If so, is that a cost I've accounted for?
- Am I comfortable defending this use of AI to a colleague, a reviewer, or an ethics board?

---

## Worked examples

These examples trace real natural science tasks through the decision flow to illustrate how the framework operates in practice.

### Example 1: Reformatting a reference list for journal submission

- **Node 1 (confidentiality):** No sensitive data — just published references. Proceed with any tool.
- **Node 2 (mode):** Professional work. Considered mode.
- **Node 3 (decomposition):** Task is already atomic — no decomposition needed.
- **Node 4a (validateability):** Very high. Clear right answer — either the citation matches the journal's style or it doesn't. Quick to check.
- **Node 4b (capability fit):** Excellent. Text-in, text-out. AI handles this well.
- **Node 5 (responsibility):** Minimal ethical concern. No perception issue — nobody expects you to format references by hand. Reproducibility is a non-issue.
- **Verdict: strong AI use case.** The flow handles this cleanly.

### Example 2: Writing a first draft of a thesis literature review

- **Node 1:** No sensitive data. Proceed.
- **Node 2:** Professional work. Considered mode.
- **Node 3 (decomposition):** This is where the framework adds real value. "Write me a literature review" decomposes into: (a) summarise individual papers, (b) identify themes and gaps across the literature, (c) write the narrative arc that positions the thesis.
- **Node 4a:** Subtask (a) — easy to validate if the student has read the papers. Subtask (b) — harder; requires domain knowledge to judge whether the right themes have been identified. Subtask (c) — hardest; the narrative is a scholarly argument, not just a summary.
- **Node 4b:** All subtasks are text-native, so capability fit is good. But AI may not have access to the most recent or niche papers — context grounding (providing the papers) helps.
- **Node 5:** Significant responsibility concerns. A thesis is the student's own scholarly contribution. If the literature review reads as AI-generated, examiners will notice. Institutional policies may apply.
- **Verdict: cautionary, but decomposition helps.** Subtask (a) with context grounding is a strong use case. Subtasks (b) and (c) require the student's own intellectual contribution. The flow correctly identifies which parts are and aren't suitable.
- **Key lesson:** A first-year student who hasn't read deeply in the field cannot validate even subtask (a) meaningfully. The validateability axis is doing real work here, and it's tied directly to the student's level of disciplinary expertise.

### Example 3: A postdoc analysing mass spectrometry data with unpublished compound identifications

- **Node 1 (confidentiality):** Pre-publication data with novel identifications. Sensitive. Can they use a local model? If yes, proceed with local model. If not, can they strip the novel identifications and work with just spectral patterns? Depends on the task. **The gate does real work here.**
- **Node 2:** Considered mode — this is for a paper.
- **Node 3 (decomposition):** Decomposes into: (a) convert proprietary instrument output to a clean table, (b) identify peaks, (c) compare against reference databases, (d) interpret novel findings.
- **Node 4a:** Subtask (a) — very easy to validate (does the table match the raw data?). Subtask (b) — harder; a wrong identification could be serious. Subtask (d) — requires deep expertise and is essentially a judgement call.
- **Node 4b:** Instrument output is often in proprietary formats. **AI capability fit may be poor unless the data has been converted first** — so subtask (a) may need to be done manually or with specialised software before AI can help with (b) or (c). Also, **node 1's tooling constraint matters here**: if routed to a local model, it may be less capable for the analytical subtasks than a frontier cloud model would be.
- **Node 5:** Reproducibility is very relevant — the analytical pipeline needs to be documented.
- **Verdict: the flow works, and it surfaces the node 1 → node 4b interaction.** The confidentiality constraint (local model) may limit AI capability, which needs to be factored in. Decomposition reveals that format conversion (a) is the enabler for everything else.

### Example 4: A grad student drafting an email to a potential PhD supervisor

- **Node 1:** No sensitive data. Proceed.
- **Node 2:** If the student is just exploring what AI produces — exploratory. If they're actually going to send it — considered.
- **Node 3 (decomposition):** Could decompose into: (a) research the supervisor's recent work and interests, (b) draft the email itself. But many students will treat this as a single task.
- **Node 4a:** The student can judge whether the email sounds right, but "right" here is subjective. Validation is possible but soft.
- **Node 4b:** Text task, AI can do it technically.
- **Node 5: This is where the flow really earns its keep.** The perception risk is high. If the supervisor perceives the email as AI-generated (and experienced academics often can), it actively undermines the student's goal. An email that's supposed to demonstrate genuine interest and initiative reads as the opposite if it feels templated. The ethical question is real even if the content is accurate.
- **Verdict: the responsibility node catches something the other nodes don't.** The task is easy for AI, moderately easy to validate, but the perception and ethical cost is the real issue. Good evidence that the responsibility node needs to be universal and taken seriously even when the technical evaluation looks favourable.

### Example 5: Writing R code for a mixed-effects model the student hasn't used before

- **Node 1:** Depends on the dataset. Ecological field data — probably fine. Clinical data — gate applies.
- **Node 2:** Considered mode if this is for a paper.
- **Node 3 (decomposition):** Could decompose into: (a) write the model code, (b) write diagnostic checks, (c) interpret the output. But many students will treat this as "write me the code."
- **Node 4a: This is the critical node.** The student can check whether the code *runs*, but they cannot validate whether the model specification is *correct* — whether the random effects structure is appropriate, whether the assumptions are met, whether the results mean what they appear to mean. **The code running is not the same as the analysis being right.** A student who has studied mixed-effects models can validate; one who hasn't, cannot.
- **Node 4b:** Good — AI is strong at code generation.
- **Node 5:** High stakes if a wrong model specification leads to incorrect conclusions in a publication. Reproducibility matters — was the prompt documented? Was the model choice justified?
- **Verdict: the flow correctly identifies this as risky despite AI being technically good at the task.** The validateability axis does the heavy lifting. This is a strong example for the disciplinary expertise argument — the student's statistical training (or lack thereof) is what determines whether this is a force multiplier or a liability.

### Example 6: A postdoc preparing slides for a conference talk

- **Node 1:** May contain unpublished results. Gate applies — check before uploading figures or data to a cloud AI.
- **Node 2:** Considered mode.
- **Node 3 (decomposition):** Decomposes into: (a) outline the talk structure, (b) draft text/bullet points for each slide, (c) build the actual slides in PowerPoint/Keynote, (d) write speaker notes, (e) create or refine figures.
- **Node 4a:** Subtasks (a) and (b) — easy to validate (does the structure make sense? do the points say what I mean?). Subtask (d) — harder to validate because it involves voice and delivery style.
- **Node 4b:** Subtasks (a) and (b) — excellent fit, text-native. Subtask (c) — poor fit, requires GUI interaction. Subtask (d) — poor fit if personal voice matters for delivery. Subtask (e) — depends on the tool; AI can draft code for figures (matplotlib, ggplot) but can't directly edit slides.
- **Node 5:** Voice matters for a talk — it should sound like the presenter. Perception risk if the talk feels generic. But the structural and content-drafting subtasks carry little risk.
- **Verdict: another strong decomposition example.** The flow reveals that the same overall task contains subtasks ranging from excellent to poor AI candidates. Without decomposition at node 3, the whole task would get a muddled evaluation.

### Example 7: A faculty member writing an NSF grant proposal (project description)

- **Node 1 (confidentiality):** Preliminary data may be unpublished. The research idea itself is sensitive — the PI may not want their unfunded proposal concept ingested into a training set that a competitor could later query. Gate applies. Local model or very careful scoping of what's shared with cloud AI.
- **Node 2:** Unambiguously considered.
- **Node 3 (decomposition):** Decomposes into: (a) literature synthesis to establish the gap, (b) articulating the central hypothesis, (c) describing the research plan and methods, (d) writing the broader impacts section, (e) drafting the budget justification, (f) ensuring compliance with the solicitation's formatting and content requirements.
- **Node 4a:** Subtask (a) — validateable if the PI knows the field well, especially with context grounding (provide the key papers). Subtask (b) — this is the PI's core intellectual contribution; even if AI drafts something, the PI must judge whether it captures their actual idea — hard to validate in a meaningful sense because it's not about correctness but about *whether this is what I mean*. Subtask (c) — the PI can validate methods they've designed, but AI might suggest approaches the PI hasn't vetted. Subtask (e) — straightforward, clear right answers. Subtask (f) — very validateable, mechanical compliance check.
- **Node 4b:** Text-native throughout, so capability fit is generally good. But solicitation requirements may be on a website the PI needs to navigate and extract manually first.
- **Node 5:** High responsibility. Grant proposals carry the PI's name and reputation. Reviewers are experienced scientists who may detect AI-generated prose. The central hypothesis and research plan must be genuinely the PI's intellectual contribution — there are real integrity questions if these are AI-generated. NSF and NIH are increasingly requiring disclosure of AI use.
- **Verdict: strongly decomposition-dependent.** Subtasks (e) and (f) are excellent AI candidates. Subtask (a) is good with context grounding. Subtasks (b) and (c) are where the PI's expertise is irreplaceable — AI might help organise thoughts, but the ideas must be the PI's. The responsibility node adds weight: reviewers notice, and disclosure requirements are tightening.

### Example 8: A faculty member peer reviewing a journal manuscript

- **Node 1 (confidentiality):** The manuscript is confidential — it's unpublished, shared under the expectation of reviewer confidentiality. **Uploading it to a cloud AI arguably violates that confidentiality.** Many journals now explicitly prohibit this. Local model or nothing.
- **Node 2:** Considered.
- **Node 3 (decomposition):** Decomposes into: (a) summarising the paper's claims and methods, (b) checking statistical and methodological rigour, (c) evaluating novelty and significance relative to the field, (d) writing constructive feedback, (e) checking references and formatting compliance.
- **Node 4a:** Subtask (a) — easy to validate if the reviewer has read the paper. Subtask (b) — depends entirely on the reviewer's statistical expertise; if they can't validate AI's assessment, it's dangerous. Subtask (c) — a pure judgement call requiring deep domain knowledge. Subtask (e) — mechanical, easy to check.
- **Node 4b:** Text-native, but constrained to local models by node 1. Less capable model may struggle with nuanced methodological assessment.
- **Node 5: This is perhaps the strongest responsibility example in the entire set.** Peer review is built on the premise that an expert human is evaluating the work. AI-generated reviews undermine the entire system, even if the content happens to be accurate. Many journals now explicitly prohibit AI use in peer review. The ethical cost is not hypothetical — it's an integrity violation with potential professional consequences.
- **Verdict: poor AI use case for most subtasks, and the responsibility node is decisive.** Subtask (e) (reference/formatting checks) is defensible. Using AI for the substantive review is ethically problematic regardless of output quality. The confidentiality gate and the responsibility node both flag serious issues independently.

### Example 9: A faculty member responding to reviewer comments on their own manuscript

- **Node 1:** The manuscript is the PI's own work, but it's under review and unpublished. Moderate sensitivity — the PI may be less concerned about confidentiality than with someone else's manuscript, but pre-publication data is still involved.
- **Node 2:** Considered.
- **Node 3 (decomposition):** Decomposes into: (a) organising and categorising reviewer comments, (b) drafting point-by-point responses, (c) revising the manuscript text, (d) running any additional analyses reviewers requested, (e) writing a cover letter to the editor.
- **Node 4a:** Subtask (a) — very easy to validate (it's just reorganisation). Subtask (b) — the PI can validate responses in their own field readily, and context grounding works well here (provide the original manuscript and the reviews, ask AI to draft responses). Subtask (c) — validateable but requires careful reading. Subtask (d) — depends on the analysis; same issues as the mixed-effects model example. Subtask (e) — straightforward.
- **Node 4b:** Text-native throughout. Good fit.
- **Node 5:** Moderate responsibility concerns. The responses represent the PI's scientific engagement with the reviewers. Editors and reviewers can sometimes detect generic or AI-flavoured responses. Transparency norms are evolving — some journals now ask about AI use in revisions.
- **Verdict: generally a reasonable AI use case, especially with decomposition and context grounding.** Subtask (a) is a clear win. Subtask (b) with the original manuscript and reviews provided as context is a strong use case — the PI can validate quickly against their own expertise. The responsibility concerns are real but less severe than peer review, because the PI is engaging with their own work.
- **Key contrast with Example 8:** Both involve reviewer–author interactions, but the ethical landscape is very different. Reviewing someone else's work with AI is a confidentiality and integrity violation. Using AI to help respond to reviews of your own work is a productivity tool — provided the responses are genuine and the PI validates them.

### Example 10: A faculty member preparing an IRB application for human subjects research

- **Node 1:** The application itself may describe sensitive populations, procedures, and data collection plans. Depending on the study, it could involve identifiable information about participants even at the protocol stage. Gate applies — consider carefully what's shared.
- **Node 2:** Considered.
- **Node 3 (decomposition):** Decomposes into: (a) drafting standard boilerplate sections (data storage plan, consent withdrawal procedures, risk mitigation for common scenarios), (b) describing the specific study procedures, (c) drafting consent forms, (d) navigating the institutional IRB submission portal and filling in forms, (e) ensuring compliance with federal regulations (Common Rule, HIPAA where applicable).
- **Node 4a:** Subtask (a) — highly validateable; boilerplate has established templates and the PI knows the standards. Subtask (b) — validateable because the PI designed the study. Subtask (c) — validateable but requires care; consent language must be precise and legally adequate. Subtask (d) — can't be validated by AI because AI can't do it (see 4b). Subtask (e) — harder to validate unless the PI has regulatory expertise; errors here have real consequences.
- **Node 4b:** Subtasks (a), (b), (c) — text-native, good fit. Subtask (d) — **poor fit**; IRB portals are typically complex institutional web interfaces with authentication, multi-step workflows, dropdown menus, and conditional logic. AI can't navigate these. Subtask (e) — text-native but may require accessing current regulatory text from government websites.
- **Node 5:** High stakes. IRB errors can delay research by months, and inadequate consent procedures have serious ethical and legal consequences. Regulatory compliance is not an area where "close enough" is acceptable.
- **Verdict: mixed, and a strong illustration of the capability fit axis.** The boilerplate drafting subtasks are good AI candidates — high validateability, text-native, and they save genuine time on tedious but important writing. The portal navigation is a clear capability barrier. The regulatory compliance subtask is where the PI needs to be honest about their own expertise: if they can validate the output, AI helps; if they can't, it's risky. This example also shows that even "boring" administrative tasks have real stakes.

### Example 11: A faculty member writing an annual progress report for a funded grant

- **Node 1:** Contains details about ongoing, unpublished research. Moderately sensitive — the funder expects a report on work in progress. The PI may not want this in a training set, but the sensitivity is lower than for a competing proposal.
- **Node 2:** Considered.
- **Node 3 (decomposition):** Decomposes into: (a) summarising research progress against the original aims, (b) listing publications, presentations, and personnel changes, (c) describing changes to the research plan, (d) writing the next-year plan, (e) compiling budget expenditure summaries.
- **Node 4a:** Subtask (a) — highly validateable if the PI provides their own notes, lab meeting slides, or manuscript drafts as context. Classic context-grounding use case. Subtask (b) — very easy to validate; factual list with clear right answers. Subtask (c) — validateable because the PI knows what changed and why. Subtask (d) — similar to grant writing, the PI must judge whether it reflects their actual plans. Subtask (e) — factual, verifiable against financial records.
- **Node 4b:** Text-native throughout. Good fit. Budget data may need to be exported from an institutional finance system first (capability barrier for extraction, but once exported, AI handles it well).
- **Node 5:** Moderate responsibility. Progress reports are official documents to funders. Accuracy matters, but the tone is more factual/administrative than a grant proposal. Less perception risk than a proposal — the funder is not evaluating the PI's intellectual creativity, they're checking whether the money was well spent.
- **Verdict: strong AI use case, especially with context grounding.** The PI has all the source material (their own notes, publications, financial data). Most subtasks are factual, validateable, and text-native. The main time savings come from subtask (a) — turning scattered notes and results into a coherent narrative, which AI does well when given the right inputs. This is the kind of task where AI demonstrably saves hours of tedious but important writing, with low risk if the PI reviews the output.

### Structural observations from the examples

**On the decision flow structure:**

1. **Task decomposition as a formal node is justified.** Multiple examples (literature review, conference talk, grant proposal, IRB application, mass spec analysis) only work well when decomposed first. Without node 3, users would evaluate composite tasks and get misleading results.

2. **Node 1 and node 4b interact.** The confidentiality gate can constrain you to a local model, which may be less capable, affecting AI capability fit downstream. The mass spec example (3) and peer review example (8) show this clearly. This interaction is worth flagging explicitly for readers.

3. **The ordering holds up well across all eleven examples.** Confidentiality first (hard constraint), mode second (sets expectations), decomposition third (reveals the real structure), axes fourth (evaluate each piece), responsibility last (universal floor). No example suggested a different ordering would be better.

**On the individual nodes:**

4. **The responsibility node is doing substantial work.** The email example (4) shows a task that passes the technical evaluation but fails on perception. The peer review example (8) is the strongest case: responsibility is decisive on its own, regardless of how the other nodes evaluate. These validate the decision to make responsibility universal.

5. **Validateability is often the decisive axis.** The mixed-effects model example (5) shows that AI capability fit can be excellent while validateability is poor — and validateability is what determines whether the use case is safe. The grant proposal example (7) shows a similar dynamic: AI can generate plausible-sounding hypothesis text, but only the PI can judge whether it captures their actual idea.

6. **Confidentiality is more complex for faculty than for students.** Grad students mostly worry about their own data. Faculty deal with other people's confidential work (peer review), sensitive proposals, regulatory documents, and institutional data. Node 1 does more work for faculty tasks, and the decisions are less straightforward.

**On pedagogy and presentation:**

7. **The contrast between Examples 8 and 9 is pedagogically valuable.** Two superficially similar tasks (both involve reviewer–author interactions on manuscripts) land in completely different places on the framework. Peer reviewing with AI is an integrity violation; responding to reviews of your own work with AI is a productivity tool. This kind of contrast helps readers develop judgement rather than just following rules.

8. **The peer review example (8) should be highlighted as a paradigmatic case** — the clearest example of a task where the responsibility node is decisive on its own and the technical evaluation is almost irrelevant.

9. **Context grounding shines for faculty tasks.** Progress reports (11), reviewer responses (9), and grant literature sections (7a) are all cases where the PI has the source material and AI's role is extraction and synthesis. This is the framework's sweet spot and a pattern worth emphasising.

10. **Even "boring" administrative tasks have real stakes.** The IRB example (10) shows that tedious boilerplate writing can still involve serious consequences if done poorly. The framework helps readers take these tasks seriously without over-complicating the ones that are straightforward.

---

## The deeper argument: disciplinary expertise as force multiplier

Disciplinary expertise is not just useful for validation — it is what makes the entire human-in-the-loop concept substantive rather than performative.

A well-trained scientist using AI for tasks they can validate is a force multiplier. A poorly-trained person using AI for tasks they can't validate is a liability — and no amount of AI capability fixes that, because the bottleneck was never the generation, it was the judgement.

This leads to a key insight: investing in disciplinary training becomes *more* important in an AI-assisted world, not less. Expertise is what separates responsible, effective AI use from the generation of plausible-looking output that no one can properly evaluate.

Furthermore, if AI handles the mechanical parts of a workflow, the time reclaimed can go towards activities where disciplinary expertise truly shines and that are irreducibly human: conversations with collaborators, mentoring students, sitting with a puzzling result and thinking it through. These person-to-person, expertise-driven activities are where science actually advances — and freeing up time for them is one of the strongest arguments for thoughtful AI adoption.

---

## Strategies for improving a task's position on the axes

The framework is not just diagnostic — it also points towards practical strategies for reshaping tasks before handing them to AI. These strategies are most relevant when a subtask scores poorly at node 4 but could be reshaped to score better.

### Context grounding: converting hard-to-validate tasks into easy ones (node 4a)

In the natural sciences, factual accuracy matters. When we care about correctness, relying on AI's training knowledge alone is risky — the output may be plausible but wrong, and validating it requires independent fact-checking, which is slow and uncertain.

The strategy: provide AI with source material that contains the correct information, and ask it to extract and process from that material rather than generate from training knowledge. Even if the correct answer is buried in the context, fragmented across documents, or embedded in a long paper, AI is good at pulling those pieces together.

This converts the validation task from *independent fact-checking* (hard, uncertain) into *comparison against source material* (fast, reliable). The human can check whether the output faithfully represents the sources they already trust, rather than having to verify claims of unknown provenance.

In exploratory mode, generating from training knowledge is fine — that's part of building intuition. But in considered mode, context grounding is a key technique for keeping tasks in the "easy to validate" zone.

### Task decomposition as a validation strategy (node 4a)

Beyond its role as a formal node in the flow, decomposition is also a strategy for improving validateability: a task that is hard to validate as a whole may consist of subtasks that are each straightforward to check individually.

### Strategies for improving AI capability fit (node 4b)

Moving a task along the AI capability fit axis means making it accessible to AI's strengths (text-in, text-out processing).

**Document conversion** — taking information locked in PDFs, proprietary formats, or web interfaces and converting it to plain text. This removes the barrier between the task and AI's sweet spot. In the natural sciences context, this includes instrument output files, lab notebooks, and institutional documents.

**API access vs. browser workflows** — if an institutional system has an API (even a clunky one), interactions can be scripted and results fed to AI as text. If not, you're stuck with the browser. Knowing whether an API exists is itself a useful skill.

**Structured data extraction** — converting messy or semi-structured data (spreadsheets with inconsistent formatting, instrument outputs, handwritten notes) into something clean and text-based. Once structured, AI can work with it readily.

**Prompt engineering and task specification** — a vague request to AI is a poor capability fit not because AI *can't* do it, but because the task hasn't been made AI-friendly. Learning to write precise, well-scoped instructions is a skill that moves tasks along this axis.

**Workflow segmentation** — separating the parts of a workflow that require human interaction with a system (logging in, clicking through menus, making selections) from the parts that are purely information processing. The human does the clicks, exports the data, and hands the processing to AI.

**The investment-to-return trap:** building tools for document conversion and workflow automation can be time-consuming and tempting — it *seems* like it would make whole categories of tasks AI-friendly. But in practice, the tools may not get reused as often as expected, or may not work as well as hoped. The person may learn lessons along the way, but this activity falls more naturally into the exploratory category than into the considered framework, where the goal is not to explore but to accomplish tasks in a professional context. Before investing heavily in conversion infrastructure, ask: will this tool genuinely pay off repeatedly, or is this exploration dressed up as productivity?

### Strategies across the flow

These strategies share a common logic: before handing a task to AI, reshape it so that the output will be easier to verify, better suited to AI's capabilities, or lower in risk. The decision flow diagnoses the problem; the strategies address it.

---

## Related work and literature context

### Task delegability frameworks

- **Lubars & Tan (2019), "Ask not what AI can do, but what AI should do: Towards a framework of task delegability"** — proposes four factors for deciding whether to delegate to AI: motivation, difficulty, risk, and trust. Trust (which maps loosely onto our validateability axis) correlated most strongly with people's actual delegation preferences. A follow-up study evaluated this framework specifically in knowledge work contexts. ([paper](https://ar5iv.labs.arxiv.org/html/1902.03245); [knowledge work evaluation](https://www.researchgate.net/publication/355146552_Task_Delegability_to_AI_Evaluation_of_a_Framework_in_a_Knowledge_Work_Context))

- **"Delegation and Verification Under AI" (2026)** — formalises a "Contract-First" principle: you should not delegate a task unless the outcome can be precisely verified. If a task is too subjective to verify, you recursively decompose it into subtasks that *can* be verified. This closely parallels our task decomposition insight and the centrality of validation. ([paper](https://arxiv.org/html/2603.02961))

### Human-in-the-loop as more than quality control

- **CUNY Graduate Center (2025), "Human in the Loop: Interpretive and Participatory AI in Research"** — argues that human input is not just quality control but the "interpretive spine" of the research process. Connects to our argument that disciplinary expertise is what makes the human-in-the-loop concept substantive rather than performative. ([post](https://gcdi.commons.gc.cuny.edu/2025/10/10/human-in-the-loop-interpretive-and-participatory-ai-in-research/))

- **"Beyond human-in-the-loop: Sensemaking between AI and human intelligence collaboration" (2025)** — moves beyond static HITL models towards conceptualising AI-human collaboration as a sociotechnical system. ([paper](https://www.sciencedirect.com/science/article/pii/S2666188825007166))

### Skill atrophy and the expertise argument

- **Google DeepMind (2026)** — proposes that AI should occasionally assign humans tasks so they don't lose the skills needed to validate AI output. Directly supports our argument that disciplinary expertise becomes *more* important with AI adoption: if researchers stop doing tasks entirely, they lose the ability to judge whether those tasks were done well. ([coverage](https://creati.ai/ai-news/2026-02-17/google-deepmind-ai-agent-delegation-framework/))

### Evaluation frameworks for AI tools in research

- **Six-tiered framework for evaluating AI models (2025)** — evaluates AI across repeatability, reproducibility, robustness, rigidity, reusability, and replaceability. More focused on model evaluation than task suitability, but the rigour dimensions are relevant. ([paper](https://www.sciencedirect.com/science/article/abs/pii/S0167779925002768))

### Data confidentiality and privacy

- **Virginia Tech guidance on AI in research** — only low-risk data may be shared with generative AI tools absent vendor agreements; data from human subjects entered into third-party AI may become publicly accessible. ([guidance](https://www.research.vt.edu/research-support/forms-guidance/sirc/guidance-using-artificial-intelligence-during-research-activities.html))
- **Re-identification risks** — even de-identified datasets can enable re-identification when combined with other data sources. ([paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC12803850/))

### Where our framework fits

Most existing frameworks are either very general (any knowledge work) or very specific (clinical decision support, software engineering). A decision-flow framework specifically grounded in natural science workflows — with concrete examples from bench science, ecology, chemistry — would fill a genuine gap. The sequential gating structure (data confidentiality → mode → decomposition → axes → responsibility), the formal decomposition step, and the disciplinary expertise argument are distinctive contributions.

---

## Still to develop

- Visual design of the decision flow for the textbook (ASCII diagram above is structural; needs a proper figure)
- How exploratory mode interacts with each node in more detail — responsibility clearly applies, but how much of validateability and capability fit should exploratory users think about?
- More examples from specific natural science disciplines (ecology fieldwork, organic synthesis, genomics pipelines)
- Guidance for readers on applying the framework to their own tasks (worksheet or template?)
- Connection to the broader course competencies (critical evaluation, integration and workflow design)
- Engage more deeply with Lubars & Tan's four-factor model — where does our framework overlap and where does it diverge?
- Consider the Contract-First recursive decomposition principle as a formalisation of our task decomposition node
- The DeepMind skill atrophy argument as supporting evidence for the disciplinary expertise thesis
- Whether the ternary diagram still works for visualising nodes 4a + 4b, or whether a different visual is needed now that they sit within a larger flow
- **Publish as a blog article** on the lab website (bustalab.org or similar) — adapt the working notes into a polished post, potentially as a subpage/blog entry. This could also serve as a public-facing companion to the textbook subchapter.
