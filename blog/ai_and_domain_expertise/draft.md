# Disciplinary Expertise in an AI-Assisted World

*Why deep training matters more now, not less*

---

One of the questions I hear most often — from students, colleagues, and audiences after talks — is some version of this: "If AI already knows my field, why do I need to?" It takes different forms. A graduate student wonders whether investing years in mastering statistical methods still makes sense when AI can generate the code in seconds. A postdoc asks whether learning to read mass spectra by eye is a useful skill or an anachronism. A department chair, planning curriculum revisions, wonders how much disciplinary content can be trimmed now that students can look anything up.

These are honest questions, and they deserve honest answers. This post is my attempt at one.

The short version: disciplinary expertise becomes more important with AI, not less. But saying that is not enough — the claim needs to be unpacked, because AI's capabilities are real and growing, and a dismissive "you still need to know your field" does not actually explain *what* expertise does that AI cannot, or *why* the gap matters. Those are the questions I want to work through here.

This is a companion to my earlier post on [when to use AI in the natural sciences](../ai_decision_flow/), which develops a decision framework for evaluating specific tasks. That framework depends on a concept it names but does not fully explore: the idea that disciplinary expertise is a "force multiplier" for AI use. This post is the fuller exploration.

<!-- MARGIN: Force multiplier — a term borrowed from military strategy, used here to mean something that amplifies the effectiveness of whatever it is combined with. Disciplinary expertise multiplied by AI produces more than either alone. AI without expertise produces output that no one can properly evaluate. -->

## What expertise actually is

To understand why expertise still matters, it helps to be precise about what it is. Expertise is not a body of facts. It is the ability to make judgements under uncertainty, recognise what is surprising or anomalous in data, and know which details matter and which can be safely ignored. An experienced organic chemist does not just know reaction mechanisms — she recognises when a yield is unexpectedly low and can generate plausible explanations. An ecologist with twenty years of fieldwork does not just know species names — he notices when a community composition suggests something has changed about the site.

<!-- MARGIN: Judgement under uncertainty — the ability to make informed decisions when the available information is incomplete, ambiguous, or contradictory. This is distinct from following a procedure or applying a rule. Most consequential scientific decisions — what to investigate, what to trust, what to publish — involve uncertainty. -->

AI has a form of knowledge, and it would be dishonest to dismiss it. Frontier models have been trained on vast bodies of scientific literature, and they can produce competent-looking summaries, generate reasonable code, and identify patterns across large datasets. In experimental tests, AI systems have demonstrated the ability to answer domain-specific questions, draft methods sections, and even propose experimental designs.

But AI's knowledge is statistical. It can tell you what is typical in a body of literature — the common approaches, the standard methods, the established conclusions. What it struggles with is telling you what is *important*: which of those common approaches has a known limitation that matters for your specific question, which established conclusion has been quietly challenged by recent work that has not yet propagated through the field, which unexpected observation in your data is an artefact of your instrument and which might be the most interesting result in the dataset.

That distinction — between typical and important — is what expertise provides. And it is not a distinction that more training data resolves, because importance is not a property of the text. It is a judgement that depends on knowing the field's open questions, its history of dead ends, and its current tensions.

## Part one: what happens to the individual

The relationship between AI and individual expertise has two sides, and they pull in opposite directions.

On one side, AI amplifies what an expert can do. A scientist with deep domain knowledge who uses AI for tasks she can verify — reformatting, first-pass literature organisation, routine code generation, data wrangling — becomes genuinely more productive. She reclaims time that can go toward the parts of the work where her expertise matters most: designing experiments, interpreting results, mentoring students. This is the force multiplier effect, and it is real.

On the other side, AI can erode or prevent the development of the very expertise it needs in order to be a force multiplier. This is what several researchers have begun calling the "expertise erosion paradox," and it deserves careful attention.

### The force multiplier

When a well-trained scientist delegates a task to AI, the delegation works because she can verify the output. She knows what correct looks like. If the AI-generated code runs but specifies the wrong model, she catches it. If the literature summary omits a key paper, she notices. If the suggested analysis is technically valid but inappropriate for the research question, she recognises the mismatch.

In the [companion post on when to use AI](../ai_decision_flow/), I call this the "validateability" axis — the question of whether you can efficiently verify what the AI produces. The axis only works if the person in the loop has genuine expertise. When they do, AI is a straightforward productivity gain: it handles the mechanical execution, the expert provides the judgement, and the combination outperforms either alone.

This is the centaur model. In centaur chess, human-AI teams consistently outperformed both humans alone and AI alone, because the human provided strategic judgement while the computer provided tactical calculation. The analogy holds well for scientific research: AI provides reach, but the researcher provides direction.

<!-- MARGIN: Centaur scientist — a concept from MIT's 2026 workshop on AI in the mathematical and physical sciences. In centaur chess, a human-AI team consistently outperformed both humans alone and AI alone. The analogy extends to research: AI provides computational and analytical reach, but the researcher provides the strategic judgement that decides where to point it. -->

### The erosion paradox

But here is the problem. The tasks that AI handles well — reading and summarising literature, writing first drafts, generating and debugging code, designing routine analyses — are not just mechanical chores. They are also the activities through which scientific expertise is *built*. A student who wrestles with a statistical model, gets it wrong, debugs it, and finally understands what the parameters mean has learned something that a student who asks AI for the code and runs it without error has not.

Lin and Sohail, writing in *Patterns* in 2026, articulate this clearly: "Unlike previous tools that augmented specific capabilities, GenAI can perform foundational intellectual activities through which scientific expertise itself develops: reading literature, designing studies, writing arguments, and generating code." Previous tools — microscopes, statistical software, sequencing machines — extended what a scientist could *do*. GenAI can substitute for what a scientist needs to *learn*.

Cotton and Scholle-Cotton, writing in *Issues in Science and Technology* in 2026, put the implication sharply: "The early-career professionals whose output benefits most from AI today may be the least prepared to lead their fields in an AI-driven future." They call this *short-circuiting the apprenticeship* — and the analogy is apt, because apprenticeship is built on productive struggle, on wrestling with material that resists you, on building understanding through effort rather than receiving answers from a tool.

<!-- MARGIN: Short-circuiting the apprenticeship — a phrase from Cotton & Scholle-Cotton (2026) describing how AI allows trainees to skip the difficult, time-consuming tasks that historically built expertise. The output is polished, but the person has not developed the judgement that would allow them to produce or evaluate such output independently. -->

The evidence extends beyond commentary. A 2025 study in *The Lancet Gastroenterology & Hepatology* found that endoscopists who routinely used AI assistance for colonoscopy performed measurably worse when the AI was removed — their precancerous lesion detection rate dropped from 28.4% to 22.4%. The skill did not simply fail to develop; it actively degraded through disuse. Skills fade not because they become unnecessary, but because they are no longer practised.

### The two sides together

This is the core tension at the individual level: AI both requires and erodes expertise. It amplifies what you know, but it can prevent you from developing that knowledge in the first place — or degrade knowledge you already had. The implication is that AI use cannot be all-or-nothing. It has to be deliberate: some tasks delegated, others done by hand specifically because the doing is what builds the capacity to judge.

Krishnan, writing in 2026, argues that PhD training should establish foundational mastery benchmarks *before* AI tools are introduced for core tasks. The logic is straightforward: you cannot validate what you have not first learned to do yourself. A student who has never manually searched and read a body of literature cannot meaningfully judge whether AI's literature summary missed something important. The apprenticeship has to come before the delegation.

<!-- MARGIN: Deliberate sequencing — the practice of structuring training so that foundational skills are built through direct practice before AI tools are used for those same tasks. The aim is not to reject AI but to ensure that students develop the expertise needed to use it responsibly. -->

A related finding offers a useful framing for how this balance plays out in practice. A 2026 study found that *fluent* AI users — those with deep expertise and comfort with the tools — actually experience more failures than novices. But their failures are visible, recoverable, and occur alongside greater success on complex tasks. The interpretation: expertise lets you know when AI has failed, which may be the single most important meta-competency in AI-assisted work. Novices do not fail less; they fail without noticing.

## Part two: what happens across a field

The individual-level picture — AI amplifies expertise but can erode it — is important, but it is not the whole story. When many researchers adopt AI simultaneously, patterns emerge at the population level that no single individual intended or controls. These are emergent effects, and several have now been documented.

### The steering problem

A large-scale study published in *Nature* in 2025, analysing 41.3 million research papers, found a striking pattern: AI-using scientists published 67% more papers and received roughly three times more citations than non-users. At the individual level, this looks like a clear win. But across 70% of sub-fields studied, AI collectively *contracted* the range of topics being investigated. Papers clustered around computationally tractable, data-rich problems. Follow-on engagement between papers dropped by 24%.

In other words, AI was making individual scientists more productive while making science as a whole less diverse. The researchers describe a "star-like structure" — many papers orbiting the same hot topics, high output but low engagement between distinct lines of work.

This is what happens when productivity is decoupled from the judgement that decides what is worth being productive about. AI makes certain kinds of work so much easier that the path of least resistance shifts — not just for individuals, but for entire fields. Without enough scientists who can recognise what *matters* as distinct from what is data-rich and statistically tractable, research drifts toward local maxima: problems that are easy to work on rather than problems that are important to solve.

<!-- MARGIN: Local maxima — a term from optimisation theory, used here metaphorically. A local maximum is a solution that looks optimal in its immediate neighbourhood but is not the best solution overall. In the research context, it describes a problem that is easy to make progress on (because it is data-rich and computationally tractable) but that may not be the most important question in the field. Expertise is what allows a scientist to look beyond the local maximum and identify where the genuinely important questions lie. -->

The steering problem is not a deskilling problem — it can happen even among experts. But it is amplified when the deskilling dynamic reduces the number of researchers with the deep judgement needed to resist the pull of tractable problems. The two levels interact: individual erosion of expertise feeds the collective narrowing of science.

### The complacency cascade

Lin and Sohail document a pattern they call the "complacency cascade": users initially scrutinise AI outputs carefully, then gradually accept them with less critical evaluation, particularly for technical content where detecting subtle errors requires concentrated attention. The trajectory is not ignorance but drift — a slow erosion of habits that were once automatic.

At the individual level, complacency is a personal risk. At the population level, it becomes a systemic one. If most researchers in a field gradually stop scrutinising AI output with full attention, the collective quality of scientific work degrades — not through any single dramatic failure, but through a thousand small acceptances of output that was almost right.

### Distorted signals

If career advancement metrics reward publication output, and if AI inflates that output without corresponding growth in expertise, then hiring committees and promotion panels are evaluating a signal that has been systematically distorted. Cotton and Scholle-Cotton warn that this could create leadership gaps in knowledge fields — researchers who rose quickly on AI-assisted productivity but lack the deep understanding to guide a research programme, mentor students through genuine difficulty, or recognise when a promising-looking line of work is headed for a dead end.

This is an institutional problem, not an individual one. Even a scientist who personally manages the balance between AI use and skill development is operating in a system where the signals used to evaluate everyone have been altered. The metrics that once served as rough proxies for expertise — publications, citations, grant productivity — become less reliable as AI decouples output from understanding.

## Part three: what AI cannot reach

The arguments above are about what AI *does* to expertise — amplifying it, eroding it, steering it. But there is a separate category of knowledge that sits outside the reach of AI altogether. Not because AI has not yet been trained on it, but because it has never been written down.

### Tacit knowledge

A large fraction of what makes a scientist expert is knowledge that is difficult to articulate in speech or writing. How to read a chromatography peak. When to trust a measurement and when to suspect contamination. How to present a negative result honestly. How to decide when an experiment is worth repeating and when to move on. How to recognise that a reaction flask looks wrong before you can say exactly why.

<!-- MARGIN: Tacit knowledge — knowledge that is difficult to articulate in writing or speech. It is learned through practice, observation, and experience rather than through instruction. In scientific contexts, tacit knowledge includes laboratory technique, pattern recognition in data, and professional judgement about how to interpret ambiguous results. It represents a large fraction of what makes an expert an expert. -->

This kind of knowledge is transmitted through years of working alongside someone more experienced: watching how they react to unexpected results, hearing them think aloud about a puzzling observation, absorbing their standards for what counts as "good enough" evidence. It is the foundation of the mentoring relationship, and it is what allows experienced scientists to make rapid, reliable judgements in situations where the formal rules do not give a clear answer.

AI cannot access this knowledge because it was never captured in the training data. It lives in the practices of working scientists — in how they move through a lab, how they read an instrument's output, how they decide what to do next when an experiment goes sideways. A 2025 preprint on AI laboratory agents acknowledged this directly: "Current AI tools have proven effective in literature analysis and code generation, but do not address the critical gap between documented knowledge and implicit lab practice."

AI can answer questions about technique. It can explain the principle behind a method, suggest troubleshooting steps, and generate protocols. But it cannot model scientific judgement in the way a mentor does — cannot show a student how to hold uncertainty wisely, how to distinguish a systematic artefact from a real phenomenon by feel, how to sit with a result that does not make sense until a better explanation arrives.

This is not a gap that will necessarily close with more data or better models. The knowledge was never written down because much of it *cannot* be written down — it is embedded in embodied practice, in the interplay between hand and eye and instrument, in the accumulation of thousands of small judgements that a scientist has made and learned from over a career. As AI's codifiable capabilities improve, this tacit dimension becomes more strategically important, not less — it is the part of expertise that remains distinctly human.

### Ethical and professional judgement in context

A related category sits at the boundary of tacit knowledge and social understanding. When is a result ready to publish? When does a reviewer's comment reveal a genuine flaw in your work versus a misunderstanding of your approach? When is a collaboration productive and when has it become extractive? When should you push back on a co-author's interpretation, and how?

These are judgements about the social and professional fabric of science. They require being embedded in a community — understanding its norms, its tensions, its unspoken expectations. AI can tell you what journal policies say about data sharing, but it cannot tell you whether sharing a particular dataset with a particular collaborator at a particular stage of the project is wise. That requires knowing the people, the power dynamics, and the stakes.

## What this means for how we train scientists

If the picture I have drawn is correct — AI amplifies expertise at the individual level, can erode it through the same mechanism, produces emergent effects at the population level, and cannot access the tacit knowledge that underpins much of scientific judgement — then how we build expertise deserves serious attention. A few implications follow.

**The case for more rigorous disciplinary training, not less.** It is tempting to conclude that if AI handles routine tasks, training can be shortened or made more superficial. The opposite is true. The routine tasks were often the vehicle through which deeper understanding was built. A student who hand-codes a statistical analysis, gets it wrong, debugs it, and finally understands what the model is doing has learned something that a student who asks AI for the code and runs it without error has not. Not every task needs to be done the hard way — but some do, and knowing which ones requires, again, the expertise of the mentor.

**Deliberate sequencing matters.** The apprenticeship has to come before the delegation. Foundational mastery benchmarks should be established before AI tools are introduced for core tasks. You cannot validate what you have not first learned to do yourself.

**Assessment needs to shift toward observable reasoning.** If AI can produce polished written outputs that are indistinguishable from expert work, then a polished output is no longer strong evidence of understanding. Assessment may need to move toward formats where reasoning is visible: oral defences, live problem-solving exercises, structured decision reviews where students explain not just what they concluded but how they got there and what alternatives they considered.

**Protect time for the irreducibly human.** If AI handles more of the mechanical work, the time reclaimed should go toward the activities where expertise is built and exercised at its highest level: mentoring conversations, collaborative thinking, sitting with a puzzling result, reading deeply rather than broadly. These activities are not luxuries that can be scheduled around the "real work." In an AI-assisted world, they *are* the real work.

**Resist the narrowing.** The population-level finding — that AI steers science toward tractable problems — is a collective action problem, but individuals and institutions can push against it. Funding bodies can specifically support data-sparse, anomaly-rich, and foundationally oriented research that AI tends to deprioritise. PIs can protect space in their research programmes for questions chosen because they are important, not because they are easy to work on with current tools.

## What this is really about

The question that prompted this post — "if AI knows my field, why do I need to?" — rests on a misunderstanding of what knowing a field means. It is not a body of information that can be stored and retrieved. It is a way of seeing: trained attention, calibrated judgement, the ability to recognise when something is wrong before you can articulate why. Much of it lives in practices that have never been written down and may not be possible to write down.

AI is not competing with that. AI is competing with the mechanical parts of scientific work — the parts that were always a means to an end. The end was always the understanding, the judgement, the ability to ask a question that opens a new line of inquiry. Those have not been automated. Whether we continue to develop them — in ourselves and in the students we train — is the question that matters.

Expertise is not what AI is replacing. Expertise is what makes AI worth using.
