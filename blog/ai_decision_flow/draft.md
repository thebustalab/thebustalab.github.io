# When Should a Scientist Use AI?

*A framework for thinking about it — not a set of rules*

---

There is a widening gap in how scientists relate to AI. On one side are researchers who have woven it into nearly every task without pausing to ask whether each use is appropriate — risking plausible text that nobody verified, analyses that run but weren't properly specified, literature summaries built on training data rather than the papers themselves. On the other side are colleagues who refuse to engage at all, forgoing genuine productivity gains on tasks where AI performs well and the risks are low.

<!-- MARGIN: Training data — the large body of text (books, articles, websites, code) that an AI model learned from during its development. When AI generates a response "from memory" rather than from documents you provide, it is drawing on this training data — which may be outdated, incomplete, or wrong. This post also uses the related terms "training knowledge" (what the model "knows" from its training data) and "training set" (the dataset used to build a model, with the added concern that your inputs to an AI service may become part of future training sets). These terms are used more or less synonymously. -->

There is also a practical reason to care about getting this right. Scientists are busy — genuinely, structurally busy. Teaching, mentoring, grant writing, committee work, and the administrative machinery of academic life consume time that could otherwise go towards the parts of science that actually advance a field: reading deeply, thinking carefully about a puzzling result, having an unhurried conversation with a collaborator or student. If AI can reliably handle some of the mechanical work — the reformatting, the boilerplate, the first-pass organisation — then the time reclaimed goes back to those irreducibly human activities. But that bargain only works if the AI use is thoughtful. Careless adoption erodes the very expertise that makes the tools worthwhile.

This post is my attempt to find a principled middle ground. Over the past year, I have been developing a framework for deciding when AI is and isn't a good fit for tasks in the natural sciences — grounded in actual research workflows, not abstract principles. It grew out of conversations with students and colleagues, and out of watching (and making) both good and bad decisions about when to reach for these tools.

A note on scope: this post focuses on *when* and *how* to use AI responsibly. It does not address the environmental costs of AI — those are real and worth taking seriously, and I plan to address them in a [separate post]().

I want to start with the ideas that kept surfacing as I built the framework.

## Five ideas I kept coming back to

**You are responsible for what AI produces.** This is the starting point, and everything else flows from it. If you put your name on a document, a dataset, or an analysis, you own it — regardless of which tool helped generate it. This means there must be a human in the loop, and that human must be meaningfully involved, not simply clicking "approve."

<!-- MARGIN: Human in the loop — a standard term in AI and engineering for a workflow where a person reviews, validates, or intervenes in an automated process before the output is used. The question this post keeps returning to is what it takes for that human involvement to be meaningful rather than ceremonial. -->

**Human in the loop has to mean something.** The phrase "human in the loop" has become a kind of incantation — invoked to justify AI use without examining what the human is actually contributing. A graduate student who has never studied mixed-effects models cannot meaningfully review AI-generated code for that analysis, even if the code runs without errors. For the loop to be more than ceremonial, the person in it needs enough expertise to catch what the AI gets wrong.

**[Disciplinary expertise]() is the force multiplier — not AI itself.** This is perhaps the most counterintuitive point: in a world where AI handles more of the mechanical work, deep training in your field becomes *more* important, not less. Expertise is what allows you to validate AI output quickly and confidently. Without it, you are generating text or analyses that you cannot properly evaluate, which is worse than not using AI at all. A well-trained scientist using AI for tasks they can verify is a force multiplier. A poorly-trained person using AI for tasks they cannot verify is a liability. This theme surfaces repeatedly throughout the examples below — I will note it briefly each time, but for the full discussion, see the [companion post on disciplinary expertise]().

**Most tasks are not one task.** When someone asks "should I use AI to write my literature review?" or "can AI help with my grant proposal?", the honest answer is almost always "it depends on which part." A literature review decomposes into summarising individual papers, identifying themes across the literature, and constructing a narrative argument. AI may be well suited to the first, moderately useful for the second, and poorly suited to the third. The habit of breaking tasks into subtasks, and evaluating each one separately, changes the question from a blunt yes/no into something much more useful.

**Personal voice and original intellectual contribution are strategic tools.** Not every piece of writing needs to sound distinctly like you, and not every task requires your own original thinking from scratch. Reformatting citations, summarising a dataset, drafting boilerplate for an IRB application — these don't require a personal voice or original ideas, and AI handles them well.

<!-- MARGIN: Boilerplate — standardised, reusable text that follows an established template. In academic contexts, this includes things like data management plan language, consent form templates, and budget justification prose that varies little between documents. --> But a recommendation letter needs to sound like the person writing it, because the recipient is reading for authentic engagement. A conference talk should sound like the presenter. And the central argument of a thesis or grant proposal must be the author's own intellectual contribution — not AI-generated prose that merely sounds plausible. Knowing when voice and originality matter — and when they don't — is itself a form of disciplinary expertise. It requires understanding your field well enough to recognise which parts of a task carry your professional identity and which are mechanical.

## From principles to a decision flow

These five ideas needed a structure — a way to apply them consistently rather than relying on ad hoc judgement each time. The result is a decision flow: a sequence of checkpoints (I call them *gates*) that a task passes through before you decide whether AI is appropriate. Rather than presenting this as an abstract diagram, let me show what it looks like applied to real tasks.

<!-- MARGIN: Cloud-hosted AI — services like ChatGPT, Claude, or Gemini, where your input is sent to the provider's servers for processing. Convenient and powerful, but your data passes through (and may be stored by) a third party. -->

<!-- MARGIN: Local model — an AI system installed on your own computer or lab workstation (e.g. Llama, Gemma). Less capable than frontier models, but your data stays entirely on your machine — nothing is transmitted externally. -->

<!-- MARGIN: Frontier model — the most capable AI models available at any given time, typically cloud-hosted and developed by major AI labs (e.g. GPT-4, Claude, Gemini). These are substantially more capable than local models, but using them means sending your data to an external service. The gap between frontier and local models is narrowing but remains significant. -->

### The easy case: reformatting a reference list

Start simple. You have a reference list and need to convert it to a different journal's citation style. No confidential data is involved. The task is atomic — it doesn't need to be broken into parts. AI is excellent at text reformatting, and the output is trivially easy to check: either each citation matches the target style or it doesn't. There is no ethical concern — nobody expects you to reformat references by hand. Not every case is this clean, but it establishes a baseline: when the task is text-native, easy to validate, and carries no ethical weight, AI is the right tool.

<!-- MARGIN: Text-native — a task that is fundamentally about reading, writing, or processing text. Summarising a paper, reformatting a citation list, or drafting an email are text-native. Navigating a web portal, running a laboratory instrument, or editing a slide deck are not. AI is strongest on text-native tasks. -->

### When the purpose is learning: exploratory mode

Not all AI use is professional output. A graduate student might paste a dataset into a chat window just to see what happens — can the model spot trends? Does it suggest a reasonable statistical approach? A postdoc might ask AI to explain a technique from a neighbouring field, not to cite the answer but to orient themselves before reading the primary literature. This is what I call *exploratory mode*: using AI to build intuition, test boundaries, and learn what the tools can and cannot do. The output is not going into a paper, a grant, or a report — the learning is the point. In exploratory mode, the evaluation can be lighter, and a lot of genuine skill-building happens here — but confidentiality obligations and the prohibition on submitting AI-generated text as your own still apply fully.

The examples that follow are all in *considered mode* — professional work where the output matters. But if you are new to AI and wondering where to start, exploratory mode is the answer: try things, see what works, develop your own sense of where AI adds value and where it doesn't.

### Where decomposition changes the answer: a thesis literature review

A first-year graduate student asks: "Can I use AI to write my literature review?" If you treat this as a single task, the answer is murky. Before even reaching the question of *how* to use AI here, there is a preliminary gate: are the papers accessible? If the full text sits behind a paywall and hasn't been made openly available by the authors, uploading it to a cloud-based AI raises the same confidentiality concerns as any other restricted material. The data confidentiality gate applies here just as it does for unpublished research data.

Assuming the papers are accessible, decompose the task and the picture sharpens.

*Summarising individual papers* — if the student provides the papers as context — a technique called *context grounding* — rather than relying on AI's training knowledge, AI can produce serviceable summaries.

<!-- MARGIN: Context grounding — the practice of providing AI with your own source documents (papers, datasets, notes) and asking it to work from those, rather than relying on what it learned during training. This converts the validation task from independent fact-checking (hard) into comparison against sources you already trust (much easier). --> The student can validate them because they have read the papers. But AI summaries tend to capture main claims while missing subtler points — methodological choices that only make sense in context, caveats buried in a discussion section. Recognising which details matter is itself a skill that comes from disciplinary expertise, and if summarisation is used extensively, the student's understanding of the literature becomes shallower as a result.

*Identifying themes and gaps across the literature* — harder. AI can suggest groupings, but whether the right themes have been identified requires domain knowledge. This is where the student's own reading and thinking must drive the work.

*Constructing the narrative arc that positions the thesis* — this is the student's original intellectual contribution. It is the core of the review: the argument about why this research matters, where the field stands, and what gap the thesis addresses. AI can help with sentence-level prose once the argument exists, but the argument itself must originate with the student. This is not a task that can be delegated and then validated — it is the thinking that the literature review exists to demonstrate.

Note, too, that all of this hinges on the student's own expertise. A student who has not read deeply in the field cannot validate even the summaries meaningfully — the same subtask is a strong AI use case for a well-read student and a poor one for a student who is just starting out.

### When the technical evaluation looks fine but something else is decisive

Consider a graduate student drafting an email to a potential PhD supervisor. The technical evaluation looks favourable — no sensitive data, text-native, easy to validate. But experienced academics often recognise AI-generated prose, and an email meant to demonstrate genuine interest reads as the opposite if it feels templated. The perception cost is decisive here: researching the supervisor's publications is a fine subtask to hand off, but the email itself needs to be the student's own words.

### The sharpest contrast: peer review versus responding to reviews

These two cases look similar but land in completely different places.

**Peer reviewing someone else's manuscript.** The manuscript is confidential. Uploading it to a cloud-based AI arguably violates that confidentiality, and many journals now explicitly prohibit this. That alone is a serious constraint. But even setting aside the data question, peer review is built on the premise that an expert human is evaluating the work. AI-generated reviews undermine the system itself, even when the content happens to be accurate. The ethical cost is not hypothetical — several journals have identified and sanctioned AI-generated reviews. This is a case where the responsibility assessment is decisive on its own, regardless of how the other considerations play out. The one defensible subtask is mechanical: checking reference formatting and compliance with journal guidelines.

**Responding to reviewer comments on your own manuscript.** The ethical picture is entirely different. You are engaging with your own work. Organising the reviewer comments into categories is a clear win — mechanical, easy to verify. Drafting point-by-point responses works well when you provide the original manuscript and the reviews as context — AI is strongest when the correct information is already present and the task is to reorganise or articulate it, and you can validate the result against your own expertise. The cover letter to the editor is straightforward boilerplate.

Two superficially similar tasks, opposite conclusions. The difference is about confidentiality, integrity, and what the task represents in scientific publishing.

### High-stakes professional work: the grant proposal

A faculty member writing an NSF proposal navigates nearly every consideration in the framework simultaneously. The research idea itself may be sensitive — many PIs would prefer their unfunded proposal concept not be ingested into a training set.

Decomposition reveals a wide spread. Budget justification and compliance checks are excellent AI candidates — mechanical, verifiable, text-native. Literature synthesis to establish the gap works well with context grounding: provide the key papers and ask AI to extract and organise, rather than generate from training knowledge. The PI can validate quickly because this is their field.

But as with the literature review, the central hypothesis and research plan are where original intellectual contribution sits. The ideas must originate with the PI. Grant agencies are increasingly requiring disclosure of AI use, and there are genuine integrity questions if these core elements are AI-generated.

This example also shows how data confidentiality interacts with capability. If the confidentiality gate routes you to a local model, that model may be less capable than a cloud-hosted alternative. A subtask that would be a strong AI candidate with a frontier model might be a poor one with a smaller local model.

### The "boring" task with real stakes: an IRB application

Drafting standard boilerplate sections of an IRB application — data storage plans, consent withdrawal procedures, risk mitigation language — is a good AI candidate: verifiable against templates, text-native, and a genuine time saver. But ensuring compliance with federal regulations is harder to validate without regulatory expertise, and navigating the institutional submission portal (authenticated web interfaces, dropdown menus, conditional logic across screens) is not a text-native task at all — current AI tools simply cannot do it. The framework applies to administrative tasks as well as scholarly ones, and "boring" does not mean "low-stakes."

## The framework underneath

Here is the structure made explicit:

1. **Data confidentiality** — does this task involve sensitive or pre-publication data? If so, can you use a local model, or de-identify sufficiently? If neither, AI is not appropriate.

<!-- MARGIN: De-identify — to remove or obscure information that could identify individuals, institutions, or specific datasets. In a research context, this might mean stripping participant names, replacing specific values with ranges, or removing location data. Even de-identified data can sometimes be re-identified when combined with other available information, so this step requires care — and disciplinary knowledge, since understanding which variables carry re-identification risk depends on knowing the data and the field. -->

2. **Exploratory versus considered mode** — are you learning and experimenting, or doing professional work? Exploratory mode is lighter on evaluation, but the responsibility check still applies.

3. **Task decomposition** — break the task into subtasks and evaluate each one independently.

4. **Validateability and AI capability fit** — for each subtask: can you efficiently verify the output, and is this the kind of task AI is practically good at? Both need to be favourable.

5. **Responsibility** — ethical implications, reproducibility, transparency, and the perception of AI involvement. This is the floor that all AI use sits on.

I have also built an [interactive version of this decision flow]() that lets you trace specific tasks through each node. It is a companion to this post, not a replacement for the reasoning above.

## What this is really about

The framework here is not a set of rules to follow mechanically. It is a way of thinking — a set of questions to ask yourself before reaching for the tool, and a structure for answering them honestly. The answers will differ across tasks, across career stages, and across disciplines. They will also change as AI capabilities evolve. What should not change is the habit of asking the questions in the first place.

I said at the outset that the time reclaimed by thoughtful AI use can go back to the parts of science that matter most — the conversations, the mentoring, the slow thinking. That remains the core argument. But it only holds if the use is genuinely thoughtful. The framework is one attempt to make "thoughtful" concrete rather than aspirational.
