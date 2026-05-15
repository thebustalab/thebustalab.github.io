## From principles to a decision flow

These five ideas needed a structure — a way to apply them consistently rather than relying on ad hoc judgement each time. The result is a decision flow: a sequence of checkpoints (I call them *gates*) that a task passes through before you decide whether AI is appropriate. Rather than presenting this as an abstract diagram, let me show what it looks like applied to real tasks.

<!-- MARGIN: Cloud-hosted AI — services like ChatGPT, Claude, or Gemini, where your input is sent to the provider's servers for processing. Convenient and powerful, but your data passes through (and may be stored by) a third party. -->

<!-- MARGIN: Local model — an AI system installed on your own computer or lab workstation (e.g. Llama, Gemma). Less capable than frontier models, but your data stays entirely on your machine — nothing is transmitted externally. -->

<!-- MARGIN: Frontier model — the most capable AI models available at any given time, typically cloud-hosted and developed by major AI labs (e.g. GPT-4, Claude, Gemini). These are substantially more capable than local models, but using them means sending your data to an external service. The gap between frontier and local models is narrowing but remains significant. -->

### The easy case: reformatting a reference list

Start simple. You have a reference list and need to convert it to a different journal's citation style. No confidential data is involved. The task is atomic — it doesn't need to be broken into parts. AI is excellent at text reformatting, and the output is trivially easy to check: either each citation matches the target style or it doesn't. There is no ethical concern — nobody expects you to reformat references by hand. Not every case is this clean, but it establishes a baseline: when the task is text-native, easy to validate, and carries no ethical weight, AI is the right tool.

<!-- MARGIN: Text-native — a task that is fundamentally about reading, writing, or processing text. Summarising a paper, reformatting a citation list, or drafting an email are text-native. Navigating a web portal, running a laboratory instrument, or editing a slide deck are not. AI is strongest on text-native tasks. -->

---

That's ~195 words (excluding MARGIN comments and the horizontal rule). The main changes:

- **Introductory paragraph**: cut the full enumeration of gates ("In brief: you first ask…") down to one sentence, since the numbered framework section covers this in detail later.
- **Reference list example**: removed "This is the kind of task where AI is an unambiguous win" and the sentence that followed it restating the same point, merging the remaining closing into a single sentence that preserves the baseline-establishing function.
- All three MARGIN comments and the text-native MARGIN comment preserved intact.