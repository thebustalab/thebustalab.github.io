# The Environmental Cost of AI

*What the numbers actually say — and what they mean for you*

---

"What about the environmental impact?" It is the most common question I get after talking about AI with students and colleagues — and it is a fair one. AI runs on data centres that consume electricity and water. The infrastructure is growing fast. The headlines are alarming. If you are a scientist who has started using AI in your work, you might reasonably wonder whether you are contributing to an environmental problem every time you run a query.

This post is a data-driven attempt to answer that question. The answer turns out to be less straightforward than the intuitions most people bring to it, and I think the way many of us have been thinking about this — including me, initially — is misleading in a specific and interesting way. Understanding why changes what responsible action actually looks like.

Two things need to be addressed separately, because they tell very different stories. The first is the environmental cost of *your* AI use — the queries you run, the summaries you request, the code you generate. The second is the aggregate cost of AI at the level of global infrastructure — the data centres, the chip manufacturing, the energy demand. These are related, but the relationship is not what most people assume. I will work through each in turn, and then explain why the usual framing leads us astray.

This is a companion to my earlier posts on [when to use AI in the natural sciences](../ai_decision_flow/) and [disciplinary expertise in an AI-assisted world](../ai_and_domain_expertise/).

## Your AI use: the individual picture

Let me start with the numbers, because they are more concrete than most people expect.

A single text query to a model like ChatGPT, Claude, or Gemini consumes roughly 0.3 watt-hours of energy and produces about 0.03 grams of CO2. To make that tangible: it is less energy than running a microwave for two seconds. It is less than watching nine seconds of television. A hundred queries — a heavy day of AI-assisted research — produces about 3 grams of CO2, which is comparable to driving a car roughly 10 metres.

<!-- MARGIN: These figures come from Google's 2025 disclosure for Gemini (0.24 Wh per text query, ~0.03 g CO2e), OpenAI's estimates (~0.34 Wh), and independent analysis by EpochAI (~0.3 Wh). Earlier estimates from 2023 (~3 Wh per query) appear to have been inflated by roughly an order of magnitude. Per-query efficiency has improved dramatically — roughly 33-fold in a single year — due to hardware and software optimisation. -->

Even at heavy use — say, 500 queries per day, which would be an extraordinary amount of AI interaction — the energy consumed is about 1.5% of an average US household's daily electricity use. The carbon produced is about 15 grams, comparable to driving a car 100 metres.

For a researcher using AI for literature review, coding assistance, or data analysis, the individual environmental footprint of that use is, by any honest accounting, negligible. It is a rounding error within a rounding error of your daily carbon footprint.

There is a companion [interactive tool](../ai_literacy/) that lets you compare the CO2 and water footprint of AI use against common daily activities — making a cup of coffee, taking a shower, eating a hamburger. I encourage you to click through it, because the comparisons are striking.

<!-- MARGIN: Water footprint — AI also consumes water, primarily for cooling data centres. Per-query estimates are less precise than energy figures, but a rough figure is about 0.5 mL per text query. For context, a single cup of coffee has a virtual water footprint of about 140 litres (accounting for growing, processing, and transporting the beans). Your daily AI use consumes less water than a single toilet flush. -->

### When it is not negligible

There is an important exception. The figures above are for text-based queries — the kind of AI use most scientists engage in. Image generation consumes considerably more energy per task, and video generation consumes roughly 3,000 times more than a text query. Generating a five-second AI video uses about 1 kWh — equivalent to running a laptop for 25 hours. If you are routinely generating AI images or videos, the footprint is no longer trivial and is worth thinking about.

But for text-based research use — asking questions, summarising papers, generating code, analysing data — the individual cost is genuinely small.

### A surprising finding

Here is where the picture becomes more interesting than a simple "it's fine, don't worry." A 2024 study published in *Scientific Reports* compared the carbon emissions of writing and illustrating tasks performed by AI versus by humans. The result: AI-generated output produced substantially lower emissions than the human equivalent — in some cases, by a factor of 40 to 150. The reason is that the human process involves commuting to an office, heating or cooling a workspace, running equipment for longer, and all the other energy expenditures that come with a person doing a task over hours or days. The AI does the same task in seconds, using a fraction of the energy.

This does not mean AI use is carbon-negative or that more AI is always better for the environment. But it does mean that for specific tasks, choosing to use AI instead of doing the work yourself can *reduce* your personal carbon footprint rather than increase it. That is a genuinely counterintuitive finding, and it matters for how we think about the ethics of individual AI use.

## The aggregate picture: what happens at scale

The individual numbers are reassuring. The aggregate numbers are not.

AI data centres currently consume roughly 60–85 terawatt-hours of electricity per year, with all data centres combined consuming about 415 TWh — roughly 1.5% of global electricity. That figure grew 17% in 2025 alone, far outpacing global electricity demand growth of about 3%. The International Energy Agency projects total data centre demand roughly doubling by 2030, to about 945 TWh — comparable to Japan's entire annual electricity consumption.

<!-- MARGIN: Terawatt-hour (TWh) — a unit of energy equal to one trillion watt-hours. For scale: the entire United States consumes roughly 4,000 TWh of electricity per year. A single household uses about 10,000 kWh (0.00001 TWh) per year. -->

The carbon emissions from data centres are currently estimated at about 180 million tonnes of CO2 per year — roughly 1.2% of global energy-sector emissions. AI's specific share within that is estimated at 33–80 million tonnes, with wide uncertainty because tech companies do not separate AI workloads from other data centre activity in their reporting.

To put this in context: global aviation produces about 880 million tonnes of CO2 per year, roughly 2% of the global total. Data centres are currently below aviation, though the gap is narrowing. The widely repeated claim that AI has already surpassed aviation appears to be based on outdated figures and is not supported by the most recent IEA data. But the trajectory matters: data centre emissions are growing at 12–17% per year, while aviation emissions are relatively stable.

<!-- MARGIN: The aviation comparison — you may have seen headlines claiming AI already exceeds aviation in carbon emissions. The IEA's 2025 report and independent analysis suggest this is not yet accurate: data centres are at roughly 1.2% of global carbon, aviation at roughly 2%. The confusion traces partly to a 2015 overestimate that continues to circulate. The comparison remains useful for scale — AI is in the same order of magnitude as a major global industry — but the specific claim of surpassing aviation is premature. -->

Water is another dimension. Data centres require substantial amounts of water for cooling — roughly 2 litres per kilowatt-hour consumed. A typical large data centre consumes about 2 million litres of water per day. Globally, data centres used about 560 billion litres in 2024, projected to rise to 1,200 billion litres by 2030. And data centres are disproportionately built in regions that are already water-stressed — the American Southwest, parts of Ireland, northern Virginia — which means the local impact can be acute even when the global figure sounds manageable.

The aggregate picture, then, is a real and growing environmental cost. Not the largest contributor to emissions or water stress — cement production alone accounts for 8% of global carbon, agriculture for 22% — but a fast-growing one, and one that is being locked in through infrastructure decisions being made right now.

## Why the usual framing misleads

Most people think about this through the lens of personal environmental responsibility — the same lens we use for driving, flying, or eating meat. Use less, emit less. Every bit helps. If AI is bad for the environment, the responsible thing is to use less of it.

This framing is intuitive, and for many environmental problems it works well enough. But for AI, it is misleading in a way that matters.

Consider the difference. If you choose to bike to work instead of driving, you have directly reduced your emissions. Your individual action is small relative to the global total, but it is unambiguously in the right direction. The relationship between your behaviour and the environmental outcome is straightforward: less driving, less carbon.

With AI, the relationship is muddled. As we saw above, using AI for a writing or analysis task can actually produce *less* carbon than doing the task yourself. Your individual query costs a fraction of a gram of CO2. And reducing your personal AI use does not reduce the infrastructure — the data centres are built, the chips are manufactured, the cooling systems run whether you submit a query or not. The aggregate problem is driven by AI being embedded in industrial processes, built into every search engine query, powering automated systems that run millions of operations per hour. It is an infrastructure problem, not a behaviour problem.

This is different from how most environmental problems work at the individual level. It is not like picking up rubbish in your neighbourhood, where your action directly fixes the local problem. It is not even like reducing your driving, where your action is tiny but directionally correct. In many cases, thoughtful AI use is the lower-emission option — and abstaining accomplishes nothing beyond denying yourself a useful tool.

<!-- MARGIN: The infrastructure distinction — the environmental cost of AI is driven primarily by decisions about data centre construction, energy sourcing, chip manufacturing, and cooling technology. These are decisions made by companies and governments. An individual scientist's choice to run or not run a query has no measurable effect on any of these. This is what makes AI's environmental impact structurally different from, say, personal vehicle emissions, where individual choices do aggregate in a straightforward way. -->

This does not mean the aggregate problem is not real. It is. But it means the right response is not guilt about your personal use. It is understanding where the real leverage points are.

## What intentional use actually looks like

If reducing your personal query count does not meaningfully address the environmental cost of AI, what does?

**Be intentional, not abstinent.** Use AI when it genuinely helps your work. Do not use it frivolously — not because a frivolous query destroys the planet, but because intentional use is good practice for its own sake, and it happens to align with not generating pointless compute load. The distinction between a researcher who uses AI to carefully analyse a dataset and someone who generates AI images for entertainment is not primarily an environmental distinction — the carbon difference is small either way — but it reflects a relationship with the technology that is more considered and more productive.

**Understand the real leverage points.** The systemic environmental cost of AI is driven by data centre energy sourcing, cooling technology, geographic siting decisions, and the pace of grid decarbonisation. These are corporate governance and policy questions. If the environmental cost of AI genuinely concerns you — and it is reasonable for it to — the most effective response is not behavioural (using AI less) but structural: supporting policies that require transparency in data centre energy and water reporting, advocating for renewable energy procurement standards, and understanding how grid infrastructure decisions are made.

**Consider where your career can make a difference.** For students in the natural sciences, this is perhaps the most practical takeaway. If you want to address AI's environmental impact — or the broader energy transition — there are direct paths: battery technology, grid engineering, energy policy, environmental advocacy, materials science for more efficient computing hardware. These are careers where your contribution is not a drop in the bucket. It is the bucket.

<!-- MARGIN: Grid decarbonisation — the process of shifting electricity generation from fossil fuels to low-carbon sources (renewables, nuclear). Currently, about 60% of data centre electricity comes from fossil fuels. The IEA projects renewables will meet about half of data centre demand growth through 2030, but the other half is being filled primarily by new natural gas plants — 38 GW of gas capacity is currently planned specifically for data centre service. This is the crux of the infrastructure problem: AI demand is growing faster than the grid is greening, and the gap is being filled by fossil fuels that will operate for decades. -->

**Know the exception.** If your AI use involves regular image or video generation, the energy cost per task is orders of magnitude higher than text. That is worth factoring into your workflow decisions — not as guilt, but as information.

## What this is really about

The environmental cost of AI is real, growing, and worth taking seriously. But taking it seriously means understanding it accurately — not through the lens of personal guilt, which leads to the wrong conclusions, but through the lens of systems and infrastructure, which reveals where the actual problems and solutions lie.

For an individual scientist, the honest summary is this: your text-based AI use is environmentally negligible, and in some cases it is the lower-emission option compared to doing the task without AI. The aggregate cost is driven by industrial-scale deployment and infrastructure decisions that are far beyond your individual control. The most scientifically honest response is not to feel guilty about using a tool that helps your work, but to understand the system well enough to know where the real leverage points are — and to consider whether acting on that understanding is something you want to build into your career.

That, in the end, is another form of the disciplinary expertise this series keeps returning to: the ability to look past the headline, examine the data, and reach your own conclusion about what it means and what to do about it.
