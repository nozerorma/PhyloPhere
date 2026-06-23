export const meta = {
  name: 'audit-module',
  description: 'Audit any pipeline module using 4 specialized subagents',
  phases: [
    { title: 'Init', detail: 'Locating files and setting up tasks' },
    { title: 'Analysis', detail: 'Running specialized subagents' },
    { title: 'Synthesis', detail: 'Compiling findings into a report' }
  ]
}

// Ensure execution is inside an async context (implicitly handled in workflow runner)
const main = async () => {
  phase('Init');
  const moduleName = args.moduleName || 'CT_ACCUMULATION';
  const modulePath = args.modulePath || 'subworkflows/CT_ACCUMULATION';
  const files = args.files || [];

  log(`Starting multi-agent audit for module: ${moduleName} under ${modulePath}`);

  const analysisPrompts = [
    {
      key: 'code-audit',
      label: 'Code Audit',
      prompt: `You are a Senior Code Auditor specializing in Nextflow workflows and Python/R script verification.
Task: Critically audit the implementation of the ${moduleName} module.
Files to look at: ${files.join(', ')}

Analyze:
1. Correctness of logic: Are files read/written properly?
2. Edge cases: Missing files, empty inputs, duplicate inputs, zero datasets.
3. Error handling: How robust is it to unexpected exceptions/errors?
4. Bug hunting: Search for indexing bugs, off-by-one errors, or wrong assertions.
Output: A detailed technical report outlining findings, bugs found, and recommendations.`
    },
    {
      key: 'methodology',
      label: 'Methodological Grounding',
      prompt: `You are a Bioinformatics Methodologist specializing in statistical genetics and comparative genomics.
Task: Audit the theoretical and methodological grounding of the ${moduleName} module.
Files to look at: ${files.join(', ')}

Analyze:
1. The Statistical Framework: Null hypothesis, test structure, assumptions, statistical soundness.
2. The Sampling/Modeling Logic: Does the modeling fit biological reality? Are there biases toward gene length, size, count?
3. Biological Caveats: Control of phylogenetic history, background rates, GC content, recombination, etc.
Output: A rigorous methodological report evaluating the design, highlighting caveats, and proposing justifications/improvements.`
    },
    {
      key: 'coherence',
      label: 'Pipeline Coherence',
      prompt: `You are an Integration Architect specializing in bioinformatics pipelines.
Task: Audit the coherence of the ${moduleName} module with upstream outputs and downstream consumers.
Files to look at: ${files.join(', ')}
You can search the repository to trace how inputs are produced upstream and how outputs are consumed downstream.

Analyze:
1. Column mappings and data types: Are columns from upstream steps parsed and used correctly? Any mismatches or contradictions?
2. Data Flow Coherence: Are all outputs processed correctly? Is there any unused/dead code or wasted computation?
3. Statistical Interpretation: Are the outputs interpreted correctly by downstream scoring/reporting steps?
Output: A detailed coherence report identifying any schema mismatches, semantic errors, or data loss.`
    },
    {
      key: 'efficiency',
      label: 'Performance & Simplification',
      prompt: `You are a Performance Engineer specializing in code optimization and design simplicity.
Task: Audit the code efficiency and simplification opportunities of the ${moduleName} module.
Files to look at: ${files.join(', ')}

Analyze:
1. Performance bottlenecks: Multiprocessing, vectorization (NumPy/Pandas), memory leaks/copies, slow I/O.
2. Simplification: Can we simplify the codebase? Any redundant/legacy features, unused arguments, or unnecessary dependencies?
Output: An efficiency and simplification report with concrete suggestions, profiling insights, and refactoring recommendations.`
    }
  ];

  phase('Analysis');
  const results = await parallel(analysisPrompts.map(p => () =>
    agent(p.prompt, { label: p.label, phase: 'Analysis' })
      .then(report => ({ key: p.key, label: p.label, report }))
  ));

  phase('Synthesis');
  log('Synthesizing report...');

  const synthesisPrompt = `You are a Lead Bioinformatics Architect.
We have conducted a multi-agent audit of the ${moduleName} module using four specialized perspectives:
1. Code Audit
2. Theoretical & Methodological Grounding
3. Pipeline Coherence
4. Performance & Simplification

Here are the reports from the subagents:
${results.map(r => `=== ${r.label} Report ===\n${r.report || 'Failed to generate'}\n`).join('\n')}

Task: Synthesize these findings into a single, cohesive, high-quality, and educational review report for the user. Highlight critical bugs, theoretical caveats, pipeline integration gaps, and performance optimizations. Provide clear, actionable recommendations.`;

  const finalReport = await agent(synthesisPrompt, { label: 'Synthesizer', phase: 'Synthesis' });
  return { finalReport, individualReports: results };
};

return await main();
