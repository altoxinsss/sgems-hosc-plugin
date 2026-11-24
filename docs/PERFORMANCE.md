# Per-run performance records

Each execution writes `performance.csv` and `performance.json` inside its run directory. Select a persistent **Record directory** to retain them. Records include the module, actual evaluator, status, source/case names, elapsed seconds, baseline working set, observed peak and peak increase.

Elapsed time covers plugin execution: input export, numerical calculation, diagnostics and result import. It excludes application startup, parameter loading and performance-record writing.

Memory is the host process working set sampled during this execution, with a requested interval of 5 ms plus phase checkpoints. `peak_working_set_increase_mib` is the observed run peak minus the pre-execution baseline. It is neither the process lifetime peak nor an exact instantaneous maximum. MiB means bytes divided by 2^20.

For comparable memory measurements, restart SGeMS for each run, load the same inputs and run one configuration. For timing summaries, use the same repetition and warm-up protocol for both evaluators. The current paper's merged records and summaries are under [paper/Data/CSV](../paper/Data/CSV), with their measurement definitions in [the paper data guide](../paper/README.md#precision-and-performance-records).

An initialization error can stop before execution starts and therefore before a performance record exists. An execution failure is recorded with its failure status when the run directory is available. Only successful records with valid measurement status should enter performance summaries.
