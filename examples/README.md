# Small SGeMS x64 examples

Load the input object first, then load the matching `.par` in the named plugin. These are synthetic software examples, not paper experiments.

| Parameter file in `Parameters/` | Input in `Inputs.prj/` | Module | Weighting |
|---|---|---|---|
| `grid_template_uniform.par` | `grid.sgems` | HOSC_TEMPLATE | Uniform |
| `grid_window_uniform.par` | `grid.sgems` | HOSC_WINDOW | Uniform |
| `points_template_uniform.par` | `points.sgems` | HOSC_TEMPLATE | Uniform |
| `points_window_uniform.par` | `points.sgems` | HOSC_WINDOW | Uniform |
| `grid_template_external.par` | `grid.sgems` | HOSC_TEMPLATE | External property `weight` |
| `grid_window_external.par` | `grid.sgems` | HOSC_WINDOW | External property `weight` |

Object names are `grid` and `points`; the statistical property is `value`. The grid has 3 × 2 × 2 cells. The point fixture includes coincident observations. Extra properties are deliberate test inputs; only `value` and, for the two external cases, `weight` are selected here.

All six examples calculate fourth-order cumulants with the AVL evaluator. Template examples use a fixed origin and three occurrences at the lagged position. Window examples cover the small source support. Set **Record directory** to a writable persistent English path if records are needed.

The examples can be checked automatically in a supplied isolated host with `python tools/run_examples.py --host <isolated-host> --out <new-output-folder>`. This executes the DLLs inside SGeMS and checks the outputs, rather than substituting a command-line calculator.

For the three declustering methods and paper settings, use the current [paper parameters](../paper/README.md). The broader regression cases live in `tests/`, not in this introductory directory.
