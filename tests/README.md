# Tests and reference fixtures

## Reviewer installation check

After installing the prebuilt package, run from the repository root:

```powershell
.\RunExamples.cmd "C:\SGeMS"
```

This runs small calculations with the supplied native Cartesian grid in a
separate SGeMS environment, exports results, and validates numerical values and
support. The output directory and validation manifest are printed by the
runner. Each run writes a separate `run_<timestamp>_<id>` directory containing
`RESULT.json`, `reported_issues.csv`, `template_result_checks.csv`, exported
`.sgems` grids, and window audit files. No compiler or Python is required.

The four shipped `examples/*.par` files are executed with their original
parameters: a fourth-order exact template, a uniform fourth-order window,
and coordinate-Gaussian and soft-high-porosity sixth-order windows. The
54 checked output nodes are recorded in `gui_example_result_checks.csv`;
`executed_parameters/` preserves the files and `RESULT.json` records their
hashes and per-case checks.

The runner deliberately submits invalid parameters before the valid cases to
check field diagnostics and same-process recovery. The four diagnostics in
`reported_issues.csv` are expected. `RESULT.json` must report `status: passed`.

For explicit paths:

```powershell
powershell.exe -NoProfile -ExecutionPolicy Bypass -File .\tools\run_examples.ps1 -SGeMSRoot "C:\SGeMS" -OutputDirectory "C:\HOSC_examples"
```

[QUICKSTART.md](../docs/QUICKSTART.md) gives equivalent interface actions and
expected values. [MANUAL_TESTING.md](../docs/MANUAL_TESTING.md) covers additional
weighting, missing-value, invalid-input, save/reload, and repeated-run cases.

## Fixtures

| File | Contents and use |
|---|---|
| `fixtures/manual/selftest_3d_cartesian.sgems` | Native Cartesian grid named `test`, dimensions `5 × 5 × 3`, unit cells, zero origin, property `phi`; ready to load |
| `fixtures/manual/selftest_2d.gslib` | 25 rows, dimensions `5 × 5 × 1`; import as a Cartesian grid |
| `fixtures/manual/selftest_3d.gslib` | 75 rows, dimensions `5 × 5 × 3`; import as a Cartesian grid |
| `fixtures/manual/expected_window.csv` | Reference values for `phi` and `spike`, three weight modes, orders 1–8, and nine 3D windows |
| `fixtures/manual/expected_template.csv` | Reference values for the documented 2D, 3D, and singleton-third-lag cases |

The native `.sgems` fixture contains only `phi`. The GSLIB files contain
`phi`, `constant`, `phi_shift`, `phi_missing0`, `locator`, and `spike`. Their
definitions and import settings are in the manual testing guide. The expected
CSVs are reference calculations, not logs of host executions.

For the basic uniform, fourth-order `3 × 3 × 3` window on the native fixture,
all nine results have raw cumulant `-6.2306722005208369e-5`, variance
`0.009114583333333334`, standardized cumulant `-0.75`, and support/ESS `27`.
The extended template CSV describes its named cases; the quick-start
`3 × 3 × 3` lag example is a separate 27-cell configuration.

## Developer checks

`vc9_core_checks.cpp` retains its historical filename and exercises the shared
weighted-statistics core: nonzero-mean low-order cumulants, independent central
moment formulas through order six, translation invariance through order eight,
standardized cumulants, effective sample size, soft weights, and relative-window
Gaussian geometry. It is a core test, not a SGeMS GUI test. Use the matching
build procedure described in [the build guide](../docs/BUILD_FROM_SOURCE.md).

`tools/verify_moving_window_contract.ps1` and
`tools/verify_template_contract.ps1` check source/interface contracts and small
layout models. They do not execute the host's grid or numerical implementation.

The independent standalone implementation has CTest checks and Python tests:

```powershell
cmake -S standalone -B standalone/build
cmake --build standalone/build --config Release
ctest --test-dir standalone/build -C Release --output-on-failure
```

See [standalone/README.md](../standalone/README.md) for toolchains, Python
requirements, and the experiment tests. Scientific result integrity is checked
separately by the tools under `Paper-data/tools/`.
