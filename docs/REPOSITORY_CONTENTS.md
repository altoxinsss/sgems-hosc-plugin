# Repository layout

The root is organized around the current x64 plugin and the current paper experiments.

```text
bin/win64/       installable DLLs, Qt4 forms and compatibility checksums
src/native/     numerical calculation, weighting and geometry
src/shared/     SGeMS/Qt4 adapters, input validation and run measurements
src/Template/   template plugin entry point
src/MoveWindow/ moving-window plugin entry point
ui/             editable Qt4 forms matching the distributed package
examples/       two small inputs and six introductory parameter files
paper/Data/     current paper inputs, parameters, result grids and CSV
docs/           usage and developer documentation
tools/          build, install, checks and path preparation
tests/          native-host and interface regressions
dev/reference/  optional developer-only comparison implementation
licenses/       retained upstream notices
```

`build/` and `output/` are ignored local products. The SDK remains external. The plugin's production build reads `src/native` directly and has no dependency on `dev/reference`.

Superseded 32-bit binaries, Qt3 scripts/screenshots, old paper data, historical experiment drivers and generated build/staging output were removed from the current distribution and preserved in a separate local archive. No calculation algorithm was changed during this reorganization. Previously tracked files remain recoverable from Git history.
