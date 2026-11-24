# Numerical contract

The production evaluators operate on the selected informed sample for each output. Moments use original values. Cumulants of order at least two use within-sample centering; first-order outputs retain the mean. All submoments use the target sample and its normalized weights.

Template source weights are converted to geometric-mean event weights over distinct records. Reusing a record in several template components counts its weight once. Point-set coincident records remain distinct observations. Moving windows use source-record weights directly and normalize locally.

Weight sources are uniform, three declustering methods, or an external observation property. Imported weights are finite and strictly positive on informed observations. The code does not clip, take absolute values or replace invalid external weights. ESS is (sum w)^2 / sum(w^2); sample count and ESS are different diagnostics.

Template lag maps and window spatial maps have different axes. Cartesian windows require complete geometric windows; default SGeMS processing also requires complete informed data. Point windows can evaluate neighborhoods at source records or on a regular grid of centers. Point neighborhoods with insufficient informed samples leave uninformed outputs.

Point regular-center grids accept box and ellipsoid windows. The axis-aligned bounding extent of a rotated window determines whether it fits completely inside the requested bounds. The sampling shape then selects observations within that window.

Supported orders depend on the evaluator and occurrence configuration. Direction/multiplicity validation, evaluator index limits and available resources apply. Standardized cumulants, when requested, use unclipped statistics and variance on identical samples and weights.

SGeMS properties store float32 values; double-precision results and resolved input/configuration are preserved in the run directory. The optional development reference is checked against the DLLs on the same exported inputs; it is not a plugin runtime dependency.
