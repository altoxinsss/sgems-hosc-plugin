# External observation weights

Load weights as a property on the same SGeMS grid or point set as the statistical values. Select **External weights**, then **Weight property**. For templates, keep **Enable selected weighting method** checked.

Every informed statistical value must have a finite, strictly positive weight on the same source record. A missing statistical value excludes that record; an invalid weight on an informed value stops the calculation. Coincident point-set records retain their separate values and weights. Imported weights are not matched by coordinates or redistributed by a declustering algorithm.

Cell, domain and kriging controls are unused in external mode. Template events use the geometric mean of weights over distinct source records; repeated occurrences do not repeat a record's weight. Moving windows use the weights of their contributing observations. Weights are normalized within the target calculation and held fixed for all its moments.

Multiplying all weights by one positive constant leaves the calculated statistics unchanged within numerical tolerance. All-one weights reproduce uniform weighting on identical support.

The two ready-to-load examples are `grid_template_external.par` and `grid_window_external.par` in [examples/Parameters](../examples/Parameters). Both use `grid.sgems`, statistical property `value`, and weight property `weight`.
