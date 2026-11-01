# Manual SGeMS integration checks

Use the small files in [tests/fixtures/manual](../tests/fixtures/manual) before
running a full reservoir model. These are synthetic software fixtures, not
experimental evidence for the paper. Run outcomes must be recorded by the
tester; this document is not a host acceptance result.

## Import and interface

Load `selftest_2d.gslib` as a **Cartesian grid** named `selftest_2d`, with
dimensions 5,5,1. Load `selftest_3d.gslib` as `selftest_3d`, with dimensions
5,5,3. Set cell dimensions to 1,1,1, origin to 0,0,0, and enable the no-data
value -999. X varies fastest. Do not import these as point sets.

Both files contain `phi`, `constant`, `phi_shift`, `phi_missing0`, `locator`,
and `spike`. Their definitions, using zero-based indices, are:

- phi = (1+x+2y+3z)/32; phi_shift = phi+1.
- constant = 0.25; phi_missing0 replaces only the first phi with no-data.
- locator = x+10y+100z; its first zero is a valid observation.
- spike = 0.125 except the central node, which is 0.5.

Check that both algorithms open. In Window, scroll to the audit-directory
input. In Template, the Exact Template fields should stay together near the
top when the panel is resized. Save and reload a set of nondefault parameters
and confirm that all values, including the final field, are restored. The
Window interface must not show inversion, normalization, legacy centering,
or coordinate-center-mode controls. Template must not show a weighting tab.

## Moving window

Select HOSC_WINDOW, Cumulant, Bell, selftest_3d/phi. Enter a unique Target
Property. Set order=4, window=3,3,3, steps=1,1,1. Disable weighting and use
uniform; input centering=true,
standardized output=true, epsilon=1e-12, diagnostics=true, clipping=false.
Enable the audit sidecar and select an existing writable test directory.

Run and inspect the new result grid and its CSV. It should have dimensions
3,3,1 and nine window centers, from (1,1,1) to (3,3,1), in X-fastest order.
Every node should have support=ESS=27, kappa2=0.009114583333333334,
kappa4=-0.00006230672200520833, lambda4=-0.75, and validity=1.

Repeat these actions with a distinct output name each time:

| Action | Expected check |
|---|---|
| Switch Bell to AVL | Same values, supports, and centers |
| Enable weighting with uniform | Same values through the shared weighted path |
| Set coordinate_gaussian, sigma=1, min_weight=0 | ESS about 22.4642796081; compare all rows with expected_window.csv |
| Set soft_high_porosity, threshold=.25, width=.02, strength=.8 | Compare all rows with expected_window.csv |
| Set soft strength=0 | Recovers uniform results |
| Select spike and run orders1–8 | Compare nonzero odd/even cumulants with the reference; disable lambda for order1 |
| Select phi_shift, uniform, order6 | Same kappa6/lambda6 as phi; order1 increases by one |
| Select constant, order6 | Raw result and kappa2 zero, numeric_valid=1, lambda uninformed and lambda_valid=0 |
| Select phi_missing0 | First result uninformed, support26, numeric_valid=0; other eight results unchanged |
| Disable standardized output | No kappa2/lambda/lambda_valid properties; audit marks lambda_not_requested |
| Disable diagnostics | No support/ess/numeric_valid/display_valid properties |
| Disable audit | No new sidecar directory, normal property output |
| Use selftest_2d with window3,3,1 and steps1,1,0 | Nine centers at Z=0, support9, kappa4=-0.00001080830891927083 |
| Use 3D steps2,2,1 | Four centers, X/Y in {1,3} |
| Select Moment, order2, lambda=false | First raw moment=0.05696614583333333 |

For a soft-weight translation check, shift the threshold by the same amount as
the values. A fixed value threshold changes the weights after a data shift.

The optional clipping output is a separate `<output>_clipped` property. With
phi, Moment order1, and percentiles25/75, it should equal
`[9,9,9,9,10,11,11,11,11]/32`; the primary values remain
`[7,8,9,9,10,11,11,12,13]/32`.

Reference CSV filters: property, weight_mode, order, node_id. Settings are
3D/window3/step1, sigma1, threshold.25, width.02, strength.8, min_weight0.
Suggested tolerances for these small fixtures are
`1e-12 + 1e-6*abs(expected)` for raw/kappa2 and
`1e-6 + 1e-6*abs(expected)` for lambda/ESS. Supports and flags must match exactly.

## Exact templates

Select HOSC_TEMPLATE, Cumulant, Bell, property phi. Disable 3D Slice.
Choose exhaustive_integer_offsets. Use order
occurrences NOD=1,1,1, LAG=1,1,1. Fill all seven remaining angle/tolerance
lists (LTOL, AZM, ATOL, BANDWH, DIP, DTOL, BANDWD) with 0,0,0.

| Reference case | Source | NLAG | Integer vectors | Display X/Y/Z | Position |
|---|---|---|---|---|---|
| T01_2d | selftest_2d | 2,3,1 | 1,0;0,1;1,1 | 1/2/0 | 1 |
| T03_3d | selftest_3d | 2,3,2 | 1,0,0;0,1,0;0,0,1 | 1/2/3 | empty |
| T04_singleton | selftest_3d | 2,3,1 | 1,0,0;0,1,0;0,0,1 | 1/2/3 | empty |

Run each case with Bell and AVL. Compare the main property and support against
`expected_template.csv`. Numeric/display validity should be one. Suggested
absolute tolerance for the float kappa4 property is 1e-9. Output axes describe
lag combinations, not reservoir node positions. Template does not emit Window
ESS or lambda properties.

Additional actions:

1. T01 with display2/1/0 gives a 3x2 transpose; support becomes
   25,20,15,20,16,12 in X-fastest order.
2. T03 with display1/2/0 and Position2 equals its second Z-lag layer.
3. T01 with constant yields valid zeros; phi_shift yields unchanged kappa4.
4. T01 with phi_missing0 gives supports24,19,19,15,14,11, without shifting
   coordinates. Values change because the common support changed.
5. T01 with LAG5,1,1 has no support in the three X-lag-index1 cells: leave
   results uninformed and set support/validity to zero.
6. T01 with NOD3,2,2 tests order8; NOD4,2,2 must be rejected as order9.

## Errors and persistence

Change only one setting per negative test and restore the working settings
afterward. Check for a readable error, preservation of source/existing results,
and the ability to run the valid case again.

- Window: order0/9/3.5; window larger than source; stepX0; 3D stepZ0;
  2D stepZ1; nonuniform weighting while disabled; sigma0; soft width0 or
  strength1.1; negative epsilon; Moment with
  lambda enabled; blank output name or required audit path.
- Template: 3D legacy_search; exact Moment; inconsistent
  list lengths; zero LAG; zero/fractional/wrong-count direction vectors;
  duplicate display directions; invalid hidden Position; order above8.
- Both algorithms: point-set input is unsupported.

For legacy-file compatibility checks, submit XML directly through a SGeMS
command file, because removed fields are no longer represented by UI widgets.
Do not rely on Load Parameters retaining an unknown XML element.

| Parameter-file action | Expected check |
|---|---|
| Run the current window examples without the four fixed fields | Same values as the reference fixtures |
| Add `normalize_weights=true`, `invert_window_output=false`, `zero_first_order_moment=false`, and `weight_center_mode=window_center` | Same values as the example with the fields absent |
| Set `normalize_weights=false`, with weighting both enabled and disabled | Rejected before calculation in both cases |
| Set `invert_window_output=true` or `zero_first_order_moment=true`, one at a time | Rejected before calculation |
| Set `weight_center_mode=parameter`, with uniform and coordinate-Gaussian modes | Rejected before calculation in both cases |
| Run a template with `Weighted_Enable` absent, then with `Disable` | Same results |
| Set `Weighted_Enable=Enable` for a valid 2D legacy template, then a valid exact template | Rejected before calculation in both cases, including with `Weight_Mode=uniform` |

After each rejected legacy-file case, submit the valid parameters again in
the same process and verify the expected numerical output.

Save parameters, alter them, reload, and rerun. Save a project containing valid
and uninformed results, close SGeMS, reopen, and check masks/values/coordinates.
Run the same case repeatedly with separate output names; previous outputs
must remain intact. Record any empty object left after a failed run.

## Slice utility

Enable 3D Slice only for this separate operation: it bypasses template
calculation. Select selftest_3d/locator in the slice source selectors.
XY/position2 gives `100+x+10y`; XZ/position3 gives `20+x+100z`;
YZ/position4 gives `3+10y+100z`. Positions are one-based.

The current slice implementation flattens XZ/YZ to a 2D display plane and
uses unit output spacing. It is not a preservation of source physical
geometry. The missing-value path also needs host verification: extracting
phi_missing0 at XY/position1 should preserve the uninformed first node; an
informed sentinel is a defect, not a valid result.

For issues, retain the saved parameter file, output CSV/GSLIB, case ID, and
error screenshot. For Window, include the associated resolved_config.txt.
