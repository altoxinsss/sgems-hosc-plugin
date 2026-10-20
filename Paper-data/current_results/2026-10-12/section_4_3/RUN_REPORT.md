# R1-11 scalability run report

## Protocol

- Stanford VI-E z=62 slice; three-direction order series with NLAG=(6,20,20), LAG=(3,3,3), and AZM=(180,90,0).
- Direction 1 is fixed at its first lag; directions 2 and 3 form the 20 x 20 output section (400 points).
- One warm-up and the configured number of recorded runs are executed for each Bell and AVL route.
- Peak memory means the maximum recorded process working set, reported in MiB.

## Fourth-order workload and moving-window cases

| Case | Bell mean (s) | AVL mean (s) | Bell/AVL | Bell/AVL peak working set (MiB) | Bell/AVL agreement |
|---|---:|---:|---:|---:|:---:|
| Template r=4 10x10 | 0.4054 | 0.2488 | 1.6298 | 17.457 / 5.391 | True |
| Template r=4 15x15 | 0.8371 | 0.5228 | 1.6012 | 17.547 / 5.480 | True |
| Template r=4 20x20 | 1.4228 | 0.8707 | 1.6341 | 17.555 / 5.531 | True |
| Template r=4 25x25 | 2.2108 | 1.4780 | 1.4958 | 17.629 / 5.691 | True |
| Moving window r=4 5x5 | 0.0644 | 0.1337 | 0.4812 | 4.375 / 4.773 | True |
| Moving window r=4 9x9 | 0.0633 | 0.2365 | 0.2678 | 4.297 / 4.473 | True |
| Moving window r=4 13x13 | 0.0655 | 0.4136 | 0.1584 | 4.277 / 4.391 | True |

## Fixed-workload order extension

| Order | Occurrences by direction | Bell mean (s) | AVL mean (s) | Bell/AVL | Bell/AVL peak working set (MiB) | Bell/AVL agreement |
|---:|:---:|---:|---:|---:|---:|:---:|
| 4 | (1,1,1) | 1.4228 | 0.8707 | 1.6341 | 17.555 / 5.531 | True |
| 5 | (2,1,1) | 2.9430 | 1.0219 | 2.8798 | 17.586 / 5.715 | True |
| 6 | (2,2,1) | 10.6654 | 1.2446 | 8.5691 | 17.578 / 6.141 | True |
| 7 | (2,2,2) | 49.1277 | 1.7725 | 27.7167 | 17.582 / 7.645 | True |
| 8 | (3,2,2) | 272.7872 | 2.1774 | 125.2834 | 17.566 / 9.055 | True |

## Sixth-order template structure

| Study | Directions | Occurrences | Bell mean (s) | AVL mean (s) | Bell/AVL | Bell/AVL peak working set (MiB) | Bell/AVL agreement |
|---|---:|:---:|---:|---:|---:|---:|:---:|
| Direction configuration | 2 | (3,2) | 9.3935 | 1.0428 | 9.0078 | 17.480 / 5.863 | True |
| Direction configuration | 3 | (2,2,1) | 10.6654 | 1.2446 | 8.5691 | 17.578 / 6.141 | True |
| Direction configuration | 4 | (2,1,1,1) | 11.7054 | 1.4059 | 8.3259 | 18.172 / 7.090 | True |
| Occurrence configuration | 3 | (3,1,1) | 11.1304 | 1.1622 | 9.5772 | 17.578 / 6.012 | True |
| Occurrence configuration | 3 | (1,3,1) | 11.2468 | 1.1795 | 9.5351 | 17.582 / 6.047 | True |
| Occurrence configuration | 3 | (2,2,1) | 10.6654 | 1.2446 | 8.5691 | 17.578 / 6.141 | True |
| Occurrence configuration | 3 | (1,1,3) | 11.1440 | 1.3276 | 8.3942 | 17.566 / 6.414 | True |

## Validity gates

- Cases: 16; routes: 32; recorded executions: 160.
- Process gate: True; numerical gate: True; overall: True.
- Largest absolute Bell/AVL discrepancy: 1.5487044528479221e-17.
- Largest relative L2 discrepancy: 5.9473251644089391e-14.

Detailed provenance is recorded in environment.json, case_matrix.csv, raw_runs.csv, 
route_summary.csv, pair_summary.csv, bell_avl_accuracy.csv, gate.json, and sha256_manifest.csv.
