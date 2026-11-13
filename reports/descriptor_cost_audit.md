# Cached electronic-descriptor job cost audit

Computed from the successful xTB logs; no new electronic calculations or model fits. One job provides all six atomic descriptor columns for every candidate site in a structure.

| Scope | Structures | Median job seconds | Mean | 95th percentile | Maximum | Sum of job seconds |
|---|---:|---:|---:|---:|---:|---:|
| conservative_bromination | 299 | 0.143 | 0.200 | 0.523 | 1.490 | 59.806 |
| all_structures | 534 | 0.141 | 0.187 | 0.504 | 1.490 | 100.043 |

Native job timings only: excludes conformer generation, Python/process-launch/file overhead, retries and encoder inference; sums are serial-equivalent job-time sums, not observed pipeline elapsed times. No new calculations run.

Both 3D encoders and the electronic descriptors used the same pre-generated geometry. No matched end-to-end timing comparison of Atom-JEPA and xTB is available. These logs cannot establish a JEPA speedup or total training/inference cost. This is semiempirical xTB, not a DFT descriptor calculation.

## Not done

Matched cold/warm pipeline benchmarks, conformer-generation costs, inference throughput on common hardware and any speed-accuracy break-even analysis.
