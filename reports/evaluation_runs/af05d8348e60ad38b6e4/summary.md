# Linear frozen-feature evaluation

Split proposal: `af05d8348e60ad38b6e4`.

| Block / label budget | Molecules | Top-1 | Top-2 | Top-2 (k>2) | MRR | NLL | Uniform NLL | Information gain |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| atom_jepa:budget50 | 299 | 0.6455 | 0.8606 | 0.8292 | 0.7955 | 0.8755 | 1.3525 | 0.4770 |
| atom_jepa:budget100 | 299 | 0.7135 | 0.9064 | 0.8852 | 0.8382 | 0.7343 | 1.3525 | 0.6182 |
| atom_jepa:budgetfull | 299 | 0.7559 | 0.9298 | 0.9139 | 0.8649 | 0.6167 | 1.3525 | 0.7358 |
| random_backbone:budget50 | 299 | 0.6143 | 0.8439 | 0.8087 | 0.7749 | 1.0289 | 1.3525 | 0.3236 |
| random_backbone:budget100 | 299 | 0.6734 | 0.9030 | 0.8811 | 0.8171 | 0.8550 | 1.3525 | 0.4975 |
| random_backbone:budgetfull | 299 | 0.7458 | 0.9465 | 0.9344 | 0.8627 | 0.6651 | 1.3525 | 0.6875 |
| morgan:budget50 | 299 | 0.4406 | 0.6983 | 0.6302 | 0.6537 | 1.2474 | 1.3525 | 0.1051 |
| morgan:budget100 | 299 | 0.5134 | 0.7897 | 0.7423 | 0.7120 | 1.1259 | 1.3525 | 0.2266 |
| morgan:budgetfull | 299 | 0.6355 | 0.9047 | 0.8832 | 0.7967 | 0.9624 | 1.3525 | 0.3902 |
| electronic:budget50 | 299 | 0.8261 | 0.9331 | 0.9180 | 0.8991 | 0.5581 | 1.3525 | 0.7944 |
| electronic:budget100 | 299 | 0.8528 | 0.9554 | 0.9454 | 0.9168 | 0.5126 | 1.3525 | 0.8400 |
| electronic:budgetfull | 299 | 0.8562 | 0.9699 | 0.9631 | 0.9211 | 0.4551 | 1.3525 | 0.8974 |
| jepa_plus_qm:budget50 | 299 | 0.8294 | 0.9521 | 0.9413 | 0.9052 | 0.5405 | 1.3525 | 0.8120 |
| jepa_plus_qm:budget100 | 299 | 0.8595 | 0.9688 | 0.9617 | 0.9236 | 0.4231 | 1.3525 | 0.9295 |
| jepa_plus_qm:budgetfull | 299 | 0.8696 | 0.9666 | 0.9590 | 0.9279 | 0.4017 | 1.3525 | 0.9508 |

## Paired top-1 differences

| Left / right | Budget | Difference | Conditional 95% interval |
|---|---|---:|---|
| atom_jepa / random_backbone | 50 | 0.0312 | [-0.0195, 0.0831] |
| atom_jepa / morgan | 50 | 0.2048 | [0.1494, 0.2590] |
| atom_jepa / electronic | 50 | -0.1806 | [-0.2347, -0.1267] |
| jepa_plus_qm / atom_jepa | 50 | 0.1839 | [0.1361, 0.2317] |
| jepa_plus_qm / electronic | 50 | 0.0033 | [-0.0206, 0.0251] |
| atom_jepa / random_backbone | 100 | 0.0401 | [-0.0100, 0.0886] |
| atom_jepa / morgan | 100 | 0.2001 | [0.1435, 0.2545] |
| atom_jepa / electronic | 100 | -0.1394 | [-0.1911, -0.0895] |
| jepa_plus_qm / atom_jepa | 100 | 0.1460 | [0.1050, 0.1889] |
| jepa_plus_qm / electronic | 100 | 0.0067 | [-0.0211, 0.0349] |
| atom_jepa / random_backbone | full | 0.0100 | [-0.0473, 0.0686] |
| atom_jepa / morgan | full | 0.1204 | [0.0561, 0.1864] |
| atom_jepa / electronic | full | -0.1003 | [-0.1611, -0.0393] |
| jepa_plus_qm / atom_jepa | full | 0.1137 | [0.0621, 0.1672] |
| jepa_plus_qm / electronic | full | 0.0134 | [-0.0237, 0.0512] |

## Regularization-grid diagnostics

| Block / budget | Final fits selecting a grid edge |
|---|---:|
| atom_jepa / 50 | 0 / 15 |
| atom_jepa / 100 | 0 / 15 |
| atom_jepa / full | 0 / 5 |
| random_backbone / 50 | 7 / 15 |
| random_backbone / 100 | 0 / 15 |
| random_backbone / full | 0 / 5 |
| morgan / 50 | 2 / 15 |
| morgan / 100 | 0 / 15 |
| morgan / full | 0 / 5 |
| electronic / 50 | 1 / 15 |
| electronic / 100 | 0 / 15 |
| electronic / full | 2 / 5 |
| jepa_plus_qm / 50 | 6 / 15 |
| jepa_plus_qm / 100 | 2 / 15 |
| jepa_plus_qm / full | 1 / 5 |

Edge selections are diagnostic; the grid is not automatically changed using outer-test performance.

## JEPA + QM score components

Relative RMS measures score amplitude, not independent information or product fractions; covariance is retained in JSON.

| Fold / budget / subset seed | JEPA RMS | QM RMS | QM relative RMS |
|---|---:|---:|---:|
| 0 / 50 / 17 | 0.6445 | 2.5123 | 0.7958 |
| 0 / 100 / 17 | 0.8488 | 2.5202 | 0.7480 |
| 0 / 50 / 29 | 0.8925 | 3.7773 | 0.8089 |
| 0 / 100 / 29 | 0.3934 | 2.9019 | 0.8806 |
| 0 / 50 / 43 | 0.7207 | 2.3421 | 0.7647 |
| 0 / 100 / 43 | 0.3356 | 2.5890 | 0.8852 |
| 0 / full / None | 1.0711 | 1.1372 | 0.5150 |
| 1 / 50 / 17 | 0.9754 | 2.1485 | 0.6878 |
| 1 / 100 / 17 | 1.0499 | 2.3282 | 0.6892 |
| 1 / 50 / 29 | 0.3822 | 1.5321 | 0.8004 |
| 1 / 100 / 29 | 1.5602 | 1.8310 | 0.5399 |
| 1 / 50 / 43 | 0.6721 | 2.5594 | 0.7920 |
| 1 / 100 / 43 | 0.6953 | 2.7481 | 0.7981 |
| 1 / full / None | 0.9808 | 2.3731 | 0.7076 |
| 2 / 50 / 17 | 0.4515 | 1.3395 | 0.7479 |
| 2 / 100 / 17 | 0.9054 | 2.3367 | 0.7207 |
| 2 / 50 / 29 | 0.4697 | 1.1949 | 0.7178 |
| 2 / 100 / 29 | 1.0413 | 2.1218 | 0.6708 |
| 2 / 50 / 43 | 0.3915 | 1.2728 | 0.7648 |
| 2 / 100 / 43 | 1.0422 | 1.8368 | 0.6380 |
| 2 / full / None | 0.9665 | 1.9894 | 0.6730 |
| 3 / 50 / 17 | 0.6950 | 6.9150 | 0.9087 |
| 3 / 100 / 17 | 1.8336 | 1.6095 | 0.4675 |
| 3 / 50 / 29 | 1.0027 | 3.5736 | 0.7809 |
| 3 / 100 / 29 | 0.8310 | 2.6108 | 0.7585 |
| 3 / 50 / 43 | 0.6207 | 2.2710 | 0.7853 |
| 3 / 100 / 43 | 0.9766 | 2.2263 | 0.6951 |
| 3 / full / None | 0.9978 | 1.9055 | 0.6563 |
| 4 / 50 / 17 | 1.3018 | 4.9943 | 0.7932 |
| 4 / 100 / 17 | 0.9097 | 2.6159 | 0.7420 |
| 4 / 50 / 29 | 0.7236 | 2.6521 | 0.7856 |
| 4 / 100 / 29 | 0.9466 | 2.7823 | 0.7461 |
| 4 / 50 / 43 | 1.4421 | 1.5386 | 0.5162 |
| 4 / 100 / 43 | 0.9205 | 2.5165 | 0.7322 |
| 4 / full / None | 0.3197 | 2.9937 | 0.9035 |

Exact ties receive expected metrics. Repetition metrics are averaged per molecule; logits are not ensembled.
Paired confidence intervals resample leakage components and are conditional on the fitted models.
Intervals are computed for the primary top-1 comparison; secondary metrics are descriptive.
All source records are retained. Bibliographic uncertainty and pretraining overlap limit generalization claims.

Full diagnostics, molecular predictions and paired intervals are in summary.json and per-fit JSON files.
