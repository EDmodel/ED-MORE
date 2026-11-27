# Fast calibration: bci

## Validation

**Failed**: the calibrated cost below the default's, and no target's RMSE more than 10% worse.

| target | default | calibrated | change |
|---|---|---|---|
| lw_up | 4.16 | 4.03 | -3% |
| le | 1.61 | 1.57 | -3% |
| h | 3.84 | 4.24 | +10% |
| gpp | 2.34 | 1.59 | -32% |
| ustar | 2.17 | 0.74 | -66% |

RMSE in units of each target's sigma. Cost: default 66635.7, calibrated 58205.

The full record with the slow tier on is the next check: run the calibrated configs (budgets closed; dry-season GPP and LE no worse than the default's).

## Keys

| key | kind | scope | MAP | 68 % | prior centre | prior z | sd ratio | note | prior source |
|---|---|---|---|---|---|---|---|---|---|
| leaf_angle_mean | trait | plant_type | 43.08 | 38.52-47.71 | 45 | -0.19 | 0.46 |  | the PFT file's value, +-10 degrees; spherical leaves average 57.3. The range is  |
| vcmax25 | trait | plant_type | 34.43 | 33.47-35.43 | 44.65 | -0.52 | 0.06 |  | EEO: the coordination of the Rubisco- and light-limited rates at the site's grow |
| stomatal_g1 | trait | plant_type | 2.99 | 2.912-3.07 | 2.78 | +0.15 | 0.05 |  | EEO: the least-cost optimum at the site's growing-season daytime climate (priors |
| stomatal_g0 | trait | plant_type | 0.0005163 | 0.0002769-0.0012 | 0.01 | -2.85 | 0.90 | beyond 2 prior sd, pushed by le: diagnose it (plan §5.1) or relabel it effective | the PFT file's value, a factor e either way: the residual conductance is poorly  |
| z0m_ratio | effective | site | 0.04948 | 0.04781-0.05124 | 0.13 | -2.85 | 0.11 |  | the base value (ED2 0.13), +-0.04: closed canopies give z0/h of about 0.06-0.15  |
| wstress_sref_stomata | effective | plant_type | 3.416 | 1.825-5.334 | 2 | +0.85 | 0.98 |  | the PFT file's value, a factor 2 either way (Sabot et al. 2022 fit 0.5-5 across  |
| kappa | observation | observation | 0.6918 | 0.6599-0.7239 | 0.65 | +0.41 | 0.31 |  | best-practice plan §3.3: the tower's respiration is about the soil chambers' alo |

## Targets

| target | chi2/n | sigma scale | model/tower | n (cal) |
|---|---|---|---|---|
| lw_up | 1.05 | 2.15 | 0.987 | 3792 |
| le | 0.98 | 1.28 | 0.944 | 6691 |
| h | 1.11 | 2.89 | 0.963 | 1782 |
| gpp | 0.95 | 1.21 | 0.992 | 1472 |
| ustar | 0.35 | 1.00 | 1.002 | 1782 |

lw_up, model/tower by local hour: 0 h 0.98, 1 h 0.98, 2 h 0.98, 3 h 0.98, 4 h 0.98, 5 h 0.98, 6 h 0.99, 7 h 1.00, 8 h 1.01, 9 h 1.01, 10 h 1.01, 11 h 1.01, 12 h 1.00, 13 h 0.99, 14 h 0.99, 15 h 0.98, 16 h 0.97, 17 h 0.97, 18 h 0.97, 19 h 0.98, 20 h 0.98, 21 h 0.98, 22 h 0.98, 23 h 0.98
lw_up, by light: SW quartile 2 0.98, SW quartile 3 0.99, SW quartile 4 1.00

le, model/tower by local hour: 6 h 1.66, 7 h 1.05, 8 h 1.04, 9 h 0.90, 10 h 0.88, 11 h 0.88, 12 h 0.87, 13 h 0.90, 14 h 0.96, 15 h 1.01, 16 h 1.13, 17 h 1.19, 18 h 0.93
le, by light: SW quartile 1 1.16, SW quartile 2 1.00, SW quartile 3 0.93, SW quartile 4 0.90

h, model/tower by local hour: 6 h 0.68, 7 h -0.01, 8 h 0.60, 9 h 0.75, 10 h 0.86, 11 h 0.97, 12 h 0.98, 13 h 1.00, 14 h 0.99, 15 h 1.02, 16 h 1.10, 17 h 2.00, 18 h -4.51
h, by light: SW quartile 1 5.28, SW quartile 2 0.99, SW quartile 3 0.94, SW quartile 4 0.91

gpp, model/tower by local hour: 6 h 0.67, 7 h 1.06, 8 h 1.00, 9 h 0.97, 10 h 0.97, 11 h 0.96, 12 h 0.94, 13 h 0.97, 14 h 1.00, 15 h 1.07, 16 h 1.09, 17 h 1.04, 18 h 0.60
gpp, by light: SW quartile 1 0.95, SW quartile 2 0.97, SW quartile 3 0.97, SW quartile 4 1.06

ustar, model/tower by local hour: 6 h 0.87, 7 h 0.85, 8 h 0.90, 9 h 0.95, 10 h 1.03, 11 h 1.06, 12 h 1.05, 13 h 1.05, 14 h 1.03, 15 h 1.02, 16 h 1.01, 17 h 1.01, 18 h 1.03
ustar, by light: SW quartile 1 1.00, SW quartile 2 1.01, SW quartile 3 1.01, SW quartile 4 1.00

## kappa

kappa = 0.692: the tower's respiration 3.28 -> 4.73 umol m-2 s-1; GPP 7.54 -> 8.41 (2.86 -> 3.19 kgC m-2 yr-1), means over the records with GPP and its respiration measured, day and night; GPP gains the scaled respiration by day only.

## Uncertainty

Laplace, from the final gradient matrix -- LOCAL ONLY (the linearity check failed).
- alternative gpp_ustar (GPP u* >= 0.325 (the fit used 0.4)): largest linear shift 0.69 sd
- alternative closure (Bowen (H and LE scaled together) against the fit's attribution shares): largest linear shift 20.38 sd -- refitted: leaf_angle_mean 43.08->35.59, vcmax25 34.43->30.84, stomatal_g1 2.99->4.825, stomatal_g0 0.0005163->0.0003521, z0m_ratio 0.04948->0.04713, wstress_sref_stomata 3.416->5.101, kappa 0.6918->0.6751
- alternative partitioning: the site TOML declares no daytime partitioning (GPP_DT, RECO_DT)
