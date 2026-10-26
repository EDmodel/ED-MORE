# Fast calibration: example04_column_biophysics

## Validation

**Failed**: the calibrated cost below the default's, and no target's RMSE more than 10% worse.

| target | default | calibrated | change |
|---|---|---|---|
| lw_up | 4.12 | 3.93 | -5% |
| le | 1.66 | 1.59 | -4% |
| h | 3.85 | 4.29 | +11% |
| gpp | 2.32 | 1.59 | -31% |
| ustar | 2.16 | 0.73 | -66% |

RMSE in units of each target's sigma. Cost: default 66454.2, calibrated 56971.4.

The full record with the slow tier on is the next check: run the calibrated configs (budgets closed; dry-season GPP and LE no worse than the default's).

## Keys

| key | kind | scope | MAP | 68 % | prior centre | prior z | sd ratio | note | prior source |
|---|---|---|---|---|---|---|---|---|---|
| leaf_angle_mean | trait | plant_type | 51.64 | 47.72-55.39 | 45 | +0.67 | 0.40 |  | the PFT file's value, +-10 degrees; spherical leaves average 57.3. The range is  |
| vcmax25 | trait | plant_type | 33.59 | 32.69-34.51 | 39.86 | -0.34 | 0.05 |  | EEO: the coordination of the Rubisco- and light-limited rates at the site's grow |
| stomatal_g1 | trait | plant_type | 2.876 | 2.807-2.946 | 2.798 | +0.05 | 0.05 |  | EEO: the least-cost optimum at the site's growing-season daytime climate (priors |
| stomatal_g0 | trait | plant_type | 0.0002994 | 0.0001872-0.0006151 | 0.01 | -3.62 | 0.98 | beyond 2 prior sd, pushed by le: diagnose it (plan §5.1) or relabel it effective | the PFT file's value, a factor e either way: the residual conductance is poorly  |
| z0m_ratio | effective | site | 0.04927 | 0.04763-0.051 | 0.13 | -2.86 | 0.11 |  | the base value (ED2 0.13), +-0.04: closed canopies give z0/h of about 0.06-0.15  |
| wstress_sref_stomata | effective | plant_type | 2.781 | 1.394-4.672 | 2 | +0.50 | 0.99 |  | the PFT file's value, a factor 2 either way (Sabot et al. 2022 fit 0.5-5 across  |
| kappa | observation | observation | 0.7137 | 0.6802-0.7469 | 0.65 | +0.62 | 0.33 |  | best-practice plan §3.3: the tower's respiration is about the soil chambers' alo |

## Targets

| target | chi2/n | sigma scale | model/tower | n (cal) |
|---|---|---|---|---|
| lw_up | 1.00 | 1.99 | 0.987 | 3792 |
| le | 0.98 | 1.28 | 0.956 | 6691 |
| h | 1.10 | 2.98 | 0.953 | 1782 |
| gpp | 0.97 | 1.19 | 0.993 | 1472 |
| ustar | 0.34 | 1.00 | 1.001 | 1782 |

lw_up, model/tower by local hour: 0 h 0.98, 1 h 0.98, 2 h 0.98, 3 h 0.98, 4 h 0.98, 5 h 0.98, 6 h 0.99, 7 h 1.00, 8 h 1.01, 9 h 1.01, 10 h 1.01, 11 h 1.00, 12 h 1.00, 13 h 0.99, 14 h 0.99, 15 h 0.98, 16 h 0.97, 17 h 0.97, 18 h 0.97, 19 h 0.98, 20 h 0.98, 21 h 0.98, 22 h 0.98, 23 h 0.98
lw_up, by light: SW quartile 2 0.98, SW quartile 3 0.99, SW quartile 4 1.00

le, model/tower by local hour: 6 h 1.73, 7 h 1.05, 8 h 1.03, 9 h 0.90, 10 h 0.88, 11 h 0.89, 12 h 0.89, 13 h 0.92, 14 h 0.98, 15 h 1.04, 16 h 1.15, 17 h 1.21, 18 h 1.00
le, by light: SW quartile 1 1.18, SW quartile 2 1.01, SW quartile 3 0.94, SW quartile 4 0.92

h, model/tower by local hour: 6 h 0.63, 7 h 0.06, 8 h 0.60, 9 h 0.74, 10 h 0.84, 11 h 0.94, 12 h 0.96, 13 h 0.98, 14 h 0.98, 15 h 1.02, 16 h 1.13, 17 h 2.09, 18 h -4.89
h, by light: SW quartile 1 5.62, SW quartile 2 1.00, SW quartile 3 0.93, SW quartile 4 0.89

gpp, model/tower by local hour: 6 h 0.68, 7 h 1.05, 8 h 0.99, 9 h 0.96, 10 h 0.98, 11 h 0.97, 12 h 0.95, 13 h 0.98, 14 h 1.01, 15 h 1.06, 16 h 1.08, 17 h 1.04, 18 h 0.62
gpp, by light: SW quartile 1 0.95, SW quartile 2 0.96, SW quartile 3 0.97, SW quartile 4 1.07

ustar, model/tower by local hour: 6 h 0.88, 7 h 0.86, 8 h 0.90, 9 h 0.95, 10 h 1.02, 11 h 1.05, 12 h 1.04, 13 h 1.04, 14 h 1.03, 15 h 1.02, 16 h 1.01, 17 h 1.01, 18 h 1.03
ustar, by light: SW quartile 1 1.00, SW quartile 2 1.01, SW quartile 3 1.00, SW quartile 4 0.99

## kappa

kappa = 0.714: the tower's respiration 3.28 -> 4.59 umol m-2 s-1; GPP 7.54 -> 8.32 (2.86 -> 3.15 kgC m-2 yr-1), means over the records with GPP and its respiration measured, day and night; GPP gains the scaled respiration by day only.

## Uncertainty

Laplace, from the final gradient matrix -- LOCAL ONLY (the linearity check failed).
- alternative gpp_ustar (GPP u* >= 0.325 (the fit used 0.4)): largest linear shift 0.66 sd
- alternative closure (Bowen (H and LE scaled together) against the fit's attribution shares): largest linear shift 21.74 sd -- refitted: leaf_angle_mean 51.64->45.8, vcmax25 33.59->30.64, stomatal_g1 2.876->4.75, stomatal_g0 0.0002994->0.0002149, z0m_ratio 0.04927->0.04593, wstress_sref_stomata 2.781->7.691, kappa 0.7137->0.6972
- alternative partitioning: the site TOML declares no daytime partitioning (GPP_DT, RECO_DT)
