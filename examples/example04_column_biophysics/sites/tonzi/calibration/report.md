# Fast calibration: tonzi

## Validation

**Passed**: the calibrated cost below the default's, and no target's RMSE more than 10% worse.

| target | default | calibrated | change |
|---|---|---|---|
| lw_up | 2.25 | 2.46 | +9% |
| le | 2.06 | 1.98 | -4% |
| h | 1.44 | 1.25 | -13% |
| gpp | 2.60 | 1.89 | -27% |
| ustar | 0.51 | 0.50 | -1% |

RMSE in units of each target's sigma. Cost: default 19213.3, calibrated 18711.9.

The full record with the slow tier on is the next check: run the calibrated configs (budgets closed; dry-season GPP and LE no worse than the default's).

## Keys

| key | kind | scope | MAP | 68 % | prior centre | prior z | sd ratio | note | prior source |
|---|---|---|---|---|---|---|---|---|---|
| vcmax25 | trait | plant_type | 46.66 | 41.47-52.35 | 113.3 | -0.95 | 0.09 |  | EEO: the coordination of the Rubisco- and light-limited rates at the site's grow |
| vcmax25_pft2 | trait | plant_type | 64.21 | 61.54-66.91 | 113.3 | -0.69 | 0.04 |  | EEO: the coordination of the Rubisco- and light-limited rates at the site's grow |
| stomatal_g1 | trait | plant_type | 0.8094 | 0.7154-0.9469 | 4.45 | -5.33 | 0.75 | beyond 2 prior sd, pushed by le: diagnose it (plan §5.1) or relabel it effective | CLM5 medlynslope for deciduous broadleaf trees (Lin et al. 2015, Nature Climate  |
| stomatal_g1_pft2 | trait | plant_type | 4.173 | 3.79-4.582 | 1.92 | +1.51 | 0.19 |  | EEO: the least-cost optimum at the site's growing-season daytime climate (priors |
| stomatal_g0 | trait | plant_type | 0.002608 | 0.002066-0.003297 | 0.01 | -1.24 | 0.21 |  | the PFT file's value, a factor e either way: the residual conductance is poorly  |
| stomatal_g0_pft2 | trait | plant_type | 0.0194 | 0.01348-0.0266 | 0.01 | +0.73 | 0.42 |  | the PFT file's value, a factor e either way: the residual conductance is poorly  |
| z0m_ratio | effective | site | 0.1446 | 0.1367-0.1525 | 0.13 | +0.36 | 0.19 |  | the base value (ED2 0.13), +-0.04: closed canopies give z0/h of about 0.06-0.15  |
| wstress_sref_stomata | effective | plant_type | 0.8656 | 0.4316-1.784 | 0.6 | +0.51 | 1.00 |  | the PFT file's value, a factor 2 either way (Sabot et al. 2022 fit 0.5-5 across  |
| kappa | observation | observation | 0.8947 | 0.8577-0.9329 | 1 | -2.14 | 0.79 |  | the gap between the provider's night-time and daytime partitionings (RECO agains |

Fixed at their priors: wstress_sref_stomata_pft2 (dead: its gradient column is zero (a key the windows never use))

## Targets

| target | chi2/n | sigma scale | model/tower | n (cal) |
|---|---|---|---|---|
| lw_up | 1.07 | 2.52 | 1.026 | 2321 |
| le | 0.90 | 1.77 | 0.912 | 4780 |
| h | 1.04 | 1.48 | 1.034 | 677 |
| gpp | 0.88 | 1.66 | 0.956 | 941 |
| ustar | 0.24 | 1.00 | 0.934 | 1208 |

lw_up, model/tower by local hour: 0 h 1.02, 1 h 1.02, 2 h 1.02, 3 h 1.02, 4 h 1.02, 5 h 1.02, 6 h 1.03, 7 h 1.05, 8 h 1.05, 9 h 1.05, 10 h 1.04, 11 h 1.04, 12 h 1.03, 13 h 1.03, 14 h 1.03, 15 h 1.02, 16 h 1.02, 17 h 1.01, 18 h 1.01, 19 h 1.01, 20 h 1.02, 21 h 1.02, 22 h 1.02, 23 h 1.02
lw_up, by light: SW quartile 2 1.02, SW quartile 3 1.02, SW quartile 4 1.04

le, model/tower by local hour: 5 h -34.17, 6 h 1.36, 7 h 0.85, 8 h 0.69, 9 h 0.71, 10 h 0.73, 11 h 0.78, 12 h 0.85, 13 h 0.91, 14 h 1.04, 15 h 1.37, 16 h 1.62, 17 h 2.11, 18 h 1.69
le, by light: SW quartile 1 1.41, SW quartile 2 0.94, SW quartile 3 0.93, SW quartile 4 0.80

h, model/tower by local hour: 8 h 0.90, 9 h 0.98, 10 h 1.09, 11 h 1.09, 12 h 1.12, 13 h 1.06, 14 h 1.03, 15 h 0.97, 16 h 0.85, 17 h 0.66, 18 h 1.09
h, by light: SW quartile 1 0.87, SW quartile 2 0.98, SW quartile 3 1.03, SW quartile 4 1.09

gpp, model/tower by local hour: 6 h 1.08, 7 h 0.83, 8 h 0.84, 9 h 0.80, 10 h 0.88, 11 h 0.90, 12 h 0.94, 13 h 0.99, 14 h 1.04, 15 h 1.05, 16 h 1.16, 17 h 1.23, 18 h 1.40
gpp, by light: SW quartile 1 1.01, SW quartile 2 0.87, SW quartile 3 0.90, SW quartile 4 1.11

ustar, model/tower by local hour: 5 h 1.08, 6 h 0.80, 7 h 0.83, 8 h 0.82, 9 h 0.91, 10 h 0.96, 11 h 0.97, 12 h 0.99, 13 h 0.96, 14 h 0.96, 15 h 0.95, 16 h 0.90, 17 h 0.90, 18 h 0.99
ustar, by light: SW quartile 1 0.92, SW quartile 2 0.87, SW quartile 3 0.96, SW quartile 4 0.97

## kappa

kappa = 0.895: the tower's respiration 2.88 -> 3.22 umol m-2 s-1; GPP 4.96 -> 5.24 (1.88 -> 1.99 kgC m-2 yr-1), means over the records with GPP and its respiration measured, day and night; GPP gains the scaled respiration by day only.

## Uncertainty

Laplace, from the final gradient matrix -- LOCAL ONLY (the linearity check failed).
- alternative gpp_ustar: the fit did not use the provider's threshold, or GPP's diagnostic finds no plateau
- alternative closure: the fit scored H and LE as measured ([closure].shares = "none"): no closure factor
- alternative partitioning (the provider's daytime partitioning (GPP_DT, RECO_DT)): largest linear shift 4.11 sd -- refitted: vcmax25 46.66->37.34, vcmax25_pft2 64.21->52.05, stomatal_g1 0.8094->1.028, stomatal_g1_pft2 4.173->4.943, stomatal_g0 0.002608->0.001921, stomatal_g0_pft2 0.0194->0.02329, z0m_ratio 0.1446->0.1452, wstress_sref_stomata 0.8656->0.4583, kappa 0.8947->0.9521
