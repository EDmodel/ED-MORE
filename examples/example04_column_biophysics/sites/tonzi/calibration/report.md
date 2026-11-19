# Fast calibration: tonzi

## Validation

**Passed**: the calibrated cost below the default's, and no target's RMSE more than 10% worse.

| target | default | calibrated | change |
|---|---|---|---|
| lw_up | 3.08 | 3.29 | +7% |
| le | 2.15 | 1.98 | -8% |
| h | 2.28 | 2.11 | -8% |
| gpp | 2.57 | 1.49 | -42% |
| ustar | 0.53 | 0.53 | +0% |

RMSE in units of each target's sigma. Cost: default 27712.7, calibrated 26073.1.

The full record with the slow tier on is the next check: run the calibrated configs (budgets closed; dry-season GPP and LE no worse than the default's).

## Keys

| key | kind | scope | MAP | 68 % | prior centre | prior z | sd ratio | note | prior source |
|---|---|---|---|---|---|---|---|---|---|
| vcmax25 | trait | plant_type | 58.22 | 52.61-64.13 | 113.3 | -0.77 | 0.08 |  | EEO: the coordination of the Rubisco- and light-limited rates at the site's grow |
| vcmax25_pft2 | trait | plant_type | 53.78 | 51.15-56.49 | 113.3 | -0.84 | 0.04 |  | EEO: the coordination of the Rubisco- and light-limited rates at the site's grow |
| stomatal_g1 | trait | plant_type | 2.057 | 1.822-2.33 | 4.45 | -2.07 | 0.33 | beyond 2 prior sd, pushed by le: diagnose it (plan §5.1) or relabel it effective | CLM5 medlynslope for deciduous broadleaf trees (Lin et al. 2015, Nature Climate  |
| stomatal_g1_pft2 | trait | plant_type | 3.077 | 2.673-3.536 | 1.92 | +0.91 | 0.27 |  | EEO: the least-cost optimum at the site's growing-season daytime climate (priors |
| stomatal_g0 | trait | plant_type | 0.006078 | 0.005328-0.006922 | 0.01 | -0.48 | 0.12 |  | the PFT file's value, a factor e either way: the residual conductance is poorly  |
| stomatal_g0_pft2 | trait | plant_type | 0.009065 | 0.003467-0.02037 | 0.01 | -0.10 | 0.88 |  | the PFT file's value, a factor e either way: the residual conductance is poorly  |
| z0m_ratio | effective | site | 0.1525 | 0.1438-0.1614 | 0.13 | +0.55 | 0.21 |  | the base value (ED2 0.13), +-0.04: closed canopies give z0/h of about 0.06-0.15  |
| wstress_sref_stomata | effective | plant_type | 6.557 | 5.412-7.513 | 0.6 | +4.29 | 0.66 |  | the PFT file's value, a factor 2 either way (Sabot et al. 2022 fit 0.5-5 across  |
| kappa | observation | observation | 0.9161 | 0.8763-0.9571 | 1 | -1.69 | 0.83 |  | the gap between the provider's night-time and daytime partitionings (RECO agains |

Fixed at their priors: wstress_sref_stomata_pft2 (dead: its gradient column is zero (a key the windows never use))

## Targets

| target | chi2/n | sigma scale | model/tower | n (cal) |
|---|---|---|---|---|
| lw_up | 1.08 | 3.00 | 1.032 | 2321 |
| le | 0.96 | 1.88 | 0.953 | 4780 |
| h | 1.03 | 2.08 | 0.947 | 677 |
| gpp | 0.71 | 1.93 | 1.048 | 941 |
| ustar | 0.28 | 1.00 | 0.932 | 1208 |

lw_up, model/tower by local hour: 0 h 1.04, 1 h 1.05, 2 h 1.04, 3 h 1.04, 4 h 1.04, 5 h 1.04, 6 h 1.03, 7 h 1.02, 8 h 1.01, 9 h 1.00, 10 h 1.00, 11 h 1.00, 12 h 1.00, 13 h 1.01, 14 h 1.02, 15 h 1.03, 16 h 1.04, 17 h 1.05, 18 h 1.05, 19 h 1.05, 20 h 1.05, 21 h 1.05, 22 h 1.05, 23 h 1.05
lw_up, by light: SW quartile 2 1.05, SW quartile 3 1.03, SW quartile 4 1.01

le, model/tower by local hour: 5 h -57.53, 6 h 1.49, 7 h 0.81, 8 h 0.66, 9 h 0.72, 10 h 0.73, 11 h 0.80, 12 h 0.88, 13 h 0.98, 14 h 1.11, 15 h 1.47, 16 h 1.78, 17 h 2.45, 18 h 2.47
le, by light: SW quartile 1 1.40, SW quartile 2 0.93, SW quartile 3 1.01, SW quartile 4 0.86

h, model/tower by local hour: 8 h 0.56, 9 h 0.69, 10 h 0.81, 11 h 0.88, 12 h 0.95, 13 h 0.97, 14 h 1.01, 15 h 1.11, 16 h 1.22, 17 h 1.54, 18 h 4.07
h, by light: SW quartile 1 1.52, SW quartile 2 0.90, SW quartile 3 0.89, SW quartile 4 0.93

gpp, model/tower by local hour: 6 h 1.08, 7 h 0.89, 8 h 0.86, 9 h 0.80, 10 h 0.93, 11 h 1.02, 12 h 1.07, 13 h 1.14, 14 h 1.18, 15 h 1.19, 16 h 1.26, 17 h 1.27, 18 h 1.60
gpp, by light: SW quartile 1 1.19, SW quartile 2 1.01, SW quartile 3 1.04, SW quartile 4 0.99

ustar, model/tower by local hour: 5 h 1.23, 6 h 0.79, 7 h 0.80, 8 h 0.76, 9 h 0.83, 10 h 0.89, 11 h 0.93, 12 h 0.96, 13 h 0.95, 14 h 0.97, 15 h 0.98, 16 h 0.97, 17 h 1.02, 18 h 1.15
ustar, by light: SW quartile 1 1.03, SW quartile 2 0.87, SW quartile 3 0.92, SW quartile 4 0.94

## kappa

kappa = 0.916: the tower's respiration 2.88 -> 3.14 umol m-2 s-1; GPP 4.96 -> 5.18 (1.88 -> 1.96 kgC m-2 yr-1), means over the records with GPP and its respiration measured, day and night; GPP gains the scaled respiration by day only.

## Uncertainty

Laplace, from the final gradient matrix.
- alternative gpp_ustar: the fit did not use the provider's threshold, or GPP's diagnostic finds no plateau
- alternative closure: the fit scored H and LE as measured ([closure].shares = "none"): no closure factor
- alternative partitioning (the provider's daytime partitioning (GPP_DT, RECO_DT)): largest linear shift 4.16 sd -- refitted: vcmax25 58.22->59.41, vcmax25_pft2 53.78->41.57, stomatal_g1 2.057->1.671, stomatal_g1_pft2 3.077->4.157, stomatal_g0 0.006078->0.006202, stomatal_g0_pft2 0.009065->0.002937, z0m_ratio 0.1525->0.1527, wstress_sref_stomata 6.557->6.295, kappa 0.9161->1.017
