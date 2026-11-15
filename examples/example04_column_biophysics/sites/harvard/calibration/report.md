# Fast calibration: harvard

## Validation

**Passed**: the calibrated cost below the default's, and no target's RMSE more than 10% worse.

| target | default | calibrated | change |
|---|---|---|---|
| lw_up | 1.79 | 1.94 | +8% |
| le | 2.90 | 2.81 | -3% |
| h | 4.61 | 4.94 | +7% |
| gpp | 2.64 | 2.29 | -13% |
| ustar | 1.26 | 0.60 | -52% |

RMSE in units of each target's sigma. Cost: default 7011.51, calibrated 6945.98.

The full record with the slow tier on is the next check: run the calibrated configs (budgets closed; dry-season GPP and LE no worse than the default's).

## Keys

| key | kind | scope | MAP | 68 % | prior centre | prior z | sd ratio | note | prior source |
|---|---|---|---|---|---|---|---|---|---|
| leaf_angle_mean | trait | plant_type | 33.83 | 26.2-42.86 | 33.9 | -0.01 | 0.85 |  | the PFT file's value, +-10 degrees; spherical leaves average 57.3. The range is  |
| vcmax25 | trait | plant_type | 56.21 | 53.31-59.19 | 73.98 | -0.50 | 0.09 |  | EEO: the coordination of the Rubisco- and light-limited rates at the site's grow |
| stomatal_g1 | trait | plant_type | 5.339 | 4.968-5.722 | 1.994 | +1.98 | 0.16 |  | EEO: the least-cost optimum at the site's growing-season daytime climate (priors |
| stomatal_g0 | trait | plant_type | 0.002484 | 0.0009821-0.006513 | 0.01 | -1.28 | 0.87 |  | the PFT file's value, a factor e either way: the residual conductance is poorly  |
| z0m_ratio | effective | site | 0.06422 | 0.05954-0.06934 | 0.13 | -2.07 | 0.22 |  | the base value (ED2 0.13), +-0.04: closed canopies give z0/h of about 0.06-0.15  |
| kappa | observation | observation | 0.9167 | 0.8798-0.9545 | 1 | -1.27 | 0.58 |  | the gap between the provider's night-time and daytime partitionings (RECO agains |

Fixed at their priors: wstress_sref_stomata (the kept data sample its process (drought) in 0 records, fewer than 48)

## Targets

| target | chi2/n | sigma scale | model/tower | n (cal) |
|---|---|---|---|---|
| lw_up | 1.05 | 2.08 | 1.002 | 1378 |
| le | 0.99 | 1.91 | 1.114 | 670 |
| h | 1.02 | 2.88 | 0.913 | 673 |
| gpp | 0.99 | 1.68 | 0.981 | 488 |
| ustar | 0.60 | 1.00 | 1.015 | 781 |

lw_up, model/tower by local hour: 0 h 0.99, 1 h 0.99, 2 h 0.99, 3 h 0.99, 4 h 0.99, 5 h 1.00, 6 h 1.02, 7 h 1.03, 8 h 1.03, 9 h 1.03, 10 h 1.03, 11 h 1.02, 12 h 1.02, 13 h 1.02, 14 h 1.01, 15 h 1.00, 16 h 0.99, 17 h 0.99, 18 h 0.98, 19 h 0.98, 20 h 0.99, 21 h 0.99, 22 h 0.99, 23 h 0.99
lw_up, by light: SW quartile 2 0.99, SW quartile 3 1.01, SW quartile 4 1.03

le, model/tower by local hour: 4 h 0.82, 5 h 0.79, 6 h 0.87, 7 h 1.01, 8 h 1.08, 9 h 1.06, 10 h 0.95, 11 h 1.15, 12 h 1.11, 13 h 1.17, 14 h 1.15, 15 h 1.19, 16 h 1.21, 17 h 1.38, 18 h 1.67
le, by light: SW quartile 1 1.30, SW quartile 2 1.10, SW quartile 3 1.08, SW quartile 4 1.12

h, model/tower by local hour: 4 h 0.33, 5 h 0.38, 6 h 0.25, 7 h 0.61, 8 h 0.64, 9 h 0.74, 10 h 0.89, 11 h 0.91, 12 h 0.90, 13 h 0.95, 14 h 1.04, 15 h 1.35, 16 h 2.22, 17 h -4.19, 18 h -0.90
h, by light: SW quartile 1 -0.66, SW quartile 2 0.99, SW quartile 3 0.90, SW quartile 4 0.81

gpp, model/tower by local hour: 5 h 0.95, 6 h 0.83, 7 h 0.92, 8 h 0.90, 9 h 0.91, 10 h 0.95, 11 h 0.96, 12 h 1.00, 13 h 1.03, 14 h 1.05, 15 h 1.07, 16 h 1.07, 17 h 1.15, 18 h 0.93
gpp, by light: SW quartile 1 0.98, SW quartile 2 0.93, SW quartile 3 0.98, SW quartile 4 1.02

ustar, model/tower by local hour: 4 h 0.82, 5 h 0.72, 6 h 0.75, 7 h 0.86, 8 h 0.95, 9 h 0.98, 10 h 1.08, 11 h 1.11, 12 h 1.05, 13 h 1.05, 14 h 1.08, 15 h 1.07, 16 h 1.11, 17 h 1.15, 18 h 1.22
ustar, by light: SW quartile 1 1.00, SW quartile 2 0.96, SW quartile 3 1.00, SW quartile 4 1.08

## kappa

kappa = 0.917: the tower's respiration 3.14 -> 3.43 umol m-2 s-1; GPP 6.05 -> 6.26 (2.29 -> 2.37 kgC m-2 yr-1), means over the records with GPP and its respiration measured, day and night; GPP gains the scaled respiration by day only.

## Uncertainty

Laplace, from the final gradient matrix.
- alternative gpp_ustar: the fit did not use the provider's threshold, or GPP's diagnostic finds no plateau
- alternative closure (Bowen (H and LE scaled together) against the fit's attribution shares): largest linear shift 10.33 sd -- refitted: leaf_angle_mean 33.83->28.37, vcmax25 56.21->58.26, stomatal_g1 5.339->10.54, stomatal_g0 0.002484->0.002783, z0m_ratio 0.06422->0.06711, kappa 0.9167->0.8867
- alternative partitioning (the provider's daytime partitioning (GPP_DT, RECO_DT)): largest linear shift 2.93 sd -- refitted: leaf_angle_mean 33.83->38.55, vcmax25 56.21->48.18, stomatal_g1 5.339->5.968, stomatal_g0 0.002484->0.002336, z0m_ratio 0.06422->0.06515, kappa 0.9167->0.9002
