# Fast calibration: harvard

## Validation

**Passed**: the calibrated cost below the default's, and no target's RMSE more than 10% worse.

| target | default | calibrated | change |
|---|---|---|---|
| lw_up | 2.12 | 2.18 | +3% |
| le | 2.41 | 2.32 | -4% |
| h | 3.72 | 3.98 | +7% |
| gpp | 2.35 | 2.05 | -13% |
| ustar | 1.24 | 0.66 | -47% |

RMSE in units of each target's sigma. Cost: default 15825.8, calibrated 15408.6.

The full record with the slow tier on is the next check: run the calibrated configs (budgets closed; dry-season GPP and LE no worse than the default's).

## Keys

| key | kind | scope | MAP | 68 % | prior centre | prior z | sd ratio | note | prior source |
|---|---|---|---|---|---|---|---|---|---|
| leaf_angle_mean | trait | plant_type | 31.19 | 23.89-40.26 | 33.9 | -0.28 | 0.88 |  | the PFT file's value, +-10 degrees; spherical leaves average 57.3. The range is  |
| vcmax25 | trait | plant_type | 59.82 | 56.56-63.17 | 73.98 | -0.39 | 0.10 |  | EEO: the coordination of the Rubisco- and light-limited rates at the site's grow |
| stomatal_g1 | trait | plant_type | 6.406 | 5.947-6.873 | 1.994 | +2.42 | 0.19 | beyond 2 prior sd, pushed by lw_up: diagnose it (plan §5.1) or relabel it effective | EEO: the least-cost optimum at the site's growing-season daytime climate (priors |
| stomatal_g1_pft2 | trait | plant_type | 0.9482 | 0.7421-1.342 | 1.994 | -1.78 | 0.93 |  | EEO: the least-cost optimum at the site's growing-season daytime climate (priors |
| stomatal_g0 | trait | plant_type | 0.003513 | 0.00148-0.008231 | 0.01 | -0.97 | 0.78 |  | the PFT file's value, a factor e either way: the residual conductance is poorly  |
| stomatal_g0_pft2 | trait | plant_type | 0.002337 | 0.0008742-0.006583 | 0.01 | -1.33 | 0.93 |  | the PFT file's value, a factor e either way: the residual conductance is poorly  |
| z0m_ratio | effective | site | 0.06502 | 0.06033-0.07015 | 0.13 | -2.03 | 0.22 |  | the base value (ED2 0.13), +-0.04: closed canopies give z0/h of about 0.06-0.15  |
| kappa | observation | observation | 0.9162 | 0.8798-0.9535 | 1 | -1.27 | 0.57 |  | the gap between the provider's night-time and daytime partitionings (RECO agains |

Fixed at their priors: wstress_sref_stomata (the kept data sample its process (drought) in 0 records, fewer than 48); wstress_sref_stomata_pft2 (the kept data sample its process (drought) in 0 records, fewer than 48); leaf_angle_mean_pft2 (uninformed: posterior/prior sd ratio 0.93 >= 0.9)

## Targets

| target | chi2/n | sigma scale | model/tower | n (cal) |
|---|---|---|---|---|
| lw_up | 1.09 | 2.00 | 1.001 | 1378 |
| le | 0.99 | 1.88 | 1.112 | 670 |
| h | 1.01 | 2.90 | 0.917 | 673 |
| gpp | 1.01 | 1.65 | 0.977 | 488 |
| ustar | 0.61 | 1.00 | 1.019 | 781 |

lw_up, model/tower by local hour: 0 h 0.99, 1 h 0.99, 2 h 0.99, 3 h 0.99, 4 h 0.99, 5 h 1.00, 6 h 1.02, 7 h 1.03, 8 h 1.03, 9 h 1.03, 10 h 1.03, 11 h 1.02, 12 h 1.02, 13 h 1.01, 14 h 1.01, 15 h 1.00, 16 h 0.99, 17 h 0.99, 18 h 0.98, 19 h 0.98, 20 h 0.99, 21 h 0.99, 22 h 0.99, 23 h 0.99
lw_up, by light: SW quartile 2 0.99, SW quartile 3 1.01, SW quartile 4 1.02

le, model/tower by local hour: 4 h 0.76, 5 h 0.78, 6 h 0.87, 7 h 1.02, 8 h 1.09, 9 h 1.07, 10 h 0.96, 11 h 1.15, 12 h 1.11, 13 h 1.17, 14 h 1.15, 15 h 1.18, 16 h 1.19, 17 h 1.33, 18 h 1.58
le, by light: SW quartile 1 1.25, SW quartile 2 1.09, SW quartile 3 1.08, SW quartile 4 1.12

h, model/tower by local hour: 4 h 0.32, 5 h 0.37, 6 h 0.27, 7 h 0.62, 8 h 0.64, 9 h 0.73, 10 h 0.89, 11 h 0.91, 12 h 0.90, 13 h 0.95, 14 h 1.05, 15 h 1.36, 16 h 2.28, 17 h -4.54, 18 h -1.04
h, by light: SW quartile 1 -0.74, SW quartile 2 1.00, SW quartile 3 0.90, SW quartile 4 0.81

gpp, model/tower by local hour: 5 h 0.97, 6 h 0.85, 7 h 0.93, 8 h 0.91, 9 h 0.92, 10 h 0.95, 11 h 0.96, 12 h 0.99, 13 h 1.02, 14 h 1.04, 15 h 1.05, 16 h 1.06, 17 h 1.14, 18 h 0.91
gpp, by light: SW quartile 1 0.98, SW quartile 2 0.94, SW quartile 3 0.98, SW quartile 4 1.01

ustar, model/tower by local hour: 4 h 0.82, 5 h 0.73, 6 h 0.75, 7 h 0.86, 8 h 0.95, 9 h 0.98, 10 h 1.08, 11 h 1.11, 12 h 1.05, 13 h 1.06, 14 h 1.09, 15 h 1.07, 16 h 1.12, 17 h 1.15, 18 h 1.23
ustar, by light: SW quartile 1 1.01, SW quartile 2 0.96, SW quartile 3 1.01, SW quartile 4 1.08

## kappa

kappa = 0.916: the tower's respiration 3.14 -> 3.43 umol m-2 s-1; GPP 6.05 -> 6.26 (2.29 -> 2.37 kgC m-2 yr-1), means over the records with GPP and its respiration measured, day and night; GPP gains the scaled respiration by day only.

## Uncertainty

Laplace, from the final gradient matrix -- LOCAL ONLY (the linearity check failed).
- alternative gpp_ustar: the fit did not use the provider's threshold, or GPP's diagnostic finds no plateau
- alternative closure (Bowen (H and LE scaled together) against the fit's attribution shares): largest linear shift 9.41 sd -- refitted: leaf_angle_mean 31.19->27.98, vcmax25 59.82->60.94, stomatal_g1 6.406->11.4, stomatal_g1_pft2 0.9482->1.288, stomatal_g0 0.003513->0.001711, stomatal_g0_pft2 0.002337->0.002794, z0m_ratio 0.06502->0.0664, kappa 0.9162->0.898
- alternative partitioning (the provider's daytime partitioning (GPP_DT, RECO_DT)): largest linear shift 3.33 sd -- refitted: leaf_angle_mean 31.19->38.34, vcmax25 59.82->49.11, stomatal_g1 6.406->6.94, stomatal_g1_pft2 0.9482->1.341, stomatal_g0 0.003513->0.003426, stomatal_g0_pft2 0.002337->0.002781, z0m_ratio 0.06502->0.06574, kappa 0.9162->0.9048
