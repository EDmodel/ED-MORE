# Fast calibration: harvard

## Validation

**Passed**: the calibrated cost below the default's, and no target's RMSE more than 10% worse.

| target | default | calibrated | change |
|---|---|---|---|
| lw_up | 1.77 | 1.87 | +6% |
| le | 2.42 | 2.30 | -5% |
| h | 3.89 | 4.15 | +7% |
| gpp | 2.48 | 2.13 | -14% |
| ustar | 1.14 | 0.58 | -49% |

RMSE in units of each target's sigma. Cost: default 14001.6, calibrated 13697.3.

The full record with the slow tier on is the next check: run the calibrated configs (budgets closed; dry-season GPP and LE no worse than the default's).

## Keys

| key | kind | scope | MAP | 68 % | prior centre | prior z | sd ratio | note | prior source |
|---|---|---|---|---|---|---|---|---|---|
| leaf_angle_mean | trait | plant_type | 32.45 | 24.93-41.59 | 33.9 | -0.15 | 0.87 |  | the PFT file's value, +-10 degrees; spherical leaves average 57.3. The range is  |
| vcmax25 | trait | plant_type | 58.28 | 55.07-61.59 | 73.98 | -0.44 | 0.10 |  | EEO: the coordination of the Rubisco- and light-limited rates at the site's grow |
| stomatal_g1 | trait | plant_type | 6.202 | 5.743-6.671 | 1.994 | +2.34 | 0.19 | beyond 2 prior sd, pushed by lw_up: diagnose it (plan §5.1) or relabel it effective | EEO: the least-cost optimum at the site's growing-season daytime climate (priors |
| stomatal_g1_pft2 | trait | plant_type | 1.105 | 0.8275-1.63 | 1.994 | -1.33 | 0.92 |  | EEO: the least-cost optimum at the site's growing-season daytime climate (priors |
| stomatal_g0 | trait | plant_type | 0.003865 | 0.001664-0.008771 | 0.01 | -0.89 | 0.76 |  | the PFT file's value, a factor e either way: the residual conductance is poorly  |
| stomatal_g0_pft2 | trait | plant_type | 0.001964 | 0.000725-0.00581 | 0.01 | -1.49 | 0.97 |  | the PFT file's value, a factor e either way: the residual conductance is poorly  |
| z0m_ratio | effective | site | 0.06296 | 0.05832-0.06805 | 0.13 | -2.12 | 0.22 |  | the base value (ED2 0.13), +-0.04: closed canopies give z0/h of about 0.06-0.15  |
| kappa | observation | observation | 0.9239 | 0.8866-0.9622 | 1 | -1.15 | 0.58 |  | the gap between the provider's night-time and daytime partitionings (RECO agains |

Fixed at their priors: wstress_sref_stomata (the kept data sample its process (drought) in 0 records, fewer than 48); wstress_sref_stomata_pft2 (the kept data sample its process (drought) in 0 records, fewer than 48); leaf_angle_mean_pft2 (uninformed: posterior/prior sd ratio 0.94 >= 0.9)

## Targets

| target | chi2/n | sigma scale | model/tower | n (cal) |
|---|---|---|---|---|
| lw_up | 1.02 | 2.03 | 1.001 | 1378 |
| le | 0.99 | 1.89 | 1.117 | 670 |
| h | 1.04 | 2.88 | 0.912 | 673 |
| gpp | 0.97 | 1.67 | 0.981 | 488 |
| ustar | 0.59 | 1.00 | 1.008 | 781 |

lw_up, model/tower by local hour: 0 h 0.99, 1 h 0.99, 2 h 0.99, 3 h 0.99, 4 h 0.99, 5 h 1.00, 6 h 1.02, 7 h 1.03, 8 h 1.03, 9 h 1.03, 10 h 1.03, 11 h 1.02, 12 h 1.02, 13 h 1.01, 14 h 1.01, 15 h 1.00, 16 h 0.99, 17 h 0.99, 18 h 0.98, 19 h 0.98, 20 h 0.99, 21 h 0.99, 22 h 0.99, 23 h 0.99
lw_up, by light: SW quartile 2 0.99, SW quartile 3 1.01, SW quartile 4 1.02

le, model/tower by local hour: 4 h 0.78, 5 h 0.79, 6 h 0.87, 7 h 1.02, 8 h 1.09, 9 h 1.07, 10 h 0.96, 11 h 1.16, 12 h 1.12, 13 h 1.18, 14 h 1.15, 15 h 1.18, 16 h 1.19, 17 h 1.34, 18 h 1.60
le, by light: SW quartile 1 1.27, SW quartile 2 1.10, SW quartile 3 1.08, SW quartile 4 1.13

h, model/tower by local hour: 4 h 0.31, 5 h 0.36, 6 h 0.27, 7 h 0.61, 8 h 0.64, 9 h 0.73, 10 h 0.89, 11 h 0.90, 12 h 0.89, 13 h 0.94, 14 h 1.04, 15 h 1.36, 16 h 2.27, 17 h -4.54, 18 h -1.04
h, by light: SW quartile 1 -0.74, SW quartile 2 1.00, SW quartile 3 0.90, SW quartile 4 0.80

gpp, model/tower by local hour: 5 h 0.97, 6 h 0.84, 7 h 0.93, 8 h 0.91, 9 h 0.92, 10 h 0.95, 11 h 0.96, 12 h 1.00, 13 h 1.03, 14 h 1.05, 15 h 1.06, 16 h 1.07, 17 h 1.15, 18 h 0.92
gpp, by light: SW quartile 1 0.98, SW quartile 2 0.94, SW quartile 3 0.98, SW quartile 4 1.01

ustar, model/tower by local hour: 4 h 0.82, 5 h 0.72, 6 h 0.75, 7 h 0.85, 8 h 0.94, 9 h 0.97, 10 h 1.07, 11 h 1.10, 12 h 1.04, 13 h 1.04, 14 h 1.07, 15 h 1.06, 16 h 1.10, 17 h 1.14, 18 h 1.22
ustar, by light: SW quartile 1 1.00, SW quartile 2 0.95, SW quartile 3 1.00, SW quartile 4 1.07

## kappa

kappa = 0.924: the tower's respiration 3.14 -> 3.40 umol m-2 s-1; GPP 6.05 -> 6.24 (2.29 -> 2.36 kgC m-2 yr-1), means over the records with GPP and its respiration measured, day and night; GPP gains the scaled respiration by day only.

## Uncertainty

Laplace, from the final gradient matrix -- LOCAL ONLY (the linearity check failed).
- alternative gpp_ustar: the fit did not use the provider's threshold, or GPP's diagnostic finds no plateau
- alternative closure (Bowen (H and LE scaled together) against the fit's attribution shares): largest linear shift 9.11 sd -- refitted: leaf_angle_mean 32.45->26.81, vcmax25 58.28->59.53, stomatal_g1 6.202->11.37, stomatal_g1_pft2 1.105->1.104, stomatal_g0 0.003865->0.001257, stomatal_g0_pft2 0.001964->0.003409, z0m_ratio 0.06296->0.0662, kappa 0.9239->0.8947
- alternative partitioning (the provider's daytime partitioning (GPP_DT, RECO_DT)): largest linear shift 3.17 sd -- refitted: leaf_angle_mean 32.45->35.45, vcmax25 58.28->49.64, stomatal_g1 6.202->6.706, stomatal_g1_pft2 1.105->1.322, stomatal_g0 0.003865->0.00243, stomatal_g0_pft2 0.001964->0.005034, z0m_ratio 0.06296->0.06497, kappa 0.9239->0.8995
