# Fast calibration: harvard

## Validation

**Passed**: the calibrated cost below the default's, and no target's RMSE more than 10% worse.

| target | default | calibrated | change |
|---|---|---|---|
| lw_up | 1.66 | 1.69 | +1% |
| le | 2.41 | 2.41 | +0% |
| h | 3.51 | 3.81 | +8% |
| gpp | 3.09 | 2.08 | -33% |
| ustar | 1.07 | 0.58 | -46% |

RMSE in units of each target's sigma. Cost: default 17731, calibrated 16555.3.

The full record with the slow tier on is the next check: run the calibrated configs (budgets closed; dry-season GPP and LE no worse than the default's).

## Keys

| key | kind | scope | MAP | 68 % | prior centre | prior z | sd ratio | note | prior source |
|---|---|---|---|---|---|---|---|---|---|
| leaf_angle_mean | trait | plant_type | 31.75 | 24.13-41.19 | 33.9 | -0.22 | 0.91 |  | the PFT file's value, +-10 degrees; spherical leaves average 57.3. The range is  |
| vcmax25 | trait | plant_type | 67.41 | 63.99-70.87 | 73.98 | -0.18 | 0.09 |  | EEO: the coordination of the Rubisco- and light-limited rates at the site's grow |
| vcmax25_pft2 | trait | plant_type | 19.64 | 17.55-22.29 | 58 | -3.03 | 0.42 | beyond 2 prior sd, pushed by gpp: diagnose it (plan §5.1) or relabel it effective | Kattge et al. 2009 (Glob. Change Biol. 15:976) for evergreen coniferous forest,  |
| capacity_low_temp_pft2 | effective | plant_type | 274.7 | 273.6-277.2 | 277.9 | -1.14 | 0.97 |  | ED2's temperate and conifer Vm_low_temp, 4.71 C (Moorcroft et al. 2001); ED2 giv |
| stomatal_g1 | trait | plant_type | 5.702 | 5.26-6.158 | 1.994 | +2.13 | 0.19 | beyond 2 prior sd, pushed by lw_up: diagnose it (plan §5.1) or relabel it effective | EEO: the least-cost optimum at the site's growing-season daytime climate (priors |
| stomatal_g1_pft2 | trait | plant_type | 2.066 | 1.445-3.034 | 1.994 | +0.07 | 0.75 |  | EEO: the least-cost optimum at the site's growing-season daytime climate (priors |
| stomatal_g0 | trait | plant_type | 0.002916 | 0.001132-0.007605 | 0.01 | -1.14 | 0.87 |  | the PFT file's value, a factor e either way: the residual conductance is poorly  |
| stomatal_g0_pft2 | trait | plant_type | 0.008103 | 0.003611-0.01649 | 0.01 | -0.21 | 0.74 |  | the PFT file's value, a factor e either way: the residual conductance is poorly  |
| z0m_ratio | effective | site | 0.06378 | 0.0598-0.06808 | 0.13 | -2.09 | 0.19 |  | the base value (ED2 0.13), +-0.04: closed canopies give z0/h of about 0.06-0.15  |
| kappa | observation | observation | 0.9243 | 0.8848-0.9647 | 1 | -1.15 | 0.62 |  | the gap between the provider's night-time and daytime partitionings (RECO agains |

Fixed at their priors: wstress_sref_stomata (the kept data sample its process (drought) in 0 records, fewer than 48); wstress_sref_stomata_pft2 (the kept data sample its process (drought) in 0 records, fewer than 48); leaf_angle_mean_pft2 (uninformed: posterior/prior sd ratio 0.92 >= 0.9)

## Targets

| target | chi2/n | sigma scale | model/tower | n (cal) |
|---|---|---|---|---|
| lw_up | 1.01 | 1.86 | 1.000 | 1852 |
| le | 0.99 | 2.06 | 1.108 | 854 |
| h | 1.03 | 2.96 | 0.862 | 857 |
| gpp | 0.95 | 1.89 | 0.983 | 618 |
| ustar | 0.53 | 1.00 | 0.994 | 1032 |

lw_up, model/tower by local hour: 0 h 0.99, 1 h 0.99, 2 h 0.99, 3 h 0.99, 4 h 0.99, 5 h 1.00, 6 h 1.01, 7 h 1.02, 8 h 1.02, 9 h 1.02, 10 h 1.02, 11 h 1.01, 12 h 1.01, 13 h 1.01, 14 h 1.00, 15 h 1.00, 16 h 0.99, 17 h 0.99, 18 h 0.99, 19 h 0.99, 20 h 0.99, 21 h 0.99, 22 h 0.99, 23 h 0.99
lw_up, by light: SW quartile 2 0.99, SW quartile 3 1.00, SW quartile 4 1.01

le, model/tower by local hour: 4 h 0.81, 5 h 0.78, 6 h 0.85, 7 h 0.98, 8 h 1.04, 9 h 1.05, 10 h 0.95, 11 h 1.14, 12 h 1.11, 13 h 1.17, 14 h 1.14, 15 h 1.18, 16 h 1.25, 17 h 1.39, 18 h 1.59
le, by light: SW quartile 1 1.29, SW quartile 2 1.09, SW quartile 3 1.07, SW quartile 4 1.12

h, model/tower by local hour: 4 h 0.32, 5 h 0.34, 6 h 0.29, 7 h 0.54, 8 h 0.59, 9 h 0.70, 10 h 0.81, 11 h 0.85, 12 h 0.87, 13 h 0.91, 14 h 1.00, 15 h 1.15, 16 h 1.63, 17 h -9.03, 18 h -0.97
h, by light: SW quartile 1 -3.13, SW quartile 2 0.95, SW quartile 3 0.83, SW quartile 4 0.77

gpp, model/tower by local hour: 5 h 0.93, 6 h 0.83, 7 h 0.92, 8 h 0.89, 9 h 0.92, 10 h 0.95, 11 h 0.97, 12 h 1.00, 13 h 1.05, 14 h 1.06, 15 h 1.07, 16 h 1.06, 17 h 1.13, 18 h 0.89
gpp, by light: SW quartile 1 0.95, SW quartile 2 0.94, SW quartile 3 0.99, SW quartile 4 1.03

ustar, model/tower by local hour: 4 h 0.82, 5 h 0.78, 6 h 0.78, 7 h 0.84, 8 h 0.92, 9 h 0.96, 10 h 1.06, 11 h 1.09, 12 h 1.05, 13 h 1.04, 14 h 1.06, 15 h 1.02, 16 h 1.05, 17 h 1.05, 18 h 1.18
ustar, by light: SW quartile 1 0.99, SW quartile 2 0.93, SW quartile 3 0.99, SW quartile 4 1.06

## kappa

kappa = 0.924: the tower's respiration 3.14 -> 3.40 umol m-2 s-1; GPP 6.05 -> 6.24 (2.29 -> 2.36 kgC m-2 yr-1), means over the records with GPP and its respiration measured, day and night; GPP gains the scaled respiration by day only.

## Uncertainty

Laplace, from the final gradient matrix -- LOCAL ONLY (the linearity check failed).
- alternative gpp_ustar: the fit did not use the provider's threshold, or GPP's diagnostic finds no plateau
- alternative closure (Bowen (H and LE scaled together) against the fit's attribution shares): largest linear shift 7.67 sd -- refitted: leaf_angle_mean 31.75->27.19, vcmax25 67.41->68.29, vcmax25_pft2 19.64->20.24, capacity_low_temp_pft2 274.7->274.4, stomatal_g1 5.702->10.08, stomatal_g1_pft2 2.066->4.53, stomatal_g0 0.002916->0.003342, stomatal_g0_pft2 0.008103->0.009513, z0m_ratio 0.06378->0.06699, kappa 0.9243->0.8862
- alternative partitioning (the provider's daytime partitioning (GPP_DT, RECO_DT)): largest linear shift 2.31 sd -- refitted: leaf_angle_mean 31.75->34.47, vcmax25 67.41->59.82, vcmax25_pft2 19.64->16.73, capacity_low_temp_pft2 274.7->276.2, stomatal_g1 5.702->6.266, stomatal_g1_pft2 2.066->2.388, stomatal_g0 0.002916->0.003002, stomatal_g0_pft2 0.008103->0.006642, z0m_ratio 0.06378->0.06561, kappa 0.9243->0.9216
