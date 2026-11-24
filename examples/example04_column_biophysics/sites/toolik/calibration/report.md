# Fast calibration: toolik

## Validation

**Passed**: the calibrated cost below the default's, and no target's RMSE more than 10% worse.

| target | default | calibrated | change |
|---|---|---|---|
| lw_up | 1.84 | 1.59 | -14% |
| le | 4.40 | 3.96 | -10% |
| h | 2.94 | 2.64 | -10% |
| gpp | 6.30 | 3.55 | -44% |
| ustar | 0.52 | 0.51 | -2% |

RMSE in units of each target's sigma. Cost: default 23044.5, calibrated 14365.9.

The full record with the slow tier on is the next check: run the calibrated configs (budgets closed; dry-season GPP and LE no worse than the default's).

## Keys

| key | kind | scope | MAP | 68 % | prior centre | prior z | sd ratio | note | prior source |
|---|---|---|---|---|---|---|---|---|---|
| vcmax25 | trait | plant_type | 13.8 | 11.24-21.8 | 88.27 | -2.91 | 0.93 | beyond 2 prior sd, pushed by gpp: diagnose it (plan §5.1) or relabel it effective | EEO: the coordination of the Rubisco- and light-limited rates at the site's grow |
| vcmax25_pft2 | trait | plant_type | 61.7 | 57.76-65.75 | 88.27 | -0.59 | 0.09 |  | EEO: the coordination of the Rubisco- and light-limited rates at the site's grow |
| vcmax25_pft3 | trait | plant_type | 16.3 | 12.33-26.99 | 88.27 | -2.49 | 0.82 | beyond 2 prior sd, pushed by gpp: diagnose it (plan §5.1) or relabel it effective | EEO: the coordination of the Rubisco- and light-limited rates at the site's grow |
| capacity_low_temp | effective | plant_type | 282.4 | 278.8-285.3 | 277.9 | +1.01 | 0.78 |  | ED2's temperate and conifer Vm_low_temp, 4.71 C (Moorcroft et al. 2001); ED2 giv |
| capacity_low_temp_pft2 | effective | plant_type | 274.3 | 273.5-276.2 | 277.9 | -1.40 | 0.92 |  | ED2's temperate and conifer Vm_low_temp, 4.71 C (Moorcroft et al. 2001); ED2 giv |
| capacity_low_temp_pft3 | effective | plant_type | 286.1 | 283.3-287.4 | 277.9 | +2.14 | 0.91 |  | ED2's temperate and conifer Vm_low_temp, 4.71 C (Moorcroft et al. 2001); ED2 giv |
| stomatal_g1_pft2 | trait | plant_type | 12.51 | 10.77-13.61 | 1.08 | +4.89 | 0.77 | beyond 2 prior sd, pushed by le: diagnose it (plan §5.1) or relabel it effective | EEO: the least-cost optimum at the site's growing-season daytime climate (priors |
| stomatal_g1_pft3 | trait | plant_type | 0.9737 | 0.6918-1.702 | 1.08 | -0.22 | 1.00 |  | EEO: the least-cost optimum at the site's growing-season daytime climate (priors |
| stomatal_g0_pft2 | trait | plant_type | 0.02959 | 0.01366-0.04955 | 0.01 | +1.30 | 0.98 |  | the PFT file's value, a factor e either way: the residual conductance is poorly  |
| kappa | observation | observation | 0.8749 | 0.8394-0.9117 | 1 | -2.55 | 0.77 |  | the gap between the provider's night-time and daytime partitionings (RECO agains |

Fixed at their priors: leaf_angle_mean (uninformed: posterior/prior sd ratio 1.00 >= 0.9); leaf_angle_mean_pft2 (uninformed: posterior/prior sd ratio 0.95 >= 0.9); leaf_angle_mean_pft3 (uninformed: posterior/prior sd ratio 1.00 >= 0.9); stomatal_g1 (uninformed: posterior/prior sd ratio 0.92 >= 0.9); stomatal_g0 (uninformed: posterior/prior sd ratio 0.97 >= 0.9); stomatal_g0_pft3 (uninformed: posterior/prior sd ratio 0.95 >= 0.9)

## Targets

| target | chi2/n | sigma scale | model/tower | n (cal) |
|---|---|---|---|---|
| lw_up | 1.01 | 1.79 | 1.005 | 2483 |
| le | 1.02 | 2.81 | 0.949 | 290 |
| h | 1.00 | 2.55 | 0.893 | 1352 |
| gpp | 1.08 | 3.00 | 0.981 | 1258 |
| ustar | 0.23 | 1.00 | 1.020 | 2069 |

lw_up, model/tower by local hour: 0 h 1.02, 1 h 1.01, 2 h 1.01, 3 h 1.00, 4 h 0.99, 5 h 0.99, 6 h 0.99, 7 h 0.98, 8 h 0.99, 9 h 0.99, 10 h 1.00, 11 h 1.00, 12 h 1.00, 13 h 1.00, 14 h 1.01, 15 h 1.01, 16 h 1.01, 17 h 1.01, 18 h 1.01, 19 h 1.01, 20 h 1.01, 21 h 1.02, 22 h 1.02, 23 h 1.02
lw_up, by light: SW quartile 1 1.01, SW quartile 2 1.01, SW quartile 3 1.00, SW quartile 4 1.00

le, model/tower by local hour: 6 h 0.57, 7 h 0.64, 8 h 0.73, 9 h 0.64, 10 h 0.68, 11 h 1.05, 12 h 0.88, 13 h 0.86, 14 h 0.97, 15 h 1.32, 16 h 1.07, 17 h 1.28, 18 h 1.16
le, by light: SW quartile 1 0.79, SW quartile 2 0.99, SW quartile 3 1.04, SW quartile 4 0.92

h, model/tower by local hour: 0 h 0.64, 1 h 0.82, 2 h 0.95, 3 h 1.06, 4 h 1.28, 5 h 1.34, 6 h -2.04, 7 h 0.01, 8 h 0.41, 9 h 0.51, 10 h 0.65, 11 h 0.76, 12 h 0.80, 13 h 0.85, 14 h 0.96, 15 h 1.09, 16 h 1.17, 17 h 1.27, 18 h 1.63, 19 h 3.35, 20 h 12.62, 21 h -2.85, 22 h -0.15, 23 h -0.13
h, by light: SW quartile 1 0.14, SW quartile 2 1.84, SW quartile 3 0.92, SW quartile 4 0.71

gpp, model/tower by local hour: 1 h 1.70, 2 h 1.19, 3 h 1.00, 4 h 0.87, 5 h 0.78, 6 h 0.78, 7 h 0.75, 8 h 0.80, 9 h 0.83, 10 h 0.86, 11 h 0.92, 12 h 0.97, 13 h 1.04, 14 h 1.02, 15 h 1.09, 16 h 1.09, 17 h 1.18, 18 h 1.21, 19 h 1.36, 20 h 1.21, 21 h 1.26, 22 h 1.23, 23 h 1.51
gpp, by light: SW quartile 1 1.28, SW quartile 2 0.96, SW quartile 3 0.89, SW quartile 4 1.01

ustar, model/tower by local hour: 0 h 1.14, 1 h 1.25, 2 h 1.23, 3 h 1.22, 4 h 1.07, 5 h 1.02, 6 h 0.89, 7 h 0.92, 8 h 0.91, 9 h 0.93, 10 h 0.94, 11 h 0.94, 12 h 0.97, 13 h 0.98, 14 h 1.00, 15 h 1.04, 16 h 1.04, 17 h 1.08, 18 h 1.05, 19 h 1.10, 20 h 1.16, 21 h 1.25, 22 h 1.21, 23 h 0.97
ustar, by light: SW quartile 1 1.11, SW quartile 2 1.09, SW quartile 3 0.97, SW quartile 4 0.96

## kappa

kappa = 0.875: the tower's respiration 0.77 -> 0.88 umol m-2 s-1; GPP 0.98 -> 1.06 (0.37 -> 0.40 kgC m-2 yr-1), means over the records with GPP and its respiration measured, day and night; GPP gains the scaled respiration by day only.

## Uncertainty

Laplace, from the final gradient matrix -- LOCAL ONLY (the linearity check failed).
- alternative gpp_ustar (GPP u* >= 0.35 (the fit used {2019: 0.17, 2020: 0.142, 2021: 0.103, 2022: 0.1, 2023: 0.105, 2024: 0.138})): largest linear shift 1.10 sd -- refitted: vcmax25 13.8->14.99, vcmax25_pft2 61.7->65.05, vcmax25_pft3 16.3->16.47, capacity_low_temp 282.4->282.9, capacity_low_temp_pft2 274.3->275.8, capacity_low_temp_pft3 286.1->284, stomatal_g1_pft2 12.51->12.62, stomatal_g1_pft3 0.9737->1.013, stomatal_g0_pft2 0.02959->0.02675, kappa 0.8749->0.915
- alternative closure: the fit already used Bowen
- alternative partitioning (the provider's daytime partitioning (GPP_DT, RECO_DT)): largest linear shift 1.39 sd -- refitted: vcmax25 13.8->17.81, vcmax25_pft2 61.7->57.97, vcmax25_pft3 16.3->16.44, capacity_low_temp 282.4->281.1, capacity_low_temp_pft2 274.3->274.5, capacity_low_temp_pft3 286.1->283.8, stomatal_g1_pft2 12.51->11.94, stomatal_g1_pft3 0.9737->1.018, stomatal_g0_pft2 0.02959->0.0297, kappa 0.8749->0.9137
