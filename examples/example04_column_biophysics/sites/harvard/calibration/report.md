# Fast calibration: harvard

## Validation

**Failed**: the calibrated cost below the default's, and no target's RMSE more than 10% worse.

| target | default | calibrated | change |
|---|---|---|---|
| lw_up | 1.49 | 1.66 | +12% |
| le | 2.24 | 2.03 | -9% |
| h | 3.32 | 3.45 | +4% |
| gpp | 2.80 | 2.13 | -24% |
| ustar | 1.12 | 0.63 | -44% |

RMSE in units of each target's sigma. Cost: default 15889.1, calibrated 15138.1.

The full record with the slow tier on is the next check: run the calibrated configs (budgets closed; dry-season GPP and LE no worse than the default's).

## Keys

| key | kind | scope | MAP | 68 % | prior centre | prior z | sd ratio | note | prior source |
|---|---|---|---|---|---|---|---|---|---|
| leaf_angle_mean | trait | plant_type | 24.37 | 18.61-32.57 | 33.9 | -1.10 | 0.96 |  | the PFT file's value, +-10 degrees; spherical leaves average 57.3. The range is  |
| vcmax25 | trait | plant_type | 92.97 | 88.33-97.46 | 73.98 | +0.52 | 0.13 |  | EEO: the coordination of the Rubisco- and light-limited rates at the site's grow |
| vcmax25_pft2 | trait | plant_type | 19.91 | 17.85-22.5 | 58 | -2.98 | 0.40 | beyond 2 prior sd, pushed by gpp: diagnose it (plan §5.1) or relabel it effective | Kattge et al. 2009 (Glob. Change Biol. 15:976) for evergreen coniferous forest,  |
| capacity_low_temp_pft2 | effective | plant_type | 274.5 | 273.6-276.8 | 277.9 | -1.27 | 0.97 |  | ED2's temperate and conifer Vm_low_temp, 4.71 C (Moorcroft et al. 2001); ED2 giv |
| stomatal_g1 | trait | plant_type | 3.438 | 3.28-3.602 | 1.994 | +1.06 | 0.09 |  | EEO: the least-cost optimum at the site's growing-season daytime climate (priors |
| stomatal_g1_pft2 | trait | plant_type | 3.619 | 3.006-4.328 | 1.994 | +1.16 | 0.36 |  | EEO: the least-cost optimum at the site's growing-season daytime climate (priors |
| stomatal_g0 | trait | plant_type | 0.001255 | 0.0005614-0.00318 | 0.01 | -1.91 | 0.85 |  | the PFT file's value, a factor e either way: the residual conductance is poorly  |
| stomatal_g0_pft2 | trait | plant_type | 0.002078 | 0.0007982-0.005801 | 0.01 | -1.44 | 0.92 |  | the PFT file's value, a factor e either way: the residual conductance is poorly  |
| z0m_ratio | effective | site | 0.05807 | 0.05452-0.06193 | 0.13 | -2.36 | 0.19 |  | the base value (ED2 0.13), +-0.04: closed canopies give z0/h of about 0.06-0.15  |
| kappa | observation | observation | 0.9053 | 0.8691-0.9426 | 1 | -1.44 | 0.58 |  | the gap between the provider's night-time and daytime partitionings (RECO agains |

Fixed at their priors: leaf_angle_mean_pft2 (uninformed: posterior/prior sd ratio 0.91 >= 0.9); wstress_sref_stomata (uninformed: posterior/prior sd ratio 1.00 >= 0.9); wstress_sref_stomata_pft2 (uninformed: posterior/prior sd ratio 1.00 >= 0.9)

## Targets

| target | chi2/n | sigma scale | model/tower | n (cal) |
|---|---|---|---|---|
| lw_up | 1.06 | 1.68 | 1.003 | 1852 |
| le | 0.97 | 1.76 | 0.981 | 2496 |
| h | 1.15 | 3.00 | 0.969 | 857 |
| gpp | 0.95 | 1.90 | 1.007 | 618 |
| ustar | 0.56 | 1.00 | 0.983 | 1032 |

lw_up, model/tower by local hour: 0 h 0.99, 1 h 0.99, 2 h 1.00, 3 h 1.00, 4 h 1.00, 5 h 1.00, 6 h 1.01, 7 h 1.02, 8 h 1.02, 9 h 1.02, 10 h 1.02, 11 h 1.01, 12 h 1.01, 13 h 1.01, 14 h 1.01, 15 h 1.00, 16 h 1.00, 17 h 1.00, 18 h 0.99, 19 h 0.99, 20 h 1.00, 21 h 1.00, 22 h 0.99, 23 h 0.99
lw_up, by light: SW quartile 2 0.99, SW quartile 3 1.01, SW quartile 4 1.01

le, model/tower by local hour: 4 h 1.63, 5 h 0.92, 6 h 0.82, 7 h 0.91, 8 h 0.92, 9 h 0.91, 10 h 0.93, 11 h 0.98, 12 h 0.97, 13 h 0.98, 14 h 1.01, 15 h 1.08, 16 h 1.09, 17 h 1.21, 18 h 1.24, 19 h 1.29
le, by light: SW quartile 1 1.07, SW quartile 2 1.01, SW quartile 3 0.95, SW quartile 4 0.98

h, model/tower by local hour: 4 h 0.25, 5 h 0.22, 6 h 0.42, 7 h 0.61, 8 h 0.64, 9 h 0.75, 10 h 0.87, 11 h 0.92, 12 h 0.95, 13 h 1.02, 14 h 1.13, 15 h 1.36, 16 h 2.07, 17 h -14.41, 18 h -2.00
h, by light: SW quartile 1 -4.34, SW quartile 2 1.12, SW quartile 3 0.92, SW quartile 4 0.84

gpp, model/tower by local hour: 5 h 0.93, 6 h 0.86, 7 h 0.97, 8 h 0.93, 9 h 0.96, 10 h 0.98, 11 h 0.99, 12 h 1.03, 13 h 1.07, 14 h 1.08, 15 h 1.08, 16 h 1.06, 17 h 1.10, 18 h 0.83
gpp, by light: SW quartile 1 0.94, SW quartile 2 0.97, SW quartile 3 1.02, SW quartile 4 1.05

ustar, model/tower by local hour: 4 h 0.81, 5 h 0.77, 6 h 0.77, 7 h 0.83, 8 h 0.91, 9 h 0.94, 10 h 1.04, 11 h 1.07, 12 h 1.03, 13 h 1.03, 14 h 1.05, 15 h 1.01, 16 h 1.04, 17 h 1.05, 18 h 1.22
ustar, by light: SW quartile 1 0.98, SW quartile 2 0.93, SW quartile 3 0.98, SW quartile 4 1.04

## kappa

kappa = 0.905: the tower's respiration 3.14 -> 3.47 umol m-2 s-1; GPP 6.05 -> 6.29 (2.29 -> 2.38 kgC m-2 yr-1), means over the records with GPP and its respiration measured, day and night; GPP gains the scaled respiration by day only.

## Uncertainty

Laplace, from the final gradient matrix -- LOCAL ONLY (the linearity check failed).
- alternative gpp_ustar: the fit did not use the provider's threshold, or GPP's diagnostic finds no plateau
- alternative closure (Bowen (H and LE scaled together) against the fit's attribution shares): largest linear shift 12.62 sd -- refitted: leaf_angle_mean 24.37->18.87, vcmax25 92.97->94.63, vcmax25_pft2 19.91->18.67, capacity_low_temp_pft2 274.5->275.1, stomatal_g1 3.438->3.505, stomatal_g1_pft2 3.619->11.09, stomatal_g0 0.001255->0.000223, stomatal_g0_pft2 0.002078->0.002101, z0m_ratio 0.05807->0.05988, kappa 0.9053->0.8811
- alternative partitioning (the provider's daytime partitioning (GPP_DT, RECO_DT)): largest linear shift 2.10 sd -- refitted: leaf_angle_mean 24.37->24.24, vcmax25 92.97->86.32, vcmax25_pft2 19.91->16.93, capacity_low_temp_pft2 274.5->275.7, stomatal_g1 3.438->3.599, stomatal_g1_pft2 3.619->4.364, stomatal_g0 0.001255->0.0006022, stomatal_g0_pft2 0.002078->0.002201, z0m_ratio 0.05807->0.05884, kappa 0.9053->0.8938
