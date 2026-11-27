# Fast calibration: harvard

## Validation

**Failed**: the calibrated cost below the default's, and no target's RMSE more than 10% worse.

| target | default | calibrated | change |
|---|---|---|---|
| lw_up | 1.57 | 1.81 | +15% |
| le | 2.16 | 1.88 | -13% |
| h | 3.11 | 3.23 | +4% |
| gpp | 2.75 | 2.13 | -22% |
| ustar | 1.14 | 0.63 | -45% |

RMSE in units of each target's sigma. Cost: default 15280.1, calibrated 14600.8.

The full record with the slow tier on is the next check: run the calibrated configs (budgets closed; dry-season GPP and LE no worse than the default's).

## Keys

| key | kind | scope | MAP | 68 % | prior centre | prior z | sd ratio | note | prior source |
|---|---|---|---|---|---|---|---|---|---|
| leaf_angle_mean | trait | plant_type | 32.71 | 24.86-42.29 | 33.9 | -0.12 | 0.91 |  | the PFT file's value, +-10 degrees; spherical leaves average 57.3. The range is  |
| vcmax25 | trait | plant_type | 90.82 | 86.03-95.47 | 73.98 | +0.46 | 0.13 |  | EEO: the coordination of the Rubisco- and light-limited rates at the site's grow |
| vcmax25_pft2 | trait | plant_type | 20.46 | 18.17-23.35 | 58 | -2.89 | 0.42 | beyond 2 prior sd, pushed by gpp: diagnose it (plan §5.1) or relabel it effective | Kattge et al. 2009 (Glob. Change Biol. 15:976) for evergreen coniferous forest,  |
| capacity_low_temp_pft2 | effective | plant_type | 277.9 | 275.8-280.7 | 277.9 | +0.01 | 0.63 |  | ED2's temperate and conifer Vm_low_temp, 4.71 C (Moorcroft et al. 2001); ED2 giv |
| stomatal_g1 | trait | plant_type | 3.304 | 3.156-3.457 | 1.994 | +0.98 | 0.09 |  | EEO: the least-cost optimum at the site's growing-season daytime climate (priors |
| stomatal_g1_pft2 | trait | plant_type | 3.253 | 2.605-4.04 | 1.994 | +0.95 | 0.43 |  | EEO: the least-cost optimum at the site's growing-season daytime climate (priors |
| stomatal_g0 | trait | plant_type | 0.001182 | 0.0005082-0.003195 | 0.01 | -1.96 | 0.91 |  | the PFT file's value, a factor e either way: the residual conductance is poorly  |
| stomatal_g0_pft2 | trait | plant_type | 0.004121 | 0.001748-0.009413 | 0.01 | -0.83 | 0.77 |  | the PFT file's value, a factor e either way: the residual conductance is poorly  |
| z0m_ratio | effective | site | 0.05806 | 0.05455-0.06187 | 0.13 | -2.36 | 0.19 |  | the base value (ED2 0.13), +-0.04: closed canopies give z0/h of about 0.06-0.15  |
| kappa | observation | observation | 0.9195 | 0.8818-0.9582 | 1 | -1.22 | 0.59 |  | the gap between the provider's night-time and daytime partitionings (RECO agains |

Fixed at their priors: leaf_angle_mean_pft2 (uninformed: posterior/prior sd ratio 0.91 >= 0.9); wstress_sref_stomata (uninformed: posterior/prior sd ratio 1.00 >= 0.9); wstress_sref_stomata_pft2 (uninformed: posterior/prior sd ratio 1.00 >= 0.9)

## Targets

| target | chi2/n | sigma scale | model/tower | n (cal) |
|---|---|---|---|---|
| lw_up | 1.11 | 1.81 | 1.001 | 1852 |
| le | 0.95 | 1.79 | 0.986 | 2496 |
| h | 1.07 | 3.00 | 1.042 | 857 |
| gpp | 0.97 | 1.92 | 1.000 | 618 |
| ustar | 0.59 | 1.00 | 0.993 | 1032 |

lw_up, model/tower by local hour: 0 h 0.99, 1 h 0.99, 2 h 0.99, 3 h 0.99, 4 h 0.99, 5 h 1.00, 6 h 1.01, 7 h 1.02, 8 h 1.02, 9 h 1.03, 10 h 1.02, 11 h 1.01, 12 h 1.01, 13 h 1.01, 14 h 1.01, 15 h 1.00, 16 h 0.99, 17 h 0.99, 18 h 0.99, 19 h 0.99, 20 h 0.99, 21 h 0.99, 22 h 0.99, 23 h 0.99
lw_up, by light: SW quartile 2 0.99, SW quartile 3 1.01, SW quartile 4 1.02

le, model/tower by local hour: 4 h 1.49, 5 h 0.82, 6 h 0.78, 7 h 0.91, 8 h 0.94, 9 h 0.92, 10 h 0.94, 11 h 0.99, 12 h 0.98, 13 h 0.98, 14 h 1.01, 15 h 1.07, 16 h 1.09, 17 h 1.23, 18 h 1.32, 19 h 1.53
le, by light: SW quartile 1 1.08, SW quartile 2 1.00, SW quartile 3 0.96, SW quartile 4 0.99

h, model/tower by local hour: 4 h 0.40, 5 h 0.62, 6 h 0.12, 7 h 0.60, 8 h 0.71, 9 h 0.85, 10 h 0.99, 11 h 1.03, 12 h 1.05, 13 h 1.10, 14 h 1.19, 15 h 1.40, 16 h 2.05, 17 h -12.76, 18 h -1.63
h, by light: SW quartile 1 -3.39, SW quartile 2 1.14, SW quartile 3 1.00, SW quartile 4 0.94

gpp, model/tower by local hour: 5 h 0.93, 6 h 0.84, 7 h 0.96, 8 h 0.93, 9 h 0.95, 10 h 0.97, 11 h 0.99, 12 h 1.02, 13 h 1.06, 14 h 1.07, 15 h 1.07, 16 h 1.05, 17 h 1.10, 18 h 0.85
gpp, by light: SW quartile 1 0.94, SW quartile 2 0.96, SW quartile 3 1.01, SW quartile 4 1.04

ustar, model/tower by local hour: 4 h 0.74, 5 h 0.71, 6 h 0.74, 7 h 0.82, 8 h 0.93, 9 h 0.97, 10 h 1.07, 11 h 1.10, 12 h 1.05, 13 h 1.05, 14 h 1.06, 15 h 1.03, 16 h 1.05, 17 h 1.05, 18 h 1.20
ustar, by light: SW quartile 1 0.96, SW quartile 2 0.93, SW quartile 3 1.00, SW quartile 4 1.07

## kappa

kappa = 0.919: the tower's respiration 3.14 -> 3.42 umol m-2 s-1; GPP 6.05 -> 6.25 (2.29 -> 2.37 kgC m-2 yr-1), means over the records with GPP and its respiration measured, day and night; GPP gains the scaled respiration by day only.

## Uncertainty

Laplace, from the final gradient matrix -- LOCAL ONLY (the linearity check failed).
- alternative gpp_ustar: the fit did not use the provider's threshold, or GPP's diagnostic finds no plateau
- alternative closure (Bowen (H and LE scaled together) against the fit's attribution shares): largest linear shift 11.37 sd -- refitted: leaf_angle_mean 32.71->20.3, vcmax25 90.82->93.16, vcmax25_pft2 20.46->19.97, capacity_low_temp_pft2 277.9->276.7, stomatal_g1 3.304->3.091, stomatal_g1_pft2 3.253->8.727, stomatal_g0 0.001182->0.0005033, stomatal_g0_pft2 0.004121->0.001561, z0m_ratio 0.05806->0.05836, kappa 0.9195->0.9069
- alternative partitioning (the provider's daytime partitioning (GPP_DT, RECO_DT)): largest linear shift 2.29 sd -- refitted: leaf_angle_mean 32.71->24.85, vcmax25 90.82->84.69, vcmax25_pft2 20.46->17.1, capacity_low_temp_pft2 277.9->276.8, stomatal_g1 3.304->3.374, stomatal_g1_pft2 3.253->4.697, stomatal_g0 0.001182->0.0007415, stomatal_g0_pft2 0.004121->0.001323, z0m_ratio 0.05806->0.0589, kappa 0.9195->0.9088
