# Demography

ED-MORE represents a site's vegetation as **cohorts**: plants of one plant functional type (PFT) and one
size, with a density $n$ [plant m⁻²], living in **patches**, pieces of the site's ground that share one
disturbance history. This page gives the laws and operators that change that structure over the slow
steps:
- plants grow, die and are recruited;
- treefalls open gaps;
- fusion, fission, the cull and patch termination keep the number of cohorts and patches bounded,
  without changing what the stand does.

The page has two halves:
- the **vital rates** (§1), the laws that decide how each cohort and patch should change;
- the **operators** that apply them (§2–§4).

They are kept apart in the code (`src/README.md` rule 4): an operator takes rates as arguments and
never evaluates a law, so the same operators apply rates computed elsewhere (the Python `apply_rates`
path). §5 lists what each operator conserves and what it hands to the slow ledger
([conservation](conservation.md)).

Where each step's operators run, and in what order, is in [order_of_processes](order_of_processes.md)
§4–§5. The keys are in `docs/configuration.md`, "Cohort and patch fusion".

## 1. The vital rates

### Growth

A cohort's growth is its carbon allocation's ([plant_carbon_allocation](plant_carbon_allocation.md)).
Demography reads it as a diameter growth rate.
- Mortality uses the **prospective** rate $g$ [cm yr⁻¹]: the diameter the step's wood growth
  $`a_{\mathrm{wood}}`$ would give, from the allometry's inverse $`D(W)`$ at the cohort's wood carbon
  $W$, less the current diameter $D$, over the step $\Delta t$:

  ```math
  g = \frac{D\big(\max(W + a_{\mathrm{wood}},\,0)\big) - D}{\Delta t} \qquad(1)
  ```

- The **growth history** holds the realized diameter growth rate of each of the last $N$ steps, where
  $N$ spans `[demography].growth_memory_days`. A new sample replaces the oldest. The mean is

  ```math
  \bar g = \frac{1}{\min(N, n_{\mathrm{s}})}\sum_{k} g_k \qquad(2)
  ```

  over the $`n_{\mathrm{s}}`$ samples held.
- The sum is taken over the buffer itself, not carried as a running total, so a window of zeros
  after a large sample cannot read slightly negative.
- A new cohort's window is as long as its life, up to $N$ steps.
- The same window bounds the absolute sink limit on wood growth
  ([plant_carbon_allocation](plant_carbon_allocation.md) eq. (3)).

### Mortality

The hazard [yr⁻¹] is Camac et al. (2018)'s, additive in a growth-independent baseline and a
low-growth term:

```math
m = \gamma + \alpha\,e^{-\beta\,g_{\mathrm{eff}}} \qquad(3)
```

- $`g_{\mathrm{eff}}`$ is the growth average $`\bar g`$ once the cohort has a sample, and the
  prospective rate $g$ before.
- The average comes from the steps **before** this one: the step's own sample enters the history
  after the rates are taken.
- $\gamma$, $\alpha$ and $\beta$ are the PFT traits `[pft].mort_gamma`, `[pft].mort_alpha` and
  `[pft].mort_beta`.

**A cohort with no leaves and no storage** cannot flush again
([plant_carbon_allocation](plant_carbon_allocation.md) §2, "A cohort with no leaves and no storage").
It dies at least at the rate that takes its density to the floor $`n_{\min}`$
(`[demography].negligible_nplant`) in one step:

```math
m_{\mathrm{floor}} = \frac{\ln\big(\max(n, n_{\min})/n_{\min}\big)}{\Delta t} \qquad(4)
```

The monthly cull then removes it (§3). This is a special case, not a starvation hazard: a cohort that
still has leaves or storage dies at the rate eq. (3) gives.

### Recruitment

A cohort's reproduction carbon $r$ [kgC plant⁻¹ per step] is its allocation's reproduction share.
- Only the share $e$ of it becomes recruits. The rest is dead seed and dead seedling, which falls with
  the canopy as litter ([soil_carbon](soil_carbon.md) §5).
- The share falls with the leaf area index $L_p$ of the cohort's patch: under a deep canopy fewer
  seeds germinate, and fewer seedlings live to the smallest tracked size.

  ```math
  e = e_{\mathrm{pft}}\,\exp(-k_{\mathrm{shade}}\,L_p) \qquad(5)
  ```

  ($`e_{\mathrm{pft}}`$ = `recruit_carbon_efficiency`, $`k_{\mathrm{shade}}`$ = `recruit_shade_decay`; both `[pft]` keys).
- The recruitment rate of a PFT in a patch [plant m⁻² yr⁻¹] sums its cohorts' contributions and the
  baseline seed rain $`s_{\mathrm{pft}}`$ (`[pft].seed_rain_recruits`):

  ```math
  R = s_{\mathrm{pft}} + \sum_{\mathrm{cohorts}} n\,\frac{r}{\Delta t}\,\frac{e}{c_{\min}} \qquad(6)
  ```

- $`c_{\min}`$ is the carbon of one plant at the smallest tracked height (`[pft].min_cohort_height`), on
  the PFT's allometry.
- A PFT with `[pft].include_pft = 0` recruits nothing, seed rain included.

### Disturbance

Treefall gaps open once a year. The share of each patch's area disturbed, from the yearly rate
$\lambda$ (`[demography].patch_disturbance_rate`) over the one-year interval, is

```math
f = 1 - e^{-\lambda \cdot 1\,\mathrm{yr}} \qquad(7)
```

A cohort at least `[demography].disturbance_survive_height` tall dies in the gap; a shorter one
survives into it (§4).

## 2. Applying the rates

**Growth and mortality are applied every slow step**, as tendencies the driver computes from the rates.

- **The carbon pools are additive.** Each of the four pools moves by its change over the step:

  ```math
  X^{n+1} = X^{n} + \Delta X \qquad(8)
  ```

  - A pool the step empties can land a few roundings below zero. A result below zero by no more
    than four epsilons of $`|X|+|\Delta X|`$ is set to 0.
    - The tendency is $`(X_{\text{new}} - X)/\Delta t`$, so the applier reaches a pool the driver
      emptied through four roundings (subtract, divide, multiply, add), each within half an epsilon.
      Four epsilons bound that with a factor of two to spare.
    - A density decline aimed at the floor takes as many: the log of a quotient, and a rate divided
      and multiplied by $`\Delta t`$.
  - Anything further below is an overdraw: an operator that took more than the pool held. It is
    left as it is, so the bounds check stops the run on that step and names the pool, rather than a
    clamp creating the carbon the overdraw took.
- **The density moves in log space**, $`n^{n+1} = n^{n}\,e^{-m\Delta t}`$, stopped at the floor.
  - A decline that reaches $`n_{\min}`$, to within the round-off of comparing the two in log space,
    ends exactly on it, so the cull finds it there.
  - A cohort already at or below the floor is not taken lower.
- **Geometry is never integrated.** Diameter, height, basal area, aboveground biomass and leaf area
  are derived from the committed carbon after every operator that moves carbon or SLA.
- **Mortality's litter is valued on the new pools.** The live store $`n\,p`$ (density times a
  per-plant pool) changes over the step as

  ```math
  n_1 p_1 - n_0 p_0 = n_0\,(p_1 - p_0) + (n_1 - n_0)\,p_1 \qquad(9)
  ```

  - That is growth at the old density, and mortality at the new pools.
  - So the carbon mortality removed is $`(n_0 - n_1)\,p_1`$, every pool of the dying plants. It goes
    to the patch's litter, split as [soil_carbon](soil_carbon.md) §5 describes.
  - Valuing it on $`p_0`$ would leave the dying plants' growth unaccounted every step: small, but
    one-signed.
  - The density drop is taken from the commit, not recomputed from the hazard, so a cohort the floor
    stopped is not credited with litter for plants still standing.
  - The dying plants' tissue and surface water goes to the ground with its enthalpy, down the same
    channel turnover shedding uses.
- **The recruit pool is credited every slow step**, $`P \leftarrow P + R\,\Delta t`$ for each PFT and
  patch.
  - Valued at $`c_{\min}`$ per plant, the credit is exactly the establishing share of the
    reproduction carbon the parents were debited that step, plus the seed rain, which arrives from
    outside the site.
  - A monthly credit from one step's rate would be a 12-point sample that depends on which days are
    month boundaries.
- **Patches age** by one step.

## 3. Cohorts

Within each patch the cohorts are kept **tallest first**: height descending, diameter descending on
ties. The overtopping LAI (§6) and the fusion scans rely on the order, and every operator that moves a
cohort re-sorts.

The monthly operators, in order, are recruitment, fusion, the cull, fission and a re-sort.

**Recruitment.** When a pool reaches `[demography].min_recruit_size` plants m⁻², one cohort is
founded.
- Its density is the pool, at `[pft].min_cohort_height` on its PFT's allometry, and the pool is reset.
- A smaller pool carries over, so rare recruiters still establish eventually.
- The recruit is born in its patch's state: at the canopy-air temperature, its tissues holding the
  water of the top soil layer's potential.
- What it holds beyond $`c_{\min}`$ per plant is drawn from outside the model. Recruitment stands in
  for everything between a seed and a sapling — germination, the seedling's photosynthesis,
  transpiration and energy balance — none of which the model represents. That is carbon, the heat of
  its tissues and their water, all declared to the ledger.
- Drawing it from the free atmosphere rather than the patch's canopy air is deliberate: crediting the
  canopy air with a seedling's uptake while representing none of its respiration, transpiration or
  shading would add one term of a missing process.

**Fusion** keeps the most crowded patch near `max_cohort` cohorts.
- It runs in passes, at least one and at most `n_cohort_fusion_iter`. The first pass uses
  `cohort_size_tol_min`; while a patch still holds more than `max_cohort` cohorts, the next pass raises
  the tolerance geometrically, toward `cohort_size_tol_max`.
- Within a patch, a taller cohort absorbs every shorter cohort of its PFT whose height differs by
  less than the tolerance times the patch's top height, unless their combined LAI would reach
  `cohort_lai_cap`.
- The tolerance scales with the stand, so a cohort that has pulled ahead in a young stand keeps its
  place in the light.
  - At a 2 % tolerance, two cohorts in a 3 m regrowth stand merge only within about 6 cm, and in a
    30 m forest within about 0.6 m.
  - Tied to the PFT's maximum height instead, the tolerance would fold each month's recruits into a
    young stand's leaders and re-merge the halves of every split. A dense stand could then never
    form the size hierarchy it thins by.
- Fusion conserves:
  - the plant number;
  - the four carbon pools (each a density-weighted mean per plant), and with them aboveground
    biomass and leaf area;
  - the cohort's water and tissue energy, blended per plant.

  The diameter and the rest of the geometry are re-derived from the merged carbon, never averaged.
- The plastic traits blend on the quantity each is per ([plant_traits](plant_traits.md) §5): SLA and
  the leaf turnover rate on leaf carbon, $V_{\mathrm{cmax25}}$ and $R_{\mathrm{d25}}$ on leaf area.
  The temperatures follow from the fused energies and heat capacities.
- The survivor keeps its own growth history and predawn water status. The two cohorts already match
  in height and LAI, so these differ little.
- Each fusion asserts that carbon, biomass and leaf area are conserved, or stops the run.

**The cull** removes the cohorts too small to track:
- those whose $`n\cdot\mathrm{AGB}`$ is below `min_cohort_agb`;
- those at or below the density floor, which mortality stops at (§2). The test includes equality, or
  a cohort that died out would be kept for ever.

A culled cohort's carbon becomes litter and its water goes to the ground.

**Fission** splits a cohort whose LAI exceeds `cohort_lai_cap` into two of half the density, in
passes that repeat until none exceeds it (at most `n_cohort_fusion_iter`).
- Their wood carbon differs by a relative $`\pm\varepsilon`$ (`split_eps`), and the other pools are
  equal.
- This conserves every pool and the biomass exactly:
  $`\tfrac12 n\,w(1+\varepsilon) + \tfrac12 n\,w(1-\varepsilon) = n\,w`$.
- The carbon is perturbed, not the diameter. The diameter is derived from the carbon (§2), so the
  split stays exact where the height is capped at `hgt_max`; a perturbed diameter would conserve
  the carbon only through the allometry's inverse.
- Both daughters keep the parent's temperatures. The second takes a new identity.

## 4. Patches

The patches are kept sorted by site, then oldest first.

The yearly operators, after the month's, are the treefall gaps, then patch fusion and termination
followed by the cohort operators they call for. Each works within one site: a patch's area is a
fraction of its site, so no operator moves area between sites.

**Treefall gaps.** A share $f$ (eq. (7)) of every patch's area becomes one new, age-0 gap per site.
- **The gap inherits the area-weighted mean of its donors' stores:** canopy air, snow, soil water and
  heat, soil carbon, the day's litter and decomposition integral, and the recruit pool. A gap does
  not sterilise the ground it opens.
- **It starts with none of its own rates**: the shed water, the CO₂ owed to the canopy air, the
  integrator's step seed and the cost estimate. It has no cohorts of its own yet.
- **Survivors** (shorter than the threshold) are copied into the gap, at the density the disturbed
  share of their patch gives over the gap's area.
- **The killed canopy's carbon becomes litter** in the gap, and its tissue and surface water goes to
  the gap's ground.
- **Each donor keeps** $`(1-f)`$ of its area with all its cohorts.
- Site area is conserved. Plant number is not: the killed canopy is the disturbance.

**Fusion** pairs patches of one site with the same disturbance history whose light profiles are
alike.
- A patch's light profile is the fraction of full sunlight left under each of `n_height_layers` height
  layers:

  ```math
  I(h) = \exp\big(-k\,L_{\ge h}\big) \qquad(10)
  ```

  with $`L_{\ge h}`$ the leaf area index at and above the layer, and $k$ =
  `light_ext_coeff_fusefiss` (0.5 by default).
- Two patches are alike when their profiles differ by at most the tolerance on average and by at most
  `patch_light_maxdev_factor` times it in any layer. Layers fully lit in both are not compared.
- The tolerance grows geometrically from `patch_light_tol` toward `patch_light_tol_max` over the
  passes. Two empty patches always fuse.
- Fusion blends every store, and the age, on area, which conserves each area-extensive amount, and
  scales the cohorts' densities by their patch's share of the fused area, which conserves plant
  number.
- The passes stop once no site holds more than `max_patch` patches. There is at least one, at
  `patch_light_tol`, and at most `n_patch_fusion_iter`.
- A pass leaves the cohort block as it is and regroups it once at its end, instead of after every
  fusion.

**Termination.** A patch below `min_patch_area` is **merged** into its site's largest surviving patch,
not deleted.
- Dropping it and renormalizing the survivors' areas would redistribute the whole site, by
  $`\tfrac{a}{1-a}\sum_i a_i X_i - a X_{\mathrm{dropped}}`$ in every conserved quantity. That vanishes
  only if the dropped patch held the survivors' mean, and a doomed patch is atypical by construction,
  usually a fresh gap.
- Merging conserves exactly; the renormalization afterwards corrects round-off only.
- The merges regroup the cohort block once, as a fusion pass does, and the cohorts are sorted
  tallest first again, the order the cohort fusion after them scans in.
- Every site keeps at least its largest patch.

## 5. What each operator conserves

| Operator | Conserves | Hands to the slow ledger |
|---|---|---|
| growth and mortality (§2) | each pool's change is its tendency's | the day's net carbon from the fast loop (in); the litter (out); the growth and storage respiration owed to the canopy air (out); the dying plants' water (out); the tissues' thermal mass |
| recruit pool (§2) | — | the seed rain (in) |
| recruitment (§3) | — | the recruit's carbon beyond $`c_{\min}`$, heat and water (in) |
| cohort fusion, fission (§3) | plant number, the four pools, biomass, leaf area: exact, asserted | the tissues' thermal mass |
| the cull (§3) | — | the culled cohort's carbon, to the litter, and its water (out) |
| treefall gaps (§4) | site area | the killed canopy's carbon, to the gap's litter, and its water (out); the thermal mass |
| patch fusion, termination (§4) | area-extensive amounts, plant number, area: exact | the thermal mass of the cohort operators they call |

"The tissues' thermal mass": a cohort's heat capacity follows its biomass, so an operator that makes,
removes or re-derives tissue at unchanged temperature changes the tissue energy with no flux
([conservation](conservation.md) §4).

## 6. The overtopping LAI

Each cohort's overtopping LAI is the leaf area index of all **taller** cohorts in its patch:

```math
L_{\uparrow,j} = \sum_{i\,:\,h_i > h_j} n_i\,A_{\mathrm{leaf},i} \qquad(11)
```

- Cohorts of equal height share one layer and do not shade each other, so the result does not depend
  on their arbitrary order. Recruits are born at one shared size and tie exactly.
- The radiation puts tied cohorts in one layer too, by the same rule.
- It is the light environment of the trait plasticity ([plant_traits](plant_traits.md) eq. (1)) and
  the competition index of the demography.
- It is refreshed after every step's re-sort and after every calendar restructuring.

## Parameters

| Key | Meaning |
|---|---|
| `[pft].mort_gamma`, `mort_alpha`, `mort_beta` | eq. (3) [yr⁻¹, yr⁻¹, yr cm⁻¹] |
| `[pft].recruit_carbon_efficiency`, `recruit_shade_decay` | eq. (5) |
| `[pft].seed_rain_recruits`, `include_pft` | eq. (6) |
| `[pft].min_cohort_height` | the smallest tracked height, at which recruits are born |
| `[demography].growth_memory_days` | the growth history's window, eq. (2) |
| `[demography].negligible_nplant`, `min_cohort_agb` | the density floor (eq. (4)) and the cull (§3) |
| `[demography].max_cohort`, `n_cohort_fusion_iter`, `cohort_size_tol_min`, `cohort_size_tol_max`, `cohort_lai_cap`, `split_eps`, `conservation_tol` | cohort fusion and fission (§3) |
| `[demography].max_patch`, `n_patch_fusion_iter`, `patch_light_tol`, `patch_light_tol_max`, `patch_light_maxdev_factor`, `n_height_layers`, `light_ext_coeff_fusefiss`, `min_patch_area` | patch fusion and termination (§4) |
| `[demography].patch_disturbance_rate`, `disturbance_survive_height` | eq. (7) and the survivors (§4) |
| `[demography].min_recruit_size` | the pool that founds a cohort (§3) |
| `[demography].demography_on`, `recruitment_on`, `cohort_fusion_on`, `cohort_fission_on`, `patch_disturbance_on`, `patch_fusion_on` | `demography_on` gates every operator below; each of the others switches one. The cull and patch termination have no switch of their own and run with `demography_on` (§3–§4) |

## Code map

| Concept | Routine |
|---|---|
| the demographic laws, eqs. (1), (3)–(7); the founding test and treefall survival | `ed_demography`: `npp_to_growth`, `camac_mortality`, `decline_to_floor_mortality`, `recruit_efficiency`, `npp_to_recruitment`, `founds_cohort`, `disturbed_fraction`, `survives_treefall` |
| $c_{\min}$ | `ed_pft_derive`: `derive_recruit_carbon` |
| the rates per cohort and patch; the month's founding and the year's treefall decisions | `ed_slow_state`: `cohort_recruit_efficiency`, `compute_vital_rates`, `decide_recruitment`, `decide_disturbance` |
| the growth history, eq. (2) | `ed_slow_update`: `update_cohort_derivatives`, `fill_cohort_deriv`; the window length `growth_window_steps` in `ed_config` |
| applying growth and mortality, eq. (8) | `ed_slow_update`: `commit_vegetation_step`, `update_cohort_states`, `advance_pool`, `advance_density` |
| mortality's litter and water, eq. (9); the cull's and the treefall's litter; the recruit pool | `ed_slow_state`: `accumulate_mortality_litter`, `shed_mortality_water`, `kill_litter`; `ed_slow_update`: `accumulate_recruit_pool` |
| sort, fusion, fission, the cull, recruitment | `ed_slow_update`: `sort_cohorts`, `fuse_cohorts`, `fuse_2_cohorts`, `split_cohorts`, `terminate_cohorts`, `apply_recruitment` |
| treefall gaps, patch fusion, termination, eq. (10) | `ed_slow_update`: `apply_patch_disturbance`, `fuse_patches`, `fuse_patch_into`, `terminate_patches`, `light_profile_of` |
| the overtopping LAI, eq. (11) | `ed_slow_update`: `update_overtopping_lai`; `ed_numerics`: `equal_height_layers` |
| the order of the operators | `ed_slow_dynamics`: `vegetation_dynamics`, `restructure_stand` |
| tests | `test/slow/test_demography.f90` (the laws), `test/slow/test_slow_update.f90` (the operators), `test/slow/test_slow_step.f90` (the step) |

## References

- **Camac, J. S., Condit, R., FitzJohn, R. G., et al. (2018).** Partitioning mortality into
  growth-dependent and growth-independent hazards across 203 tropical tree species. *PNAS* 115(49):
  12459–12464 — the mortality hazard, eq. (3).
- **Moorcroft, P. R., Hurtt, G. C. & Pacala, S. W. (2001).** A method for scaling vegetation
  dynamics: the ecosystem demography model (ED). *Ecol. Monogr.* 71:557 — cohorts and patches.
- **Longo, M., et al. (2019).** The biophysics, ecology, and biogeochemistry of functionally diverse,
  vertically and horizontally heterogeneous ecosystems: the Ecosystem Demography model, version 2.2.
  *GMD* 12:4309 — ED2's fusion, fission and treefall disturbance.
