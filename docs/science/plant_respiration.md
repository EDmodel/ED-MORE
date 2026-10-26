# Plant respiration

Autotrophic respiration returns a large fraction of gross photosynthesis to the atmosphere. ED-MORE charges
it in **three places on two timescales**: leaf dark respiration inside the gas-exchange solver and stem
and fine-root **maintenance** respiration in `ed_respiration`, both on the fast step (`dt_fast`), and
**growth (construction)** respiration once per slow step (`[run].dt_slow`, a day by default) inside the
carbon allocator. This page collects the three limbs, states the accounting rule that keeps them from
being double-counted, and gives the equations for the two non-leaf maintenance terms: ED2 forms (Chambers
et al. 2004 for stem), re-referenced to the 25 °C datum. Their inputs come from the fast loop's
state: the cohort's wood temperature, the soil temperature by layer, and the cohort's size and pools.

## 1. The autotrophic budget

```math
R_a \;=\; \underbrace{R_d}_{\text{leaf, fast step}} \;+\; \underbrace{R_{\mathrm{stem}} + R_{\mathrm{root}}}_{\text{maintenance, fast step}} \;+\; \underbrace{R_g}_{\text{growth, slow step}} \qquad(1)
```

- **Leaf dark respiration** $R_d$ is part of the photosynthesis solution, not a separate flux: the leaf kernel
  scales a per-PFT (and plastic) $`R_{d25}`$ to leaf temperature and returns $`A = A_g - R_d`$. It is derived
  in [`leaf_gas_exchange.md`](leaf_gas_exchange.md) and only *placed* here.
- **Stem and fine-root maintenance** are computed every fast step from the two kernels in §2 and §3.
- **Growth respiration** $R_g$ is charged once per slow step by the allocator, on realized growth (§5).

The fast loop integrates the three fast-step limbs onto the cohort as `leaf_resp_accum`,
`stem_resp_accum` and `root_resp_accum` $`[\mathrm{kgC\,plant^{-1}}]`$ (converting with `umol_2_kgc` at
each fast step), and the slow step sums exactly those three into the allocator's `resp_maint`, so

```math
\mathrm{NPP} \;=\; G \;-\; \big(R_d + R_{\mathrm{stem}} + R_{\mathrm{root}}\big) \qquad(2)
```

**NPP as reported is GPP net of MAINTENANCE respiration only.** Growth respiration is deliberately not
subtracted here: it is charged *inside* the allocator, on the growth actually built, so subtracting it
again would double-count it. That one sentence is the whole accounting rule for the autotrophic budget.

## 2. Stem maintenance respiration

Woody tissue respires on two bases, each with its own per-PFT rate at 25 °C: the **stem surface
area** (the ED2 / Chambers et al. 2004 form) and the **sapwood volume**, the living wood (Ryan 1990).
A PFT sets one rate to 0 to use the other alone; `stem_resp_sapwood25` defaults to 0, so a PFT file
that does not name it respires by surface alone:

```math
R_{\mathrm{stem}} \;=\; s(T_w)\,\Big[\,r^{\mathrm{stem}}_{25}\,A_{\mathrm{stem}}
\;+\; r^{\mathrm{sap}}_{25}\,V_{\mathrm{sap}}\Big],
\qquad
A_{\mathrm{stem}} \;=\; \frac{\pi\,(10^{-2}\,\mathrm{dbh})\,h \;+\; \pi\,\mathrm{WAI}/n}{f_{\mathrm{ag}}},
\qquad
V_{\mathrm{sap}} \;=\; \frac{2\,C_{\mathrm{sap}}}{1000\,\rho} \qquad(3)
```

with $`r^{\mathrm{stem}}_{25}`$ = `stem_resp_factor25` the per-PFT surface rate
$`[\mathrm{\mu mol\,CO_2\,m^{-2}\,stem\,s^{-1}}]`$ and $`r^{\mathrm{sap}}_{25}`$ = `stem_resp_sapwood25`
the per-PFT sapwood rate $`[\mathrm{\mu mol\,CO_2\,m^{-3}\,sapwood\,s^{-1}}]`$, both at 25 °C;
$`C_{\mathrm{sap}}`$ the cohort's sapwood carbon [kgC/plant] (the sapwood fraction of its wood carbon,
which counts the coarse roots as $`A_{\mathrm{stem}}`$ does through $`f_{\mathrm{ag}}`$), as dry mass
over the PFT's wood density $\rho$ [g cm⁻³]; $h$ the cohort height [m], $n$ = `nplant`
$`[\mathrm{plant\,m^{-2}}]`$, and WAI the **per-ground** wood area index (hence the $`/n`$ returning the
branch term to a per-plant area). The result is per plant
$`[\mathrm{\mu mol\,CO_2\,plant^{-1}\,s^{-1}}]`$; the caller multiplies by `nplant` for the patch total.
Four details are load-bearing:

- **Temperature driver.** $`T_w`$ is the cohort's own **wood temperature** (`fast_dynamic_state_t%wood_temp`),
  the prognostic store solved by the vegetation energy balance — not canopy air, not leaf temperature. Its
  larger heat capacity damps the diurnal cycle ([`vegetation_energy_dynamics.md`](vegetation_energy_dynamics.md)).
- **Aboveground fraction.** The cylinder-plus-branch area is the *aboveground* surface; dividing by
  $`f_{\mathrm{ag}}`$ = `aboveground_frac` scales it to the whole woody plant, so belowground structure
  respires in proportion. It is a per-PFT trait gathered per cohort — the same number demography and
  cohort fusion read — never a run-uniform constant.
- **The woody gate.** `is_woody = .false.` (a grass) returns exactly zero, whatever the rates.
- **Surface or sapwood.** Maintenance respiration of stems is linear in sapwood (live-cell) volume in
  conifers (Ryan 1990; Ryan et al. 1995), while tropical stem efflux follows surface and volume both
  (Cavaleri et al. 2006; Rowland et al. 2018). The two bases spread a stand's stem respiration over
  sizes differently: per m² of leaf, the surface term rises about 2.6-fold from a 1 cm to a 120 cm stem
  on the tropical allometry, the sapwood term about 40-fold, so the sapwood basis charges saplings far
  less and canopy trees more. Volume rather than carbon: biomass grows as $`\rho^{0.976}`$, so
  $`V_{\mathrm{sap}}`$ barely depends on the wood density, and denser sapwood respires less per gram
  (Westerband et al. 2022). A volume basis is only as good as the sapwood fraction behind it (Rowland et
  al. 2018), which is why a sapwood rate is best set against a stand's total stem respiration.

The surface rate has no size factor of its own. ED2's optional DBH exponent on it (~0.0041 cm⁻¹,
$`10^{\,0.0041\,\mathrm{dbh}}`$) is noisy and pulls against the shrinking-live-fraction argument for
large stems; how stem respiration grows with size is set instead by the choice of basis above, and a
PFT that wants it to grow faster than the surface uses the sapwood term.

## 3. Fine-root maintenance respiration

Fine roots respire **per unit fine-root carbon**, the ED2 per-`broot` form, with the temperature response
summed over the soil layers the roots occupy:

```math
R_{\mathrm{root}} \;=\; r^{\mathrm{root}}_{25}\;\; B_{\mathrm{root}} \sum_k \phi_k\, s(T_k) \qquad(4)
```

with $`r^{\mathrm{root}}_{25}`$ = `root_resp_factor25` $`[\mathrm{\mu mol\,CO_2\,kgC^{-1}\,s^{-1}}]`$ at
25 °C, $`B_{\mathrm{root}}`$ the cohort's fine-root carbon pool `fineroot_carbon`
$`[\mathrm{kgC\,plant^{-1}}]`$, $`\phi_k`$ the patch's fine-root profile over the active layers
(`root_frac`, summing to 1) and $`T_k`$ the layer's soil temperature. Unlike the stem limb, this one
tracks a **pool** directly: a cohort that sheds fine roots stops paying for them the same day.

**Temperature driver.** The response is summed over the layers rather than taken at the
root-weighted mean temperature: $s$ is nonlinear, so $`\sum_k \phi_k\, s(T_k) \ne s\big(\sum_k \phi_k T_k\big)`$.
The peaked form is convex below its optimum and concave near it, so the error a mean temperature
would make flips sign with the season rather than biasing one way. At Ithaca on a 2 m column the
layer sum differs from the response at the mean by +4.8 % in December and −8.3 % in June, but only
−0.17 % in the annual mean: a seasonal correction to root respiration, not an annual-budget one. The
sum is patch-uniform, so `root_zone_temp_scale` evaluates it once per patch per fast step and the
kernel broadcasts it over the cohorts. A cohort with no fine roots respires nothing.

## 4. The shared temperature response

Both maintenance limbs, and leaf $R_d$, scale from a 25 °C baseline through the **same two functions**
in `ed_therm_lib`:

```math
s_{\mathrm{arr}}(T) = \exp\!\Big[\frac{E_a}{R\,T_{25}}\Big(1 - \frac{T_{25}}{T}\Big)\Big],
\qquad
s(T) = s_{\mathrm{arr}}(T)\,\frac{f_H(T_{25})}{f_H(T)}, \quad f_H(T) = 1 + \exp\!\Big[\frac{\Delta S\,T - H_d}{R\,T}\Big] \qquad(5)
```

The plain Arrhenius form rises without bound; the **peaked** form multiplies it by a high-temperature
deactivation envelope (Medlyn et al. 2002), normalised so that $`s(T_{25}) = 1`$ — a thermal optimum with a
roll-off above it. $`T_{25}`$ = `t_ref_photo` = 298.15 K is the model's single reference, so
`*_resp_factor25` **is** the rate at 25 °C and the kernels perform no reference conversion; a parameter
quoted at another datum (ED2 references stem and root respiration at 15 °C) is converted **once, where the
PFT file is written**, never at runtime.

Leaf $R_d$ reaches these functions through the dispatcher `temp_response(form, …)`, so its form follows the
run's configured selection; stem and root call `peaked_arrhenius_scale` directly and are always peaked.
Both routes end in the same code, and that is the point: the temperature response (`ed_therm_lib`) sits
in the shared library at the root of the dependency DAG precisely so leaf, stem and root cannot drift
into three near-copies of an Arrhenius — a change to the deactivation normalisation or to the clamped
`safe_exp` reaches every tissue at once. Stem and root take the leaf's $`E_a, H_d, \Delta S`$ too,
`[leaf_physiology]`'s `ea_rd`, `hd_rd` and `ds_rd` (46390 / 200000 / 490 in the reference config), set in
one place. Stem and root always use them in the peaked form, while leaf $`R_d`$ takes
`temp_response_form`, so the three tissues share one thermal shape under the default `peaked`.

## 5. Growth respiration on realized growth

Building one unit of growth tissue costs $`(1+g)`$ carbon, of which the fraction $g$ =
`growth_resp_factor` is respired. The allocator applies that cost **inside** each funding step, per pool,
so the charge is

```math
R_g \;=\; g\,\big(a_{\mathrm{leaf}} + a_{\mathrm{root}} + a_{\mathrm{wood}} + a_{\mathrm{repro}} + E\big) \qquad(6)
```

over the growth $`a_p`$ the cohort actually built this step, and the root exudate $E$ of the sink limit on
wood growth, which costs the same to make and export ([`plant_carbon_allocation.md`](plant_carbon_allocation.md) §2; zero without a
limit). **Storage refill is exempt** — moving carbon
into a nonstructural pool is a transfer of sugar, not construction.

Charging on *realized* growth rather than on *available* carbon is the correct choice for three reasons.
It resolves, exactly and without iteration, the circularity that construction cost reduces the carbon
available for growth, which changes the growth. It gets the degenerate cases right: a dormant or
carbon-starved cohort that funds no growth pays no construction cost, whereas a charge on the
pre-allocation balance would bill it for tissue it never built. And it keeps carbon bound for storage out
of the charge, which a balance-based form cannot. This is a real term in the ecosystem carbon balance,
not a rounding detail: it is $`g/(1+g)`$ of all carbon entering growth, 23 % at the shipped $g = 0.3$.

## 6. Where the respired carbon goes

All four limbs reach the canopy air space as part of the biotic CO₂ source driving the prognostic
`can_co2`, assembled at the start of each fast step (`build_frozen`) as

```math
\mathrm{NEE}_{\mathrm{biotic}} \;=\; R_d + R_{\mathrm{stem}} + R_{\mathrm{root}} + R_h \;-\; G \;+\; \dot{c}_{\mathrm{slow}} \qquad(7)
```

$`[\mathrm{\mu mol\,m^{-2}\,ground\,s^{-1}}]`$, positive toward the canopy air. The three fast-step limbs
are exhaled in the very fast step that computes them. Growth respiration cannot be, because it is known
only once per slow step, so the slow step's allocation driver commits the per-patch total as
`slow_co2_rate`, a frozen rate the *next* fast window emits. That rate carries an opposite-signed starvation
correction (`growth_resp - deficit`): when storage could not fund the maintenance the fast loop had already
reported, the over-report is netted back here rather than by rewriting the previous step. Real maintenance
respiration is limited by the substrate the plant has, so a maintenance flux no pool could pay is read as
an over-report of the fast loop, which respires without that limit (§7); an operator-split model corrects
it on the next step, in one channel with one sign convention.

$`R_h`$ is the heterotrophic limb from the CENTURY pools ([`soil_carbon.md`](soil_carbon.md)); the
canopy-air budget these fluxes enter is in [`canopy_air_space_biophysics.md`](canopy_air_space_biophysics.md).

## 7. What is not here

- **No thermal acclimation of respiration** — no running-mean-temperature reference shift for any
  tissue: the 25 °C baselines are fixed traits, and the kernels carry no acclimation-temperature
  argument (photosynthetic capacity does acclimate, [`leaf_gas_exchange.md`](leaf_gas_exchange.md), "Thermal acclimation").
- **No reproductive-tissue or storage maintenance** — seed carbon and the non-structural pool respire
  nothing.
- **No per-layer root biomass** — the temperature response is resolved by layer (§3), but the fine-root
  pool is one per-plant pool spread over the layers by the root profile, not a prognostic pool per layer.
- **No substrate or moisture limitation of maintenance respiration** — it runs at the temperature-set rate
  whatever the carbon status; §6's starvation correction is after-the-fact bookkeeping, not a rate law.
- **No size factor on the stem surface rate** (ED2's optional DBH exponent, §2). Stem respiration
  still grows with size through its basis: the example 05 PFTs scale it through the sapwood term, which
  charges saplings less and canopy trees more than the surface term (§2).

The open items are in [`ROADMAP.md`](../ROADMAP.md) §7.

## Parameters

Per-PFT traits the respiration consumes (PFT file, `[pft]`).

| Symbol | Config key | Meaning |
|---|---|---|
| — | `is_woody` | woody gate; `.false.` ⇒ stem respiration is identically 0 |
| $`r^{\mathrm{stem}}_{25}`$ | `stem_resp_factor25` | stem rate per surface $`[\mathrm{\mu mol\,CO_2\,m^{-2}\,stem\,s^{-1}}]`$ at 25 °C (shipped 0.06) |
| $`r^{\mathrm{sap}}_{25}`$ | `stem_resp_sapwood25` | stem rate per sapwood volume $`[\mathrm{\mu mol\,CO_2\,m^{-3}\,sapwood\,s^{-1}}]`$ at 25 °C (optional, default 0) |
| $`r^{\mathrm{root}}_{25}`$ | `root_resp_factor25` | fine-root baseline $`[\mathrm{\mu mol\,CO_2\,kgC^{-1}\,s^{-1}}]`$ at 25 °C (shipped 0.30) |
| $`f_{\mathrm{ag}}`$ | `aboveground_frac` | aboveground fraction of woody carbon; scales the stem area (eq. (3)) |
| $g$ | `growth_resp_factor` | construction-cost fraction charged on realized growth (shipped 0.3) |
| $`R_{d25}`$ | `rd25` (plastic per cohort) | leaf dark-respiration capacity — see [`plant_traits.md`](plant_traits.md), [`leaf_gas_exchange.md`](leaf_gas_exchange.md) |
| $`E_a, H_d, \Delta S`$ | `[leaf_physiology]` `ea_rd`, `hd_rd`, `ds_rd` | peaked-Arrhenius parameters, the leaf $R_d$'s, shared by stem and root |

## Code map

| Concept | Routine |
|---|---|
| stem maintenance respiration (eq. (3)) | `ed_respiration`: `stem_maintenance_respiration` (`elemental pure`: scalar = one cohort, array = one patch) |
| fine-root maintenance respiration (eq. (4)) | `ed_respiration`: `fine_root_maintenance_respiration`, with the layer-summed temperature scale `root_zone_temp_scale` |
| shared temperature response (eq. (5)) | `ed_therm_lib`: `arrhenius_scale`, `peaked_arrhenius_scale`, `temp_response` |
| leaf dark respiration $R_d$ | `ed_leaf_gas_exchange` (inside the Ci solve; see [`leaf_gas_exchange.md`](leaf_gas_exchange.md)) |
| per-fast-step call site + patch totals | `ed_fast_integ_frozen`: `canopy_maintenance_respiration`, `build_frozen` (eq. (7)) |
| fast step → slow step integration | `ed_fast_state`: `accumulate_fast_to_slow` → `leaf_resp_accum` / `stem_resp_accum` / `root_resp_accum` on the cohort SoA |
| slow-step maintenance sum handed to the allocator | `ed_cohort_state_vars`: `cohort_carbon_supply` (`resp_maint`) |
| growth respiration (eq. (6)) | `ed_plant_carbon_allocation`: `growth_respiration`, charged inside `fill_carbon_demand` / `plant_carbon_allocation` |
| growth respiration → canopy air | `ed_slow_state`: `compute_carbon_allocation` → `patch%slow_co2_rate` → `ed_fast_integ_frozen` |
| per-PFT traits + TOML loading | `ed_pft_params` (`pft_table_t`), `ed_config_io` (`[pft]` block; every key is required except `stem_resp_sapwood25`, §2) |
| tests | `test/fast/test_respiration.f90` (the 25 °C identity, the layered root response, the sapwood term, the grass gate, `aboveground_frac` inverse scaling; the soil's heterotrophic respiration there too, [soil_carbon.md](soil_carbon.md)); growth respiration in `test/slow/test_plant_carbon_allocation.f90` |

## References

- **Chambers et al. (2004)**, *Ecol. Appl.* 14:S72 — respiration of a central Amazonian forest, the stem
  maintenance form ED2 uses.
- **Ryan (1990)**, *Can. J. For. Res.* 20:48 — stem maintenance respiration scales with sapwood volume.
- **Medlyn et al. (2002)**, *Plant Cell Environ.* 25:1167 — the peaked temperature response of Rubisco
  and electron-transport rates.
- **Ryan et al. (1995)**, *Oecologia* 101:133 — woody-tissue maintenance respiration of four conifers
  tracks sapwood volume.
- **Cavaleri et al. (2006)**, *Glob. Change Biol.* 12:2442 — wood CO₂ efflux in a primary tropical rain
  forest, by surface area, volume and biomass.
- **Rowland et al. (2018)**, *New Phytol.* 218:1393 — drought stress and tree size determine stem CO₂
  efflux in a tropical forest.
- **Westerband et al. (2022)**, *Ann. Bot.* 129:633 — nitrogen and tissue density as drivers of
  woody-tissue respiration.
