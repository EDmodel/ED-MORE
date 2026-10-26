# Plant carbon allocation

This page explains how ED-MORE turns each cohort's carbon budget into tissue growth, leaf and fine-root
turnover and litter. One **stateless, elemental** kernel, `plant_carbon_allocation`
(`ed_plant_carbon_allocation`), does it. It covers the science of ED2's `growth_balive` and
`structural_growth`, **unified into one step**, and it follows FATES PARTEH Hypothesis 1 ("Allometrically
Guided, Carbon Only"): every pool has an allometric target, the net carbon gain fills toward the
targets in **priority order**, and the residual advances stature (wood). All quantities are carbon
$`[\mathrm{kgC\,plant^{-1}}]`$ over one slow step; every biomass↔carbon conversion is folded into the PFT
(plant functional type) traits once at initialization, so the kernel never converts.

**Timescale and inputs.** The allocation runs once per slow step (`[run].dt_slow`, a day by default),
in the vegetation phase of the slow step: after leaf phenology and (if on) trait acclimation, before the
vital rates and the commit of growth and mortality (`ed_slow_dynamics`: `vegetation_dynamics`;
[`order_of_processes.md`](order_of_processes.md)). Its carbon comes from the fast loop, which
accumulates each cohort's GPP and maintenance respiration over the step's fast steps; its leaf
targets come from the cohort's size, its plastic SLA and the phenology tendencies.

## 1. The carbon budget

The step's gross primary production $`G`$ and maintenance respiration $`R_m`$ (leaf, stem and fine root;
both accumulated by the fast loop, or a stub when `fast_biophysics_on` is false) give the **net carbon**:

```math
C_{\mathrm{net}} = G - R_m \qquad(1)
```

which may be negative. Growth (construction) respiration is **not** subtracted here: it is charged inside
the ladder, on realized growth only (§3). The kernel also receives the current storage (nonstructural
carbon) $`C_{store}`$, three demands (leaf, fine root and storage, each a deficit toward its target; the
leaf and fine-root demands are formed after this step's turnover, §4) and the reproduction fraction
`repro_frac`. Wood takes the residual.

### Leaf and fine-root resorption

A fraction $`f_r`$ = `retained_carbon_fraction` of **all** leaf and fine-root loss returns to the
non-structural pool instead of entering litter. With $B$ the background leaf turnover and $S$ the
senescence of `leaf_turnover_step` ([`plant_phenology.md`](plant_phenology.md) §5, eq. (9)), and $R$ the
fine-root turnover:

```math
\Delta C_{leaf} = -(B + S), \qquad
\Delta C_{root} = -R, \qquad
\Delta C_{store} \mathrel{+}= f_r (B + S + R), \qquad
\text{litter} = (1-f_r)(B + S + R) \qquad(2)
```

**Each pool loses the full loss either way.** Crediting storage while removing only the litter
share would *create* carbon. The split is exact by construction.

Background turnover is resorbed like senescence: `leaf_turnover_rate_toc` and `fineroot_turnover_rate`
are tissue loss rates, not litterfall, so the litter is the unresorbed share of the loss.

**Default $`f_r = 0.3`$.** Most measured resorption is of N and P rather than C, so the carbon
fraction is uncertain and is a per-PFT key; examples 03–05 set it to 0. Resorption moves
carbon from the litter path to the plant reserve, and the reserve then funds growth, so a larger
$`f_r`$ raises storage, GPP and LAI and lowers soil carbon.

## 2. The priority ladder

Targets come from allometry: leaf $`L^{*}=`$ `size_to_leaf_carbon`, fine root $`q\,L^{*}`$ (with $q$ the
root:leaf ratio) and storage $`c_s\,(1+q)\,L^{*}`$ (`storage_carbon_target`), the carbon of $`c_s`$ full
flushes. A **demand** is the deficit `target − pool`; the leaf and fine-root demands arrive already
**flush-capped** ($`\Delta_{\mathrm{fl}} = r_{\mathrm{fl}} L^{*} dt`$), so a dormant canopy
($`r_{\mathrm{fl}}=0`$) presents zero growth demand. When $`C_{\mathrm{net}}\ge0`$ the kernel funds, in
order:

1. **leaf + fine-root growth** toward the flush-capped demand, funded from NPP, then storage;
2. **storage refill** toward $`c_s (1+q) L^{*}`$, from NPP only and with no construction cost;
3. **reproduction**, a fraction $`f_{\mathrm{repro}}`$ of the post-storage residual (zero below the maturity height);
4. **wood**, the residual sink: everything left becomes structural growth.

Wood needs no explicit demand (it absorbs the remainder), and any cohort below its PFT's maturity
height, `min_repro_height` [m], sets $`f_{\mathrm{repro}}=0`$. The height is set outside ED-MORE, one per
PFT. Thomas (1996)'s relative size at onset of maturity is a guide: his 37 Malaysian species matured at
0.20–0.75 of their asymptotic height, 0.42 on average, so a PFT's `min_repro_height` can be that share of
the height it reaches on its own curve. The height is the same at every site: where Chave et al.
(2014)'s stress index E makes trees short, a cohort reaches it at a larger diameter, or not at all.

When $`C_{\mathrm{net}}<0`$ (for example a leafless canopy at bud-break, where GPP is near zero but stem
and root maintenance still runs), the maintenance debt is paid **from storage first**, and then leaf and
fine-root growth may still draw the **remaining** reserves, so spring leaf-out is storage-funded and
not deadlocked. A plant whose storage cannot cover its maintenance has a **deficit**: the
maintenance respiration no pool paid for. The driver nets the deficit out of the next day's CO₂ release
([`plant_respiration.md`](plant_respiration.md) §6); no tissue is destroyed, and the kernel never mutates
a pool. The full priority order is therefore **maintenance debt → leaf/fine-root growth (NPP then
storage) → storage refill (NPP only) → reproduction → wood**.

### The size allometry

The targets and every cohort's geometry come from its PFT's size allometry (`ed_allometry_lib`, one
record per PFT), whose coefficients are the PFT file's `[allometry]` block
([`configuration.md`](../configuration.md), "The PFT file"). Each key takes one value for every PFT or
one per PFT, so an understorey specialist and a canopy emergent can each have their own curve, cap and
leaf display. With $D$ the diameter [cm], $H$ the height [m] and $\rho$ the wood density [g cm⁻³]:

- **Height** follows one of two curves, each capped at the PFT's `hgt_max`: the pan-tropical power
  law $`\ln H = b_1 + b_2 \ln D`$ (`b1Ht`, `b2Ht`), or a log-quadratic,
  $`\ln H = a - E + b \ln D + c (\ln D)^2 + s^2/2`$ (`height_allometry = "chave2014"`), with $E$ the
  site's stress index (`chave_e`, the same for every PFT) and the PFT's coefficients defaulting to
  Chave et al. (2014) eq. 6a: $a$ = 0.893, $b$ = 0.760, $c$ = −0.0340 and the fit's residual standard
  error $s$ = 0.243. The $`s^2/2`$ turns the fit's log-scale median into its mean, as BIOMASS's
  `retrieveH` does. A site's own heights give a local fit of $a$, $b$ and $c$: Barro Colorado's, fitted
  uniformly in $`\ln D`$ over 1–300 cm to Cano et al. (2019)'s curve for 9884 trees of the Barro Colorado
  Nature Monument, is $a$ = 0.8898, $b$ = 0.7908, $c$ = −0.04923 (within 4.3 % of that curve). Height to
  diameter inverts the uncapped curve; on Chave's curve that is the smaller root of the quadratic in
  $`\ln D`$, the rising branch, which holds up to a diameter of about 7×10⁴ cm.
- **Biomass** is Chave et al. (2014) eq. 4 in carbon, $`\mathrm{AGB} = c_1\,\rho^{c_2}\,(D^2 H)^{c_2}`$
  [kgC plant⁻¹] (`agb_c1`, `agb_c2`), with $`c_2 = 0.976`$ and $`c_1 = 0.0673/2`$: Chave's coefficient for
  dry mass, halved at the carbon fraction 0.5. Without the halving, the biomass would be dry mass
  labelled as carbon.
- **Leaf area** of one stem, one-sided, is $`A_{leaf} = l_1\,(D^2 H)^{l_2}`$ [m²] (`lai_b1`, `lai_b2`), so
  a cohort's LAI is $`n_{plant}\,A_{leaf}`$. The shipped $`l_1 = 0.4677/2`$ and $`l_2 = 0.641`$ are the
  pan-tropical BAAD fit, halved as ED2 applies it ([`ed2_comparison.md`](../ed2_comparison.md)).
- **Pool targets** restate the same allometry in carbon, so the two descriptions of a plant cannot
  disagree. The leaf target is $`A_{leaf}/\mathrm{SLA}`$ (`size_to_leaf_carbon`), so a cohort at its
  target shows exactly the allometric leaf area ([`plant_traits.md`](plant_traits.md) §4). The wood
  carbon is $`\mathrm{AGB}/f_{ag}`$ (`size_to_wood_carbon`, $`f_{ag}`$ = `aboveground_frac`): the total
  woody carbon, coarse roots included, whose aboveground share is exactly the Chave biomass.
- **Carbon to diameter.** Wood carbon is the prognostic size, and the diameter follows from it
  (`wood_to_dbh`); `carbon_to_structure` derives a cohort's diameter, height, basal area, AGB and leaf
  area from its wood and leaf carbon. Under the power law the inverse is closed form,
  $`\mathrm{AGB} \propto D^{P}`$ with $`P = c_2(2 + b_2)`$ below the height cap and $`P = 2c_2`$ above it,
  the two meeting at the cap; under the log-quadratic it is the smaller root of a quadratic in $`\ln D`$.

### Sink limitation of wood growth

Carbon supply is not the only limit on growth: cambial cell division and expansion run at a bounded
rate however much carbon the leaves fix (sink limitation; Körner 2015). Wood growth is therefore capped
by `growth_sink_limitation` (`ed_plant_carbon_allocation`), which returns the most wood the step can
build. It has two forms, each a per-PFT trait that is off at zero; with both set the tighter applies:

```math
a_{\mathrm{wood}}^{\max} = \min\!\Big(r_{\max}\,C_{\mathrm{wood}}\,\Delta t,\;\;
W\big(D + \max(g_{\max}D^{c}\,T - \Delta D_T,\,0)\big) - W(D)\Big) \qquad(3)
```

The first is **relative**: $`r_{\max}`$ is the PFT's `max_relative_growth_rate` [yr⁻¹] and
$`C_{\mathrm{wood}}`$ the cohort's wood carbon at the step's start. The second is **absolute**:
$`g_{\max}`$ is the PFT's `max_absolute_growth_rate` [cm yr⁻¹], the most a 1 cm stem's diameter can grow,
$c$ its `max_absolute_growth_exponent` (default 0) with $D$ in cm, and $`W(D)`$ the wood carbon of a stem of
diameter $D$ on the model's allometry (`size_to_wood_carbon`). The absolute cap is a **rate over the growth
window** $T$ (`[demography].growth_memory_days`, 365 days in the shipped config, the window mortality
averages growth over), not over each step: $`\Delta D_T`$ is the diameter the stem has already grown in
the window, read from the cohort's ring buffer of realized growth, and $T$ includes this step (for a
cohort younger than the window, its age). A step the limit binds therefore brings the window's growth to
exactly $`g_{\max}D^{c}T`$, and a stem that grew nothing through a dry season may make up the year's
diameter in the wet season. A per-step cap, $`g_{\max}D^{c}\Delta t`$, would instead throw away every
day without spare carbon, so a sapling with carbon to spare all wet season could never reach the cap's
annual rate.

The cambium lays down a bounded width of wood, so the absolute form scales with the stem's surface
rather than its mass: per unit wood it is far more generous to a sapling than to a large tree. At BCI
the upper quantiles of census diameter growth at a given light rise about as $`D^{0.5}`$ (example 05),
between a size-independent width ($`c=0`$) and the relative form ($`c\approx1`$).

When the residual would build more than $`a_{\mathrm{wood}}^{\max}`$, the kernel builds
$`a_{\mathrm{wood}}^{\max}`$ and turns the rest of the residual into **root exudate**,
$`E = C_{\mathrm{resid}}/(1+g) - a_{\mathrm{wood}}^{\max}`$. Making and exporting exudate costs the same
$`(1+g)`$ per unit as tissue, so growth respiration (eq. (4)) is the same whether or not the limit binds;
the limit only changes where the carbon ends up. The driver puts $E$ in the patch's below-ground labile
litter (`labile_soil`), where soil respiration returns it to the air. A cohort below the limit is
unaffected, so the limit only bites where carbon is plentiful (open-grown trees), and with both limits
at 0 (the default) there is no limit. The form is deliberately simple; a mechanistic sink (temperature-
or water-limited cambial activity) replaces the body of `growth_sink_limitation` without touching the
allocator.

### A cohort with no leaves and no storage

When $`C_{\mathrm{net}} < 0`$ the maintenance debt is paid from storage first; what storage cannot
cover is the **deficit**. No pool pays it and no tissue is destroyed: the fast loop has already
respired it, and the driver corrects that on the next step ([`plant_respiration.md`](plant_respiration.md)
§6).

One case is different: a cohort whose deficit is larger than **all the carbon in its leaves**. Its
storage is spent, and burning every leaf it has could not pay one step of its upkeep, so nothing is
left to pay for a flush and it cannot recover. It dies. Its mortality for the step
(`decline_to_floor_mortality`) takes its density to the density floor, `[demography].negligible_nplant`,
down the path every death takes: its carbon goes to litter and its water to the ground. The next
monthly cull removes the cohort, so it is gone within a month.

The test needs no threshold: its scale is the cohort's own maintenance respiration. A canopy holds far
more leaf carbon than one step's upkeep, and a cohort with storage pays it, so a leafless cohort living
on its reserves through a dry season or a winter is not affected. Every other cohort dies at the
growth-driven hazard; there is no separate starvation hazard. A growth- or storage-based starvation
mortality is not modelled.

## 3. Growth respiration on realized growth

Building one unit of **growth tissue** (leaf, fine root, wood or reproduction), or exporting one unit of
root exudate (§2), consumes $`(1+g)`$ carbon, where $g$ = `growth_resp_factor`: the fraction $g$ is
respired as **growth (construction) respiration**. Storage refill is a 1:1 carbon transfer
(nonstructural sugar has no construction cost), so it is exempt. Charging $`(1+g)`$ *inside* the funding
step resolves, exactly and without iteration, the circularity that growth respiration reduces the carbon
available for growth, which changes the growth. If the pools built are
$`a_{\mathrm{leaf}},a_{\mathrm{root}},a_{\mathrm{wood}},a_{\mathrm{repro}}`$ and the exudate is $E$, then

```math
R_g = g\,(a_{\mathrm{leaf}} + a_{\mathrm{root}} + a_{\mathrm{wood}} + a_{\mathrm{repro}} + E) \qquad(4)
```

A charge on the whole pre-allocation balance would also bill the carbon bound for storage, which builds
nothing ([`plant_respiration.md`](plant_respiration.md) §5).

## 4. Leaf display as emergent replaceability

Turnover is taken before the ladder. The leaf loss (background turnover plus senescence) and the flush
cap come from one routine, `leaf_turnover_step` ([`plant_phenology.md`](plant_phenology.md) §5), not a
replaceable/non-replaceable split. Whether a loss is replaced is **emergent**: the loss opens a leaf
deficit, and the flush-capped growth step (step 1 of the ladder, §2) refills it *if* flushing is active.
A flushing canopy holds full cover because its turnover is continuously refilled; a deciduous canopy in
dormancy (flush tendency $`\to0`$, so flush cap $`\to0`$) is not refilled and senesces to near bare; an
evergreen's senescence stops at its leaf-cover floor. The losses decay the current pool, so a
partially-built canopy is never over-shed.

## 5. Carbon closure and litter

The kernel is **growth-only**: it returns the per-pool growth $`a_{\mathrm{leaf}},a_{\mathrm{root}},
a_{\mathrm{wood}},a_{\mathrm{repro}}\ge0`$ and the net storage change $`\mathrm{npp}_{\mathrm{store}}=`$
refill − drawdown. The **turnover and shed are applied upstream**, before the ladder
(`cohort_carbon_demand`): the pools handed to the kernel are already net of this step's shed, and the
driver forms the net leaf and root change $`\mathrm{npp}_{\mathrm{leaf}}=a_{\mathrm{leaf}}-\ell_{\mathrm{leaf}}`$
(and likewise for root). The **growth-side** budget the kernel closes on every call is

```math
\big(a_{\mathrm{leaf}} + a_{\mathrm{root}} + a_{\mathrm{wood}} + a_{\mathrm{repro}} + \mathrm{npp}_{\mathrm{store}} + E\big) \;-\; \mathrm{deficit} \;=\; (G - R_m) \;-\; R_g \qquad(5)
```

with `deficit` the maintenance no pool paid for (zero unless storage ran out; the pools' change then
exceeds $`G - R_m - R_g`$ by that amount, which the driver nets out of the next day's CO₂) and $E$ the
root exudate of the sink limit (zero without one). Adding the driver's terms closes the whole plant: the
**plant-pool change = GPP − (maintenance + growth respiration) − litter − exudate − reproduction**, where
the litter is the unresorbed share of the shed, $`(1-f_r)(\ell_{\mathrm{leaf}}+\ell_{\mathrm{root}})`$
(eq. (2)). The driver puts it in the patch's litter, from which the soil-carbon pools and $`R_h`$ follow
([`soil_carbon.md`](soil_carbon.md)). The seed carbon that does not establish is litter too
([`demography.md`](demography.md) §1).

## Interface with other modules

`ed_plant_carbon_allocation` is a **pure kernel library**: it `use`s only the kinds, the constants, the
allometry and the phenology's leaf turnover (no `polygon_t`, no config, no PFT table). Everything is
passed as plain scalars (and the sink limit its PFT's allometry), so the **driver**
(`ed_slow_state`: `compute_carbon_allocation`) is the single place that assembles the inputs and
disposes of the outputs, calling the kernel once per cohort as an `elemental` sweep over the cohort
Structure-of-Arrays. This is the only plant-flux call on the slow carbon path.

**Inputs the driver gathers per cohort:**

| input | source |
|---|---|
| `gpp`, `resp_maint` | fast-loop accumulators on the cohort (`gpp_accum`, `leaf_resp_accum`, `stem_resp_accum`, `root_resp_accum`) when `fast_biophysics_on`; otherwise `gpp_ref` × leaf area × step, and no maintenance (`cohort_carbon_supply`) |
| `leaf_demand`, `fineroot_demand` | allometric deficit (`size_to_leaf_carbon`, using the **plastic** `cohort%sla`) after this step's shed, **flush-capped** by the phenology flush rate |
| `storage_demand`, `storage`, `repro_frac` | cohort state (`nonstructural_carbon`) + PFT traits (`storage_cushion`, `reproduction_investment_fraction`, zero below `min_repro_height`) |
| `growth_resp_frac` | PFT trait `growth_resp_factor` |
| this step's **shed** (litter) | `cohort_carbon_demand` calls `leaf_turnover_step` (background turnover at the cohort's `leaf_turnover_rate` while flushing and no leaves drop, plus senescence down to `min_leaf_cover`) and turns fine roots over on the same `turnover_gate`, for the carbon amounts it pre-subtracts from the pools |
| `wood_growth_max` | `growth_sink_limitation`, with the window's growth read from the cohort's ring buffer (`growth_window_so_far`) |

**Outputs the driver disposes of:**

| output | destination |
|---|---|
| per-pool growth `growth_{leaf,fineroot,wood,repro}`, `npp_store` | the driver forms the net per-pool change `growth − shed` in the step record `slow_step_data_t` (`ed_slow_state`; the net pool change, not the gross NPP the `CO_NPP_*` diagnostics report), which `commit_vegetation_step` (`ed_slow_update`) hands to `update_cohort_derivatives`, the `wood_carbon → dbh` flip (`carbon_to_structure`) and `update_cohort_states` |
| `growth_resp` | autotrophic-respiration accounting (with maintenance respiration) |
| `deficit` | netted, with the opposite sign, into the patch's `slow_co2_rate` ([`plant_respiration.md`](plant_respiration.md) §6); no tissue is destroyed. A deficit larger than the cohort's leaf carbon means it has no leaves and no storage, and it dies (§2) |
| leaf + fine-root **litter** (`shed`) | the resorbed share to storage, the rest to the patch's litter (`add_necromass`); the exudate to `labile_soil` |

`growth_respiration` lives in `ed_plant_carbon_allocation` beside the allocation kernel. The kernel is
orthogonal to the core engine (`ed_slow_update`), which only *applies* the tendency arrays the driver
backs out.

## Parameters

Per-PFT traits the allocation consumes (PFT file, `[pft]` unless noted).

| Symbol | Config key | Meaning |
|---|---|---|
| $g$ | `growth_resp_factor` | construction-cost fraction charged on growth |
| $`c_s`$ | `storage_cushion` | storage target as a multiple of the leaf + fine-root target |
| $q$ | `root_to_leaf_ratio` | fine-root : leaf target ratio |
| $`f_{\mathrm{repro}}`$ | `reproduction_investment_fraction` | fraction of the residual → reproduction (above maturity) |
| — | `min_repro_height` [m] | the height a cohort reproduces from |
| $`f_r`$ | `retained_carbon_fraction` (optional, default 0.3) | share of leaf and fine-root loss resorbed to storage, eq. (2) |
| $`r_{\max}`$ | `max_relative_growth_rate` (optional, default 0 = off) | sink limit on wood growth per unit wood carbon [yr⁻¹]; the residual above it is exuded |
| $`g_{\max}`$ | `max_absolute_growth_rate` (optional, default 0 = off) | sink limit on the diameter growth of a 1 cm stem [cm yr⁻¹], a budget over the growth window; the residual above it is exuded |
| $T$ | `[demography].growth_memory_days` | the window the absolute limit (and mortality's growth average) runs over |
| $c$ | `max_absolute_growth_exponent` (optional, default 0) | how that limit scales with diameter, $`g_{\max}D^{c}`$ |
| $`L^{*}`$ inputs | `sla`, `[allometry]` keys (`hgt_max`, …) | full-canopy leaf carbon via `size_to_leaf_carbon` |

The flush and senescence rates, `min_leaf_cover` and `bare_leaf_cover` are phenology traits
([`plant_phenology.md`](plant_phenology.md), Parameters); the allocation consumes their *carbon amounts*.

## Code map

| Concept | Routine |
|---|---|
| allocation ladder (elemental) | `ed_plant_carbon_allocation`: `plant_carbon_allocation` (+ private `fill_carbon_demand`, the $`(1+g)`$ funder) |
| growth (construction) respiration | `ed_plant_carbon_allocation`: `growth_respiration(npp_growth, g)` |
| sink limit on wood growth | `ed_plant_carbon_allocation`: `growth_sink_limitation` (relative and absolute limits, the tighter applies), with its window from `growth_window_so_far` |
| leaf loss (turnover + senescence) + flush cap | `ed_plant_phenology`: `leaf_turnover_step` |
| turnover-first carbon demand per cohort | `ed_plant_carbon_allocation`: `cohort_carbon_demand` |
| the step's GPP and maintenance respiration | `ed_cohort_state_vars`: `cohort_carbon_supply` |
| size allometry, pool targets | `ed_allometry_lib`: `size_to_leaf_carbon`, `size_to_wood_carbon`, `storage_carbon_target`, `wood_to_dbh` |
| per-cohort orchestration (slow step) | `ed_slow_state`: `compute_carbon_allocation` (shed first, then the flush-capped demands, one elemental call, then the net change and litter) |
| geometry flip (wood_carbon → dbh) | `ed_slow_update`: `update_cohort_derivatives`; `ed_allometry_lib`: `carbon_to_structure` |

## References

- **Chave et al. (2014)**, *Glob. Change Biol.* 20:3177 — pan-tropical tree allometry (biomass eq. 4,
  height eq. 6a).
- **Körner (2015)**, *Curr. Opin. Plant Biol.* 25:107 — sink limitation of growth.
- **Thomas (1996)**, *Oikos* 76:145 — relative size at the onset of maturity in 37 Malaysian rain
  forest tree species.
- **Cano et al. (2019)**, *Biogeosciences* 16:847 — tree height allometry of Barro Colorado Island.
- **Falster et al. (2015)**, *Ecology* 96:1445 — BAAD, the plant allometry database behind the leaf-area fit.
