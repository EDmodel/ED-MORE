# Leaf gas exchange

How ED-MORE computes the CO₂ uptake and water loss of a leaf: the photosynthetic demand for CO₂, the
stomatal conductance that supplies it, the effect of drought on both, and the solver that makes them
consistent. The code is
[`src/fast_dynamics/processes/ed_leaf_gas_exchange.f90`](../../src/fast_dynamics/processes/ed_leaf_gas_exchange.f90);
the public entry points are `solve_leaf_gas_exchange` and `leaf_gas_exchange_batch`. Each PFT's
(plant functional type's) parameters are flattened from the configuration once, at load, into the
leaf table of the PFT table (`ed_pft_derive`).

The kernel takes one leaf's environment as given (incident PAR, leaf temperature, vapour pressure
deficit, CO₂, pressure, water potentials, boundary-layer conductance) and does no radiation, energy
balance or hydraulics itself. The radiation page gives how the cohort's PAR is formed
([canopy_radiation_transfer](canopy_radiation_transfer.md), "PAR passed to photosynthesis"); the cohort
code divides it by the leaf area index (`cohort_leaf_par`).

Symbols: $A$ = net CO₂ assimilation, $A_g$ = gross assimilation, $R_d$ = leaf (dark) respiration
($A = A_g - R_d$); $C_a, C_s, C_i$ = CO₂ mole fraction (µmol mol⁻¹) in the ambient air, at the leaf
surface and in the intercellular space; $g_s$ = stomatal conductance **to water vapour**; $g_b$ = leaf
boundary-layer conductance to water vapour; $D$ = leaf-to-air vapour deficit; $\psi$ = water
potential. Two fixed diffusion constants appear throughout: the H₂O:CO₂ diffusivity ratios $1.6$
through the stomata (`gsw_2_gsc`) and $1.4$ through the boundary layer (`gbw_2_gbc`).

---

## 1. The coupled A–Cᵢ–gₛ system

Assimilation is fixed by the intersection of two curves in $C_i$:

- **Biochemical demand** — the Farquhar–von-Caemmerer–Berry (C3) or Collatz (C4) rate $A(C_i)$
  (§2), which *rises* with $C_i$.
- **Diffusive supply** — CO₂ drawn down through the stomata:

```math
A = \frac{g_s}{1.6}\,(C_s - C_i), \qquad C_s = C_a - 1.4\,\frac{A}{g_b} \qquad(1)
```

  (the $C_s$ term applies only when the boundary layer is resolved, `use_boundary_layer`;
  otherwise $C_s = C_a$). This *falls* with $C_i$.

The stomatal models set $g_s$ as a function of $A$, $C_s$ and $D$, so the whole system reduces to
**one nonlinear equation in $C_i$**, solved on a bracket by the Illinois method (§5).

---

## 2. Photosynthetic demand $A(C_i)$

**C3 (Farquhar et al. 1980).** Three potential rates, co-limited:

```math
A_c = V_{cmax}\frac{C_i-\Gamma^*}{C_i + K_c\left(1 + O/K_o\right)},\quad
A_j = J\frac{C_i-\Gamma^*}{4C_i + 8\Gamma^*},\quad
A_p = 3\,\mathrm{TPU} \qquad(2)
```

where $A_c$ is Rubisco-limited, $A_j$ RuBP-regeneration (light) limited and $A_p$
triose-phosphate-limited; $`\Gamma^*`$ is the CO₂ compensation point without respiration, and $K_c, K_o$
are the Rubisco Michaelis constants for CO₂ and O₂ ($O$ is the O₂ concentration). The
electron-transport rate $J$ is the smaller root of the non-rectangular hyperbola
$`\theta_J J^2 - (I_2 + J_{max})J + I_2 J_{max}=0`$, with
$`I_2 = \tfrac12\,\phi_{\text{PSII}}\,\alpha_{\text{leaf}}\,\mathrm{PAR}`$ (`electron_transport_j`). The shipped
$`\phi_{\text{PSII}}`$ = 0.74 gives $J$ an initial slope of 0.37 electrons per absorbed photon. At four
electrons per CO₂ (the $4C_i + 8\Gamma^*$ in $A_j$) that is 0.093 CO₂ per absorbed photon, the measured
maximum quantum yield of C3 leaves (Long et al. 1993). The maximum PSII efficiency Fv/Fm (0.85) in
its place would give the O₂-evolution yield, 0.106 (Björkman & Demmig 1987).

**$`\Gamma^*`$ is proportional to the oxygen partial pressure**, because it is set by Rubisco's
specificity for CO₂ over O₂: $`\Gamma^* = \tfrac{1}{2} O / S_{c/o}`$. The O₂ partial pressure is
the O₂ mole fraction times the air's pressure, as the CO₂ one is, so $`\Gamma^*`$ as a mole fraction
depends on the O₂ mole fraction and **not on pressure**. ED-MORE computes it as a mole fraction from the
O₂ mole fraction and the pressure its reference value is stated at:

```math
\Gamma^*(T, O) = \frac{\Gamma^*_{25}}{P_\mathrm{ref}}\,f_{\mathrm{Arr}}(T)\;\frac{O}{O_\mathrm{ref}},
\qquad O_\mathrm{ref} = 0.209,\quad P_\mathrm{ref} = 100\ \mathrm{kPa} \qquad(3)
```

The shipped `gstar25` = 4.275 Pa is Bernacchi et al. (2001)'s 42.75 µmol mol⁻¹ at 21 % O₂, written
as a partial pressure at 100 kPa, so $`\Gamma^*`$ at 25 °C is 42.75 µmol mol⁻¹ at any elevation.
$`K_c`$ and $`K_o`$, the Michaelis constants of dissolved CO₂ and O₂, are partial pressures and are
divided by the local pressure, so $`O/K_o`$ is the ratio of the two partial pressures
(`rubisco_kinetics` in `ed_leaf_gas_exchange`; `ed.plant.leaf.rubisco_kinetics` in Python).
Dividing `gstar25` by the local pressure instead would raise $`\Gamma^*`$ with elevation and cut
gross photosynthesis in dim light, where photorespiration sets the quantum yield.

The O₂ factor makes `[leaf_physiology].o2_mol_frac` reach both places oxygen enters the C3 demand,
the Michaelis term $`K_c(1 + O/K_o)`$ and the compensation point. It equals 1 in any default run. It
matters for palaeo-atmosphere work (Carboniferous and Permian O₂ is reconstructed near 30–35 %) and for
the low-O₂ sensitivity experiment, the standard way to isolate photorespiration.

**C4 (Collatz et al. 1992):** $A_c = V_{cmax}$; $A_j = \phi_{C4}\,\alpha_{\text{leaf}}\,\mathrm{PAR}$ is a
light-limited slope with the quantum yield $\phi_{C4}$ (`quantum_yield_c4`); and $A_p = k_p C_i$ is the
PEP-carboxylase CO₂ limitation. $`\Gamma^* \approx 0`$ because the CO₂-concentrating mechanism suppresses
photorespiration. `kp25` is a slope on the CO₂ partial pressure stated at standard pressure, so the
code scales $`k_p`$ by $`p/p_{std}`$ and $`k_p C_i`$ follows the partial pressure at any altitude.

The three rates are combined into $A_g$ either as a sharp $\min(A_c,A_j,A_p)$
(`[leaf_physiology].colimitation = "min"`) or as two nested smoothing quadratics
(`"quadratic"`, the same `quadratic_smaller_root` used for $J$). Then $A = A_g - R_d$.

**Temperature.** $`K_c`$, $`K_o`$ and $`\Gamma^*`$ scale from 25 °C by the plain Arrhenius form
$`f_{\mathrm{Arr}}`$ of eq. (3), with the activation energies of Bernacchi et al. (2001).
$`V_{cmax}`$, $`J_{max}`$, TPU, $`R_d`$ and the C4 $`k_p`$ take the form that
`[leaf_physiology].temp_response_form` selects, plain or peaked (`ed_therm_lib`;
[plant_respiration](plant_respiration.md) eq. (5)). TPU and $`k_p`$ have no temperature parameters of
their own: both take $`V_{cmax}`$'s, and with them its high-temperature deactivation under the peaked
form.

### The co-limitation curvatures are not $`\theta_J`$

When two limiting rates are combined smoothly, a curvature says how sharply the leaf switches from one
to the other. Each pathway has its own pair of these co-limitation curvatures, and neither is
`theta_j`:

| | first transition | second transition |
|---|---|---|
| C3 | `theta_cj_c3` — $A_c$ vs $A_j$ | `theta_ip_c3` — $\min(A_c,A_j)$ vs $A_p$ |
| C4 | `theta_cj_c4` — $A_c$ vs $A_j$ | `theta_ic_c4` — $\min(A_c,A_j)$ vs $A_p$ |

They share units and a functional form, and nothing else. `theta_j` is the curvature of the
electron-transport hyperbola and describes the light saturation of $J$. Implementations generally put the
co-limitation curvature much closer to 1 (CLM and FATES use about 0.98 for the C3 $A_c/A_j$ transition)
than the 0.7–0.9 that fits the $J$ hyperbola.

Using a $J$ curvature for the co-limitation would cost assimilation exactly where the model spends most
of its time, because $A_c$ and $A_j$ sit close together at ambient CO₂. The table gives the loss against
the sharp $\min(A_c,A_j,A_p)$ at 25 °C, 101.325 kPa and $C_i = 280$ µmol mol⁻¹, with the ratios of the
first shipped PFT ($`J_{max} = 1.8\,V_{cmax}`$, $`\mathrm{TPU} = 0.167\,V_{cmax}`$,
$`R_d = 0.015\,V_{cmax}`$); the net-rate column spans PAR 1500–2000 µmol m⁻² s⁻¹ and
$V_{cmax,25}$ from 60 to 150, with $J$ from the hyperbola:

| `theta_cj_c3` / `theta_ip_c3` | gross $A$, saturating light ($J = J_{max}$) | net $A$, PAR 1500–2000 µmol m⁻² s⁻¹ |
|---|---|---|
| 0.85 / 0.85 (a $J$ curvature on both transitions) | about −27 % | −31 to −34 % |
| 0.85 / 0.95 (on $A_c$/$A_j$ only) | about −24 % | −27 to −31 % |
| 0.98 / 0.95 (shipped) | about −10 % | −12 to −15 % |

At saturating light the shortfall is the same at every $V_{cmax,25}$, because all three rates scale
with $V_{cmax}$. The net rate loses more than the gross because $R_d$ does not shrink with it. The loss
is a systematic offset and not a regime effect, so no measured $V_{cmax}$ would reproduce a measured
assimilation rate if a $J$ curvature were used here. The shipped $V_{cmax,25}$ presets sit mid-range for
their PFT descriptions and are not tuned to offset a smoothing loss.

Kernels: `rubisco_kinetics`, `assimilation_demand_c3`, `assimilation_demand_c4`, `electron_transport_j`,
`combine_limits`.

---

### Thermal acclimation (opt-in)

Photosynthetic capacity adapts to the temperature a stand has grown in. With
`[leaf_physiology].thermal_acclimation = true` (default `false`), the peaked form's entropy term and
the capacity ratio follow the growth temperature $T_g$ (Kattge & Knorr 2007):

```math
\Delta S_{V} = 668.39 - 1.07\,T_g, \qquad
\Delta S_{J} = 659.70 - 0.75\,T_g, \qquad
\frac{J_{max,25}}{V_{cmax,25}} = 2.59 - 0.035\,T_g \qquad(4)
```

with $T_g$ in °C. The coefficients are `acclim_ds_vcmax_a`, `acclim_ds_vcmax_b`, `acclim_ds_jmax_a`,
`acclim_ds_jmax_b`, `acclim_jv_a` and `acclim_jv_b`, at the values shown; the capacity ratio is
floored at 0.1 far outside the fitted range. $`\Delta S`$ sets where the response peaks, so a
warm-grown stand runs a higher optimum: the $V_{cmax}$ optimum is 28.2 °C for $T_g$ = 10 °C and 38.2 °C
for $T_g$ = 30 °C.

Three choices here are deliberate:

- **The ratio acclimates too.** Shifting $`\Delta S`$ alone would change the shape of each response
  but pin the two branches in a fixed proportion, which is not what the study measured: a warm-grown
  plant invests relatively less in electron transport.
- **$T_g$ is the growth temperature, not the leaf temperature.** The fit is calibrated against the
  mean air temperature of the preceding weeks, so ED-MORE tracks an exponential running mean of the
  daily-mean air temperature (window `acclim_window_days`, default 30 d) and feeds that. An
  instantaneous leaf temperature would use the relation far outside what it was fitted to.
- **It is one value per site, not per cohort**, for the same reason: each site keeps a running mean of
  its own air temperature, after its terrain correction ([forcing](forcing.md) §8), and its own leaf
  table. Per-cohort acclimation would need a relation calibrated on leaf temperature
  ([ROADMAP](../ROADMAP.md)).

The running mean is prognostic and slow, so it is written to the state file; rebuilding a month of
memory after every restart would make the acclimated optimum jump. A mean not yet seeded takes the
first day's value outright, so the optimum starts from the site's own climate. It is applied by
refreshing the per-PFT leaf table once per slow step from the reference coefficients in the
configuration, which are never overwritten, so the refresh is idempotent (`acclimate_site_leaf_tables`
in `ed_slow_state`, `acclimate_leaf_photo_table` in `ed_plant_traits`). `thermal_acclimation`
requires the peaked form: with the Arrhenius form there is no $`\Delta S`$ to shift, and the loader
refuses the combination.

## 3. Stomatal-conductance models

The model is selected by `[leaf_physiology].stomatal_model` (`"leuning"`, `"medlyn"` or `"katul"`;
`SM_LEUNING`, `SM_MEDLYN`, `SM_KATUL` in the code). The first two give an **explicit** $g_s$ law;
Katul defines $g_s$ **implicitly** through an optimization.

### 3.1 Leuning (1995)

```math
g_s = g_0 + \frac{g_1\,A}{(C_s - \Gamma^*)\,(1 + D/D_0)} \qquad(5)
```

A semi-empirical law: conductance tracks assimilation, discounted by the CO₂ drawdown
$`(C_s-\Gamma^*)`$ and by a hyperbolic humidity response $(1 + D/D_0)$. Parameters: `stomatal_g0`
($g_0$), `stomatal_g1` ($g_1$), `stomatal_d0` ($D_0$, in Pa). Code: `stomata_gs_leuning`.

### 3.2 Medlyn et al. (2011) — unified stomatal optimization (USO)

```math
g_s = g_0 + 1.6\left(1 + \frac{g_1}{\sqrt{D}}\right)\frac{A}{C_s}, \qquad D \text{ in kPa} \qquad(6)
```

with $D$ floored at `[leaf_physiology].vpd_min` (0.05 kPa), since $g_1/\sqrt{D}$ is undefined at 0.
It comes from the same $\max(A-\lambda E)$ optimization as Katul's model, solved analytically under a
linearized $A(C_i)$, which gives a closed form. Here $g_1 \propto \lambda^{-1/2}$, a relation the
water-stress limbs use (§4). Parameters: `stomatal_g0`, `stomatal_g1`. Code: `stomata_gs_medlyn`.

### 3.3 Katul et al. (2010) — marginal-WUE optimality

Stomata are assumed to **maximize carbon gain net of a water cost**:

```math
\max_{g_s}\bigl(A - \lambda E\bigr) \qquad(7)
```

where $E = g_s D$ is transpiration and $\lambda$ is the **marginal water-use efficiency**, the
"price" of water (µmol CO₂ per mol H₂O, the unit `katul_lambda25` is given in). The optimum satisfies
$\partial A/\partial E = \lambda$.
Writing both fluxes in $C_i$-space (with $A = \tfrac{g_s}{1.6}(C_s - C_i)$, so
$`g_s = 1.6\,A/(C_s-C_i)`$ and $`E = 1.6\,D\,A/(C_s-C_i)`$) and setting
$`dA/dC_i = \lambda\,dE/dC_i`$ rearranges to the residual root-found on $C_i$:

```math
A'\,(C_s - C_i)^2 \;=\; 1.6\,D\,\lambda\,\bigl(A'(C_s - C_i) + A\bigr), \qquad A' \equiv \frac{dA}{dC_i} \qquad(8)
```

$A'$ is taken by central difference, because the co-limited $A(C_i)$ has no clean analytic slope. Unlike
the explicit models, $g_s$ is recovered after solving, as $`g_s = 1.6\,A/(C_s-C_i)`$; if that would
fall below the cuticular floor $g_0$, the solver solves again with $g_s$ pinned at $g_0$ (§5).
Parameters: `katul_lambda25` ($\lambda_{25}$) and the water-stress terms below. Code:
`residual_optimality`, `katul_lambda`.

---

## 4. Water stress — a fast and a slow control

Drought reduces carbon gain in two ways, and ED-MORE treats them separately, following the two-factor
decomposition of Sabot et al. (2022) and Zhou et al. (2013): it lowers the photosynthetic **capacity**
(biochemistry) and it narrows the **stomatal** aperture. The two act on different time scales, as in
Xu et al. (2016)'s ED2 scheme: capacity follows the leaf's own current water potential, the aperture
follows the previous day's predawn potential.

### 4.1 The fast control: capacity, on the current leaf potential — $\beta_{ns}$

A linear ramp in the **current leaf** water potential (the step-start $`\psi_{leaf}`$ that the
plant-water solver holds) scales the potential rates. It applies to all stomatal models, because
capacity downregulation is a biochemistry effect that does not depend on the aperture scheme. It scales
$V_{cmax}$, $J_{max}$ and TPU for C3, and $V_{cmax}$, the light-limited slope $A_j$ and $k_p$ for C4.
$R_d$ is not scaled.

```math
\beta_{ns} = \mathrm{clamp}\!\left(\frac{\psi_{leaf} - \psi_{close}}{\psi_{open} - \psi_{close}},\,0,\,1\right),
\qquad \{V_{cmax}, J_{max}, \mathrm{TPU}\}_{C3},\ \{V_{cmax}, A_j, k_p\}_{C4} \mathrel{*}= \beta_{ns} \qquad(9)
```

$`\psi_{open}`$ is the PFT's leaf turgor-loss point $`\psi_{tlp}=\pi_0\varepsilon/(\pi_0+\varepsilon)`$
and $`\psi_{close}`$ is twice it, both from the pressure–volume curve the hydraulics solver stores leaf
water on (each PFT's own `[hydraulics]` `leaf_pi0` and `leaf_elastic_mod`; `pft_wstress_psi_open`,
`pft_wstress_psi_close`). A PFT file may set `wstress_psi_open` and `wstress_psi_close` instead. A leaf
keeps its capacity while it holds turgor and loses it as it wilts.

Every stomatal model ties $`g_s`$ to $`A`$, so this ramp closes the stomata within the day: as the leaf
wilts, $`A`$ and with it $`g_s`$ fall continuously to the cuticular $`g_0`$. At
$`\psi_{leaf}\le\psi_{close}`$ the leaf by day has $`A_g=0`$, $`A=-R_d`$ and $`g_s=g_0`$, and it loses only
the cuticular water $`g_0 D`$. That makes $`g_0`$ the leaf's minimum conductance, a physical trait and
not a fitted intercept: a large $`g_0`$ keeps transpiration running from a wilted leaf.

### 4.2 The slow control: stomatal aperture, on the predawn potential — $\beta_s$ (Sabot 2022)

An exponential in the **predawn** water potential narrows the stomatal aperture:

```math
\beta_s = \min\!\bigl(1,\ \exp\bigl(s_{ref}\,(\psi - \psi_{onset})\bigr)\bigr) \qquad(10)
```

- **Leuning and Medlyn:** scale the slope, $`g_{1,\text{eff}} = g_1\,\beta_s`$.
- **Katul:** raise the marginal price, $`\lambda = \lambda_{25}\,\beta_s^{-e}`$ (`katul_lambda`).
  Because Medlyn gives $g_1 \propto \lambda^{-1/2}$, the exponent $e = 2$ makes the Katul and Medlyn
  responses identical for the same $\beta_s$. Leuning's $g_1$ has no such link to $\lambda$.

Parameters: `wstress_sref_stomata` ($s_{ref}$, 2 MPa⁻¹ in the shipped PFT file), `wstress_lambda_exp`
($e$) and `stomata_psi_onset` ($`\psi_{onset}`$, default half the PFT's turgor-loss point, so the gravity
head a tall tree has in wet soil does not read as drought). Above the onset $`\beta_s = 1`$: Sabot et
al. (2022, eq. 5) apply no stress while the soil is at field capacity, and $`\psi_{onset}`$ stands in
for that threshold because a tall tree's gravity head is present even in wet soil.

**What drives $`\beta_s`$.** The stomatal driver is a supply-side potential $\psi$ (`leaf_env_t%psi`;
the default 0 is well-watered), while the capacity driver is the current $\psi_{leaf}$, the
demand-side tension. Sabot writes $\psi$ as the soil or predawn potential. ED-MORE supplies the cohort's
own predawn potential, taken as the **previous day's daily-maximum leaf water potential**. That is the
quantity the field measurements behind $`s_{ref}`$ describe. It is also the plant's own integration over
its rooted profile, so it carries rooting depth, per-PFT vulnerability and the gravity head, which no
single soil layer's $`\psi`$ does. Cohorts with no value yet are seeded from the surface-layer soil
potential, the one place a soil $`\psi`$ enters this path. The kernel field is named plainly `psi`
because the supply potential equals the soil's only in wet soil: the wood–soil relaxation time runs
from seconds in wet soil to days in dry soil, so a droughted cohort does not equilibrate overnight.
This mirrors ED2's Katul `stoma_beta` ($`s_{ref}\,e`$ is ED2's `stoma_beta`, and $\psi$ is ED2's
`dmax_leaf_psi`); the differences are in §4.3.

The value is carried on the cohort block as two fields, `dmax_psi_leaf`, the published value the
kernel reads, and `dmax_psi_leaf_accum`, the running accumulator. They are a double buffer, not a
redundant pair. Within a day one is read-only and the other write-only. A single field cannot do both:
reset at the start of the day, it leaves nothing to read; never reset, its running maximum ratchets to
the least negative $`\psi`$ the run has ever seen and disables the closure. Merging them would also
change the physics, because the kernel would read the partial maximum of the current day, which early in
the morning barely differs from the instantaneous $`\psi_{leaf}`$, the tight sub-step feedback the slow
control exists to avoid. A one-field design would have to sample $`\psi_{leaf}`$ once a day at true
predawn, which adds a clock dependency to the fast loop and a sensitivity to which sub-step is caught.

### 4.3 Where this diverges from ED2

ED2's `farq_katul.f90` (Xu et al. 2016's scheme, `H2O_PLANT_LIM = 4`) has the same two time scales: a
turgor-loss capacity term on the leaf potential at the start of each step, and a Manzoni-style
`stoma_beta` on the day's running maximum potential. The forms differ:

| | ED-MORE | ED2 (`farq_katul.f90`) |
|---|---|---|
| **capacity shape** | linear ramp in $\psi_{leaf}$ from $\psi_{tlp}$ to $2\psi_{tlp}$ | $`1/\bigl(1+0.1(\psi_{leaf}/\psi_{tlp})^{6}\bigr)`$, a 6th-power turgor-loss curve |
| **capacity targets** | $V_{cmax}$, $J_{max}$ **and** TPU | $J_{max}$ and quantum yield $\alpha$ only ($V_{cmax}$ and TPU are commented out) |
| **stomatal scope** | Leuning/Medlyn $g_1$ **and** Katul $\lambda$ | Katul $\lambda$ only |
| **stomatal driver** | predawn (previous-day daily-max) leaf $\psi$ | `dmax_leaf_psi`, the current day's running maximum |
| **parameter split** | $s_{ref}$ (shared with the $g_1$ limb) $\times$ $e$ (Katul-only) | one `stoma_beta` |

**Capacity shape and targets.** Both anchor the capacity term on $\psi_{tlp}$. ED2 applies it to
$J_{max}$ and $\alpha$ and leaves $V_{cmax}$ untouched, so it acts only on the light-limited branch.
ED-MORE downregulates all three potential rates together, which is what the underlying observations
describe, and so the effect does not depend on which limitation a leaf is under. ED-MORE also applies it
for every stomatal model, not just Katul, because capacity downregulation is biochemistry. The ramp is
linear, so its response to an error in $\psi_{leaf}$ is linear too; ED2's sixth power would amplify the
same error more.

**Stomatal scope.** Restricting $\beta_s$ to Katul, as ED2 does, would leave the two explicit schemes
with no stomatal drought response. ED-MORE applies it to all three; for Medlyn the exponent $e = 2$ makes
it equivalent to Katul's, because Medlyn gives $`g_1 \propto \lambda^{-1/2}`$ (§4.2).

**Parameter split.** ED-MORE's $`s_{ref}\,e`$ is ED2's `stoma_beta`, the same exponential in the same
potential, so a value calibrated for one transfers to the other as
$`\texttt{stoma\_beta} = -s_{ref}\,e`$. For Katul alone only the product matters, so the split is
degenerate there; it is useful because $s_{ref}$ alone drives the Leuning and Medlyn limb, which ED2
has no equivalent of.

### 4.4 No separate shut-off

$`\beta_s`$ scales $`g_1`$ only, so on its own it leaves a drought-stressed leaf the conductance
$`g_0 + 1.6\,A/C_s`$ (Medlyn's "1 +" term) for as long as it still photosynthesises. The fast control
removes the $`A`$-driven part as the leaf wilts, so no further factor on $`g_s`$ is needed. What is left
is $`g_0`$, the cuticular loss, which a real leaf has too. At a VPD of 3.5 kPa the 0.01 mol m⁻² s⁻¹ of
the shipped PFT file loses 0.35 mmol m⁻² s⁻¹, whereas a measured minimum conductance of
0.3 mmol m⁻² s⁻¹ loses about 0.01 mmol m⁻² s⁻¹. That is why $`g_0`$ should be set to a measured minimum
conductance.

**The substomatal air is taken as saturated.** It is in equilibrium with leaf water at $`\psi_{leaf}`$,
so strictly its relative humidity is the Kelvin value $`\mathrm{RH}=\exp\!\big(\psi/(\rho_w R_v T)\big)`$
and not 1. Transpiration on $`e_i=\mathrm{RH}\cdot e_{sat}`$ would shrink the vapour gradient
continuously as the leaf dries and reverse it once $`e_i<e_a`$ (foliar water uptake), which would
arrest transpiration with no threshold parameter, at $`\psi=\rho_w R_v T\ln(\mathrm{RH}_{air})`$. ED-MORE
applies this relation to soil evaporation (`ground_evaporation_rate`, `alpha_soil`) but not to leaves.

The leaf humidity cannot simply scale the latent flux. A tissue's latent flux is linearised about the
canopy-air temperature as $`\mathrm{ref}+\mathrm{slope}\,(T-T_{cas})`$ (`latent_pathway_t`,
[vegetation_energy_dynamics](vegetation_energy_dynamics.md)), so scaling both terms by RH moves the
base point and the slope inconsistently, and on a nearly turgid leaf it gives condensation at midday.
The fault is the linearisation, not the magnitude. Doing it properly means rebuilding the flux and its
temperature derivative from $`\mathrm{RH}\cdot q_{sat}(T_{leaf})`$ inside the leaf energy balance. Its
likely value is as a route to foliar water uptake, not as a drought response ([ROADMAP](../ROADMAP.md)).

The Kelvin correction is kept out of the stomatal model on purpose: the Medlyn and Leuning $`g_1`$
were calibrated against a VPD based on $`e_{sat}`$, so feeding them a corrected VPD would silently
re-tune $`g_1`$.

---

## 5. The coupled solver

`solve_leaf_gas_exchange` brackets $`C_i \in (\Gamma^*, C_{i,open}]`$ and finds the root of the residual
by the Illinois method (the module's private `root_ci`, on the Illinois steps of `ed_numerics`):
regula falsi on the bracket, which halves the kept end's residual whenever the same end is kept twice,
so both ends close in. It needs far fewer residual evaluations than bisection and stops when the
bracket is narrower than the tolerance `CI_TOL_PPM` (10⁻³ µmol mol⁻¹):

- **Night / closed branch** — if the best-case net rate $A(C_a) \le 0$, stomata sit at $g_0$ and
  the solve returns immediately.
- **The bracket's top** — $`C_{i,open}`$ is $C_i$ with the stomata fully open, where $`C_i = C_s`$:
  $C_a$ without a boundary layer, and with one the root of $`C_i = C_a - 1.4\,A(C_i)/g_b`$, found by
  the same root finder (`GS_OPEN`). Katul's residual, eq. (8), is $`-1.6\,D\lambda A < 0`$ there, so it
  brackets; at $C_a$ the boundary layer can leave it positive, which would send a well-lit leaf to the
  cuticular fallback. The explicit models bracket to $C_a$, which holds their root as well.
- **Two residual providers**, chosen by `ci_residual` — `residual_explicit_gs` (Leuning/Medlyn
  diffusion identity, eq. (1), and the g0-pinned fallback with $g_s := g_0$) and
  `residual_optimality` (Katul, eq. (8)). Both are `pure`.
- **Cuticular fallback** (`GS_CUTICULAR`) — if the model yields no sign change (for example Katul
  under strong stress wanting $g_s < g_0$), the solver retries with $g_s$ pinned at $g_0$, which
  always brackets when $A(C_a) > 0$.
- **Guards** — $`g_s = g_0`$ where $`A \le 0`$; Leuning floors $`C_s-\Gamma^*`$ and Medlyn $`C_s`$ just
  above zero, and Medlyn's and Katul's VPD at `vpd_min` (at $D = 0$ water costs nothing, so Katul's
  optimum would open the stomata without bound); Katul's $\lambda$ takes $`\beta_s`$ no lower than
  `BETA_FLOOR` (10⁻⁴), so it stays finite; the lower bracket sits `LO_EPS_PPM` (10⁻³ µmol mol⁻¹) above
  $`\Gamma^*`$; Katul's derivative in $C_i$ is a finite difference of 10⁻³ $C_i$, at least 10⁻² µmol mol⁻¹.

Then $g_s$, $C_s$, $C_i$ and the transpiration $E$ are assembled consistently. $E$ uses the vapour
deficit over the pressure and, with the boundary layer resolved, $g_s$ and $g_b$ in series.

---

## Parameters

Per-PFT traits are in the PFT file (`[pft]`, one value per PFT); the shared biochemistry and the model
selectors are in `[leaf_physiology]`.

| Symbol | Config key | Meaning |
|---|---|---|
| $V_{cmax,25}$, $J_{max}/V_{cmax}$, $\mathrm{TPU}/V_{cmax}$, $R_{d,25}/V_{cmax}$ | `[pft]` `vcmax25`, `jmax_vcmax_ratio`, `tpu_vcmax_ratio`, `rd_vcmax_ratio` | capacities at 25 °C and their ratios (a cohort's plastic $V_{cmax,25}$ and $R_{d,25}$ replace the table's) |
| $k_{p,25}$, $\phi_{C4}$ | `[pft]` `kp25`, `quantum_yield_c4` | C4 PEP-carboxylase slope and light-limited quantum yield |
| $\theta_J$ | `[pft]` `theta_j` | curvature of the electron-transport hyperbola only |
| co-limitation curvatures | `[pft]` `theta_cj_c3`, `theta_ip_c3`, `theta_cj_c4`, `theta_ic_c4` | the two transitions of each pathway (§2) |
| $g_0$ | `[pft]` `stomatal_g0` | cuticular (minimum) conductance |
| $g_1$ | `[pft]` `stomatal_g1` | stomatal slope (Leuning dimensionless, Medlyn kPa⁰·⁵) |
| $D_0$ | `[pft]` `stomatal_d0` | Leuning humidity sensitivity [Pa] |
| $\lambda_{25}$ | `[pft]` `katul_lambda25` | Katul marginal water-use efficiency (well-watered) |
| $\psi_{open},\psi_{close}$ | `[pft]` `wstress_psi_open`, `wstress_psi_close` | fast-limb ramp bounds (optional; default $`\psi_{tlp}`$, $`2\psi_{tlp}`$) |
| $\psi_{onset}$ | `[pft]` `stomata_psi_onset` | where the stomatal limb starts (optional; default $`\psi_{tlp}/2`$) |
| $s_{ref}$ | `[pft]` `wstress_sref_stomata` | stomatal-limb sensitivity (Sabot) |
| $e$ | `[pft]` `wstress_lambda_exp` | exponent of $\beta_s$ in the Katul $\lambda$ |
| none | `[leaf_physiology]` `stomatal_model`, `temp_response_form`, `colimitation`, `use_boundary_layer` | model selectors |
| $K_c,K_o,\Gamma^*$ | `[leaf_physiology]` `kc25`, `ko25`, `gstar25`, `ea_kc`, `ea_ko`, `ea_gstar` | Rubisco kinetics at 25 °C [Pa] and their activation energies |
| none | `[leaf_physiology]` `ea_vcmax`, `ea_jmax`, `ea_rd`, `hd_*`, `ds_*` | activation, deactivation and entropy terms of the temperature responses |
| $O$ | `[leaf_physiology]` `o2_mol_frac` | atmospheric O₂ mole fraction (0.209) |
| $\alpha_{leaf}$, $\phi_{PSII}$ | `[leaf_physiology]` `leaf_absorptance`, `phi_psii` | leaf PAR absorptance (0.85); electrons per photon in low light, times 2 (0.74) |
| none | `[leaf_physiology]` `vpd_min` | least VPD the Medlyn and Katul models use (0.05 kPa) |
| none | `[leaf_physiology]` `thermal_acclimation`, `acclim_*` | thermal acclimation (§2) |

## Code map

| Concept | Routine (`ed_leaf_gas_exchange` unless noted) |
|---|---|
| Rubisco kinetics $K_c$, $K_o$, $`\Gamma^*`$ in µmol mol⁻¹, eq. (3) | `rubisco_kinetics` |
| C3 and C4 demand, eq. (2) | `assimilation_demand_c3`, `assimilation_demand_c4`, `combine_limits` |
| electron transport $J$ | `electron_transport_j` |
| temperature responses and acclimation, eq. (4) | `ed_therm_lib`: `temp_response`, `kattge_knorr_entropy`, `kattge_knorr_jv_ratio`; `ed_plant_traits`: `acclimate_leaf_photo_table` |
| Leuning, Medlyn, Katul $\lambda$, eqs. (5)–(7) | `stomata_gs_leuning`, `stomata_gs_medlyn`, `katul_lambda` |
| the two water-stress controls, eqs. (9)–(10) | `capacity_stress_factor` ($`\beta_{ns}`$), `stomatal_stress_factor` ($`\beta_s`$); `ed_pft_derive`: `pft_wstress_psi_open`, `pft_wstress_psi_close`, `pft_stomata_psi_onset` |
| coupled solve | `solve_leaf_gas_exchange`; over n leaves, `leaf_gas_exchange_batch` |
| a leaf's table entry with the cohort's capacities | `leaf_params_at_capacity` |
| a cohort's incident PAR per unit of its own leaf area | `cohort_leaf_par` |
| residuals, eqs. (1), (8) | `ci_residual`, then `residual_explicit_gs` or `residual_optimality` |
| root finder | `root_ci` (private), on `ed_numerics`' `root_bracket_t`, `regula_falsi_point` and `illinois_narrow` |

## References

- Farquhar, von Caemmerer & Berry (1980). *Planta* 149:78–90: FvCB C3 photosynthesis.
- Collatz, Ribas-Carbo & Berry (1992). *Aust. J. Plant Physiol.* 19:519–538: C4.
- Bernacchi et al. (2001). *Plant Cell Environ.* 24:253–259: Rubisco kinetics and their temperature response.
- Kattge & Knorr (2007). *Plant Cell Environ.* 30:1176–1190: thermal acclimation of photosynthetic capacity.
- Long et al. (1993) and Björkman & Demmig (1987): maximum quantum yield of C3 leaves.
- Leuning (1995). *Plant Cell Environ.* 18:339–355: semi-empirical $g_s$.
- Medlyn et al. (2011). *Glob. Change Biol.* 17:2134–2144: USO $g_s$.
- Katul, Manzoni et al. (2010). *Ann. Bot.* 105:431–442: optimization theory ($\lambda$).
- Manzoni et al. (2011). *Funct. Ecol.* 25:456–467: $\lambda$ under water stress.
- Vico, Manzoni et al. (2013). *Agric. For. Meteorol.* 182–183:191–199.
- Sabot et al. (2022). *JAMES* 14(4):e2021MS002761: two-limb ($\beta_{stomata}$ and $\beta_{nonstomata}$).
- Zhou et al. (2013). *Agric. For. Meteorol.* 182–183:204–214: combined stomatal and non-stomatal.
- Xu et al. (2016). *New Phytol.* 212:80–95: the ED2 plant-hydraulics and Katul implementation.
