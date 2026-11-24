# Ground surface layers

The ground surface layers are what lies on the soil, top down: the **water layer** (snow, ice and ponded water, one store), the **litter** and the **moss**. Each layer has its own temperature. Where a layer exists it is part of the ground's surface: it absorbs the radiation the ground receives, exchanges heat and vapour with the canopy air, holds and passes water, and insulates what lies below it. Without any layer the ground is the bare soil.

This page gives each layer's state and laws, the stage that prepares them before each fast step's radiation, their energy balance (one coupled stack over the soil top), their water (nodes of the soil-plant water network), how they come and go, how a snow pack buries plants, and what reads them. It is one of the per-store pages under [coupled_biophysics](coupled_biophysics.md); the water network they belong to is [soil_plant_water](soil_plant_water.md).

---

## 1. The layers and their state

| Layer | Type | Prognostic | Diagnosed |
|---|---|---|---|
| water: snow, ice and ponded water | `water_layer_state_t` | water $W$ [kg m⁻²], energy $E$ [J m⁻²], the ice matrix's depth $`d_m`$ [m] | temperature $T$, liquid fraction $`f_{liq}`$ |
| litter, then moss (the porous layers) | `porous_layer_state_t` | dry matter $M$ [kg m⁻²], water $W$ [kg m⁻²], energy $E$ [J m⁻²] | $T$, $`f_{liq}`$ |

Both extend `ground_store_t` (water, energy, temperature, liquid fraction), so one pair of moves serves every hand-off between them: `store_move` moves water with its energy, and `store_take_exchange` takes an exchange's energy and vapour, stopping the run if a store would give more vapour than it holds. The layers are fixed slots, water over litter over moss, held in `ground_surface_state_t`; an absent layer holds nothing.

**One threshold.** Matter is a layer from `[ground_surface_layer].min_depth` (0.1 mm) deep. A porous layer's matter is registered whole from that depth and not at all below it (`porous_layer_matter`); its carbon stays in its pool either way. A water layer thinner than that is a **film** (§2): it holds its water at the temperature of what it lies on and exchanges nothing.

**The water layer.** Its phase is read off its energy as pure water's, on the plateau at the triple point $`T_{3ple}`$ (`internal_energy_to_temp`): ice below it, liquid above it, a mixture on it. Its ice $`I = (1-f_{liq})W`$ forms a matrix of depth $`d_m`$, and its liquid $L$ fills the matrix's pores first. The ice holds liquid up to the share `liquid_holding_frac` $h$ of the layer's mass, LEAF-3's 1:9; the liquid beyond that is the **free water** $F$, the pond, which the water network moves (`water_layer_free_water`):

```math
F = \max\!\left(0,\ L - \frac{h}{1-h}\,I\right) \qquad(1)
```

The layer's depth is its matrix's, or its water's own volume where that is larger (`water_layer_depth`), so a pond, which has no matrix, is as deep as its water:

```math
d = \max\!\left(d_m,\ \frac{I}{\rho_{ice}} + \frac{L}{\rho_w}\right) \qquad(2)
```

The matrix follows its ice (`matrix_depth_update`). Snowfall $S$ adds depth at `fresh_density` $`\rho_{fresh}`$. Ice that grows beyond the old ice plus the snowfall, rain or free water that freezes, fills the pores first and deepens the matrix only once it is solid ice; ice that is lost, to melt or sublimation, takes its share of the depth, so the matrix keeps its density:

```math
d_m' = \begin{cases}
\max\!\left(d_m + S/\rho_{fresh},\ I'/\rho_{ice}\right) & I' \ge I + S \\
\left(d_m + S/\rho_{fresh}\right) I'/(I + S) & I' < I + S
\end{cases} \qquad(3)
```

with $I$ the ice before and $I'$ after. Rain therefore adds mass and no depth. There is no compaction.

**The porous layers.** A porous layer's depth is its matter over its material's bulk density, $`d = M/\rho_b`$, and it holds $`W_{max} = \theta_s\,\rho_w\,d`$ saturated, $`\theta_s`$ its retention curve's water at saturation; its fill is $`W/W_{max}`$ (`porous_layer_depth`, `porous_layer_capacity`, `porous_layer_fill`). Its water freezes on its material's freezing curve, as a soil layer's does ([soil_biophysics](soil_biophysics.md), "Frozen soil"; `porous_layer_temp`, `porous_energy_state`): at temperature $T$ its liquid is the curve's water at the head $`L_f(T - T_{3ple})/(gT)`$, and its energy is

```math
E = C_d T + W\left[f_{liq}\,u_{liq}(T) + (1 - f_{liq})\,c_{ice} T\right], \qquad C_d = c_{org} M \qquad(4)
```

with $`c_{org}`$ = 1926 J kg⁻¹ K⁻¹ for dry organic matter (Farouki 1986). Water above saturation freezes at $`T_{3ple}`$, on a plateau.

## 2. The stage before the radiation

At each fast step's start, before the radiation, `ground_stage` (called first by `build_frozen`) prepares the layers, so that the two-stream's ground and the ground's balance see one cover of the water layer:

1. The ledgers take the layers' water and energy before anything moves (`ground_surface_store`).
2. The canopy intercepts rain only (`intercept_rain`); snowfall passes to the ground.
3. Each porous layer's water above saturation joins the water layer, with the energy that keeps the porous layer at its temperature (`porous_layers_spill`).
4. The water layer takes the throughfall's rain and snow and the water the tissues shed (`water_layer_stage`), each with its enthalpy: snow as ice at $`\min(T_{3ple}, T_{air})`$, rain as liquid at $`T_{air}`$, the air temperature at the canopy-air top (`precip_enthalpy`). Its matrix follows, eq. (3).
5. A water layer thinner than `min_depth` is a film. On a porous layer, the film and that layer come to one temperature, the layer on its freezing curve (`film_on_porous_layer`); on the soil, the film takes the top soil node's temperature and phase. The heat that moves is taken from that layer or that node. A film holds no liquid in its ice: all its liquid is free water.
6. The free water, eq. (1), is set apart at the layer's temperature: it becomes the water network's inflow to the pond for the step. A layer all liquid, a pond, sets all of its water apart.
7. The cover, the water layer and its free water together (`water_layer_exposed`, `water_layer_cover_of`), is the step's.

The radiation then runs through the canopy onto the ground the stage leaves (`canopy_ground_radiation`).

**Cover.** The water layer covers the fraction of Niu & Yang (2007), over its bulk density $`\rho = W/d`$ (`water_layer_cover`):

```math
f = \tanh\!\left(\frac{d}{s}\right), \qquad s = a\,z_{0g}\left(\frac{\rho}{\rho_{ref}}\right)^{m} \qquad(5)
```

with $a$ = `cover_scale_factor` (2.5), $`z_{0g}`$ = `[aerodynamics].soil_rough` (0.01 m), $`\rho_{ref}`$ = 100 kg m⁻³ and $m$ = `cover_density_exponent` (1). Fresh snow's scale is 6.25 cm: 1 cm covers 16 % of the ground and 10 cm covers 92 %. A film covers nothing. The same law covers a pond, at the density of water: a 5 mm pond covers 2 %.

**Optics.** The layer's albedo in each band weighs three surfaces (`water_layer_albedo`, `water_layer_optics_shares`): the water standing above the matrix, the share $`s = \min(1, \max(0, I/\rho_{ice} + L/\rho_w - d_m)/d)`$ of the depth, at standing water's albedo; and the matrix, at snow's albedo over its porous share $`p = (\rho_{ice} - \rho_m)/(\rho_{ice} - \rho_{fresh})`$, clipped to [0, 1], and at solid ice's over the rest, $`\rho_m = I/d_m`$ the matrix's density:

```math
\alpha = s\,\alpha_{liq} + (1 - s)\left[p\,\alpha_{snow} + (1 - p)\,\alpha_{ice}\right], \qquad
\alpha_{snow} = (1 - f_{liq})\,\alpha_{fresh} + f_{liq}\,\alpha_{aged} \qquad(6)
```

Fresh snow has $s = 0$ and $p = 1$; a pond has $s = 1$; a pond frozen solid has $p = 0$. Snow's albedo ages through its liquid fraction, not with time.

**Conduction.** The matrix conducts at Calonne et al. (2011)'s regression of snow's effective conductivity on density (as Brondex et al. 2023 give it, their eq. A3, fitted at 271 K over about 100–550 kg m⁻³), used as a shape scaled by `conductivity_at_fresh` $`k_{fresh}`$, so fresh snow keeps $`k_{fresh}`$ (`snow_matrix_conductivity`); the water standing above the matrix conducts as still liquid water, $`k_w`$ = 0.57 W m⁻¹ K⁻¹ (`water_layer_resistance`):

```math
k_m = k_{fresh}\,\frac{k_C(\rho_m)}{k_C(\rho_{fresh})}, \quad
k_C(\rho) = 2.5\times10^{-6}\rho^2 - 1.23\times10^{-4}\rho + 0.024, \quad
R_w = \frac{d_m}{k_m} + \frac{d - d_m}{k_w} \qquad(7)
```

Above its fitted range, towards solid ice, the law is used as that shape beyond its range: solid ice conducts about 13 times more than fresh snow.

A porous layer holds its water in its lower part and dries from its top (§4), so its heat crosses the dried top at the material's dry conductivity, in series with the wet part at its wet one (`porous_layer_resistance`):

```math
R_p = d\left[\frac{1 - W/W_{max}}{k_{dry}} + \frac{W/W_{max}}{k_{wet}}\right] \qquad(8)
```

3 cm of dried moss over 3 cm of wet moss conduct 0.11 W m⁻¹ K⁻¹.

## 3. Energy: one stack over the soil top

Every evaluation of the march solves the layers' energy together (`ground_surface_exchange`, at every stage of either integrator), against the canopy air and the top soil node as the march holds them there. The **nodes** are the water layer, where it holds water and covers some ground, and each porous layer present, top down; the top soil node is the stack's held lower boundary.

**The top-face convention.** A node's one temperature is its top face's: it emits, and exchanges heat and vapour with the canopy air, at that temperature, and its heat to the node below crosses its whole depth. A pack's cold content is that of a pack as cold as its surface; only a layered pack would separate the two. The conductance between node $i$ and the node below it is $`c_i/R_i`$, and between the bottom node and the top soil node $`c_n/(R_n + r_{soil})`$, with $R$ the layer's resistance, eq. (7) or (8), $`r_{soil} = |z_1|/k_{soil,1}`$ the top soil node's half layer at its conductivity at the step's start, and $c$ the contact: the water layer's cover for the water layer, 1 for a porous layer. The water layer thus passes heat only over its cover, and a deepening pack decouples the soil.

**Each node's balance.** Over the step, node $i$'s energy $`e_i`$ satisfies a backward-Euler balance (`stack_residual`):

```math
e_i - e_i^0 = \Delta t\left[A_i\left(R_{sw,i} + \varepsilon_i\left(L^{\downarrow} - \sigma T_i^4\right) - H_i - \lambda_i E_i\right) - \Lambda_i + Q_i\right] \qquad(9)
```

- $`A_i`$ is the node's open area: the water layer's cover; the top porous layer's uncovered share, $1-f$; 0 for a porous layer under another.
- $`R_{sw,i}`$ is the shortwave it absorbs at its own optics, and its longwave is formed at its own temperature (§3, "The ground's radiation").
- $`H_i = g_a\rho c_p(T_i - T_{cas})`$ is its sensible heat through the ground's one aerodynamic conductance $`g_a`$, the leg its vapour takes too (§4).
- $`\lambda_i E_i`$ is the water layer's vapour, live in the march: $`E = g_a\rho\,[q_{sat}(T, f_{liq}) - q_{cas}]`$ over its cover, saturated over its own mixture of ice and liquid, eq. (11), and leaving at the layer's state. A vapour that would take more than the layer holds is held at what empties it, and the balance solved again, so a pack that sublimates whole leaves with all its energy.
- $`\Lambda_i`$ is a porous layer's vapour, held at the water network's rate (§4), which takes its whole enthalpy $`h_{vap}(T_i)`$ from the layer's energy at the layer's temperature. The layer's phase law then reads what is left, so the balance does not depend on which phase gave the vapour.
- $`Q_i`$ is conduction to and from its neighbours, and the heat of the water that moves through the stack (`link_heat`, `advect_through`): water passing down leaves a node at that node's temperature and arrives at the next at the giver's, and water rising from the soil arrives at the top soil node's.

Every term is linearised about the node's start temperature, as each tissue's is. The water layer's energy is on the pure-water plateau (`water_phase`), and a porous node's on its freezing curve, eq. (4) (`porous_energy_state`). **Newton** on the node energies, with the tridiagonal Jacobian of eq. (9), solves the stack (`surface_stack_balance`); should it not converge, Gauss–Seidel sweeps solve each node in turn with its neighbours held (`surface_node_balance`: on a plateau, the ice line or the liquid line in closed form; on a freezing curve by a bracketed Newton, `curve_node_root`). Each balance is implicit in its own temperature, so a node of little heat capacity relaxes towards its balance and cannot overshoot it. A run's log reports the stack's Newton iterations and fallbacks at its end (`surface stack`).

**The soil's top face.** The top soil node receives the bottom node's conduction, the heat of the water crossing into it at the temperature it leaves the stack at, less what rises from it, and less its own evaporation's latent heat. Where no porous layer lies, the bare soil exchanges over the uncovered ground itself: its net radiation at the step's start and its sensible heat.

**Commit.** After the march each layer takes the time integral of its exchange (`water_layer_commit`, `store_take_exchange`): its energy and the vapour it lost, the water layer also the free water the network returns and the canopy air's condensate. The layers, the canopy air and the soil close together in the whole-column ledgers.

### The ground's radiation

The two-stream sees one ground: the uncovered ground (bare soil, or the top porous layer, at the soil's optics, `[soil].ground_albedo_*` and `ground_emissivity`) and the water layer, blended by its cover $f$ (`snow_blend_ground_optics`). Its albedo and emissivity are the cover-weighted means, and it emits at the temperature $`T_{eff}`$ at which the blended emissivity $`\bar\varepsilon`$ emits what the two emit at their own temperatures, the uncovered ground at its top's $`T_g`$ (`ground_top_temp`: the top porous layer's, else the top soil node's):

```math
\bar\varepsilon\,\sigma T_{eff}^4 = f\,\varepsilon_w\,\sigma T_w^4 + (1 - f)\,\varepsilon_g\,\sigma T_g^4 \qquad(10)
```

What reaches the ground in each band is then shared at each surface's own optics (`ground_tile_radiation`), so the canopy and the ground agree on the ground's radiation. The shortwave is split by albedo, not by area: on half-covered ground the snow's half absorbs far less than the uncovered half.

### Vapour saturates over ice

The water layer's saturation humidity is taken over its own mixture of ice and liquid, through the saturation vapour pressure

```math
e_{sat}(T,f_{liq}) = f_{liq}\,e_l(T) + (1-f_{liq})\,e_i(T), \qquad
e_p(T) = e_{3ple}\left(\frac{T_{3ple}}{T}\right)^{(c_p-c_{p,v})/R_v}
         \exp\left[\frac{L_p(T_{3ple})}{R_v\,T_{3ple}} - \frac{L_p(T)}{R_v\,T}\right] \qquad(11)
```

with $p$ the phase: liquid ($`c_p=c_{p,l}`$, $`L_p=L_v`$) or ice ($`c_p=c_{p,i}`$, $`L_p=L_s`$), $`L_v(T)`$ and $`L_s(T)`$ the latent heats of [canopy air space biophysics](canopy_air_space_biophysics.md) (`latent_heat_vap`), $`R_v`$ the vapour's gas constant and $`e_{3ple}`$ the saturation vapour pressure at the triple point (`es_3ple`, 611.655 Pa; Ambaum 2020, eq. 14). This is Ambaum's (2020) closed form of the Clausius–Clapeyron equation for a latent heat that is linear in temperature, so the slope is **by definition**

```math
\frac{\partial e_p}{\partial T} = \frac{L_p(T)\,e_p(T)}{R_v\,T^2} \qquad(12)
```

(`sat_vapor_pressure_temp_deriv`), with the same $`L_p`$ the energy budget uses: the second law holds between the vapour pressure and the heat a flux of vapour carries. An empirical fit such as Bolton's or Magnus's has no such link. Ice's latent heat at $`T_{3ple}`$ is the liquid's plus the heat of fusion, so both branches equal $`e_{3ple}`$ **exactly** at $`T_{3ple}`$, where a melting pack sits, and the blend is continuous in temperature and in $`f_{liq}`$: a pack that freezes or melts slides between the curves rather than stepping. The slope of the blend is taken at fixed $`f_{liq}`$.

$`f_{liq}`$ is the layer's own liquid fraction, not a temperature threshold. Its energy is ice-aware (vapour that leaves an ice-referenced store takes vaporisation plus fusion), so its vapour must be too. Over ice $`e_{sat}`$ is 9 % lower at −10 °C, 18 % at −20 °C and 26 % at −30 °C, so the liquid curve would overstate it, and with it the driving gradient, by 10 %, 22 % and 34 %. With the model's own heat capacities and latent heat, the ice form is within 0.01 % of Murphy & Koop (2005) at −10 °C and 0.1 % at −30 °C, and the liquid form within 0.04 % at −10 °C, 0.09 % at 20 °C and 0.24 % at 40 °C.

A frozen porous layer or soil layer does not use the blend: its liquid sits on its freezing curve, whose potential puts the liquid's vapour pressure at the ice's, so Philip's humidity on that potential over liquid saturation already carries the ice ([soil_biophysics](soil_biophysics.md), "Frozen soil"). The blend is not used for dewpoint conversion or a diagnostic vapour-pressure deficit: dewpoint is defined over liquid.

## 4. Water

**Nodes of the water network.** The water layer's free water is the network's top node, the pond, and each porous layer present is a node between the pond and the top soil layer, top down ([soil_plant_water](soil_plant_water.md) §1). The network moves water only; the heat it carries is the stack's (§3). A porous layer is a soil layer to the solver, on its material's retention curve and saturated conductivity at its depth, rebuilt each step (`porous_layer_curve`), so the litter and the moss fill, drain to the soil and are drawn on by a drying soil as a soil layer is. Its ice is held over the step as a soil layer's is (`soil_ice_state`): it adds its water, impedes the flow by $`10^{-6\,\theta_{ice}/\theta_s}`$ (CLM5) and leaves its pores to the liquid. A layer holding less than its residual carries what it lacks as an offset over the step, one rule for porous and soil layers. The roots take nothing from the porous layers. A porous layer pushed past saturation stores through the specific storage and pushes water up to the pond, and what it holds past saturation joins the water layer at the next stage (§2).

**Faces.** Every face below the pond takes the soil's thickness-weighted geometric mean of its two conductivities. The pond's face onto a porous layer wets: water flowing down takes the pond's saturated conductivity, water flowing up the layer's own ([soil_plant_water](soil_plant_water.md) §1).

**Evaporation: one law, one leg.** Each porous layer and the top soil layer evaporate by the soil's law: Philip's pore humidity $`\alpha(\psi) = \exp(\psi g/R_vT)`$ at the store's potential, implicit in the network, with its supply limit at residual. Each porous layer holds its water in its lower part and dries from its top, as the soil's top layer does under its dry surface layer (Swenson & Lawrence 2014), and by the same law (`dry_layer_resistance`, `slot_evaporation_setup`): a dried depth that is the whole layer when it holds nothing, crossed at water vapour's diffusivity in air $`D_v`$ times the Millington–Quirk tortuosity of its open pores, $`\theta_{s}^{4/3}`$. Layer $i$'s path to the ground's surface node is

```math
g_i = \frac{1}{r_{dry,i} + \sum_{j<i} r_{cover,j}}, \qquad
r_{dry,i} = \frac{d_i\,(\theta_{s,i} - \theta_i)/\theta_{s,i}}{D_v\,\theta_{s,i}^{4/3}}, \qquad
r_{cover,j} = \frac{d_j\,\theta_{s,j}^2}{D_v\,\epsilon_j^{10/3}} \qquad(13)
```

with $`r_{cover,j}`$ the open pores $`\epsilon_j = \theta_{s,j} - \theta_j`$ of each layer above it, so the moss evaporates through the litter's open pores; the top soil layer's vapour crosses every porous layer's open pores, in series with its own dry surface layer. A saturated layer has no open pores and passes no vapour. The paths of every store meet at the ground's surface node in parallel, and the node passes their sum to the canopy air through the ground's one aerodynamic leg $`1/g_a`$, the leg its sensible heat takes (`ground_vapour`). The node's humidity is in closed form, $`q_s = (g_a q_{cas} + \sum_k g_k q_k)/(g_a + \sum_k g_k)`$ over the stores that evaporate into it, and store $i$ evaporates

```math
E_i = (1 - f)\,\rho\,g_i\,\left[q_i - q_s\right]^+ \qquad(14)
```

over the uncovered ground, where the water layer leaves it. A store drier than the node takes nothing, from the air or from the other stores: there is no dew into a porous layer. Its supply limit tapers its path to nothing at residual (`evaporation_supply`). The layer's dried depth is its own, set by its water, so it is implicit in the network. 3 cm of dry moss add about 1,700 s m⁻¹ at 280 K, many times the ground's aerodynamic resistance. A drying litter comes to equilibrium with the canopy air: at 30 % relative humidity near −1.7×10⁴ m and θ ≈ 0.012.

**Condensate.** Vapour the canopy air condenses onto the ground joins the water layer, the ground's top water store.

## 5. How the layers come and go

**The porous layers follow their matter** (`follow_porous_sources`, `follow_porous_layers`), at a run's start, at a restart and at every slow step, by each material's `mass_source`:

- `metabolic_litter`: $`M = 2\,C_{fast,grnd}`$, the above-ground metabolic litter pool's carbon over half of dry organic matter being carbon (`c2b`). The structural pool holds the coarse woody debris too, which spread over a gap patch would make a layer tens of centimetres deep, so it is left out.
- `site_depth`: a site's `[soil].moss_depth` at the material's bulk density. The moss is prescribed: it neither grows nor photosynthesises.
- `none`: no layer.

Matter that arrives or leaves moves its dry heat at the layer's temperature (the top soil layer's, for a layer it starts), a declared term of the slow ledger. Arriving matter brings the water its own volume holds at the matric head it and the top soil layer share once that water has left the soil layer, drawn from the top soil layer with that water's energy there, at the soil's temperature and liquid fraction: a matched move, which the slow ledger needs no declaration for. A large fall of litter on dry soil so takes what the soil can give at a common head, and never draws the soil below its air-dry head (CLM5's −10⁴ m, `PSI_AIR_DRY`). At a run's start the water comes from outside, at the top soil layer's head, temperature and liquid fraction. Water above what a layer's new depth holds saturated, and all of a layer's water and energy when its matter is gone, joins the water layer. A restart rebuilds a layer its state file lacks and removes one its `mass_source` does not give.

**The water layer** starts a run at `initial_mass` of ice at `initial_temp`, its matrix at `fresh_density` (`seed_snow`). A state file that holds the pack and the pond apart, or the litter and moss as one record, is read into the layers: the pond joins the water layer, whose temperature is read off the sum.

**Mixing two patches.** Fusion and gap formation area-weight each layer's water, energy, matter and matrix depth (`blend_ground_surface`), so the site's totals are kept; the temperatures are read again from the blended energies.

## 6. Burial by snow

A snow pack buries the crowns it reaches (`gather_patch`, `crown_exposed_fraction`). A cohort's crown spans from its base, `[allometry].crown_base_fraction` of its top, to its top, its height bent down by the pack by the share `snow_bending`; the crown's share above the matrix's depth $`d_m`$ is exposed:

```math
x = \begin{cases} 1 & d_m \le z_b \\ (z_t - d_m)/(z_t - z_b) & z_b < d_m < z_t \\ 0 & d_m \ge z_t \end{cases}, \qquad
z_t = H\,(1 - b), \quad z_b = c\,z_t \qquad(15)
```

with $H$ the height, $b$ = `snow_bending` and $c$ = `crown_base_fraction`. A pond, which has no matrix, buries nothing. The exposure holds over each slow step's fast loop.

- **What the canopy sees.** Each cohort's leaf and wood area take its exposed share: its radiation, gas exchange and transpiration, its canopy-air exchange and its interception. While a pack lies, the canopy stands at its bent tops, $`\max_j H_j(1 - b_j)`$, for the roughness ([canopy_aerodynamics](canopy_aerodynamics.md) §1).
- **Respiration.** Leaf and stem respiration run on the whole tissue, buried or not.
- **Buried tissue's heat.** A cohort's buried share joins the water layer's node: its heat capacity and energy join the node's, and after the step it takes the node's temperature, the heat it gains or loses a matched pair with the water layer's energy (`buried_tissue_store`, `commit_tissue_energy`). The heat of the water that reaches a cohort splits between its exposed and buried shares the same way (`tissue_water_heat`).
- **Phenology.** A wholly buried cohort's cues read no light, and its buried share reads the ground top's temperature ([plant_phenology](plant_phenology.md)).

## 7. What reads the layers

- **The aerodynamics:** the under-canopy stability reads the ground top's temperature, the top porous layer's, else the top soil node's (`ground_top_temp`).
- **Decomposition:** the above-ground litter pools (metabolic and structural, both lying on the ground) take the ground top's temperature and, as their moisture, the porous layers' saturation above residual, $`\sum(W - W_{res})/\sum(W_{max} - W_{res})`$ over the layers, as the soil's moisture is read, instead of the top soil layer's ([soil_carbon](soil_carbon.md)).
- **The ledgers:** the layers' water and energy are stores of the whole-column budgets and of the slow ledger (`ground_surface_store`). The water they exchange with the soil and each other is a transfer between stores with its enthalpy; their vapour leaves through the canopy air.
- **The outputs:**

| Output | Meaning |
|---|---|
| `mean_water_layer_mass_py`, `_pa` | the water layer's water: snow, ice and the liquid they hold, the pond included |
| `mean_water_layer_ice_py` | its ice |
| `mean_water_layer_depth_py` | its depth, eq. (2) |
| `mean_water_layer_cover_py` | its cover, eq. (5) |
| `mean_free_water_py`, `_pa` | the pond the water network ends the step with |
| `mean_litter_water_py`, `mean_moss_water_py` (and `_pa`) | each porous layer's water |
| `mean_ground_top_temp_py`, `_pa` | the ground top's temperature: the top porous layer's, else the top soil node's |
| `mean_ice_melt_py` | the water layer's net ice to liquid |
| `mean_water_layer_evap_py` | the water layer's vapour: sublimation and the pond's evaporation, deposition negative |
| `mean_infiltration_py` | the pond's flow into the ground: into the top porous layer where there is one |
| `mean_ground_evap_py` | the vapour from the water layer, the litter, the moss and the soil |

## 8. Limitations

- **One bulk water layer.** Snow is a single layer, with no compaction and no ageing of its albedo with time: snowfall adds depth at a settled pack's density, not new snow's. It cannot hold an ice lid over liquid, only their mixture on the plateau, and still water's conductivity understates a convecting pond's heat transport.
- **The pond's contact area.** The free water infiltrates and runs off over the whole patch, while its energy exchanges and conducts over its cover, eq. (5): a 5 mm pond covers 2 % of the ground.
- **Litter as a full sheet.** A litter layer covers the whole ground from 0.1 mm (4.5 g m⁻²) and switches on and off at that depth. Real thin litter is patchy; a cover law by mass would follow its area.
- **The litter's hydraulics.** The litter is a continuum porous layer on the retention and conductivity curves the config gives (by default Letts et al. 2000's surface moss, as CLASS takes it), so a thin litter over a wet layer below draws water up and evaporates near the rate the air allows. Forest-floor measurements support continuum flow through the forest floor as a whole (Kosugi et al. 2001). Models of the litter layer itself instead exclude capillary rise into it (Napoly et al. 2017), and the evidence is mixed (Matthews 2014).
- **Each porous layer dries from its top.** A Sphagnum moss draws water up to its tips and keeps its top wetter, so its dried depth would be thinner.
- **The porous layers' optics are the soil's**, albedo and emissivity alike; the soil's emissivity is CLM5's 0.96, and no measured emissivity of litter or moss is used.
- **What reads the ground's temperature** (the stability, decomposition, a buried cohort's phenology) reads the top porous layer's or the soil's, never the water layer's.
- **Held rates.** The porous layers' and the soil's evaporation are the water network's, held over the march at the step's start humidity; the next step's network restores the response. Burial holds over each slow step's fast loop.

## Parameters

The `[ground_surface_layer]` block, every key optional ([configuration](../configuration.md)): `min_depth` (0.1 mm), the depth from which matter is a layer, and one table per material.

`[ground_surface_layer.water]`:

| Symbol | Key | Default | Meaning | Source |
|---|---|---|---|---|
| $`\rho_{fresh}`$ | `fresh_density` | 250 kg m⁻³ | the bulk density snowfall adds depth at | a settled pack's, since the layer does not compact: CLM's bulk density of snow, its initial estimate of a pack's (Lawrence et al. 2018, §2.3.2; Oleson et al. 2013, §2.2.2); new snow is lighter, 50–169 kg m⁻³ (Oleson et al. 2013, eq. 7.21, after Anderson 1976) |
| $`k_{fresh}`$ | `conductivity_at_fresh` | 0.15 W m⁻¹ K⁻¹ | the matrix's conductivity at `fresh_density`, scaled by eq. (7) | Calonne et al. (2011)'s regression at `fresh_density` (0.1495) |
| | `albedo_vis_fresh`, `albedo_nir_fresh` | 0.95, 0.65 | fresh snow's albedo | BATS's new snow at a zenith angle below 60° (Dickinson et al. 1993, §2e, eq. 5) |
| | `albedo_vis_aged`, `albedo_nir_aged` | 0.518, 0.435 | wet (aged) snow's albedo | ED2's (`ed_params`), after Roesch et al. (2002) |
| $`\alpha_{ice}`$ | `albedo_vis_ice`, `albedo_nir_ice` | 0.60, 0.40 | solid ice's albedo | CLM5's lake ice (Lawrence et al. 2018, eqs. 2.3.61–2.3.62) |
| $`\alpha_{liq}`$ | `albedo_vis_liquid`, `albedo_nir_liquid` | 0.10, 0.10 | standing water's albedo | CLM5's diffuse open lake water (the same, eq. 2.3.60) |
| $`\varepsilon_w`$ | `emissivity` | 0.97 | longwave emissivity | CLM5's snow (the same, §2.4), as ED2's |
| $h$ | `liquid_holding_frac` | 0.10 | the layer's mass share its ice holds as liquid | LEAF-3's 1:9 (Walko et al. 2000), as ED2 percolates |
| $a$, $m$ | `cover_scale_factor`, `cover_density_exponent` | 2.5, 1.0 | the cover's scale over the ground's roughness, and its growth with density, eq. (5) | Niu & Yang (2007), eq. 4, with the exponent they give for CLM (Oleson et al. 2010, eq. 3.50) |
| | `pond_capacity` | 5 kg m⁻² | the free water above which it runs off | no verified source: CLM's surface water has no fixed capacity, its runoff threshold set by the slope (Oleson et al. 2013, §7.3.2), and CLM4.5's 10 kg m⁻² is the top soil layer's allowance above saturation before drainage (the same, §7.6), another store |
| | `initial_mass`, `initial_temp` | 0 kg m⁻², 270 K | the layer at a fresh start | |

`[ground_surface_layer.litter]` and `[ground_surface_layer.moss]`. Each material's hydraulics are a parameter set the config gives: a retention curve, van Genuchten's or Campbell's, and a saturated conductivity, as a soil layer's.

| Key | Litter | Moss | Meaning | Source |
|---|---|---|---|---|
| `mass_source` | `metabolic_litter` | `site_depth` | where the dry matter comes from: `none`, `metabolic_litter` or `site_depth` | |
| `bulk_density` | 45 kg m⁻³ | 43 kg m⁻³ | dry mass per volume | litter: Napoly et al. (2017); moss: moist acidic tussock tundra, the organic layer's top 5 cm (McLaren 2021) |
| `k_dry`, `k_wet` | 0.10, 0.1294 W m⁻¹ K⁻¹ | 0.06, 0.50 W m⁻¹ K⁻¹ | conductivity dry and saturated, eq. (8) | litter: de Vries (1963)'s 0.1 + 0.03 θ, via Napoly et al. (2017), at θ = 0 and `theta_sat`; moss: Chadburn et al. (2015) |
| `retention` | `campbell` | `van_genuchten` | the retention curve's form | litter: Letts et al. (2000)'s surface moss, as CLASS takes it for litter and moss (Wu et al. 2016); moss: Zhang et al. (2010) |
| `theta_sat`, `theta_res` | 0.98, 0.01 | 0.85, 0.18 | water at saturation and residual [m³ m⁻³] | litter: Wu et al. (2016) Table 1, surface moss (0–10 cm); moss: Zhang et al. (2010) Table 4, surface moss and lichen (0–0.1 m) |
| `curve_par_a`, `curve_par_n` | −0.0103 m, 2.3 | 50 m⁻¹, 1.45 | Campbell's air-entry head and $b$; van Genuchten's α and $n$ | the same |
| `ksat` | 1.83×10⁻³ m s⁻¹ | 1.3×10⁻³ m s⁻¹ | saturated conductivity | litter: the same; moss: O'Connor et al. (2020), the Arctic Foothills' acrotelm |

`[soil].moss_depth` [m], which a `[[sites]]` table may set, is 0 or at least `min_depth`, and needs a table whose `mass_source` is `site_depth`; one table at most takes the metabolic litter. The checks refuse a non-positive bulk density, conductivity or saturated conductivity, a residual outside $`[0, \theta_s)`$, a saturation outside (0, 1), curve parameters of the wrong sign (van Genuchten's α > 0 and n > 1; Campbell's air-entry head < 0 and b > 0), a `fresh_density` not below solid ice's, an albedo outside [0, 1), an emissivity outside (0, 1], and a `liquid_holding_frac` outside [0, 1). The keys of the `[snow]` and `[surface_layer]` blocks stop a run, naming the key that holds their value.

Fixed in the code: $`\rho_{ref}`$ = 100 kg m⁻³ in eq. (5) (Niu & Yang 2007's fresh snow), $`k_w`$ = 0.57 W m⁻¹ K⁻¹ for still water, $`c_{org}`$ = 1926 J kg⁻¹ K⁻¹ (Farouki 1986), and the stack solver's iteration limits and tolerance (`STACK_*` in `ed_surface`), which are numerical.

## Code map

| Concept | Routine |
|---|---|
| the `[ground_surface_layer]` block | `ed_ground_surface_params`: `ground_surface_layer_params_t`, `water_material_t`, `porous_material_t`, `POROUS_LITTER`, `POROUS_MOSS`, `MASS_FROM_*`; `ed_config_io`: `load_ground_surface_params`; `ed_config_check`: `check_ground_surface_layer`, `check_water_material`, `check_porous_material`, `MOVED_KEYS` |
| the state and the moves between stores | `ed_patch_state_vars`: `ground_store_t`, `water_layer_state_t`, `porous_layer_state_t`, `ground_surface_state_t`, `store_move`, `store_take_exchange`, `ground_surface_store`, `ground_top_temp`, `porous_layer_curve`, `porous_layers_diagnose`, `porous_layers_spill`, `blend_ground_surface` |
| the layers' laws: phase, depth, matrix, cover, optics, conduction, the film | `ed_ground_surface_lib`: `porous_layer_temp`, `porous_energy_state`, `porous_energy_at_temp`, `film_on_porous_layer`, `porous_layer_depth`, `porous_layer_capacity`, `porous_layer_fill`, `porous_layer_resistance`, `porous_layer_matter`, `water_layer_free_water`, `water_layer_depth`, `matrix_depth_update`, `water_layer_cover`, `water_layer_albedo`, `water_layer_optics_shares`, `snow_matrix_conductivity`, `water_layer_resistance` |
| the stage, and the radiation after it | `ed_fast_integ_frozen`: `ground_stage`, `intercept_rain`; `ed_surface`: `water_layer_stage`, `water_layer_exposed`, `water_layer_cover_of`; `ed_fast_integ_radiation`: `canopy_ground_radiation`; `ed_surface`: `snow_blend_ground_optics`, `ground_tile_radiation` |
| the stack's energy balance and its commit | `ed_surface`: `ground_surface_exchange`, `surface_stack_balance`, `surface_node_balance`, `stack_residual`, `link_heat`, `advect_through`, `curve_node_root`, `water_layer_commit`; `ed_fast_integ_update`: `update_dynamic_state` |
| saturation vapour pressure, eqs. (11)–(12) | `ed_therm_lib`: `sat_vapor_pressure`, `sat_vapor_pressure_temp_deriv`, `sat_specific_humidity` |
| the layers in the water network | `ed_fast_integ_frozen`: `porous_network_drive`, `ground_inflow`, `porous_nodes`; the network: [soil_plant_water](soil_plant_water.md) |
| evaporation through the dried tops and the shared leg | `ed_soil_water`: `ground_vapour`, `store_vapour_path`, `evaporation_supply`, `slot_evaporation_setup`, `dry_layer_resistance`, `porous_cover_path` |
| the layers follow their matter; the start | `ed_slow_dynamics`: `follow_porous_sources`; `ed_slow_update`: `follow_porous_layers`; `ed_init`: `seed_porous_layers`, `seed_snow`; `ed_init_restart`: `read_older_surface_layers` |
| burial | `ed_allometry_lib`: `crown_exposed_fraction`; `ed_fast_state`: `gather_patch`; `ed_fast_integ_frozen`: `buried_tissue_store`; `ed_plant_energy`: `commit_tissue_energy`, `tissue_water_heat` |

## References
- Ambaum, M. H. P. (2020), *Q. J. R. Meteorol. Soc.* 146:4252, doi:10.1002/qj.3899 — saturation vapour pressure over water and ice.
- Brondex, J., et al. (2023), *Geosci. Model Dev.* 16:7075, doi:10.5194/gmd-16-7075-2023 — gives Calonne et al. (2011)'s regression and its fitted range (eq. A3).
- Calonne, N., Flin, F., Morin, S., Lesaffre, B., Rolland du Roscoat, S., and Geindreau, C. (2011), *Geophys. Res. Lett.* 38:L23501, doi:10.1029/2011GL049234 — snow's effective thermal conductivity against its density.
- Chadburn, S. et al. (2015), *Geosci. Model Dev.* 8:1493 — the moss layer's thermal conductivity in JULES.
- de Vries, D. A. (1963), in *Physics of Plant Environment*, North-Holland — thermal properties of soils.
- Dickinson, R. E., Henderson-Sellers, A. & Kennedy, P. J. (1993), *Biosphere-Atmosphere Transfer Scheme (BATS) Version 1e as Coupled to the NCAR Community Climate Model*, NCAR/TN-387+STR, doi:10.5065/D67W6959 — new snow's albedo (§2e, eq. 5).
- Farouki, O. T. (1986), *Thermal Properties of Soils*, Trans Tech — the heat capacity of organic matter.
- Kosugi, K., Mori, K. & Yasuda, H. (2001), *J. Hydrol.* 246:96–108, doi:10.1016/S0022-1694(01)00366-3 — the hydraulic properties of a forest floor.
- Lawrence, D. et al. (2018), *Technical Description of version 5.0 of the Community Land Model (CLM)*, NCAR — lake ice, open-water and snow optics.
- Letts, M. G. et al. (2000), *Atmosphere-Ocean* 38:141 — the hydraulic parameters of organic soils.
- Matthews, S. (2014), *Int. J. Wildland Fire* 23:78–92, doi:10.1071/WF13005 — a review of the water dynamics of litter and fuels.
- McLaren, J. (2021), *Soil biogeochemical variables collected on the Arctic LTER experimental plots in moist acidic, moist non-acidic, wet shrub and shrub tundra, Toolik Field Station, Alaska 2015*, Arctic LTER, knb-lter-arc.20129.1 — the bulk density of the top 5 cm of the organic soil.
- Murphy, D. M. & Koop, T. (2005), *Q. J. R. Meteorol. Soc.* 131:1539 — vapour pressure over ice.
- Napoly, A. et al. (2017), *Geosci. Model Dev.* 10:1621–1644 — the ISBA-MEB forest-floor litter layer.
- Niu, G.-Y. & Yang, Z.-L. (2007), *J. Geophys. Res.* 112:D21101, doi:10.1029/2007JD008674 — snow-cover fraction (eq. 4).
- O'Connor, M. T. et al. (2020), *Geophys. Res. Lett.* 47:e2020GL087646 — the saturated conductivity of the Arctic Foothills' acrotelm.
- Oleson, K. W. et al. (2010), *Technical Description of version 4.0 of the Community Land Model (CLM)*, NCAR/TN-478+STR, doi:10.5065/D6FB50WZ — the snow-cover fraction's density exponent for CLM (eq. 3.50).
- Oleson, K. W. et al. (2013), *Technical Description of version 4.5 of the Community Land Model (CLM)*, NCAR/TN-503+STR, doi:10.5065/D6RR1W7M — the bulk density of snow (§2.2.2) and of new snow (eq. 7.21); surface water (§7.3.2) and soil ponding (§7.6).
- Roesch, A., Wild, M., Pinker, R. & Ohmura, A. (2002), *J. Geophys. Res.* 107(D14):4221, doi:10.1029/2001JD000809 — snow albedo on land, the source of ED2's aged-snow albedo.
- Swenson, S. C. & Lawrence, D. M. (2014), *J. Geophys. Res. Atmos.* 119:10299 — the dry surface layer's resistance to evaporation.
- Walko, R. L. et al. (2000), *J. Appl. Meteorol.* 39(6):931–944, doi:10.1175/1520-0450(2000)039<0931:CABHMF>2.0.CO;2 — LEAF-3, its snow's liquid holding.
- Wu, Y., Verseghy, D. L. & Melton, J. R. (2016), *Geosci. Model Dev.* 9:2639 — CLASS's organic soil parameters, Table 1.
- Zhang, Y., Carey, S. K., Quinton, W. L. et al. (2010), *Hydrol. Earth Syst. Sci.* 14:729 — the retention of surface moss at Scotty Creek, Table 4.
