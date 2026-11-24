# The surface organic layer

The organic layer is the dead and living matter on the soil: leaf litter, and moss under it. It sits between the snow pack (or the canopy air) and the top soil layer. Where it exists it is the ground's surface: it absorbs the radiation the ground receives, exchanges heat and vapour with the canopy air, holds and passes water as part of the soil-plant water network, and insulates the soil below it. Without it the ground is the bare soil, exactly as before.

This is one of the per-store pages under [coupled_biophysics](coupled_biophysics.md). The organic layer is built on the snow pack's template ([snow_biophysics](snow_biophysics.md)): one node with its own energy and water, its balance relaxed over the step with its own heat store, coupled to the soil top.

---

## The layer and its slots

The layer has two slots, top down: **litter** over **moss**. Each slot carries four numbers:

| | Unit | Litter | Moss |
|---|---|---|---|
| dry matter $`M`$ | kg m⁻² | the above-ground metabolic litter pool's | a site's `[soil].moss_depth` × its bulk density |
| water $`W`$ (liquid and ice) | kg m⁻² | prognostic | prognostic |
| internal energy $`E`$ | J m⁻² | prognostic | prognostic |
| bulk density $`\rho_b`$ | kg m⁻³ | its material's | its material's |

The depth is $`d = M/\rho_b`$, and the water a slot holds saturated is $`W_{max} = \theta_{s}\,\rho_w\,d`$, with $`\theta_{s}`$ its material's water at saturation, its retention curve's porosity (§ Water). Its fill is $`W/W_{max}`$.

A slot is a layer only from 0.1 mm deep (`[ground_surface_layer].min_depth`, `porous_layer_matter`). Thinner, a film of a few micrometres such as a newly opened gap's first litter, it would hold a fraction of a joule per kelvin and a trace of water, which no balance carries: in the solver's tests a 1 µm litter (0.1 J K⁻¹) evaporated its last water and fell to 96 K, and 19 µm under a trace pond forced sub-steps. Its matter stays in its pool, and the ground's surface is what lies below. The litter becomes a layer when its pool reaches the depth, bringing water at the top soil's head as new litter does, and stops being one when its pool falls below it, its water and that water's energy joining the pond.

**One thermal node.** The two slots share one temperature, read off their summed energy, water and dry heat capacity on the freezing plateau, as the snow's is (`internal_energy_to_temp`):

```math
E = \sum_i E_i = C_d T + W\left[f\,c_{liq}(T - T_{sc}) + (1 - f)\,c_{ice} T\right], \qquad C_d = c_{org} \sum_i M_i
```

with $`c_{org}`$ = 1926 J kg⁻¹ K⁻¹ for dry organic matter (Farouki 1986). The water is free water: it freezes at $`T_{3ple}`$ with no depression. Each slot keeps its share of the energy at the common temperature, so a slot can later be given its own node without a change of state.

## Where the matter comes from

- **Litter** follows the above-ground metabolic litter pool each slow step, after the soil-carbon step: $`M = 2\,C_{fast,grnd}`$, half of dry organic matter being carbon (`c2b`). The structural pool is left out: it holds the coarse woody debris too, and a big tree's dead wood spread over its gap patch would make a layer tens of centimetres deep. It joins once the soil carbon splits fine litter from woody debris. Matter that arrives or leaves does so at the layer's temperature (the top soil layer's, in a layer it starts), and the dry heat it brings or takes is a declared term of the slow ledger. Matter that arrives brings the water its own volume holds at the matric head it and the top soil layer share once that water has left the layer (the layer's head read from its liquid, on its freezing curve where it is frozen), drawn from the top soil layer with the energy that water has there, at the layer's temperature and liquid fraction: a matched move, which the slow ledger needs no declaration for. So a large fall of litter on dry soil takes what the soil can give at a common head, and never draws the soil below its air-dry head (CLM5's −10⁴ m, `PSI_AIR_DRY`): drawn to its residual, the soil's head would have no bound, and the roots in it none either. Matter arriving on soil drier than that comes dry, and a layer it leaves below its residual carries its deficit as an offset (Water, above). Water above what the slot's new depth holds saturated, or all its water when the litter is gone, passes to the water layer with its energy. `[ground_surface_layer.litter].mass_source = "none"` keeps the litter slot empty.
- **Moss** is prescribed: its depth is a site key, `[soil].moss_depth` (0 by default, so no moss), at its material's bulk density. The checks refuse a moss thinner than the 0.1 mm from which a slot is a layer (`[ground_surface_layer].min_depth`), which would otherwise be dropped.
- **One path for both layers** (`follow_porous_sources`, `follow_porous_layers`): the litter and the moss follow the matter their `mass_source` gives at a run's start, at a restart and at every slow step. A fresh start builds them from nothing, a restart from what its state file held, so a layer the file lacks (an older file without the moss) is made, and one its `mass_source` no longer gives is removed. At a start the water the matter brings arrives from outside, at the top soil layer's head, temperature and liquid fraction.
- **At a start** the slots are filled with this matter, each holding the water its curve holds at the top soil layer's matric head, as the water network would bring it to, at the top soil layer's temperature.

## Water

- **Nodes of the water network.** Each slot is a store of the soil-plant water network ([soil_plant_water](soil_plant_water.md) §1), between the pond and the top soil layer: the pond over the litter, the litter over the moss, the moss over the soil. Water moves between them down its potential, gravity and capillarity, by the soil's own face law, so the litter and the moss fill, drain to the soil, and are drawn on by a drying soil as a soil layer is. The pond's water reaches the ground through the top slot; a slot pushed past saturation pushes water back up to the pond. A slot is a soil layer to the solver, on its material's retention curve and saturated conductivity at its depth, rebuilt each step: the litter on Campbell's curve, the moss on van Genuchten's (§ Parameters). At −10 kPa the moss holds 0.29 m³ m⁻³, within the 0.26–0.49 measured at 2–6 cm in Arctic tussock tundra in summer, and the litter 0.14, near the 0.12 drained litter holds (Putuhena & Cordery 1996). The roots take nothing from the slots.
- **Frozen slots.** The layer's water freezes at $`T_{3ple}`$ with no freezing curve in its energy. For the water network a slot's ice is held over the step as a soil layer's is, its liquid on its freezing curve at the layer's temperature (`soil_liquid_capacity`, as a frozen soil layer's): it adds its water, slows the flow by its impedance $`10^{-6\,\theta_{ice}/\theta_s}`$ (CLM5), and leaves its pores to the liquid. That liquid is the solver's start only: the layer's energy still freezes all its water at the melting point; a freezing curve in the layer's energy belongs with the unfrozen water still to come. On the melt plateau the two part: the curve holds the slot all liquid at $`T_{3ple}`$, while its energy may hold most of it as ice. There the slot's impedance takes the energy's ice, $`(1 - f_{liq})\,\theta`$, so a mostly frozen litter is throttled as any ice throttles; and the water it gives up leaves as its energy holds it (below, "The water's heat"). A frozen saturated moss passes almost nothing, so rain on it ponds; dry frozen litter takes it. A slot holding less than its residual (litter just arrived, dry) carries the water it lacks as an offset over the step, so its potential is its curve's dry end and its mass exact.
- **The snow's meltwater** reaches the pond, as on bare ground, and the slots through it.
- **Above saturation.** What a slot holds above saturation, which a saturated or frozen slot can neither take nor pass on, joins the pond at the next step's start, with its energy at the layer's temperature (`organic_spill`).
- **Evaporation** from each porous layer, the litter and the moss beneath it alike, leaves its water through its own dried top and the open pores of the layers above it, on the snow-free ground. Each layer holds its water in its lower part and dries from its top, as the soil's top layer does under its dry surface layer (Swenson & Lawrence 2014), and by the same law (`dry_layer_resistance`, `dry_layer_setup`): a dried depth that forms below the layer's saturation and is the whole layer when it holds nothing, $`d_i\,(\theta_{s,i} - \theta_i)/\theta_{s,i}`$, crossed at water vapour's diffusivity in air times the Millington–Quirk tortuosity of its open pores, $`\theta_{s,i}^{4/3}`$, so

```math
g_i = \frac{1}{r_{dry,i}(\theta_i) + \sum_{j<i} r_{cover,j}}, \qquad
r_{dry,i} = \frac{d_i\,(\theta_{s,i} - \theta_i)/\theta_{s,i}}{D_v\,\theta_{s,i}^{4/3}}
```

is its path's conductance up to the ground's surface node, with $`r_{cover,j}`$ the open pores of each layer above (as the soil's vapour crosses them, below), held over the step. The paths of every layer and of the soil's top meet at that node in parallel, and the node passes their sum to the canopy air through the ground's one aerodynamic leg $`1/g_{ground}`$, the leg its sensible heat takes, so the leg is counted once however many layers evaporate. Layer $`i`$ evaporates $`(1 - f_{w})\,\rho\,g_i\,[\alpha(\psi_i)\,q_{sat}(T_i) - q_s]^+`$, the node's humidity $`q_s`$ in closed form ([soil_biophysics](soil_biophysics.md), "The soil water column", `ground_vapour`). The dried top is the layer's own and set by its own water, so it is implicit in the water network, as the soil's dry surface layer is; the moss under a litter evaporates through the litter's open pores. A layer full to its top exchanges through no dried depth. A drying one exchanges through it: 3 cm of dry moss add about 1,700 s/m at 280 K, many times the ground's own aerodynamic resistance. The water network takes each layer's evaporation by the soil's own law (`store_evaporation`, with the layer's boundary from `slot_evaporation_setup`): Philip's pore humidity $`\alpha(\psi) = \exp(\psi g/R_vT)`$ at the layer's potential, implicit in the network, with no dew: a layer drier than the node takes nothing, from the air or from the other layers; and the soil's supply limit at the residual. A drying layer so comes to its equilibrium with the canopy air's humidity: the litter at 30 % relative humidity near −1.7×10⁴ m and θ ≈ 0.012, at 90 % near −1.5×10³ m and θ ≈ 0.016. The energy stage then holds the network's rate, as it holds the soil's evaporation. The diffusivity is the layer's temperature's, as the soil's.
- **The soil beneath** evaporates through the layer: its vapour diffuses up through the slots' air-filled pores, a resistance $`r_{cover} = \sum_i d_i\,\theta_{s,i}^2 / (D_v\,\varepsilon_i^{10/3})`$ (Millington–Quirk) in series with the dry surface layer's up to the ground's surface node, and from there through the shared leg, where $`\varepsilon_i = \theta_{s,i} - W_i/(\rho_w d_i)`$: a saturated slot has no open pores and passes no vapour.

## Heat

The layer's balance over the step is the pack's own (`surface_node_balance`): a linearized backward-Euler balance of its energy, with its dry heat capacity on each of the ice line, the melt plateau and the liquid line, and one root in closed form. Its terms are:

- on the snow-free ground $`(1 - f_{snow})`$: its net radiation at the ground's albedo and emissivity (`[soil].ground_albedo_*`, `ground_emissivity`) and its own temperature, less the sensible heat to the canopy air;
- its evaporation, the water network's rate per unit ground area held over the step, wherever the pack leaves it. Its water is apart from what the layer keeps over the step, the network's end water, and the vapour takes its whole enthalpy at the layer's temperature from the layer's energy; the layer's phase law then reads what is left, so the balance does not depend on which phase gave the vapour. Valuing the vapour as its latent heat at the step's start plus its water's own energy at the end, the ice's below its freezing curve's onset and the liquid's above, made the latent heat jump by the latent heat of fusion where the layer's temperature crossed the onset: 9.5×10⁻⁶ kg m⁻² s⁻¹ from a saturated litter at the triple point made 3.2 W m⁻², and the stage's canopy-air balance over that jump had no root (a tundra's late-summer snowfall on wet litter);
- under the pack: the pack's base flux, which conducts into the layer through the whole snow depth and the layer's upper half;
- over the whole patch: conduction into the top soil node through the layer's lower half and the top soil node's half layer, in series:

```math
G = \frac{T - T_{soil,1}}{R/2 + |z_1|/k_{soil,1}}, \qquad R = \sum_i d_i\left[\frac{1 - W_i/W_{max,i}}{k_{dry,i}} + \frac{W_i/W_{max,i}}{k_{wet,i}}\right]
```

Each slot holds its water in its lower part and dries from its top (§ Water), so its heat crosses the dried top at the slot's dry conductivity, in series with the wet part at its wet one: the dried depth that resists the layer's vapour resists its heat too. 3 cm of dried moss over 3 cm of wet moss conduct 0.11 W m⁻¹ K⁻¹, where a conductivity linear in the fill would give 0.28.

The balance is implicit in the layer's own temperature, so a thin layer of little heat capacity relaxes towards its balance without overshooting it. The soil receives $`G`$ as its top flux, less its own evaporation's latent heat, as it receives the pack's. The layer takes the exchange's time integral after the march, as the pack does, with the vapour the water network took from each slot.

**The water's heat.** The water the network moves carries its enthalpy, operator-split from the march: after the network's step each slot takes its new water, and the layer the enthalpy its flows carried, shared among the slots at the layer's new temperature. Water comes in from the pond at the pond's mixture and up from the soil at the top soil layer's temperature; water passing through the layer in a sub-step leaves, to the soil or back to the pond, at the temperature it came in at, as far as what came in covers what leaves, and the rest as the layer's energy holds its water, its liquid share at its temperature and the rest ice, as its spill does. Valued as liquid, the water a nearly frozen litter on the melt plateau gives up would take latent heat the litter does not have (0.05 kg m⁻² drained from frozen litter at Harvard took it to 187 K); the soil refreezes the ice share by its freezing curve, the pond in the snow stage (`pond_freezes`). Valued at the layer's start temperature, the water a thin litter passes on would heat or cool it by its whole difference from that temperature: 1.7 kg m⁻² through 0.6 mm of litter, a store of 700 J K⁻¹, took one to 524 K. The march then balances the layer as it now is; the soil's top face takes the water the bottom slot passes down at the temperature it leaves at ([soil_plant_water](soil_plant_water.md) §5). Rain on frozen litter melts some of its ice and warms it, as it should.

## What reads the layer

- **The two-stream and the aerodynamics:** the ground emits at the layer's temperature, and the under-canopy stability reads it.
- **Decomposition:** the above-ground litter pools (metabolic and structural; both lie on the ground, though only the metabolic one makes the layer) take the layer's temperature and, as their moisture, its saturation above residual, the water the porous layers hold beyond their residual over what they hold beyond it saturated, as the soil's moisture is read, instead of the top soil layer's ([soil_carbon](soil_carbon.md)).
- **Outputs:** `mean_ground_top_temp_py` (the layer's temperature, the soil top's where there is no layer), `mean_litter_water_py` and `mean_moss_water_py`, and their per-patch forms.

## Ledgers

The layer's water and energy are stores of the whole-column water and energy budgets, and of the slow ledger. The water it exchanges with the pond and the soil is a transfer between stores, with its enthalpy. Its evaporation leaves through the canopy air's vapour. The slow step's change of dry matter carries a declared energy. `mean_infiltration_py` is the pond's infiltration into the top slot where there is a layer.

## Parameters

The `[ground_surface_layer]` block's `min_depth` (0.1 mm, the depth from which a slot's matter is a layer) and its `litter` and `moss` tables, every key optional ([configuration](../configuration.md)):

| Key | Default | Meaning | Source |
|---|---|---|---|
| `litter.mass_source` | metabolic_litter | where the litter's matter comes from: `none`, `metabolic_litter` or `site_depth` | |
| `litter.bulk_density` | 45 kg m⁻³ | dry mass per volume | Napoly et al. (2017) |
| `litter.k_dry`, `litter.k_wet` | 0.10, 0.1294 W m⁻¹ K⁻¹ | conductivity dry and saturated | de Vries (1963): 0.1 + 0.03 θ, via Napoly et al. (2017), at θ = 0 and `litter.theta_sat` |
| `litter.retention` | campbell | its retention curve's form | Letts et al. (2000), as CLASS takes its surface moss for litter and moss (Wu et al. 2016) |
| `litter.theta_sat`, `litter.theta_res` | 0.98, 0.01 m³ m⁻³ | water at saturation and residual | the same: the surface moss (0–10 cm) row, Wu et al. (2016) Table 1 |
| `litter.curve_par_a`, `litter.curve_par_n` | −0.0103 m, 2.3 | Campbell's air-entry head and $`b`$ | the same |
| `litter.ksat` | 1.83×10⁻³ m s⁻¹ | saturated conductivity | the same |
| `moss.mass_source` | site_depth | where the moss's matter comes from: a site's `[soil].moss_depth` | |
| `moss.bulk_density` | 43 kg m⁻³ | dry mass per volume | moist acidic tussock tundra, top 5 cm of the organic layer (ARC LTER, knb-lter-arc.20129.1) |
| `moss.k_dry`, `moss.k_wet` | 0.06, 0.50 W m⁻¹ K⁻¹ | conductivity dry and saturated | Chadburn et al. (2015) |
| `moss.retention` | van_genuchten | its retention curve's form | Zhang et al. (2010) |
| `moss.theta_sat`, `moss.theta_res` | 0.85, 0.18 m³ m⁻³ | water at saturation and residual | Zhang et al. (2010) Table 4, surface moss and lichen (0–0.1 m) at Scotty Creek |
| `moss.curve_par_a`, `moss.curve_par_n` | 50 m⁻¹, 1.45 | van Genuchten's α and n | the same |
| `moss.ksat` | 1.3×10⁻³ m s⁻¹ | saturated conductivity | O'Connor et al. (2020), the Arctic Foothills' acrotelm |

No retention curve of leaf litter has been measured; Letts et al.'s surface moss is the curve CLASS gives both its litter and its moss, and at −10 kPa it holds 0.14 m³ m⁻³, near the water drained litter holds, 0.12 (Putuhena & Cordery 1996, via Napoly et al. 2017).

and `[soil].moss_depth` [m], which a `[[sites]]` table may set; it is 0 or at least `min_depth`, and needs a table whose `mass_source` is `site_depth`. One table at most takes the metabolic litter. The checks refuse a non-positive bulk density, conductivity or saturated conductivity, a residual outside $`[0, \theta_s)`$, a saturation outside (0, 1), and curve parameters of the wrong sign (van Genuchten's α > 0 and n > 1; Campbell's air-entry head < 0 and b > 0).

Every slot that is a layer, from 0.1 mm, is a node of the water network; the solver's node classes ([soil_plant_water](soil_plant_water.md) §3) eliminate a slot whose time scale is far below the sub-step and detach a group of nodes it cannot see, so neither a third of a millimetre of litter nor a frozen slot needs a floor of its own in the solver.

## Limits

- The two slots are one thermal node: moss under litter does not hold a temperature of its own.
- Its optics are the ground's.
- The moss is prescribed: it neither grows nor photosynthesises.
- The litter is the metabolic pool alone: the structural leaf litter is not in the layer until the soil carbon splits it from the woody debris.
- Each slot dries from its top. A Sphagnum moss draws water up to its tips and keeps its top wetter, so its dried depth would be thinner.
- The slots have no freezing curve: their water freezes at the melting point, all of it.
- Within a step the slots' evaporation does not follow the canopy air's changing humidity; the next step's surface solve restores the response, as for the soil's.

## Code map

| Concept | Routine |
|---|---|
| the `[ground_surface_layer]` block | `ed_ground_surface_params`: `ground_surface_layer_params_t`, `porous_material_t`, `POROUS_LITTER`, `POROUS_MOSS`, `MASS_FROM_*`; `ed_config_io`: `load_ground_surface_params`; `ed_config_check`: `check_ground_surface_layer` |
| the porous layers' state, their sums, each one's curve, their energy shared among them, their spill and their commit, the ground top's temperature, the blend and the total the ledgers read | `ed_patch_state_vars`: `porous_layer_state_t` (in `ground_surface_state_t`), `porous_layers_sums`, `porous_layers_diagnose`, `porous_layers_split_energy`, `porous_layer_curve`, `porous_layers_spill`, `porous_layers_take_exchange`, `ground_top_temp`, `blend_ground_surface`, `ground_surface_store` |
| a porous layer's laws: its depth, capacity, fill, resistance and the matter that makes a layer | `ed_ground_surface_lib`: `porous_layer_depth`, `porous_layer_capacity`, `porous_layer_fill`, `porous_layer_resistance`, `porous_layer_matter` |
| the balance and the ground surface's blend | `ed_surface`: `surface_node_balance`, `ground_surface_exchange` |
| the slots in the water network, and their water and enthalpy after it | `ed_fast_integ_frozen`: `organic_network_drive`, `commit_organic_water`; the network: [soil_plant_water](soil_plant_water.md) |
| the soil's vapour through it, and its own through its dried top | `ed_soil_water`: `organic_cover_path`, `vapour_diffusivity`, `ground_evaporation_rate`; `ed_fast_integ_frozen`: `organic_ground_sums`, `ground_inflow` |
| the layers follow their matter, at the start, a restart and every slow step | `ed_slow_update`: `follow_porous_layers`; `ed_slow_dynamics`: `follow_porous_sources`; `ed_init`: `seed_porous_layers` |

## References
- Chadburn, S. et al. (2015), *Geosci. Model Dev.* 8:1493 — the moss layer's thermal conductivity in JULES.
- de Vries, D. A. (1963), in *Physics of Plant Environment*, North-Holland — thermal properties of soils.
- Farouki, O. T. (1986), *Thermal Properties of Soils*, Trans Tech — the heat capacity of organic matter.
- Letts, M. G. et al. (2000), *Atmosphere-Ocean* 38:141 — the hydraulic parameters of organic soils, the surface moss row the litter's.
- Napoly, A. et al. (2017), *Geosci. Model Dev.* 10:1621 — the ISBA-MEB forest-floor litter layer.
- O'Connor, M. T. et al. (2020), *Geophys. Res. Lett.* 47:e2020GL087646 — the saturated conductivity of the Arctic Foothills' acrotelm.
- Putuhena, W. M. & Cordery, I. (1996), *J. Hydrol.* 180:283 — the forest floor's interception capacity.
- Swenson, S. C. & Lawrence, D. M. (2014), *J. Geophys. Res. Atmos.* 119:10299 — the dry surface layer's resistance to evaporation.
- Wu, Y., Verseghy, D. L. & Melton, J. R. (2016), *Geosci. Model Dev.* 9:2639 — CLASS's organic soil parameters, Table 1.
- Zhang, Y., Carey, S. K., Quinton, W. L. et al. (2010), *Hydrol. Earth Syst. Sci.* 14:729 — the retention of surface moss at Scotty Creek, Table 4.
- McLaren, J. (2021), *Soil biogeochemical variables collected on the Arctic LTER experimental plots in
  moist acidic, moist non-acidic, wet shrub and shrub tundra, Toolik Field Station, Alaska 2015*, Arctic LTER,
  knb-lter-arc.20129.1 — the bulk density of the top 5 cm of the organic layer.
