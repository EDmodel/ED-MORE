# The surface organic layer

The organic layer is the dead and living matter on the soil: leaf litter, and moss under it. It sits between the snow pack (or the canopy air) and the top soil layer. Where it exists it is the ground's surface: it absorbs the radiation the ground receives, exchanges heat and vapour with the canopy air, holds water in its pores, and insulates the soil below it. Without it the ground is the bare soil, exactly as before.

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

The depth is $`d = M/\rho_b`$, and the water a slot holds at capacity is $`W_{max} = \theta_{cap}\,\rho_w\,d`$, with $`\theta_{cap}`$ its material's water capacity per volume. Its fill is $`W/W_{max}`$.

**One thermal node.** The two slots share one temperature, read off their summed energy, water and dry heat capacity on the freezing plateau, as the snow's is (`internal_energy_to_temp`):

```math
E = \sum_i E_i = C_d T + W\left[f\,c_{liq}(T - T_{sc}) + (1 - f)\,c_{ice} T\right], \qquad C_d = c_{org} \sum_i M_i
```

with $`c_{org}`$ = 1926 J kg⁻¹ K⁻¹ for dry organic matter (Farouki 1986). The water is free water: it freezes at $`T_{3ple}`$ with no depression. Each slot keeps its share of the energy at the common temperature, so a slot can later be given its own node without a change of state.

## Where the matter comes from

- **Litter** follows the above-ground metabolic litter pool each slow step, after the soil-carbon step: $`M = 2\,C_{fast,grnd}`$, half of dry organic matter being carbon (`c2b_organic`). The structural pool is left out: it holds the coarse woody debris too, and a big tree's dead wood spread over its gap patch would make a layer tens of centimetres deep. It joins once the soil carbon splits fine litter from woody debris. Matter that arrives or leaves does so at the layer's temperature, and the energy it brings or takes is a declared term of the slow ledger. Water above the slot's new capacity, or all its water when the litter is gone, passes to the pond with its energy. `[surface_layer].litter_on = false` keeps the litter slot empty.
- **Moss** is prescribed: its depth is a site key, `[soil].moss_depth` (0 by default, so no moss), at its material's bulk density.
- **At a start** the slots are filled with this matter, each holding the top soil layer's share of saturation of its capacity, at the top soil layer's temperature.

## Water

- **Filling.** At each fast step's start the pond, and then the water reaching the ground (throughfall, or the snow's meltwater under a pack, with the tissues' shed water), fill the slots top down, each to its capacity. The water brings its energy. What the slots cannot hold reaches the pond, which infiltrates or runs off as before. A pond therefore forms only on a full layer.
- **Evaporation and dew** from and onto the layer are the wet-surface rate times its fill, on the snow-free ground:

```math
E_{org} = (1 - f_{snow})\,\frac{W}{W_{max}}\,\rho\,g_{ground}\,\left[q_{sat}(T) - q_{cas}\right]
```

A wet layer exchanges as open water does, both ways, so the exchange has no corner where humid canopy air sits near the layer's saturation; a dry one does not exchange. The evaporation takes no more water than the layer holds, and dew beyond its capacity spills to the pond at the next step's start.
- **The soil beneath** evaporates through the layer: its vapour diffuses up through the slots' air-filled pores, a resistance $`r_{cover} = \sum_i d_i\,\phi^2 / (D_v\,\varepsilon_i^{10/3})`$ (Millington–Quirk) in series with the dry surface layer's and the aerodynamic one. Here $`\phi`$ = 0.9 is organic matter's porosity (Letts et al. 2000) and $`\varepsilon_i = \phi - W_i/(\rho_w d_i)`$.

## Heat

The layer's balance over the step is solved like the pack's (`organic_surface_balance`): a linearized backward-Euler balance of its energy, with its dry heat capacity on each of the ice line, the melt plateau and the liquid line, and one root in closed form. Its terms are:

- on the snow-free ground $`(1 - f_{snow})`$: its net radiation at the ground's albedo and emissivity (`[soil].ground_albedo_*`, `ground_emissivity`) and its own temperature, less the sensible heat and the evaporation to the canopy air;
- under the pack: the pack's base flux, which conducts into the layer through the half snow depth and the layer's upper half;
- over the whole patch: conduction into the top soil node through the layer's lower half and the top soil node's half layer, in series:

```math
G = \frac{T - T_{soil,1}}{R/2 + |z_1|/k_{soil,1}}, \qquad R = \sum_i \frac{d_i}{k_i}, \qquad k_i = k_{dry} + (k_{wet} - k_{dry})\,\frac{W_i}{W_{max,i}}
```

The balance is implicit in the layer's own temperature, so a thin layer of little heat capacity relaxes towards its balance without overshooting it. The soil receives $`G`$ as its top flux, less its own evaporation's latent heat, as it receives the pack's. The layer takes the exchange's time integral after the march, as the pack does.

## What reads the layer

- **The two-stream and the aerodynamics:** the ground emits at the layer's temperature, and the under-canopy stability reads it.
- **Decomposition:** the above-ground litter pools (metabolic and structural; both lie on the ground, though only the metabolic one makes the layer) take the layer's temperature and its fill as their moisture, instead of the top soil layer's ([soil_carbon](soil_carbon.md)).
- **Outputs:** `mean_organic_temp_py` (the soil top's where there is no layer) and `mean_organic_water_py`, and their per-patch forms.

## Ledgers

The layer's water and energy are stores of the whole-column water and energy budgets, and of the slow ledger. The water it takes from the pond is a transfer between two stores. The water it takes from the ground's inflow arrives with the inflow, which the budgets already count. Its evaporation leaves through the canopy air's vapour. The slow step's change of dry matter carries a declared energy.

## Parameters

The `[surface_layer]` block, every key optional ([configuration](../configuration.md)):

| Key | Default | Meaning | Source |
|---|---|---|---|
| `litter_on` | true | the litter slot | |
| `litter_bulk_density` | 45 kg m⁻³ | dry mass per volume | Napoly et al. (2017) |
| `litter_water_capacity` | 0.12 m³ m⁻³ | water it holds | Putuhena & Cordery (1996), via Napoly et al. (2017) |
| `litter_k_dry`, `litter_k_wet` | 0.10, 0.1036 W m⁻¹ K⁻¹ | conductivity dry and at capacity | de Vries (1963): 0.1 + 0.03 θ, via Napoly et al. (2017) |
| `moss_bulk_density` | 43 kg m⁻³ | dry mass per volume | moist acidic tussock tundra, top 5 cm of the organic layer (ARC LTER, knb-lter-arc.20129.1) |
| `moss_water_capacity` | 0.90 m³ m⁻³ | water it holds | Chadburn et al. (2015) |
| `moss_k_dry`, `moss_k_wet` | 0.06, 0.50 W m⁻¹ K⁻¹ | conductivity dry and saturated | Chadburn et al. (2015) |

and `[soil].moss_depth` [m], which a `[[sites]]` table may set. The checks refuse a non-positive bulk density or conductivity, and a water capacity outside (0, 0.9], the pores' share.

## Limits

- The two slots are one thermal node: moss under litter does not hold a temperature of its own.
- Its optics are the ground's.
- The moss is prescribed: it neither grows nor photosynthesises.
- The litter is the metabolic pool alone: the structural leaf litter is not in the layer until the soil carbon splits it from the woody debris.

## Code map

| Concept | Routine |
|---|---|
| the layer's state, its sums, its water, its spill and its commit | `ed_organic_lib`: `organic_state_t`, `organic_layer_sums`, `organic_diagnose`, `organic_split_energy`, `organic_take_water`, `organic_spill`, `organic_take_exchange` |
| the balance and the ground surface's blend | `ed_surface`: `organic_surface_balance`, `ground_surface_exchange` |
| the stage: the pond and the inflow fill it | `ed_fast_integ_frozen`: `build_frozen`, `ground_inflow` |
| the soil's vapour through it | `ed_soil_water`: `organic_cover_path`, `ground_evaporation_rate` |
| the litter follows the pools | `ed_slow_update`: `follow_litter_pools`; `ed_slow_dynamics`: `litter_follows_pools` |
| the start | `ed_init`: `seed_organic_layer` |

## References
- Chadburn, S. et al. (2015), *Geosci. Model Dev.* 8:1493 — the moss layer's porosity and conductivity in JULES.
- de Vries, D. A. (1963), in *Physics of Plant Environment*, North-Holland — thermal properties of soils.
- Farouki, O. T. (1986), *Thermal Properties of Soils*, Trans Tech — the heat capacity of organic matter.
- Letts, M. G. et al. (2000), *Atmosphere-Ocean* 38:141 — the porosity of organic soils.
- Napoly, A. et al. (2017), *Geosci. Model Dev.* 10:1621 — the ISBA-MEB forest-floor litter layer.
- Putuhena, W. M. & Cordery, I. (1996), *J. Hydrol.* 180:283 — the forest floor's interception capacity.
- McLaren, J. (2021), *Soil biogeochemical variables collected on the Arctic LTER experimental plots in
  moist acidic, moist non-acidic, wet shrub and shrub tundra, Toolik Field Station, Alaska 2015*, Arctic LTER,
  knb-lter-arc.20129.1 — the bulk density of the top 5 cm of the organic layer.
