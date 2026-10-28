# Snow biophysics

The snow pack is a single bulk layer of water-equivalent mass and internal energy, stacked between the canopy air space and the top soil layer. This page gives how it accumulates, drains, exchanges energy and vapour with the air, and covers the ground. Where the pack exists, its energy balance replaces the bare-ground skin and it becomes the soil's top thermal boundary, blended by a cover fraction so the transition between snow-on and snow-off is continuous.

This is one of the per-store pages under [coupled_biophysics](coupled_biophysics.md), which describes the two principles the stores share (prognostic internal energy; closed-budget discipline) and how they are woven each `dt_fast`.

---

## The snow store

The snow kernels live in `ed_surface`; the pack is a single bulk layer. The prognostic state is the water-equivalent mass `swe` [kg m⁻²], the **extensive** internal energy [J m⁻²] and the depth (`swe` over `rho_snow`). Temperature and liquid fraction are read-offs of `internal_energy_to_temp` (with no dry heat capacity), so **melt and refreeze are the internal-energy plateau**, exactly as for soil. Each `dt_fast` runs **accumulate**, then **drain**, before the fast loop's march; the **exchange** happens inside the march and the **commit** after it:

- **Accumulate** (`snow_accumulate`). Snowfall lands as ice at $`\min(T_{3ple},T_{air})`$ and rain-on-snow as liquid at $`T_{air}`$ (it refreezes later through the temperature inversion). A pack is created only from `min_new_snow_mass` of new snow in one `dt_fast`; sub-threshold snow joins the rain as liquid (the canopy films, then the pond).
- **Drain** (`snow_drain_meltwater`). The free liquid the pack holds above its holding capacity (`liquid_holding_frac`), from the last step's melt and this step's rain, percolates out as a paired (mass, enthalpy) hand-off to the pond. It reaches the pond at the melt temperature and infiltrates from there in this step's water network. Full melt-out dumps the residual and reverts to bare ground. Meltwater therefore leaves the pack the step after it melts.
- **Exchange** (`snow_surface_balance`, called through `ground_surface_exchange` at every evaluation of the march). A linearized backward-Euler balance of the pack's enthalpy over the step, against the canopy air and the top soil node as the march holds them at that evaluation, as each tissue's is. Its terms are the shortwave the snow absorbs at its own albedo, plus the longwave it absorbs less its own emission, less sensible heat, less sublimation or evaporation, less conduction at the base. Each is linearized about the pack's temperature $`T^n`$, so the net flux $`R(T)`$ is linear and falls as $`T`$ rises. The pack's enthalpy $`E(T)`$ rises with $`T`$ along the ice line, across the melt plateau and along the liquid line, so the balance $`E(T) = E^n + \Delta t\,f_{snow}\,R(T)`$ has one root. It is found in closed form on the plateau when the plateau's enthalpy range holds it, else on the ice or the liquid line, and the fluxes are taken at that root. A thin pack under an open patch's bare-ground conductance therefore relaxes toward its balance and cannot overshoot it, and a melting pack under cold air freezes through the plateau and stops where the balance does.
  - Sublimation takes its latent heat from the balance, and the mass that leaves takes the pack's specific energy with it, so the rest keeps its temperature. A pack that would sublimate whole does so, the vapour taking all the energy left.
  - The canopy air's temperature and humidity are the march's own at each evaluation, so the pack's sensible, latent and vapour exchange is as implicit in them as the canopy air's other sources: a deposition from warm saturated air cannot take more vapour than the canopy air holds.
  - The **snow-base to soil-top conductance** is the series resistance of the half snow layer and the top soil node. It falls as the pack deepens ($`k_{snow}\ll k_{soil}`$), which decouples the soil and caps the winter soil surface, and it is the soil's top boundary against the march's own top-node temperature. The top node's $`k_{soil}`$ is the soil heat solver's conductivity of layer 1 at the step's start (its water, ice fraction and texture), so a frozen or a dry top layer conducts as the soil says it does.
- **Commit** (`snow_take_exchange`). The pack takes the march's time integral of the exchange, the energy and vapour the canopy air and the soil received from it, so the pack, the canopy air and the soil close together in the whole-column ledgers. Its temperature is read off the new energy.

### The ground's radiation

The two-stream sees one ground, the bare soil and the pack blended by the cover fraction $`f`$ (`snow_blend_ground_optics`). Its albedo and emissivity are the cover-weighted means, and it emits at the temperature $`T_{eff}`$ at which the blended emissivity $`\bar\varepsilon`$ gives the two tiles' emission, $`\bar\varepsilon\,\sigma T_{eff}^4 = f\,\varepsilon_{snow}\sigma T_{snow}^4 + (1-f)\,\varepsilon_{soil}\sigma T_{soil}^4`$.

What reaches the ground in each band, $`S_b`$ and $`L^\downarrow`$, is then shared at each tile's own optics (`ground_tile_radiation`). Per unit of its own area the bare soil absorbs $`\sum_b S_b(1-\alpha_{soil,b}) + \varepsilon_{soil}(L^\downarrow-\sigma T_{soil}^4)`$, at the top soil node's temperature at the step's start. The snow absorbs $`\sum_b S_b(1-\alpha_{snow,b})`$ of shortwave, with $`\varepsilon_{snow}(L^\downarrow-\sigma T^4)`$ of longwave formed in its balance at its own temperature. Weighted by the cover, the two give the blended ground's absorption, so the canopy and the ground agree on the ground's radiation. The shortwave is split by albedo, not by area: on half-covered ground the snow's half absorbs far less than the soil's half, where an even split would give each the same.

### Sublimation saturates over ice

The vapour flux is $`E=g_{net}\rho\,[q_{sat}(T_s,f_{liq})-q_{CAS}]`$, with $`g_{net}`$ the ground-to-canopy-air conductance, $`\rho`$ the air density and $`q_{CAS}`$ the canopy air's specific humidity. The saturation humidity $`q_{sat}`$ is taken on the **ice** curve when the pack is frozen, through the saturation vapour pressure

```math
e_{sat}(T,f_{liq}) = f_{liq}\,e_l(T) + (1-f_{liq})\,e_i(T), \qquad
e_p(T) = e_{3ple}\left(\frac{T_{3ple}}{T}\right)^{(c_p-c_{p,v})/R_v}
         \exp\left[\frac{L_p(T_{3ple})}{R_v\,T_{3ple}} - \frac{L_p(T)}{R_v\,T}\right] \qquad(1)
```

with $`p`$ the phase: liquid ($`c_p=c_{p,l}`$, $`L_p=L_v`$) or ice ($`c_p=c_{p,i}`$, $`L_p=L_s`$), $`L_v(T)`$ and $`L_s(T)`$ the latent heats of [canopy air space biophysics](canopy_air_space_biophysics.md) (`latent_heat_vap`), $`R_v`$ the vapour's gas constant and $`e_{3ple}`$ the saturation vapour pressure at the triple point (`es_3ple`, 611.655 Pa; Ambaum 2020, eq. 14). This is Ambaum's (2020) closed form of the Clausius–Clapeyron equation for a latent heat that is linear in temperature, so the slope is **by definition**

```math
\frac{\partial e_p}{\partial T} = \frac{L_p(T)\,e_p(T)}{R_v\,T^2} \qquad(1b)
```

(`sat_vapor_pressure_temp_deriv`), with the same $`L_p`$ the energy budget uses: the second law holds between the vapour pressure and the heat a flux of vapour carries. The empirical Bolton and Magnus forms it replaces had no such link. Ice's latent heat at $`T_{3ple}`$ is the liquid's plus the heat of fusion, so both branches equal $`e_{3ple}`$ **exactly** at $`T_{3ple}`$, which is where a melting pack sits, and the blend is continuous in temperature *and* in $`f_{liq}`$. A pack that freezes or melts slides between the curves rather than stepping, which matters because the flux lands in the right-hand side an adaptive integrator steps through. The slope of the blend is taken at fixed $`f_{liq}`$.

$`f_{liq}`$ is the pack's own prognostic liquid fraction, not a temperature threshold. The energy side is ice-aware (removing `enthalpy_vapor` from an ice-referenced layer debits sublimation, vaporization plus fusion, automatically), so the vapour side must be too: the liquid curve alone would treat the surface as ice for *enthalpy* and as liquid for *vapour pressure*. Over ice $`e_{sat}`$ is **9 % lower at −10 °C, 18 % at −20 °C and 26 % at −30 °C**, so the liquid curve would overstate it, and with it the driving gradient, by 10 %, 22 % and 34 %. With the model's own heat capacities and latent heat, the ice form is within 0.01 % of Murphy & Koop (2005) at −10 °C and 0.1 % at −30 °C, and the liquid form within 0.04 % at −10 °C, 0.09 % at 20 °C and 0.24 % at 40 °C.

Ground evaporation from a frozen soil does not use the blend: its liquid sits on the soil freezing curve, whose potential puts the liquid's vapour pressure at the ice's, so Philip's humidity on that potential over liquid saturation already carries the ice ([soil_biophysics](soil_biophysics.md), "Frozen soil"). The blend is deliberately **not** used for dewpoint conversion or diagnostic vapour-pressure deficit: dewpoint is *defined* over liquid, so an ice branch there would mis-convert the forcing.

### Snow-cover fraction

The **Niu–Yang (2007) snow-cover fraction** is $`\mathrm{snowfac}=\tanh(\text{depth}/\text{scale})`$ (`snow_cover_fraction`), with scale = `ny07_a` × `[aerodynamics].soil_rough` × ($`\rho_{snow}`$/100 kg m⁻³)^`ny07_m`. The scale is Niu & Yang's ground roughness, the micro-relief a thin pack has to fill: 6.25 cm at the defaults, so 1 cm of snow covers 16 % of the ground and 10 cm covers 92 %. The cover fraction ramps the ground's optics and the blend of the soil's top boundary, and it area-weights all boundary exchange, so a thin patchy pack barely exchanges (stable and continuous as $`\mathrm{snowfac}\to0`$). With no snow the result is bit-identical to a run without the store.

Ground evaporation is such a tile flux: only the snow-free share evaporates, at the bare soil's own resistances, $`E=(1-\mathrm{snowfac})\,\rho\,(q_g-q_{air})/(r_{aero}+r_{soil})`$. Inflating $`r_{aero}`$ by $`1/(1-\mathrm{snowfac})`$ instead would differ whenever the dry-surface-layer resistance $`r_{soil}`$ is positive: it divides only the aerodynamic leg, so it over-predicts bare-soil evaporation at intermediate cover, most where a dry surface makes $`r_{soil}`$ dominant. The two agree only at no cover and full cover. The canopy's roughness sees the pack through its own burial fraction, the pack's depth over the canopy height, not through this cover fraction ([canopy_aerodynamics](canopy_aerodynamics.md) §1).

---

## Prognostic state

| Store | Type | Prognostic variable(s) | Diagnosed |
|---|---|---|---|
| Snow / surface water | `snow_state_t` | `swe`, `snow_energy` [J m⁻²], `snow_depth` | `snow_temp`, `snow_fliq` |

## Parameters

The `[snow]` block; every key is optional ([configuration](../configuration.md)).

| Symbol | Key | Default | Meaning |
|---|---|---|---|
| $`\rho_{snow}`$ | `rho_snow` | 250 kg m⁻³ | bulk density, fixed (depth = `swe` / `rho_snow`) |
| $`k_{snow}`$ | `k_snow` | 0.15 W m⁻¹ K⁻¹ | thermal conductivity |
| | `albedo_vis_fresh`, `albedo_nir_fresh` | 0.85, 0.65 | fresh-snow albedo, visible and near infrared |
| | `albedo_vis_aged`, `albedo_nir_aged` | 0.518, 0.435 | aged-snow albedo; the pack's own liquid fraction blends fresh to aged |
| $`\varepsilon_{snow}`$ | `snow_emiss` | 0.97 | longwave emissivity |
| | `liquid_holding_frac` | 0.10 | liquid fraction the pack holds before it percolates |
| | `min_new_snow_mass` | 10⁻³ kg m⁻² | snowfall in one `dt_fast` that starts a pack |
| | `tiny_snow_mass` | 10⁻⁶ kg m⁻² | pack mass below which it melts out into the soil |
| | `ny07_a`, `ny07_m` | 2.5, 1.0 | shape factor and density exponent of the snow-cover scale |

## Code map

| Concept | Routine |
|---|---|
| snow energy, base conductance | `ed_surface`: `snow_surface_balance`, `snow_take_exchange`, `snow_base_conductance` |
| the stage before the march, the exchange in it, the commit | `ed_surface`: `prepare_snow_stage` (called by `ed_fast_integ_frozen`: `build_frozen`), `ground_surface_exchange` (called by `ed_fast_integ_derivs`: `surface_derivs`); `ed_fast_integ_update`: `update_dynamic_state` |
| snow mass, cover, melt | `ed_surface`: `snow_accumulate`, `snow_cover_fraction`, `snow_drain_meltwater` |
| the ground's optics and radiation | `ed_surface`: `snow_blend_ground_optics`, `snow_albedos`, `ground_tile_radiation` |
| saturation vapour pressure, eq. (1) | `ed_therm_lib`: `sat_vapor_pressure`, `sat_specific_humidity` |
| prognostic snow type | `ed_patch_state_vars`: `snow_state_t` |

## References
- Ambaum (2020), *Q. J. R. Meteorol. Soc.* 146:4252, doi:10.1002/qj.3899 — saturation vapour pressure over water and ice.
- Murphy & Koop (2005), *Q. J. R. Meteorol. Soc.* 131:1539 — vapour pressure over ice.
- Niu & Yang (2007), *JGR* 112:D21101 — snow-cover fraction.
