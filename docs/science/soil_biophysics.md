# Soil biophysics

The soil is two coupled stores on one grid of layers, with $z$ negative below the surface as in ED2: a **water** store (mixed-form Richards' equation) and a **thermal** store (heat diffusion in internal energy). This page gives their equations, the layers, the hydraulic and thermal properties the soil maps provide, frozen soil, and the coupling of water and heat. The water is solved together with the pond and the plants' water as one network ([soil_plant_water](soil_plant_water.md)); this page is the soil's part.

This is one of the per-store pages under [coupled_biophysics](coupled_biophysics.md), which describes the two principles the stores share (**prognostic internal energy**, not temperature; **closed-budget discipline**) and how they are woven each `dt_fast`. Here the first makes freeze and thaw a read-off of each layer's freezing curve (§ Frozen soil), and the second makes both stores close a machine-precision budget every step.

---

## The soil water column

A patch's soil moisture $`\theta_k`$ and its ponded surface water are advanced over `dt_fast` by the water network's step (`water_step`: implicit Euler on adaptive sub-steps, Newton on every store at once; [soil_plant_water](soil_plant_water.md)). For the soil layers the network solves the mixed-form Richards equation. Per layer,

```math
C(\psi_k)\,\frac{\partial\psi_k}{\partial t}
= \frac{\partial}{\partial z}\!\left[K(\theta)\Big(\frac{\partial\psi}{\partial z}+1\Big)\right] - S_k
\qquad(1)
```

with matric potential $`\psi(\theta)`$, hydraulic conductivity $`K(\theta)`$ and specific moisture capacity $`C=d\theta/d\psi`$ from the **van Genuchten–Mualem** (default) or the **Campbell / Clapp–Hornberger** retention curve (`ed_soil_hydr_lib`). In effective saturation $`S_e=(\theta-\theta_{res})/(\theta_{sat}-\theta_{res})`$ and matric head $h<0$ [m], van Genuchten's curve is $`S_e=[1+(\alpha|h|)^n]^{-m}`$ with $`m=1-1/n`$. Campbell's is $`h=\psi_{sat}\,S_e^{-b}`$, saturated between its air-entry head $`\psi_{sat}<0`$ and zero, with $`K=K_{sat}\,S_e^{2b+3}`$. Both families keep their parameters in one pair, `[soil].curve_par_a` and `[soil].curve_par_n` (α [m⁻¹] and $n$, or $`\psi_{sat}`$ [m] and $b$), so every reader of a profile (the water network, the predawn seed, the diagnostics) evaluates the curve the profile was built with. The van Genuchten conductivity is Mualem's,

```math
K = K_{sat}\,S_e^{L}\,\big[1-(1-S_e^{1/m})^m\big]^2, \qquad L = -1 \qquad(2)
```

with pore connectivity $L=-1$ (Schaap & Leij 2000, `VG_PORE_CONNECTIVITY`), as ED2 uses with this curve. With Mualem's 0.5 a soil of $n$ near 1.4 would conduct less than ED2's floor of 1.16×10⁻¹³ m s⁻¹ from about −0.7 MPa, far wetter than measured (UNSODA); $L=-1$ keeps it above that floor to about −2 MPa. There is no film-flow term for drier soil. The conductivity at a layer face follows ED2's rule (`rk4_derivs`): $\ln K$ interpolated linearly from the two nodes to the face, a thickness-weighted geometric mean,

```math
K_{k+1/2} = K_k^{\,1-w}\,K_{k+1}^{\,w}, \qquad w = \frac{\Delta z_k}{\Delta z_k+\Delta z_{k+1}} \qquad(3)
```

since the face lies $`\Delta z_k/2`$ below node $k$ and $`\Delta z_{k+1}/2`$ above node $k+1$; equal layers give the plain geometric mean. Across a wet layer over a dry one the face conductivity is small, so a wetting front, or capillary rise into a dry profile, advances more slowly than it would under an upstream pick.

Each layer is solved in its effective saturation on the dry side and its potential on the wet side. It stores water above saturation through a specific storage instead of being clipped, and as it dries its potential falls without bound towards $`\theta_{res}`$ and the conductivity of its faces falls with it, so no flow can empty it ([soil_plant_water](soil_plant_water.md) §2).

The **root sink** is the network's own: each cohort's uptake from each layer it roots in, solved with the layer's water. The soil has no wilting point of its own: the soil-to-root conductance falls with $`K(\theta)`$ and a layer drier than the wood gives nothing, so a threshold here would cut the same water twice.

**Ground evaporation** (`ground_evaporation_rate`, at the end-of-sub-step state) is Philip's pore-space relative humidity of the top layer through a Swenson–Lawrence dry-surface-layer (DSL) resistance, in series with the aerodynamic resistance $`r_{aero}=1/g_{g,net}`$ from the aerodynamics kernel. As the top layer dries, Philip's humidity goes to zero and evaporation stops by itself. The DSL resistance is $`r_{soil} = (d_{max}/D_{dsl})\,\max(\theta_{init}-\theta_1,0)/(\theta_{init}-\theta_{air})`$ with the diffusivity through the layer $`D_{dsl} = D_v\,\phi_{air}^{10/3}/\phi^2`$, where $`d_{max}`$ is `[soil].dsl_dmax`, $`\theta_{init}`$ is `[soil].dsl_theta_init` times the saturated water content, $`\theta_{air}`$ is the air-dry water content at a matric head of −10⁴ m (CLM5's), $`\phi`$ is the porosity, $`\phi_{air}=\phi-\theta_{air}`$, and $`D_v`$ is the vapour diffusivity in free air. The layer's thickness and tortuosity are referenced to the air-dry water content: referenced to the bulk $`\theta_1`$ instead, the tortuosity would fall faster than the layer thickens and a wetter soil would evaporate less. The tortuosity is Millington–Quirk's rather than CLM5's form in the Clapp–Hornberger $b$; with $`\phi_{air}`$ fixed the two differ by a constant factor, so the moisture response has the same shape. This flux alone sets the ground latent flux: it drives both the canopy air's humidity and the ground energy balance's latent heat, so no water is counted twice.

### Boundary conditions, and what the aquifer option actually means

The **top boundary** is the pond, a store of its own. Rain enters it, and it exchanges water with the top layer by the Darcy flux

```math
q_{pond} = \sqrt{K(\theta_1)\,K_{sat}}\;\Big(1 + \frac{h-\psi_1}{-z_1}\Big) \qquad(4)
```

from the pond, treated as a saturated layer as thick as the top one, to the top node at elevation $`z_1<0`$. This is eq. (3) with $w = 1/2$; $h$ is the pond's head and $`\psi_1`$ is in metres here. The flux is solved implicitly and is two-way: a top layer pushed past saturation exfiltrates to the pond. The pond holds up to `[soil].w_pond_max`; its first millimetre is a depression store, where its head grows as the square root of what it holds ([soil_plant_water](soil_plant_water.md) §2). Above `w_pond_max`, water runs off at the pond's temperature.

The geometric mean with $`K_{sat}`$ is what lets the pond re-wet a dried surface. Near residual water content the top layer's own $`K(\theta_1)`$ is about $10^{-15}$ of $`K_{sat}`$ for a van Genuchten loam, so a face conductivity of $`K(\theta_1)`$ alone would seal the surface: after a dry season the pond would overflow and the soil would never re-wet.

The **bottom** either drains or sits on an aquifer. Drainage is the unit-gradient flux scaled by a drainage factor $f \in [0, 1]$, $`q = f\,K(\theta_n)`$, always **downward**. Free drainage, a leaky base and a sealed one are values of $f$ (`drain_factor` in the code), not separate cases:

| option | flux at the base | use for |
|---|---|---|
| `free_drain` (default) | drainage, $f = 1$: $`q = K(\theta_n)`$ | a deep, unseen water table |
| leaky (soil maps only) | drainage, $`f = f_{leak}`$, the map's | soil on weathered rock that drains slowly |
| `bedrock` | drainage, $f = 0$: $`q = 0`$ | an impermeable base |
| `aquifer` | $`q = \sqrt{K(\theta_n)\,K_{sat}}\,(\psi_n/\Delta + 1)`$, $`\Delta = \Delta z_n/2`$, **two-way** | a shallow water table at the column base |

`[soil].water_bottom_bc` names `free_drain`, `bedrock` or `aquifer`. A leaky base comes only from the soil maps (§ The column from the soil maps), which give each cell its own $`f_{leak}`$.

**`aquifer` is a boundary condition, not a store.** The water table is *defined* to sit at the column base, so the modelled soil column **is** the unsaturated zone above it. There is no aquifer bucket, no baseflow and no water-table depth, so there is no Dunne saturation-excess runoff and no Zeng–Decker equilibrium correction. The face between the deepest node and the table is eq. (3) with the table as a saturated node half a layer below: the geometric mean of $`K(\theta_n)`$ and $`K_{sat}`$, as at the pond face and as ED2 joins its aquifer, in both directions. Water rising from the table crosses a half-layer whose conductivity runs from $`K(\theta_n)`$ to $`K_{sat}`$, so taking $`K_{sat}`$ alone would be the upper bound: a 0.5 m bottom layer at −2 m head would draw 7 $`K_{sat}`$ however dry it is, where the mean gives 0.2–0.4 $`K_{sat}`$. The mean is also smooth through zero gradient, so the network's Newton meets no corner there.

Two consequences are easy to get wrong:

- The flux **reverses upward** once $`|\psi_n| > \Delta`$, that is once the deepest layer is dry enough to pull against the saturated zone. So `aquifer` is the **wet-site** boundary, right for riparian, floodplain and wetland columns and wrong for an upland one, where it supplies water a real hillslope would not.
- The head term is what makes `aquifer` differ from `free_drain`: the deep-water-table limit $`q = K(\theta_n)`$ is free drainage, whatever store is put behind it.

Every step closes a machine-precision water budget, whose residual is reported in the water ledger (output `mean_resid_water_py` for the whole column):

```math
\mathrm{resid}_W = \Delta W_{stores}
 - \Delta t\,(P - E_{soil} - \text{drainage} - \text{uptake} - \text{runoff}) \approx 0 \qquad(5)
```

The step also exports the time-mean Darcy flux at each face and each layer's uptake, which carry water enthalpy into the soil-energy solve (§ Water–thermal coupling).

### The soil layers

Water and heat share one grid. Every profile cuts one list of layer bottoms at its own depth (`soil_layer_grid`): SoilGrids 2.0's depths to 1 m (5, 15, 30, 60 and 100 cm), 150 and 200 cm, then 3, 4.5, 6.5, 9, 12, 16, 21 and 27 m, each layer at most twice as thick as the one above. The last layer ends at the profile depth, and a remainder thinner than half of the list layer it falls in joins the layer above, so no profile ends in a sliver.

- **The top layer is 5 cm in every profile.** That is thin enough for the diurnal heat wave, whose damping depth is about 0.1 m. The drying surface needs no thinner layer: the dry-surface-layer resistance stands in for it.
- **A shallow soil gets a few layers rather than thin ones:** 0.3 m is 5, 10 and 15 cm.
- **Soil maps line up with the layers,** because every profile is cut from the same list.

`[soil].layer_bottoms` replaces the list for a site with a measured profile, and the profile then ends at its last bottom unless `[soil].depth` cuts it.

### Hydraulic parameters from soil texture

`soil_texture_to_hydraulics` (`ed_soil_hydr_lib`) gives a layer's van Genuchten curve and saturated conductivity from the quantities SoilGrids 2.0 maps per layer: sand $Sa$, silt $Si$ and clay $Cl$ (% of the fine earth, < 2 mm), organic carbon $OC$ (%), cation exchange capacity CEC (cmol(c) kg⁻¹), pH in water, the fine earth's bulk density $\rho_b$ (g cm⁻³) and the volume share of coarse fragments $f_v$. The curve is Hodnett & Tomasella's (2002, Table 8), fitted to tropical soils and kept by ED2 as `soil_hydro_scheme` 2; each polynomial is in percent:

```math
\begin{aligned}
\ln\alpha &= \big(-2.294 - 3.526\,Si + 2.440\,OC - 0.076\,\mathrm{CEC} - 11.331\,\mathrm{pH} + 0.019\,Si^2\big)/100\\
\ln n &= \big(62.986 - 0.833\,Cl - 0.529\,OC + 0.593\,\mathrm{pH} + 0.007\,Cl^2 - 0.014\,Sa\,Si\big)/100\\
\theta_{sat} &= \big(81.799 + 0.099\,Cl - 31.42\,\rho_b + 0.018\,\mathrm{CEC} + 0.451\,\mathrm{pH} - 0.0005\,Sa\,Cl\big)/100\\
\theta_{res} &= \max\Big[0,\ \big(22.733 - 0.164\,Sa + 0.235\,\mathrm{CEC} - 0.831\,\mathrm{pH} + 0.0018\,Cl^2 + 0.0026\,Sa\,Cl\big)/100\Big]
\end{aligned}
\qquad(6)
```

with $\alpha$ per kPa, converted to per metre of water head. The saturated conductivity of the fine earth follows Ottoni et al. (2019) from the water that drains between saturation and −33 kPa:

```math
K_{sat} = 19.31\ \mathrm{m\,d^{-1}}\,\big(\theta_{sat} - \theta(-33\ \mathrm{kPa})\big)^{1.948} \qquad(7)
```

Coarse fragments hold no water, so the layer's water contents are the fine earth's scaled by $1-f_v$. They slow the flow by Saxton & Rawls (2006, eqs. 19 and 22), written for their share of the mass $f_m$, with $\rho_p = 2.65$ g cm⁻³ the density of the grains and of the fragments:

```math
\theta_{sat,res} \leftarrow (1-f_v)\,\theta_{sat,res}, \qquad
K_{sat} \leftarrow K_{sat}\,\frac{1-f_m}{1-f_m\,(1-\tfrac{3}{2}\rho_b/\rho_p)}, \qquad
f_m = \frac{f_v}{(\rho_b/\rho_p)(1-f_v)+f_v} \qquad(8)
```

**Why this function.** Against 12,297 measured retention curves (GSHP; Gupta et al. 2022) it gives the smallest error in plant-available water ($\theta$ from −33 to −1500 kPa) of the functions tried (Cosby 1984, Saxton & Rawls 2006, Tomasella & Hodnett 1998 and Rosetta3): a median of 0.026 m³ m⁻³ from measured texture and 0.035 from SoilGrids at 0.1°, and 0.023 in tropical clays, where Rosetta3 overstates it by 0.07–0.09 and Cosby by 0.04. Its $K_{sat}$ is among the least biased against SoilKsatDB (Gupta et al. 2021) in the tropics, and at Barro Colorado Island (BCI) its curve follows the flat measured dry branch (Kupers et al. 2019). It is closed-form, so it costs nothing at run time.

**Where it is weak.** In boreal, organic-rich and acid soils, high CEC and low pH raise $`\theta_{res}`$, so the available water comes out about 0.04 m³ m⁻³ low (Hyytiälä's topsoil: $\theta$ at −1500 kPa 0.39 against 0.07 measured). Like every texture-based function, $K_{sat}$ scatters by about a factor of ten and runs low where soil structure adds large pores.

### The column from the soil maps

With the `[maps]` soil maps (ED-SoilMaps v2; [configuration](../configuration.md), "Input maps"), `apply_map_inputs` (`ed_config`) builds each cell's profile from what `ed_maps` reads. A region's polygons and a site on the ED_ERA5land grid build it the same way.

- **The layers.** The profile ends at the map's `soil_depth`, on the layers of § The soil layers.
- **Each layer's soil.** A layer takes the texture and chemistry of the map layers it overlaps (SoilGrids' six, to 2 m), each weighted by the overlap. Below 2 m it takes the deepest map layer's.
- **Hydraulics and heat, per layer.** Each layer gets its own curve and $K_{sat}$ by eqs. (6)–(8). Its thermal properties come from the same texture by CLM5's mineral-soil formulas (Lawrence & Slater 2008), with $\theta_{sat}$ the layer's porosity and silt left out as in CLM5:

  ```math
  k_{solid} = \frac{8.80\,Sa + 2.92\,Cl}{Sa + Cl}, \qquad
  k_{dry} = \frac{0.135\,\rho_d + 64.7}{2700 - 0.947\,\rho_d}, \qquad
  c_{dry} = \frac{2.128\,Sa + 2.385\,Cl}{Sa + Cl}\,10^6\,(1-\theta_{sat}) \qquad(9)
  ```

  Here $`\rho_d = 2700\,(1-\theta_{sat})`$ kg m⁻³ is the dry bulk density; $k$ is in W m⁻¹ K⁻¹ and $c_{dry}$, the solids' heat capacity per m³ of layer, is in J m⁻³ K⁻¹ (the water's own heat capacity is counted separately).
- **The start.** Each layer starts at its own field capacity, $\theta(\psi_{fc})$.
- **The bottom water boundary** is the map's `bc_class`:
  - free drain, drainage at $f = 1$;
  - leaky, drainage at the map's `f_leak` (a leaky cell whose `f_leak` is outside (0, 1] stops the run);
  - sealed, drainage at $f = 0$ (bedrock);
  - aquifer.
- **The thermal anchor.** With `[soil].thermal_bottom_bc = "dirichlet"`, the bottom conducts to the map's `deep_temp` (in ED-SoilMaps v2, the cell's mean 2 m air temperature plus 1 K), at $d/\sqrt2$ below the bottom node, $d = 1.973$ m (§ The bottom thermal boundary).

A cell whose maps hold no value stops the run, naming the cell. The map products are made complete in their own pipelines; ED-MORE does not fill them.

Each layer has its own curve rather than the profile one fitted curve: the top layers hold the surface's organic, high-CEC soil and the deep layers the subsoil.

---

## The soil thermal column

`soil_energy_step_implicit` (`ed_soil_heat`) advances each layer's prognostic volumetric internal energy $`E_k`$ [J m⁻³] by backward Euler, on the same layers as the water. The balance of layer $k$ over the step is

```math
\Delta z_k\,\frac{E_k^{n+1}-E_k^n}{\Delta t}
= (hf_k-hf_{k-1}) + (qwf_k-qwf_{k-1}) + \Delta z_k\,q_{src,k} \qquad(10)
```

with the faces at the temperatures $`T_k(E_k^{n+1})`$ that the new energies give on the layers' freezing curves (§ Frozen soil):
- the conduction $`hf_k=-\kappa_{k+1/2}\,(T_k-T_{k+1})/\Delta z_{k+1/2}`$ (positive up), whose face conductivity is the thickness-weighted harmonic mean of the layers' ice-aware $`\kappa(\theta,f_{liq})`$ (eq. (11)), held at the step's start;
- the water's enthalpy $`qwf_k`$, upwind on the giving layer's temperature (§ Water–thermal coupling);
- a top Neumann flux $`G_{top}`$, and the bottom boundary of § The bottom thermal boundary;
- the sources $`q_{src,k}`$, the enthalpy that leaves each layer with root uptake, exfiltration and drainage.

**Newton on the energies.** The step's Jacobian is tridiagonal: the faces' slopes in the temperatures, divided by each layer's apparent heat capacity $`dE/dT`$, plus $`\Delta z_k/\Delta t`$ on the diagonal. The apparent heat capacity is the sensible heat plus the latent heat of the liquid the freezing curve frees per kelvin, so the latent heat is inside each step's balance. A step that does not lower the worst residual is halved. The iteration stops when the worst residual, as the temperature change it would need, is under 10⁻⁸ K. A thawed column's balance is linear, so it takes one iteration; a column freezing from the top takes a few.

The energies are then **committed** from the faces at the solution, $`E_k^{n+1} = E_k^n + \Delta t\,[(hf_k-hf_{k-1}) + (qwf_k-qwf_{k-1})]/\Delta z_k + \Delta t\,q_{src,k}`$, so the budget closes to round-off whatever Newton's residual. Temperature and liquid fraction are re-diagnosed from the committed energy. The closed residual is `energy_resid` $`=\Delta E - \Delta t\,(G_{top}-\text{bottom}-\sum\text{water\_heat\_out}+qwf_n-qwf_0)\approx0`$, where `water_heat_out` is the enthalpy leaving each layer with root uptake, exfiltration and drainage, and the infiltrating water's enthalpy $`qwf_0`$ is among the inputs.

A sibling, `soil_energy_time_deriv`, exposes the same flux divergence as an explicit tendency (faces at $`T^n`$) for the Cash–Karp RK45 integrator; ESDIRK2's implicit stages take the backward-Euler step.

**Why the energies.** A solve for $`T^{n+1}`$ with the sensible heat capacity at the step's start would take no latent heat into the step: a layer reaching freezing would overshoot in temperature, and the conservative commit would put the latent heat back only afterwards, so the freezing front would depend on $`\Delta t`$. With Newton on the energies, a five-layer column frozen from the top for two days gives the same profile with 30-minute and 1-minute steps (`test_soil_heat` requires agreement within 0.1 K).

The conductivity is Johansen's, made ice-aware (`soil_thermal_conductivity`, `ed_therm_lib`). With the relative saturation $`S_r=\theta/\theta_{sat}`$, floored at 0.05 where the logarithm leaves its range of validity, the Kersten number blends Johansen's unfrozen form with the linear frozen one by the liquid fraction $`f_{liq}`$, and the saturated conductivity is the geometric mean of the solids, the liquid and the ice:

```math
K_e = f_{liq}\,(\log_{10}S_r + 1) + (1-f_{liq})\,S_r,\quad 0\le K_e\le 1, \qquad
\kappa_{sat} = k_{solid}^{\,1-\theta_{sat}}\,k_w^{\,\theta_{sat}f_{liq}}\,k_i^{\,\theta_{sat}(1-f_{liq})}, \qquad
\kappa = K_e\,\kappa_{sat} + (1-K_e)\,k_{dry} \qquad(11)
```

$`k_{solid}`$ and $`k_{dry}`$ come from the texture (eq. (9)); $`k_w`$ and $`k_i`$ are the conductivities of liquid water and ice.

### The bottom thermal boundary

Two boundary conditions are selectable, through `[soil].thermal_bottom_bc`.

**`geothermal` (Neumann, the default).** A prescribed flux at the base, held at exactly zero, so the base is an **adiabatic wall**. The flux term is plumbed through the kernel (the soil heat step's `geothermal` input enters $`hf_n`$ and the ledger's `bottom_heat`), but nothing assigns it a non-zero value. Real continental geothermal flux is about 0.05–0.09 W m⁻², two to three orders below the diurnal $`G_{top}`$ signal, so neglecting that flux is defensible; the wall itself is the consequential part. A zero-flux base **reflects** the downward thermal wave instead of transmitting it, and with a 2.0 m column against an annual damping depth near 2 m the annual cycle has barely attenuated by the time it arrives.

**`dirichlet` (a deep temperature anchor).** The bottom node conducts to a plane held at
`deep_temp`, a distance $`\ell = \texttt{deep\_depth} - |z_{node,n}|`$ below it, through the bottom
layer's own conductivity:

```math
hf_n = q_{geo} - g_{deep}\,(T_n - T_{deep}), \qquad g_{deep} = \frac{\kappa_n}{\ell} \qquad(12)
```

with $`q_{geo}`$ the prescribed flux (0) and $`g_{deep} = 0`$ under the `geothermal` boundary, so one expression serves both. The anchor enters the backward-Euler **matrix** (it adds $`g_{deep}`$ to the bottom diagonal and $`g_{deep}T_{deep}`$ to the residual), not the right-hand side alone, so the solve stays unconditionally stable. The same `bottom_heat_face` expression serves the implicit step at $`T^{n+1}`$ and the explicit `soil_energy_time_deriv` at $`T^n`$, so the two paths cannot drift apart.

#### Why the anchor sits *below* the column, and how deep

Pinning the base **face** to a fixed temperature is not the answer: an adiabatic base reflects the annual wave with coefficient $`+1`$, and a pinned face reflects it just as hard with coefficient $`-1`$. Neither is transparent. A semi-infinite continuation of the soil presents an impedance $`\kappa(1+i)/d`$ to the column at its base, where $`d=\sqrt{2\alpha/\omega}`$ is the annual damping depth. A resistive link of length $`\ell`$ presents $`\kappa/\ell`$, and the reflection coefficient is smallest when the two magnitudes match:

```math
\frac{\kappa}{\ell} = \left|\frac{\kappa(1+i)}{d}\right| = \frac{\sqrt{2}\,\kappa}{d}
\quad\Longrightarrow\quad \ell = \frac{d}{\sqrt{2}} \qquad(13)
```

So the anchor **depth** is a physical choice, and the default is derived rather than fitted. For the default column (2 m on the default layers) in a mid-latitude mineral soil at $`\theta = 0.3`$, $`d = 1.973`$ m and the bottom node is at $`|z_{node,7}| = 1.75`$ m, giving `[soil].deep_depth` $`= 1.75 + 1.973/\sqrt{2} = 3.145`$ m.

The anchor is checked against the analytic profile in `test_soil_annual_damping`: a homogeneous 2 m column of ten layers thickening with depth (bottom node at 1.727 m, so the anchor at 3.12 m), driven by a harmonic $`G_{top}`$ to periodic steady state. This is the constant-$`\kappa`$, constant-$`C`$ case, where $`\exp(-z/d)`$ holds exactly. The test asserts that the adiabatic base over-amplifies the bottom node by more than 1.6 times and that the anchor is within 10 % of the analytic amplitude and at least five times closer over the whole profile; the values from that setup are:

| bottom BC | amplitude at $`-1.73`$ m, relative to the surface | vs analytic | RMS error over the profile |
|---|---|---|---|
| analytic (semi-infinite) | 0.421 | — | — |
| `geothermal` (adiabatic) | 0.764 | **+82 %** | 0.145 |
| `dirichlet`, `deep_depth = 3.12` | 0.412 | **−2 %** | 0.015 |

A sweep of `deep_depth` puts the measured optimum at 3.0–3.1 m, so the derivation is right to within one sweep step, and the minimum is flat: anything from 2.8 to 4.0 m is still three to seven times better than the adiabatic wall. The harness itself is sound: a 12 m column with the adiabatic base reproduces $`\exp(-z/d)`$ to 0.1 % over the top three damping depths, because there the boundary is too far away to matter.

**What this does not fix.** A purely resistive termination cannot reflect less than 0.41 in amplitude, whatever $`\ell`$ is: matching a complex impedance with a real one leaves the phase wrong by 45°. Closing that last gap needs heat *capacity* below the column, that is passive deep thermal layers under the hydrologically active column. The anchor buys roughly a factor of ten in this metric, not exactness.

**`deep_temp` is the user's choice, and the loader requires it.** It is the mean annual soil temperature below the damping depth, close to the mean annual air temperature of the driving forcing, and a site property like latitude. An error in it is a steady flux $`g_{deep}\,\Delta T`$ into the column base, so it biases deep-soil temperature in the annual mean. There is no default: a silent one would bring back, in a new place, the bias this boundary condition exists to remove.

## Frozen soil

Below the triple point a layer's water is shared between liquid and ice by the **soil freezing curve** (CLM5's supercooled liquid; Niu & Yang 2006). Liquid water and ice are in equilibrium in the pores when the liquid's potential equals the freezing-point-depression head (Clausius–Clapeyron):

```math
h_{frz}(T) = \frac{L_f\,(T - T_f)}{g\,T},\qquad T_f = t_{3ple} \qquad(14)
```

The layer's own retention curve at that head is the liquid it can hold,
$`\theta_{liq,max}(T)=\theta(h_{frz}(T))`$ (`soil_liquid_capacity`). A layer is:
- all liquid at and above $`T_f`$, and below it while $`\theta\le\theta_{liq,max}(T)`$;
- otherwise liquid $`\theta_{liq,max}(T)`$ and ice $`\theta-\theta_{liq,max}(T)`$.

At $`T_f`$ the curve holds $`\theta_{sat}`$. Water a layer stores beyond saturation (through the specific storage, [soil_plant_water](soil_plant_water.md) §2) therefore freezes at $`T_f`$ itself, on a plateau as narrow as that excess.

The curve treats ice as it would air: the soil "dries" as it freezes. A van Genuchten loam ($`\theta_{res}`$ 0.078, $`\alpha`$ 3.6 m⁻¹, $`n`$ 1.56) holds 0.090 m³ m⁻³ of liquid at −1 °C and 0.081 at −10 °C, a little more than its residual water, so a frozen soil is never quite dry. Freezing happens over a range of temperatures just below $`T_f`$, so a wet layer lingers there as the latent heat comes out (the zero curtain) instead of sitting exactly at $`T_f`$.

**Energy and temperature.** A layer's energy at temperature $`T`$ is

```math
E(T) = C_{dry}\,T + \rho_w\theta\,\big[f_{liq}\,u_{liq}(T) + (1-f_{liq})\,c_{ice}T\big],\qquad
f_{liq}=\min\!\big(1,\theta_{liq,max}(T)/\theta\big) \qquad(15)
```

with $`C_{dry}`$ the solids' heat capacity per m³, $`\rho_w`$ the density of water, $`u_{liq}(T)`$ the liquid's internal energy (`internal_energy_liquid`) and $`c_{ice}`$ the specific heat of ice. $`E`$ rises with $`T`$. Its slope, the apparent heat capacity, adds to the sensible heat the latent heat of the liquid freed per kelvin, $`\rho_w\,(u_{liq}-c_{ice}T)\,d\theta_{liq,max}/dT`$.

`soil_layer_temp` inverts it:
- the all-liquid temperature when the curve holds all the water;
- otherwise Newton, kept inside the bracket between the all-liquid temperature and the lower of the
  all-ice temperature and $`T_f`$, with bisection when a step leaves it.

The thermal step (§ The soil thermal column) and every diagnosis use this one inversion.

**Water.** Over a water step each layer's ice is held at its start, $`\theta_{ice}=(1-f_{liq})\,\theta`$
(`soil_ice_t`, built by `soil_ice_state`). The water network then reads the layer as follows:
- **The liquid's potential is the curve's at the liquid's content**, $`\psi(\theta-\theta_{ice})`$, on the full pore space. A frozen layer's liquid therefore sits at $`h_{frz}(T)`$, and it draws water from the wet soil beside it (cryosuction).
- **The ice impedes the flow:** each node's conductivity is the liquid's times $`10^{-6\,\theta_{ice}/\theta_{sat}}`$ (CLM5; Swenson et al. 2012). The face conductivity is the weighted geometric mean of its two nodes' (eq. (3)), so a face takes 10 to the mean of its nodes' ice fractions, which is CLM5's face impedance. The pond face and the aquifer face lie inside the top and bottom layers, so each takes its layer's impedance in full: the pond conducts at $`K_{sat,1}`$ times the top layer's impedance, and the table at $`K_{sat,n}`$ times the bottom layer's.
- **Pores the ice leaves.** Once the liquid fills the pores the ice leaves, $`\theta_{liq}=\theta_{sat}-\theta_{ice}`$, the layer holds no more water. The liquid's head rises from the head the pores filled at (pressure in a gas-free frozen soil), and above zero the layer stores through the specific storage, as a thawed layer does. Cryosuction therefore fills a frozen layer to saturation and no further.
- **Evaporation** reads the liquid. Philip's humidity is taken on the liquid's potential, with saturation over liquid: on the freezing curve the liquid's vapour pressure is the ice's, so blending in the ice's saturation as well would count the ice twice. The dry surface layer reads the liquid plus `[soil].dsl_theta_init` of the ice, so it forms below `dsl_theta_init` of the pores the ice leaves, which is CLM5's effective porosity.

The energy step that follows sees the new water through $`\theta`$, so water drawn into a frozen layer freezes there at the next diagnosis, and its latent heat enters the layer's energy balance.

## Water–thermal coupling

The water and thermal stores are coupled by construction: the thermal step reads the just-updated soil moisture (for $`\kappa`$ and the freezing curve), and the water the network moves carries liquid enthalpy.
- Across each face, water carries the enthalpy of the node it leaves.
- Root uptake takes enthalpy out of each layer at that layer's own temperature (`root_heat`, per layer). Exfiltration to the pond leaves the top layer in the same way (`water_heat_out`).
- Infiltration crosses the top face at the pond's temperature.
- Drainage has no face term: its enthalpy leaves with the bottom layer's sink, `water_heat_out`, beside the root uptake.
- Soil evaporation is not a face either: its water leaves as vapour, and its enthalpy (liquid plus latent) has already left in the ground's latent flux inside $`G_{top}`$.

Every face, the top one included, is applied by one upwind rule at the time level of the conduction it sits beside ($`T^{n+1}`$ in the step, $`T^n`$ in the tendency). In the supercooled-liquid datum each face term is large, about 1 MJ per kg of water moved, so the heat the moving water delivers is a small difference of large numbers, and a face valued at another temperature would swamp it.

## Canopy interception

The canopy's surface water, whose throughfall and drip are the rain that reaches the pond, is described with the tissue it wets in [vegetation_energy_dynamics](vegetation_energy_dynamics.md) ("The canopy's surface water").

## Soil heterotrophic respiration

The heterotrophic soil CO₂ flux $`R_h`$ comes from one formula, the CENTURY transfer matrix (`ed_respiration`: `heterotrophic_respiration`, reached through `patch_heterotrophic_respiration`). It is a soil-carbon decomposition process (the pools' slow step is `ed_soil_carbon` in `src/slow_dynamics/processes/`; [soil_carbon](soil_carbon.md)), but its flux enters the fast loop as a source of the canopy air's CO₂ ([canopy_air_space_biophysics](canopy_air_space_biophysics.md)). It is computed once per fast step.

---

## Prognostic state

| Store | Type | Prognostic variable(s) | Diagnosed |
|---|---|---|---|
| Soil water | `soil_water_state_t` | `theta(:)`, `w_surface`, `w_surface_energy` | none |
| Soil thermal | `soil_energy_state_t` | `soil_energy(:)` [J m⁻³] | `soil_temp`, `soil_fliq` |

## Parameters

The `[soil]` block ([configuration](../configuration.md) lists every key, and the soil maps replace the column keys). Every key is optional.

| Symbol | Key | Default | Meaning |
|---|---|---|---|
| | `depth` | 2.0 m | profile depth |
| | `layer_bottoms` | the template | measured layer bottoms [m, positive down] |
| | `retention` | `"van_genuchten"` | `"van_genuchten"` or `"campbell"` |
| $\theta_{sat}$, $\theta_{res}$ | `theta_sat`, `theta_res` | 0.43, 0.078 m³ m⁻³ | porosity, residual water content |
| $K_{sat}$ | `ksat` | 2.89×10⁻⁶ m s⁻¹ | saturated conductivity |
| α or $\psi_{sat}$; $n$ or $b$ | `curve_par_a`, `curve_par_n` | 3.6 m⁻¹, 1.56 | van Genuchten α and $n$, or Campbell $\psi_{sat}$ (< 0, m) and $b$ |
| $\psi_{fc}$ | `psi_fc` | −3.37 m | field-capacity matric head; sets each layer's field capacity |
| $k_{solid}$, $k_{dry}$, $C_{dry}$ | `solid_conductivity`, `dry_conductivity`, `dry_heat_capacity` | 3.0, 0.15 W m⁻¹ K⁻¹; 2.0×10⁶ J m⁻³ K⁻¹ | thermal texture, eqs. (9)–(11) |
| | `w_pond_max` | 5.0 kg m⁻² | ponding capacity before runoff |
| $d_{max}$ | `dsl_dmax` | 0.015 m | greatest thickness of the dry surface layer |
| | `dsl_theta_init` | 0.8 | top-layer $\theta/\theta_{sat}$ below which the dry surface layer forms |
| | `water_bottom_bc` | `"free_drain"` | `"free_drain"`, `"bedrock"` or `"aquifer"` |
| $f$ | | 1.0 | drainage factor: 1 for `"free_drain"`, 0 for `"bedrock"`, the map's `f_leak` for a leaky base; no key of its own |
| | `thermal_bottom_bc` | `"geothermal"` | `"geothermal"` or `"dirichlet"` |
| $T_{deep}$ | `deep_temp` | none | anchor temperature [K]; required by `"dirichlet"` |
| | `deep_depth` | 3.145 m | anchor depth, eq. (13) |

Fixed in the code: the Mualem pore connectivity $L=-1$, the specific storage 10⁻⁴ m⁻¹, the ice-impedance exponent 6, the air-dry head −10⁴ m, the coefficients of eqs. (6)–(8) and the pond's 1 mm depression store.

## Code map

| Concept | Routine |
|---|---|
| soil water (the network) | `ed_fast_integ_backward_euler`: `water_step` ([soil_plant_water](soil_plant_water.md)); a layer's water relations: `ed_soil_hydr_lib`: `soil_water_relations` |
| soil water flows, eqs. (3)–(4) | `ed_soil_water`: `soil_water_flux`, `bottom_flux`, `ground_evaporation_rate` (its per-step terms from `ground_evaporation_setup`); the water reaching the pond, `pond_inflow` |
| soil layers and profile | `ed_soil_derive`: `SOIL_LAYER_TEMPLATE`, `soil_layer_grid`, `soil_layers`, `build_soil_profile`; the types `soil_params_t` (`[soil]`) and `soil_profile_t` (what the solvers read), in `ed_soil_params` |
| retention curves (van Genuchten, Campbell), eq. (2) | `ed_soil_hydr_lib`, with the matric potential `soil_matric_potential` [MPa] |
| hydraulic parameters from texture, eqs. (6)–(8) | `ed_soil_hydr_lib`: `soil_texture_to_hydraulics` |
| thermal properties from texture, eqs. (9) and (11) | `ed_therm_lib`: `soil_texture_to_thermal`, `soil_thermal_conductivity` |
| soil thermal step (backward Euler, Newton on the energies), eq. (10) | `ed_soil_heat`: `soil_energy_step_implicit` |
| the soil's thermal forcing from the water network's flows | `ed_soil_heat`: `set_soil_heat_env` |
| a layer's temperature from its energy, and back, eq. (15) | `ed_soil_heat`: `soil_layer_temp`, `soil_layer_energy`, `soil_layer_energy_slope` |
| the soil freezing curve, eq. (14) | `ed_soil_hydr_lib`: `soil_liquid_capacity` |
| a layer's ice in the water network | `ed_soil_hydr_lib`: `soil_ice_t`, `soil_ice_state` |
| explicit thermal tendency (RK45) | `ed_soil_heat`: `soil_energy_time_deriv` |
| bottom thermal boundary, eqs. (12)–(13) | `ed_soil_heat`: `bottom_heat_face`, `deep_conductance` |
| the profile from the soil maps | `ed_maps`: `read_map_inputs`; `ed_config`: `apply_map_inputs`; `ed_soil_derive`: `map_layer_texture` |
| soil heterotrophic respiration | `ed_respiration`: `heterotrophic_respiration`; `ed_fast_integ_frozen`: `patch_heterotrophic_respiration` |
| prognostic soil types | `ed_patch_state_vars`: `soil_water_state_t`, `soil_energy_state_t` |

## References
- ED2 (`ED/src/dynamics/`): `lsm_hyd.f90`, `soil_respiration.f90` — the soil stores in the reference model.
- Celia, Bouloutas & Zarba (1990), *Water Resour. Res.* 26:1483 — mixed-form Richards.
- Brenner & Cancès (2017), *Adv. Water Resour.* — variable switching for Richards' equation.
- Mualem (1976), *Water Resour. Res.* 12:513; Schaap & Leij (2000), *Soil Sci. Soc. Am. J.* 64:843 — the conductivity and its pore connectivity.
- van Genuchten (1980); Clapp & Hornberger (1978) — soil retention curves.
- Hodnett & Tomasella (2002), *Geoderma* 108:155 — van Genuchten parameters from texture for tropical soils.
- Ottoni et al. (2019), *J. Hydrol.* 575:1345 — saturated conductivity from the drainable porosity.
- Saxton & Rawls (2006), *Soil Sci. Soc. Am. J.* 70:1569 — coarse fragments.
- Gupta et al. (2021), *Earth Syst. Sci. Data* 13:1593 (SoilKsatDB) and (2022), *Sci. Data* 9:444 (GSHP) — the measurements the function was chosen on.
- Swenson & Lawrence (2014), *JGR-Atmos.* 119:10299 — dry-surface-layer soil evaporation.
- Lawrence & Slater (2008), *Clim. Dyn.* 30:145 — the mineral and dry-soil thermal properties CLM uses, eq. (9).
- Niu & Yang (2006), *J. Hydrometeorol.* 7:937 — the supercooled liquid of a frozen soil, eq. (14).
- Swenson, Lawrence & Lee (2012), *J. Adv. Model. Earth Syst.* 4:M08002 — the ice impedance on the conductivity.
- Lawrence et al. (2018), *Technical Description of version 5.0 of the Community Land Model (CLM)*, NCAR — the supercooled soil water and frozen-soil hydrology.
