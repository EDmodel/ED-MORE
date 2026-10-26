# Canopy aerodynamics

How ED-MORE computes the turbulent transfer that couples the canopy air space (CAS) to the free
atmosphere above and to the leaf, wood and ground surfaces inside it. The page gives the roughness
geometry, the Monin–Obukhov surface layer, the conductances of the canopy air's three variables, the
ground conductance, the in-canopy wind and the leaf and wood boundary layers.

The kernel is a pure function of the free-atmosphere forcing, the CAS state and the canopy geometry. It
produces every turbulence quantity the other fast-loop kernels take as input: the friction velocity
$`u_*`$; the scalar profile factors that set the atmosphere↔CAS conductances of enthalpy, humidity and
CO₂; the per-cohort in-canopy wind and leaf and wood boundary-layer conductances; and the ground↔CAS
conductance. Like ED2's `canopy_turbulence8`, it carries no integrated state and is recomputed from the
current CAS and vegetation state. The physics combines ED2 and CLM: CLM5's four-range Monin–Obukhov
surface layer, ED2's forced-plus-free-convection Nusselt boundary layers, ED2's per-cohort wind
extinction, and a CLM4-like ground conductance.

Symbols: $`u_*`$ = friction velocity; $\zeta = z/L$ = Monin–Obukhov stability parameter ($L$ = Obukhov
length); $`z_{0m}, z_{0h}, z_{0q}`$ = roughness lengths for momentum, heat and vapour; $d$ = displacement
height; $`\kappa`$ = von Karman constant; $`\Psi_m,\Psi_h`$ = integrated stability functions;
$`\mathrm{Rib}`$ = bulk Richardson number; temp1 and temp2 = the scalar transfer factors for heat and for
vapour and CO₂ (CLM's names; `heat_transfer_factor` and `vapour_transfer_factor` in `aero_out_t`).
Everything is per unit ground area. All routines are `pure` and write per-cohort outputs into
caller-owned arrays.

## 1. Roughness, displacement, and the reference geometry

Roughness and displacement scale with the canopy top height $h$, and a snow pack that buries the canopy
blends the roughness toward that of snow, by the burial fraction `snowfac`:

```math
z_{0m} = (1-\mathrm{snowfac})\,c_{z0}\,h + \mathrm{snowfac}\,z_{0,snow}, \qquad
d = c_d\,h\,(1-\mathrm{snowfac}) \qquad(1)
```

with $`c_{z0}=`$ `z0m_ratio` (0.13, ED2 `vh2vr`), $`c_d=`$ `d_ratio` (0.63, ED2 `vh2dh`) and
$`z_{0,snow}=`$ `snow_rough`, which also floors $`z_{0m}`$ (`canopy_roughness`). The burial fraction is
ED2's: the pack's depth over the canopy height, capped at 0.99 (`canopy_burial_fraction`). It is taken
from the pack at the start of the forcing window, and the roughness is fixed for the window. It is
distinct from the ground's snow-cover fraction ([snow_biophysics](snow_biophysics.md)). The height
above the displacement, $`z_{ldis}=\max(z_{ref}-d,\,2z_{0m})`$, is the surface-layer thickness the
Monin–Obukhov solve integrates over. The CAS depth is
$`\max(\text{min\_canopy\_depth},\,h+\text{canopy\_freeboard})`$, the tallest cohort plus a freeboard
([canopy_air_space_biophysics](canopy_air_space_biophysics.md)); it sets the CAS air mass per ground
area and hence the storage capacity of each of its three variables.

## 2. The CLM5 Monin-Obukhov surface layer

`mo_surface_layer` returns $`u_*`$ and the scalar profile factor from the gradients between the
atmosphere and the CAS. It starts from the analytic **MoninObukIni** bulk-Richardson guess and then
sweeps until the coupling $`u_*\cdot\mathrm{temp1}`$ moves by less than $`10^{-4}`$ of itself from one
sweep to the next, with at most `n_iter_mo` (20) sweeps. Near neutral two sweeps suffice. Under hot
canopy air and a weak wind the sweeps swing about the root, because $\zeta$ and the convective wind
feed each other, so the loop runs to convergence instead of a fixed count. With virtual potential
temperatures
$`\theta_v=\theta\,(1+0.60783\,q)`$ (the factor is $`M_{dry}/M_{h2o}-1`$, `virtual_temp_factor`) and the
CAS↔atmosphere differences $`\Delta\theta_v`$:

```math
\mathrm{Rib} = \frac{g\,z_{ldis}\,\Delta\theta_v}{\theta_{v,atm}\,u_m^2}, \qquad
u_m = \begin{cases}\bar u & \Delta\theta_v\ge0\ \text{(stable)}\\[2pt]
\sqrt{\bar u^2+w_c^2} & \Delta\theta_v<0\ \text{(convective)}\end{cases}, \qquad
\bar u = \max(u_{ref},\,u_{b,min}) \qquad(2)
```

($w_c$ = `wc`, the first guess's convective velocity; $`u_{b,min}`$ = `ubmin`, the reference-wind
floor). Each sweep evaluates the CLM5 momentum and scalar **profile denominators** $D_m(\zeta)$ and
$D_h(\zeta)$ over four stability ranges (very-unstable / unstable / stable / very-stable, `d_mom` and
`d_heat`), then updates $`u_*`$, the scalar factor, $\zeta$ and, as CLM5 does, the wind $u_m$ with a
convective velocity from the sweep's own buoyancy flux:

```math
u_* = \max\!\Big(u_{*,min},\,\frac{\kappa\,u_m}{D_m(\zeta)}\Big), \qquad
\mathrm{temp1} = \frac{\kappa}{D_h(\zeta)}, \qquad
\zeta = \frac{z_{ldis}\,\kappa\,g\,\theta_{v*}}{u_*^2\,\theta_{v,atm}} \qquad(3)
```

```math
u_m = \begin{cases}\bar u & \zeta\ge0\\[2pt]
\sqrt{\bar u^2+w_c^2},\ \ w_c = \beta\Big(\dfrac{-g\,u_*\,\theta_{v*}\,z_i}{\theta_{v,atm}}\Big)^{1/3} & \zeta<0\end{cases}
\qquad(4)
```

with $\beta = 1$ and the mixed-layer depth $`z_i`$ = 1000 m (CLM5's constants). As in CLM5, each sweep's
$\zeta$ is held in $`[-100,\ \zeta_{max,stable}]`$ (`zeta_max_stable`, 0.5). Without the upper bound a
stable night runs past $\zeta = 1$, where $`u_*`$ falls to $`u_{*,min}`$ and the coupling no longer
depends on the wind. CLM5 also keeps $|\zeta| \ge 0.01$ at each sweep; ED-MORE keeps it only in the first
guess, because nothing divides by $\zeta$ near neutral and the bound would put a jump in $`u_*`$ there.
The buoyancy flux $`-u_*\theta_{v*}`$ vanishes as the air turns neutral, so $w_c$ does too, and the
wind, with it $`u_*`$ and the conductances, passes through neutral without a jump. Only the first
guess, eq. (2), uses the fixed `wc`.

The unstable-range stability integrals are the standard Businger-Dyer / Paulson forms with
$`x=(1-16\zeta)^{1/4}`$:

```math
\Psi_m(\zeta) = 2\ln\frac{1+x}{2} + \ln\frac{1+x^2}{2} - 2\arctan x + \frac{\pi}{2}, \qquad
\Psi_h(\zeta) = 2\ln\frac{1+x^2}{2} \qquad(5)
```

(`stabfunc1` and `stabfunc2`; the stable branch is $`-5\zeta`$). The roughness term
$`\zeta_{0m}=z_{0m}\zeta/z_{ldis}=z_{0m}/L`$ is formed without dividing by $\zeta$, and $\zeta$ divides
only in the very unstable branches ($\zeta<-\zeta_m$ for momentum, $\zeta<-\zeta_t$ for scalars) and
the very stable branch ($\zeta>1$), so the solve is safe under floating-point trapping (`-fpe0`). The
neutral limit is $`u_*\to\kappa u/\ln(z_{ldis}/z_{0m})`$ and
$`\mathrm{temp1}\to\kappa/\ln(z_{ldis}/z_{0m})`$.

### The scalar transfer factors that couple the canopy air

The canopy uses $`z_{0h}=z_{0m}`$ and $`z_{0q}=z_{0h}`$, so the heat and vapour profile factors are equal:
$`\mathrm{temp2}=\mathrm{temp1}`$. These dimensionless Monin–Obukhov coefficients set the
**atmosphere↔CAS conductances** of the canopy air's three variables, the most important output of this
kernel:

```math
g_{ah} = \rho\,u_*\,\mathrm{temp1}, \qquad
g_{aw} = \rho\,u_*\,\mathrm{temp2}, \qquad
g_{ac} = \rho_{mol}\,u_*\,\mathrm{temp2} \qquad(6)
```

$g_{ah}$ [kg m⁻² s⁻¹] drives the enthalpy, $g_{aw}$ the humidity, and $g_{ac}$ [mol m⁻² s⁻¹]
the CO₂ (molar, using the dry-air molar density). The profile factor matters: dropping it
($`\mathrm{temp1}=1`$) would over-couple the CAS to the free atmosphere by about $`1/\mathrm{temp1}`$
and damp the nocturnal CO₂ build-up under the canopy. The kernel also reports the temperature scale
$`t_*=\mathrm{temp1}\,\Delta\theta`$, and $`\mathrm{Rib}`$, $\zeta$ and $L$ as diagnostics.

**The reference level.** $`z_{ref}`$ is the canopy-air top, where the forcing is moved to
([forcing](forcing.md) §8), and both potential temperatures are referenced to it: $`\theta_{atm}`$ is
the air temperature there and $`\theta_{CAS}`$ the canopy-air temperature, so
$`\Delta\theta = T_{atm} - T_{CAS}`$. That is the difference the canopy air's energy budget
exchanges, $`g_{ah}\,(h_{CAS}-h_{atm})`$ on actual temperatures, so the reported sensible heat
$`H = g_{ah}\,c_p\,(T_{CAS}-T_{atm})`$ is the flux the budget books. The two temperatures must share the
reference. Referencing $`\theta_{atm}`$ to the ground, $`T_{atm}+(g/c_p)\,z_{ref}`$, while $`\theta_{CAS}`$
is not, would add $`(g/c_p)\,z_{ref}`$ (0.4 K at a 40 m canopy) of spurious stable stratification to every
solve and lower the reported H by $`g_{ah}\,g\,z_{ref}`$.

### Re-solved at every integrator stage

Most of what the fast loop computes per cohort (photosynthesis, stomatal conductance, the leaf boundary
layers, radiation, hydraulics) is held fixed for one `dt_fast` (see [numerical_scheme](numerical_scheme.md)
§2). The conductances of eq. (6) are the exception. `mo_surface_layer` is solved again at each stage of
the time integrator: under RK45 (Cash–Karp) at the stage's own canopy-air temperature and humidity, under
ESDIRK2 at the state the stage starts from, held through its Newton solve. Only $`g_{ah}`$, $`g_{aw}`$ and
$`g_{ac}`$ are updated (`canopy_air_conductances`, from the step's reference level and canopy geometry);
the per-cohort boundary layers stay fixed.

The reason is a feedback that runs through this kernel. A warmer canopy air makes the surface layer more
unstable, which raises $`u_*`$ and $`\mathrm{temp1}`$, which vents the canopy air harder and cools it.
The gain is strong: a canopy air a fraction of a kelvin warmer multiplies $`g_{ah}`$ severalfold.
Evaluated one whole `dt_fast` behind the state it responds to, this negative feedback becomes a numerical
oscillation of the canopy-air temperature between consecutive steps, even though every conservation
budget closes. Re-solving per stage removes it. The measurements, and the check that the other fixed
coefficients do not matter, are in [numerical_scheme](numerical_scheme.md) §2. The re-solve is cheap,
because the aerodynamics is a small part of the per-step work next to leaf gas exchange.

There is no switch for this. ED2 does the same (`update_diagnostic_vars` calls `canopy_turbulence8` at
every Runge–Kutta stage), and the fixed alternative is unstable at production step sizes on most stand
heights.

## 3. Ground ↔ CAS conductance

The ground-to-CAS scalar conductance blends CLM5's bare-soil term with CLM4's dense-canopy form
(`ggbare`/`ggveg`/`ggnet`), weighted by the patch's open fraction $W$ as CLM weighs them
(Zeng et al. 2005; the stability factor is Sakaguchi and Zeng 2009's):

```math
g_{g,bare} = \frac{\kappa\,u_*}{a\,(z_{0g}\,u_*/\nu)^{b}}, \qquad
g_{g,veg} = \frac{c_s\,u_*}{1+\gamma_g\,\mathrm{stab}}, \qquad
g_{g,net} = W\,g_{g,bare} + (1-W)\,g_{g,veg}, \qquad
W = e^{-(\mathrm{LAI}+\mathrm{WAI})} \qquad(7)
```

The bare term is CLM5's soil coefficient (`soil_conductance`; Zeng & Dickinson 1998, $a = 0.13$,
$b = 0.45$), with $`z_{0g}=`$ `soil_rough` (0.01 m) and $\nu$ the air's kinematic viscosity at the
canopy air's temperature. It spans only the thin layer between the soil and its roughness height. The
turbulent layer above it is the canopy air's own exchange with the atmosphere, eq. (6), so on bare
ground the two in series are CLM5's bare-soil resistance from the soil to $`z_{ref}`$. A patch with no
cohorts, or only shorter ones, takes the canopy height whose roughness is the soil's,
`soil_rough`/`z0m_ratio` (7.7 cm), so its $`z_{0m}`$ is $`z_{0g}`$. ED2 uses $`u_*\,\mathrm{temp1}`$
here, the canopy air's own factor, which counts the surface layer twice
([ED2 comparison](../ed2_comparison.md) §4.9).

In the dense term, `stab` is a buoyancy factor built from the CAS↔ground temperature difference and
bounded to $[0,10]$, $`c_s=`$ `cs_dense`, $`\gamma_g=`$ `gamma_g`, and LAI and WAI are the patch's leaf and
wood area indices (`open_canopy_fraction`; the kernel takes $W$ as `opencan_frac`). A patch with no
cohorts ($W=1$) couples its ground as bare soil; under a closed canopy the dense form sets the
conductance ($W=0.0025$ at LAI + WAI = 6). CLM blends the transfer coefficients that multiply one
velocity scale; with both terms on $`u_*`$, that is this linear blend of the conductances. Snow does not
bury the canopy in $W$ (CLM's burial would raise $W$ through the exposed leaf and wood area). The **soil-evaporation aerodynamic resistance** is $`r_{aero}=1/g_{g,net}`$, and
$`g_{g,net}`$ is the ground sensible and latent conductance in the surface energy balance.

## 4. In-canopy wind and the leaf/wood boundary layers

The in-canopy wind sets the forced convection at each cohort's leaves and wood.

**Canopy-top wind** comes from the log profile with the surface-layer stability correction
(`reduced_wind`): the momentum denominator of eq. (3), taken from $`z_{0m}`$ up to the canopy top at
the surface layer's Obukhov length, so the profile and $`u_*`$ share CLM5's four ranges:

```math
u_h = \frac{u_*}{\kappa}\,D_m\big(\zeta_h;\ h-d,\ z_{0m}\big), \qquad
\zeta_h = \zeta\,\frac{h-d}{z_{ldis}} \qquad(8)
```

floored at `ugbmin`. When $`u_*`$ is not on its floor the profile passes through $`u_m`$ at
$`z_{ldis}`$, so $`u_h<u_m`$; on a calm night with $`u_*=u_{*,min}`$ it lies above that profile.

From the top cohort down, the wind attenuates through each cohort's crowns by ED2's per-cohort
extinction, here with **closed crowns**: every cohort's crowns cover its patch (crown area 1, as in the
two-stream radiation, [canopy_radiation_transfer](canopy_radiation_transfer.md) §3), so the wind meets
all of the cohort's leaf area. Cohorts of one height are one layer, as they are in the radiation
([canopy_radiation_transfer](canopy_radiation_transfer.md) §4): LAI below is the layer's summed leaf
area, every member meets the layer's wind, and a cohort split into two of half the density meets the
same wind as the whole one. The wind at a layer's mid-crown (also floored at `ugbmin`) and the value
handed to the layer below are

```math
u_{ico} = u_h\,e^{-0.25\,\mathrm{LAI}}, \qquad
u_h \leftarrow u_h\,e^{-0.50\,\mathrm{LAI}} \qquad(9)
```

so the cascade walks from the top down. The kernel's cohort arrays run bottom to top (the radiation
convention) while the cohort block is gathered tallest first, so `canopy_aerodynamics` takes the
gather's ascending order (`static%order`), puts the cohorts bottom to top, and puts its outputs back;
the permutation is the identity for $n\le1$.

**Boundary-layer conductance** of each cohort's leaves (flat plate) and wood (cylinder) is ED2's Nusselt
kernel `boundary_layer_conductance`: the sum of forced convection (Reynolds number) and free convection
(Grashof number), where the free term matters at low wind. The Reynolds number is
$`\mathrm{Re}=u\,\ell/\nu`$ and the Grashof number
$`\mathrm{Gr}=g\,|T_{elem}-T_{CAS}|\,\ell^3/(T_{CAS}\,\nu^2)`$, over a characteristic length $\ell$
(leaf width or branch diameter). $`T_{elem}`$ is the leaf temperature for the leaf layer and the wood
temperature for the wood layer. For the cylinder each forced term carries an offset,
$`o+a\,\mathrm{Re}^{\,n}`$ (`ocyli_*`). Each Nusselt number is the larger of its laminar and turbulent
forms:

```math
\mathrm{Nu}_{forced}=\max\!\big(a_{lam}\mathrm{Re}^{\,n_{lam}},\,a_{turb}\mathrm{Re}^{\,n_{turb}}\big),
\qquad
\mathrm{Nu}_{free}=\max\!\big(b_{lam}\mathrm{Gr}^{\,m_{lam}},\,b_{turb}\mathrm{Gr}^{\,m_{turb}}\big)
\qquad(10)
```

```math
g_{bh} = \max\!\Big(g_{bh,min},\ \big(\mathrm{Nu}_{forced}+\mathrm{Nu}_{free}\big)\,
\frac{\alpha_{th}}{\ell}\Big) \qquad(11)
```

The kinematic viscosity $\nu$ and thermal diffusivity $\alpha_{th}$ are linear-in-temperature about a
0 °C reference (`kin_visc0`, `dkin_visc`, `th_diff0`, `dth_diff`, `t_ref_air`): 1.52×10⁻⁵ and
2.15×10⁻⁵ m² s⁻¹ at 20 °C. ED2 forms the Reynolds number with $\alpha_{th}$; ED-MORE uses $\nu$, its
definition ([ED2 comparison](../ed2_comparison.md) §4.9). The floor is
$`g_{bh,min}=`$ `gbhmos_min`. The vapour conductance is $`g_{bw}=1.075\,g_{bh}`$ (`gbh_2_gbw`).
$g_{bh}$ is a velocity [m s⁻¹]; the energy kernel uses it directly, and the leaf gas exchange forms the
molar $`g_{b,mol}=g_{bw}\,\rho_{mol}`$ with $`\rho_{mol}=p/(R\,T_{leaf})`$, the molar density of air at
the leaf temperature [mol m⁻³].

## Parameters

The `[aerodynamics]` block, held in `aero_params_t`; every key is optional.

| Symbol | Config key | Meaning |
|---|---|---|
| $`\kappa`$ | none (`von_karman` in `ed_constants`) | von Karman constant (0.4) |
| $`c_{z0}`$ | `z0m_ratio` | roughness length over canopy height (0.13) |
| $`c_d`$ | `d_ratio` | displacement height over canopy height (0.63) |
| $`z_{0,snow}`$ | `snow_rough` | snow-surface roughness, a buried canopy's $`z_{0m}`$ and the $`z_{0m}`$ floor (0.001 m) |
| $`z_{0g}`$ | `soil_rough` | bare soil's roughness length (0.01 m, CLM's `zlnd`): the soil coefficient of eq. (7), a bare patch's $`z_{0m}`$ and the length scale of the snow-cover fraction ([snow_biophysics](snow_biophysics.md)) |
| none | `n_iter_mo` | the most Monin–Obukhov sweeps (20); they stop when $`u_*\cdot\mathrm{temp1}`$ settles to $`10^{-4}`$ |
| $`\zeta_m,\zeta_t`$ | `zeta_m`, `zeta_t` | momentum and scalar range transitions of CLM5 (1.574, 0.465) |
| none | `zeta_max_stable` | cap on the stable $\zeta$, of the first guess and of each sweep (0.5) |
| $w_c$ | `wc` | convective velocity of the first guess, eq. (2) (0.5 m s⁻¹); the sweeps take $w_c$ from eq. (4) |
| $`u_{*,min}`$ | `ustmin` | friction-velocity floor (0.10 m s⁻¹) |
| $`u_{b,min}`$ | `ubmin` | reference-wind floor in eq. (2) (0.65 m s⁻¹) |
| none | `ugbmin` | floor on the canopy-top and in-canopy wind (0.25 m s⁻¹) |
| $`c_s,\gamma_g`$ | `cs_dense`, `gamma_g` | CLM4 dense-canopy ground conductance (0.004, 0.5); its stability parameter is capped at 10 |
| none | `gbh_2_gbw` | boundary-layer heat-to-vapour conductance ratio (1.075) |
| $`g_{bh,min}`$ | `gbhmos_min` | boundary-layer conductance floor (10⁻⁹ m s⁻¹) |
| none | `kin_visc0`, `dkin_visc`, `th_diff0`, `dth_diff`, `t_ref_air` | the air's kinematic viscosity and thermal diffusivity, linear in temperature about `t_ref_air` (273.15 K) |
| none | `min_canopy_depth` | CAS depth floor (5 m) |
| none | `canopy_freeboard` | CAS depth above the tallest cohort (5 m) |
| none | `aflat_*`, `bflat_*` | flat-plate (leaf) Nusselt coefficients (`aero_params_t` only, no TOML key) |
| none | `ocyli_*`, `acyli_*`, `bcyli_*` | cylinder (wood) Nusselt coefficients (`aero_params_t` only, no TOML key) |

The leaf width and branch diameter that set $\ell$ are PFT traits (`[pft]`).

## Code map

| Concept | Routine (`ed_canopy_aerodynamics` unless noted) |
|---|---|
| the per-patch entry point | `canopy_aerodynamics` (cohorts reordered, then `aerodynamics_bottom_to_top`) |
| roughness and displacement, eq. (1) | `canopy_roughness`; burial fraction `canopy_burial_fraction` |
| Monin–Obukhov solve ($`u_*`$, temp1, $\zeta$), eqs. (2)–(4) | `mo_surface_layer` |
| profile denominators $D_m, D_h$ | `d_mom`, `d_heat` |
| stability functions $`\Psi_m,\Psi_h`$, eq. (5) | `stabfunc1`, `stabfunc2` |
| atmosphere↔CAS conductances, eq. (6) | `cas_atm_conductances` |
| per-stage re-solve of eq. (6) | `canopy_air_conductances`, called by `refresh_cas_conductances` in `ed_fast_integ_derivs` |
| ground conductance and open fraction $W$, eq. (7) | `soil_conductance`, `open_canopy_fraction`; blended inside `aerodynamics_bottom_to_top` |
| canopy-top wind, eq. (8) | `reduced_wind` |
| in-canopy wind, eq. (9) | `aerodynamics_bottom_to_top` (inline, closed-crown extinction) |
| leaf and wood Nusselt boundary layer, eqs. (10)–(11) | `boundary_layer_conductance` |
| types | `ed_biophysics_params`: `aero_params_t`; `aero_env_t`, `aero_args_t`, `aero_out_t` |
| filling the environment | `set_aero_env_atm` (the reference level), `set_aero_env_canopy` (the canopy air and ground) |
| cohorts of one height are one layer | `ed_numerics`: `equal_height_layers` |

## References

- ED2 `../ED2/ED/src/dynamics/canopy_struct_dynamics.f90` (`canopy_turbulence8`): boundary layers and
  in-canopy wind.
- CLM5 Technical Note (Oleson et al. 2018): Monin–Obukhov surface layer, stability functions.
- Businger et al. (1971); Paulson (1970): surface-layer flux–profile relations.
- Zeng & Dickinson (1998), *J. Climate* 11:2628–2644: the bare-soil coefficient of eq. (7).
- Zeng, X., et al. (2005), Treatment of undercanopy turbulence in land models, *J. Climate* 18,
  5086–5094: the open-fraction blend of eq. (7); Sakaguchi, K., and X. Zeng (2009), *J. Geophys.
  Res.* 114, D01107: its under-canopy stability factor.
