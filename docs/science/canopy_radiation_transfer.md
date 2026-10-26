# Canopy radiative transfer

How MEDS shares the incident shortwave and longwave radiation among the cohorts of a patch and the
ground. The page gives the leaf-angle distribution, the optical properties of leaves and wood, the
multi-layer **two-stream** solution, the ground boundary, and the fluxes handed to photosynthesis, the
energy balance and the output.

The scheme is ED2's two-stream multi-layer canopy optics (`icanrad=2` in
`../ED2/ED/src/dynamics/twostream_rad.f90`; Liou 2002 ch. 6; Longo et al. 2019), solved by the adding
method of §3. One solver runs every band; the default bands are visible (VIS), near-infrared (NIR) and
thermal longwave (LW). Its products are, per cohort, the absorbed leaf VIS, which the leaf-photosynthesis
kernel takes as PAR (photosynthetically active radiation), and the absorbed leaf and wood shortwave and
longwave, which drive the leaf and wood energy balance. The code takes plain arrays and value types and
never a `polygon_t`.

Symbols: $\mu = \cos\theta_{sun}$ (cosine of the solar zenith angle; the beam extinction $k$ guards
$\mu\to0$); $\rho, \tau$ = leaf or wood reflectance and transmittance per band; $\omega = \rho+\tau$ =
single-scatter albedo; $g$ = leaf-angle scattering asymmetry; $\beta$ = diffuse **backscatter**
fraction; $\beta_0$ = direct-beam **upscatter** fraction; $k$ = direct-beam extinction per unit area
index; ELAI, EWAI = the clumping-corrected leaf and wood area indices. All fluxes are absolute
`W m⁻²`; nothing is normalized to unit incidence.

The pipeline has four stages: (1) the leaf-angle distribution and its integrals (`meds_optics_lib`);
(2) the per-PFT (plant functional type) scattering table, independent of $\mu$, built once at load
(`build_rad_optics`, `derive_rad_optics` in `meds_pft_derive`); (3) the per-cohort blend of leaf and wood
optics with the solar geometry (`blend_cohort_optics`); and (4) the single-band adding solver. Stages 3
and 4 run inside `canopy_radiation` (`meds_canopy_radiation`). The optical-property kernels are in the
shared library `meds_optics_lib` (`src/share/functions/`); the assembly, the solver and the per-patch
entry point are together in `meds_canopy_radiation`.

## 1. Leaf-angle distribution and the G-function

How leaves are tilted decides how much light a canopy intercepts and how it scatters it. The leaf
inclination distribution (LIDF) is a two-parameter **Beta** density on the normalized inclination
$`t=\theta/(\pi/2)\in[0,1]`$ (Goel & Strebel 1984; the SCOPE / 4SAIL family, Verhoef 1984,
van der Tol et al. 2009). One family covers the Verhoef archetypes (planophile, erectophile,
plagiophile, extremophile, and uniform $=\mathrm{Beta}(1,1)$). MEDS uses it in place of ED2/CLM's
$`\phi_1/\phi_2/\bar\mu`$ and $`(1+\chi_L)^2`$ backscatter approximation (design record:
[`radiative_transfer_design.md`](../dev_plans/archive/radiative_transfer_design.md)).

The class weights are integrals of the unnormalized Beta kernel over the SCOPE 13-class inclination
grid, renormalized to sum to one (the normalizing constant cancels):

```math
w_c \;\propto\; \int_{t_{c-1}}^{t_c} t^{\,p-1}(1-t)^{\,q-1}\,dt , \qquad \sum_c w_c = 1 \qquad(1)
```

evaluated by a 64-panel midpoint rule (`beta_lidf` over `beta_log_kernel`). It needs no special
function and never evaluates the integrable endpoint singularities at $t=0,1$. The kernel is summed
relative to its largest value on the grid, so a narrow distribution (a spread of a few degrees, or 0),
whose unnormalized density would underflow, still sums to one. The mean and standard deviation of the
leaf angle, $`(\overline{\theta},\,\sigma_\theta)`$, map to $(p,q)$ by moment matching
(`moments_to_beta_params`).

Two quantities feed the two-stream optics. The **second moment**, independent of $\mu$, sets the split
between back and forward scattering:

```math
\mathrm{bf} \;=\; \langle\cos^2\theta_{leaf}\rangle \;=\; \sum_c w_c\,\cos^2\theta_c \qquad(2)
```

($`\mathrm{bf}=1`$ for horizontal leaves, $0$ for vertical, $1/3$ for spherical; `leaf_bf`). The **Ross
G-function**, the mean leaf-area projection toward the sun, depends on $\mu$ and uses the SCOPE
`volscatt` projection $`\chi_s(\theta_c,\mu)`$:

```math
G(\mu) \;=\; \sum_c w_c\,\chi_s(\theta_c,\mu), \qquad
k \;=\; \frac{G(\mu)}{\mu} \qquad(3)
```

`gfun_direct` returns $G(\mu)$; the direct-beam extinction coefficient is $`k=G(\mu)/\mu`$. $G$ does
not depend on the band, so `canopy_radiation` evaluates it once per PFT present and reuses it for every
cohort of that PFT and every beam band.

## 2. Single-scatter albedo, asymmetry, and the scattering coefficients

For each PFT and band the $\mu$-independent optics are built once (`derive_rad_optics` over
`scatter_pair`):

```math
\omega = \rho + \tau, \qquad
g = \mathrm{bf}\,\frac{\rho-\tau}{\omega} \qquad(4)
```

($g=0$ for a perfect absorber, $\omega\to0$). For the **thermal band** the optics table takes
$`\rho = 1-\varepsilon`$ and $\tau=0$ (ED2/SCOPE convention), so $`\omega = 1-\varepsilon`$ and
$`g=\mathrm{bf}`$, with no special case.

For each cohort and band, `blend_cohort_optics` weights the leaf and wood optics by their
clumping-corrected areas $`\mathrm{ELAI}=\Omega_{leaf}\,\mathrm{LAI}`$,
$`\mathrm{EWAI}=\Omega_{wood}\,\mathrm{WAI}`$ (weights $`w_l=\mathrm{ELAI}/\mathrm{ETAI}`$,
$`w_w=\mathrm{EWAI}/\mathrm{ETAI}`$, $`\mathrm{ETAI}=\mathrm{ELAI}+\mathrm{EWAI}`$), then adds the
solar geometry:

```math
\omega = w_l\,\omega_{leaf} + w_w\,\omega_{wood}, \qquad
g = w_l\,g_{leaf} + w_w\,g_{wood} \qquad(5)
```

```math
\beta = \tfrac12\,(1 + g), \qquad
\beta_0 = \tfrac12\!\left(1 + \frac{g}{k}\right) \qquad(6)
```

$\beta$ is the diffuse backscatter fraction and $\beta_0$ the direct-beam upscatter fraction
(computed for beam bands only). The routine also returns the **leaf absorption share**, the fraction of
a cohort's absorbed radiation that lands on leaves rather than wood, weighted by absorptivity $1-\omega$
and clumped area:

```math
\mathrm{leaf\_frac} = \frac{(1-\omega_{leaf})\,\mathrm{ELAI}}
{(1-\omega_{leaf})\,\mathrm{ELAI} + (1-\omega_{wood})\,\mathrm{EWAI}} \qquad(7)
```

## 3. The single-band two-stream (adding) solve

`solve_band` (`meds_canopy_radiation`) is the one single-band solver. It works on the vertical stack of
cohort layers, ordered from the bottom (index 1) to the top (index n), with the sky above. Shortwave
bands carry a direct beam and no emission; the thermal band carries emission and no beam. The diffuse
multiple-scattering operator is the same in every band; only the source terms and the boundary
conditions differ. The solver uses the **adding** (layer-recursion) method: two sweeps whose cost grows
linearly with the number of layers, bottom-up for the reflectance and upward source of the stack below
each layer, then top-down for the fluxes. Each layer covers the whole ground (crown area index 1). A
layer's diffuse operator, the beam it scatters and its emission are each the exact solution for a
homogeneous layer, so the adding method is exact for the stack, and a uniform canopy's light does not
depend on how it is cut into cohorts.

**Diffuse reflectance and transmittance of a layer.** For a homogeneous layer of effective area index
ETAI (ELAI + EWAI), `layer_rt` gives the classical Stokes two-stream coefficients, with backward and
forward diffuse scattering $`\sigma_b=\omega\beta`$, $`\sigma_f=\omega(1-\beta)`$, attenuation
$`\mathrm{att}=1-\sigma_f`$:

```math
m=\sqrt{\mathrm{att}^2-\sigma_b^2},\quad
r_\infty=\frac{\mathrm{att}-m}{\sigma_b},\quad
e=e^{-m\,\mathrm{ETAI}} \qquad(8)
```

```math
t_{dd}=\frac{(1-r_\infty^2)\,e}{1-r_\infty^2 e^2}, \qquad
r_{dd}=\frac{r_\infty(1-e^2)}{1-r_\infty^2 e^2} \qquad(9)
```

($\omega$ is capped at $1-10^{-6}$ to avoid the conservative-scattering singularity at $\omega\to1$;
$`\sigma_b\to0`$ collapses to pure Beer attenuation with $r_{dd}=0$.)

**Direct beam.** The collimated beam attenuates by Beer's law from the top down,

```math
I^\downarrow_0(i) = I^\downarrow_0(i{+}1)\,e^{-k_i\,\mathrm{ETAI}_i} \qquad(10)
```

and inside a layer, at depth $x$ in area index, it scatters $`\omega k S_0 e^{-kx}`$ into the diffuse
field, a fraction $\beta_0$ of it upward ($`S_0 = I^\downarrow_0(i{+}1)`$, the beam at the layer's top).
The layer's up- and down-sources are the part of this light that leaves its top and its bottom. They
come from the two-stream equations solved through the layer's depth with no diffuse entering either
face (`layer_beam_sources`). With the attenuation att and backscatter $`\sigma_b`$ of eq. (8), the
particular solution is

```math
u_p = a\,e^{-kx},\quad d_p = b\,e^{-kx},\qquad
a = \frac{\omega k S_0\,[\beta_0(\mathrm{att}-k) + \sigma_b(1-\beta_0)]}{m^2-k^2},\quad
b = \frac{\omega k S_0\,[(1-\beta_0)(\mathrm{att}+k) + \sigma_b\beta_0]}{m^2-k^2} \qquad(11)
```

The diffuse modes $`e^{-mx}`$ and $`e^{m(x-\mathrm{ETAI})}`$, in the ratio $`r_\infty`$ of up to
down, take the boundary conditions. $k = m$ is a removable singularity, avoided by moving $k$ off $m$ by
a relative $`10^{-6}`$.

A thin layer's sources tend to the single-scatter split, $`\beta_0\,\omega\,\Delta I`$ up and
$`(1-\beta_0)\,\omega\,\Delta I`$ down of the beam it intercepts. A thick layer's do not, because light
scattered deep in the layer is absorbed on its way out. Injecting the single-scatter split at the faces
instead would make the result depend on how the canopy is cut into cohorts; the exact source gives the
same albedo for any cut.

For the thermal band each layer adds a blackbody source $`(1-r_{dd}-t_{dd})\,\sigma T^4`$ to both
faces. Each tissue emits at its absorptivity (Kirchhoff's law), so a cohort's $`\sigma T^4`$ is the
mean of its tissues' at their start-of-step temperatures, weighted as their absorption is, by the
thermal band's `leaf_frac` (eq. (7)):
$`\sigma T^4 = \mathrm{leaf\_frac}\,\sigma T_{leaf}^4 + (1-\mathrm{leaf\_frac})\,\sigma T_{wood}^4`$. The
longwave field therefore sees both tissues.

**The recursion.** The bottom-up sweep builds the reflectance $`r_{bel}`$ and upward source $`s_{bel}`$
of the growing stack {ground + layers $1..i$}. It starts from the ground reflectance and the ground's
upward source (thermal emission plus the reflected direct beam reaching the ground):

```math
r_{bel}(i) = r_{dd}(i) + \frac{t_{dd}(i)^2\,r_{bel}(i{-}1)}{1-r_{bel}(i{-}1)\,r_{dd}(i)},\qquad
s_{bel}(i) = s_{up}(i) + t_{dd}(i)\,\frac{s_{bel}(i{-}1)+r_{bel}(i{-}1)\,s_{dn}(i)}{1-r_{bel}(i{-}1)\,r_{dd}(i)}
\qquad(12)
```

The top-down sweep, which starts from the incident diffuse radiation at the canopy top, then recovers
the upward and downward diffuse fluxes at every interface. The radiation **absorbed by a layer** is the
beam plus diffuse flux divergence (in minus out), so energy is conserved by construction:

```math
A_i = \big(I^\downarrow_0(i{+}1)-I^\downarrow_0(i)\big)
      + \big(d^\downarrow_{in}(i)+u^\uparrow_{up}(i)\big)
      - \big(u^\uparrow_{top}(i)+d^\downarrow_{bot}(i)\big) \qquad(13)
```

The solver also returns the downwelling flux below the canopy $`d_{ground}`$, the upwelling flux from
the ground $`u_{ground}`$, and the patch **albedo**, the upward flux leaving the top over the total
incident flux.

## 4. Ground optics and the returned fluxes

This section gives the lower boundary of the two-stream and the fluxes `canopy_radiation` returns.

`ground_optics` closes the two-stream at the bottom. For a shortwave band the ground reflectance is the
per-band soil albedo and its emission is zero. For the thermal band the reflectance is
$1-\varepsilon_{soil}$ and the emission is $`\varepsilon_{soil}\,\sigma\,T_{soil}^4`$. Its input
(`ground_optics_env_t`) is **bare soil plus snow**: the per-band soil albedo, the emissivity and the
skin temperature, which the driver blends toward the snowpack's albedo and emissivity by the Niu–Yang
snow-cover fraction, so a partial pack changes the ground continuously. The blended temperature is the
one at which the blended emissivity emits what the soil and the pack emit at their own temperatures.
What reaches the ground in each band is then shared between soil and snow at their own optics
([snow_biophysics](snow_biophysics.md), "The ground's radiation"). Litter, standing water and soil
moisture do not change the ground optics ([ROADMAP](../ROADMAP.md)).

`canopy_radiation` (the per-patch entry point) loops over the bands, blends the optics per cohort and
passes each band to `solve_band`.

**Co-dominant cohorts share a layer.** Cohorts of exactly equal height form one radiative layer. Its
ETAI is the sum of the members', its $\omega$, $g$ (hence $\beta$, $\beta_0$) and $k$ are ETAI-weighted
means of theirs, and its $`\sigma T^4`$ is weighted, as within a cohort, by the absorptivity-weighted
area $`(1-\omega_i)\,\mathrm{ETAI}_i`$. The solve's absorbed flux $`A_\ell`$ is net of the layer's
emission. Every member of the layer meets the same radiation, so what the layer takes in, $`A_\ell`$
plus the emission $`E_\ell=\sum_{j\in\ell}E_j`$ of its members, goes back to the members by
absorptivity-weighted area (the weight `leaf_frac` is built from), and each is charged its own
emission:

```math
A_{b,i} = \big(A_{b,\ell}+E_{b,\ell}\big)\,\frac{(1-\omega_i)\,\mathrm{ETAI}_i}{\sum_{j\in\ell}(1-\omega_j)\,\mathrm{ETAI}_j}
          - E_{b,i}
\qquad(14)
```

A cohort's emission is its tissues', each from both faces of its area at its absorptivity,
$`E^{leaf}_{b,i} = 2\,(1-\omega_{leaf})\,\mathrm{ELAI}_i\,\sigma T_{leaf}^4`$ and
$`E^{wood}_{b,i} = 2\,(1-\omega_{wood})\,\mathrm{EWAI}_i\,\sigma T_{wood}^4`$, zero in a shortwave
band. Every member sees the light reaching the layer's top. A one-member layer is exactly the per-cohort
solve. Stacking tied cohorts instead would put whichever sorts first on top, so the light split between
PFTs would follow the PFT index rather than the ecology; the canopy total hardly changes, so only the
split between the tied cohorts would show it. Ties arise at recruitment, when every PFT recruits at
the one smallest tracked height, `min_cohort_height`.

`canopy_radiation` then splits each cohort's absorbed radiation into leaves and wood the same way:
they share what the cohort takes in by the `leaf_frac` weight the optics blend produced, and each is
charged its own emission:

```math
A^{leaf}_{b,i} = \big(A_{b,i}+E_{b,i}\big)\,\mathrm{leaf\_frac}_i - E^{leaf}_{b,i}, \qquad
A^{wood}_{b,i} = \big(A_{b,i}+E_{b,i}\big)\,(1-\mathrm{leaf\_frac}_i) - E^{wood}_{b,i} \qquad(15)
```

In a shortwave band this is the split by `leaf_frac`. In the thermal band it is exact for the
two-stream's well-mixed layer: leaves and wood meet the same radiation through its depth and absorb it
in the proportion `leaf_frac`, and the solve's net is what the layer absorbs less what it emits. Wood
warmer than its leaves therefore warms them: with the sky, the ground and the leaves at one temperature
and the wood warmer, the leaves gain longwave and the wood loses it. Splitting the cohort's net by
`leaf_frac` alone would charge the cool leaves a share of the warm wood's emission.

A patch with no resolvable cohorts short-circuits: all radiation reaches the ground, and the albedo is
the ground's alone.

**PAR passed to photosynthesis.** The absorbed leaf VIS $`A^{leaf}_{\mathrm{VIS},i}`$ is an absorbed
flux, but the leaf gas-exchange kernel takes incident PAR, because it applies the leaf absorptance
itself when it forms the electron-transport light $`I_2`$ ([leaf_gas_exchange](leaf_gas_exchange.md)
§2). The fast-loop driver (`apply_rt_forcing` in `meds_fast_forcing`) therefore passes
$`A^{leaf}_{\mathrm{VIS}} / \alpha_{leaf}`$ as the incident-equivalent PAR. Passing the absorbed flux
would apply the absorptance twice. The absorbed leaf and wood shortwave and longwave go to the energy
balance as absorbed, with no division.

## Surface radiative fluxes on the diagnostic path

To compare a run with a radiometer or a satellite product, the output carries the top-of-canopy
shortwave and longwave fluxes. The two-stream forms the upward flux every step as
`albedo(b) x incid(b)` for each band $b$. Five variables carry it:

| variable | meaning |
|---|---|
| `mean_sw_in_vis_py`, `mean_sw_in_nir_py` | incident shortwave at canopy top, per band |
| `mean_sw_up_vis_py`, `mean_sw_up_nir_py` | upwelling shortwave leaving the canopy top |
| `mean_lw_up_py` | upwelling longwave, surface emission included, plus the change in emission within the step that the tissues' and the ground's balances booked, so the reported budget closes on the booked net radiation ([vegetation_energy_dynamics](vegetation_energy_dynamics.md), "What the output reports") |

The output holds fluxes and not a time-averaged albedo. A period-mean albedo is the mean of a ratio,
which is not the ratio of the means, and at night the shortwave ratio is 0/0. The albedo a reader
wants, and what a satellite product reports, is

```math
\alpha_b = \frac{\sum_t \mathrm{SW}^\uparrow_b}{\sum_t \mathrm{SW}^\downarrow_b} \qquad(16)
```

over the period of interest. The contrast between the low VIS and high NIR albedo is the vegetation
red edge: chlorophyll absorbs the visible and leaves scatter the near-infrared, the signal that
vegetation indices are built on. A snow-covered ground raises both.

## Parameters

| Symbol | Source | Meaning |
|---|---|---|
| $`\overline{\theta}, \sigma_\theta`$ | `[pft].leaf_angle_mean`, `leaf_angle_std` | mean and standard deviation of the leaf inclination [deg] |
| $p, q$ | derived from $`\overline{\theta}, \sigma_\theta`$ (`moments_to_beta_params`) | Beta-LIDF shape parameters |
| $`\rho_{leaf},\tau_{leaf}`$ | `[pft].leaf_reflect_vis`, `leaf_reflect_nir`, `leaf_transmit_vis`, `leaf_transmit_nir` | leaf reflectance and transmittance |
| $`\rho_{wood},\tau_{wood}`$ | `[pft].wood_reflect_vis`, `wood_reflect_nir`, `wood_transmit_vis`, `wood_transmit_nir` | wood reflectance and transmittance |
| $`\varepsilon_{leaf},\varepsilon_{wood}`$ | `[pft].leaf_emissivity`, `wood_emissivity` | longwave: $`\rho=1-\varepsilon,\ \tau=0`$ |
| $`\Omega_{leaf},\Omega_{wood}`$ | `[pft].leaf_clumping`, `wood_clumping` | clumping factors in $(0,1]$ |
| $`\varepsilon_{soil}`$ | `[soil].ground_emissivity` (0.95) | ground thermal emissivity |
| — | `[soil].ground_albedo_vis`, `[soil].ground_albedo_nir` (0.15, 0.30) | bare-soil shortwave albedo per band |
| $`\alpha_{leaf}`$ | `[leaf_physiology].leaf_absorptance` | leaf VIS absorptance (the PAR conversion above) |

## Code map

| Concept | Routine |
|---|---|
| Beta LIDF and its parameters | `meds_optics_lib`: `beta_lidf`, `beta_log_kernel`; `meds_pft_derive`: `moments_to_beta_params` |
| moments $`\mathrm{bf}`$, $G(\mu)$ | `meds_optics_lib`: `leaf_bf`, `gfun_direct` |
| $\omega, g$ single-scatter pair | `meds_optics_lib`: `scatter_pair` |
| per-PFT optics table ($\omega, g$), eq. (4) | `meds_pft_derive`: `build_rad_optics`, `derive_rad_optics`, once at load onto the PFT table (`cfg%pft%optics`) |
| per-cohort blend ($\beta, \beta_0, k$, `leaf_frac`), eqs. (5)–(7) | `meds_canopy_radiation`: `blend_cohort_optics` |
| ground reflectance and emission | `meds_canopy_radiation`: `ground_optics` (`ground_optics_env_t`) |
| single-band adding solve, eqs. (8)–(13) | `meds_canopy_radiation`: `solve_band`, `layer_rt`, `layer_beam_sources` |
| the per-patch entry point, tied-cohort layers, each tissue's emission and the leaf/wood split, eqs. (14)–(15) | `meds_canopy_radiation`: `canopy_radiation` |
| cohorts of one height are one layer (as in the wind and the overtopping LAI) | `meds_numerics`: `equal_height_layers` |
| types | `meds_pft_params`: `rad_pft_optics_t` and the bands `RAD_VIS`, `RAD_NIR`, `RAD_LW`; `meds_canopy_radiation`: `rad_env_t`, `rad_out_t`, `ground_optics_env_t` |
| fast-loop join (PAR conversion, cohorts bottom to top, top-of-canopy fluxes) | `meds_fast_forcing`: `apply_rt_forcing` |
| the ground's optics under snow, and each tile's radiation | `meds_surface`: `snow_blend_ground_optics`, `ground_tile_radiation` |

## References

- Longo et al. (2019), *GMD* 12:4309: ED-2.2 technical description (the `icanrad=2` two-stream).
- Liou (2002), *An Introduction to Atmospheric Radiation*, ch. 6: two-stream theory.
- Verhoef (1984), *Remote Sens. Environ.* 16:125–141: SAIL, the diffuse two-stream layer operator.
- van der Tol et al. (2009), *Biogeosciences* 6:3109: SCOPE (leaf angle and `volscatt`).
- Goel & Strebel (1984), *Agron. J.* 76:800: Beta leaf-inclination distribution.
- ED2 `../ED2/ED/src/dynamics/twostream_rad.f90`; design records
  [`radiative_transfer_design.md`](../dev_plans/archive/radiative_transfer_design.md) and
  [`MEDS_ENERGY_BALANCE_DESIGN.md`](../dev_plans/archive/MEDS_ENERGY_BALANCE_DESIGN.md).
