# Plant hydraulics

ED-MORE resolves water transport per individual as two water stores per cohort, the leaves and the wood (sapwood and fine roots). The xylem joins the two, the roots connect the wood to the soil layers, transpiration drives the top and soil water the bottom. The prognostic state is each tissue's water mass, carried in the cohort state; its potential is read off the tissue's pressure–volume curve. Everything is in **MPa** (never metres of head), matching the leaf gas-exchange input `psi_leaf` and the `[pft].wstress_psi_*` traits of each plant functional type (PFT). The physical reference is Xu et al. (2016, X16), revised for coupled stores, a nonlinear pressure–volume curve and a Kirchhoff-integrated conductance.

This page gives the plant's curves and flow laws. Each fast step solves them together with the soil's, as one network: [soil_plant_water](soil_plant_water.md).

Total head is $`\Psi_i = \psi_i + \rho g\,z_i`$ with $`\rho g = 9.80665\times10^{-3}\ \mathrm{MPa\,m^{-1}}`$ (`grav_head`, water's density times standard gravity; the soil uses the same constant). Signs: $\psi \le 0$ (tension), $z$ measured upward.

## 1. The water balance

Per-store water mass balance (kg H₂O per plant):

```math
C_i(\psi_i)\,\frac{d\psi_i}{dt} \;=\; \sum_{j\sim i} K_{ij}(\psi)\,\big(\psi_j - \psi_i + g_{ij}\big) \;-\; S_i \qquad(1)
```

with $C_i$ the **capacitance** [kg MPa⁻¹] (§2), $K_{ij}$ the **edge conductance** [kg s⁻¹ MPa⁻¹] (§3),
$`g_{ij}=\rho g\,(z_j-z_i)`$ the gravity offset, and $S_i$ a sink. Each cohort has two stores, its
leaves L and its lumped wood W:

```math
C_L\,\dot\psi_L = K_{LW}\big(\psi_W-\psi_L+g_{WL}\big) - E, \qquad
C_W\,\dot\psi_W = K_{LW}\big(\psi_L-\psi_W-g_{WL}\big) + Q_{\text{root}} \qquad(2)
```

where $`g_{WL}=-\rho g\,H`$ (leaf sits a height $H$ above the wood datum), $E$ is transpiration
[kg s⁻¹], and $Q_{\text{root}}$ is root water uptake (§4), the sum of the uptake from every layer the
cohort roots in.

## 2. Pressure–volume curves and capacitance

Tissue water storage follows the nonlinear **Bartlett / Tyree–Hammel** pressure–volume relation, per tissue (leaf, wood). Its traits are the osmotic potential at full turgor $\pi_0<0$, the bulk elastic modulus $\varepsilon>0$, the apoplastic fraction $a_f\in[0,1)$ and the saturated water content $w_{sat}$ [kg H₂O / kgC].

**Turgor loss point** and the symplastic RWC there:

```math
\psi_{tlp}=\frac{\pi_0\,\varepsilon}{\pi_0+\varepsilon}, \qquad R_{tlp}=\frac{\pi_0+\varepsilon}{\varepsilon} \qquad(3)
```

**Potential ↔ symplastic RWC** ($R$):

```math
\psi(R)=\begin{cases}\varepsilon\,(R-R_{tlp})+\pi_0/R & R\ge R_{tlp}\ \text{(turgid)}\\[2pt]
\pi_0/R & R< R_{tlp}\ \text{(flaccid)}\end{cases}
\qquad
R(\psi)=\begin{cases}\dfrac{b+\sqrt{b^2-4\varepsilon\pi_0}}{2\varepsilon},\ b=\psi+\varepsilon+\pi_0 & \psi\ge\psi_{tlp}\\[6pt]
\pi_0/\psi & \psi<\psi_{tlp}\end{cases} \qquad(4)
```

**Tissue water and capacitance** ($`W_{sat}=w_{sat}\cdot\text{biomass}`$; the apoplast is a constant
reservoir):

```math
W(\psi)=(1-a_f)\,W_{sat}\,R(\psi)+a_f\,W_{sat}, \qquad
C(\psi)=\frac{dW}{d\psi}=(1-a_f)\,W_{sat}\,\frac{dR}{d\psi} \qquad(5)
```

```math
\frac{dR}{d\psi}=\begin{cases}\big(\varepsilon-\pi_0/R^2\big)^{-1} & \text{turgid}\\[2pt]
-R^2/\pi_0 & \text{flaccid}\end{cases} \qquad(6)
```

At full turgor ($R=1$) the turgid branch gives $`dR/d\psi = 1/(\varepsilon+|\pi_0|)`$, so there the
capacitance equals the pure-elastic (X16-style) linear value
$`C_{\text{lin}}=(1-a_f)\,W_{sat}/(\varepsilon+|\pi_0|)`$.

**A leafless cohort's leaf node.** A cohort with no leaves has an empty leaf store and no curve to read: eq. (4) would put its leaves at the floor of the relative water content, about −15,000 MPa. Its leaf node sits instead at the wood potential less the gravity head to the crown, $`\psi_L = \psi_W - \rho g\,H`$ (`leaf_node_psi`). This is the value the leaf node relaxes to with no transpiration and the one Xu et al. (2016) give a leafless cohort. A bare cohort's predawn potential, which the stomata and the phenology water cue read, therefore follows its stem and roots.

## 3. Xylem vulnerability and the Kirchhoff conductance law

The retained conductance fraction (1 − PLC) is a Weibull-like curve in $r=\psi/\psi_{50}\ge 0$
($\psi_{50}<0$ is the potential at 50 % loss; $a=$ `wood_kexp` the shape):

```math
k(\psi)=\frac{1}{1+(\psi/\psi_{50})^{a}}=\frac{1}{1+r^{a}} \qquad(7)
```

**Vulnerability enters through the Kirchhoff (matric flux) potential — never pointwise.** Define

```math
\Phi(\psi)=\int_0^{\psi} k(s)\,ds=\psi_{50}\!\int_0^{r}\frac{du}{1+u^{a}}
=\begin{cases}\psi_{50}\ln(1+r) & a=1\\ \psi_{50}\arctan(r) & a=2\\
\psi_{50}\,G(r)\ \text{(lookup table)} & \text{else}\end{cases} \qquad(8)
```

The **edge conductance** over a finite drop is the $\psi$-averaged retained fraction. Without gravity it gives the exact steady flux through the path. With gravity the flow is this conductance times the gravity-inclusive drive, which overstates the steady flux where the xylem is embolised and the head is a large share of the drop. The exact flux $Q$ (per $`k_{cond}`$) solves $`\int_{\psi_L}^{\psi_W} d\psi/(Q + \rho gH + Q\,(\psi/\psi_{50})^{a}) = 1`$, an implicit equation, not a function of $\Phi$ at the two ends; near equilibrium its conductance is the harmonic mean of $k$ along the path, where the law takes the arithmetic mean. At $`\psi_W=-1`$, $`\psi_L=-3`$, $`\rho gH = 0.5`$ and $`\psi_{50}=-1.5`$ MPa the law overstates the flux by 3.4 % at $a = 2$ and 15.5 % at $a = 4$; on realistic ranges (at Barro Colorado Island, BCI, with $a = 2$: drops up to 0.5 MPa and trees up to 36 m) it is under 0.6 %. The exact flux would need a root solve inside every evaluation of the water network's Newton, so the law stays:

```math
K_{ij}=k_{\text{cond}}\,\frac{\Phi(\psi_{up})-\Phi(\psi_{down})}{\psi_{up}-\psi_{down}}
\;\xrightarrow[\Delta\psi\to0]{}\; k_{\text{cond}}\,k(\psi) \qquad(9)
```

$k_{\text{cond}}$ is the maximum (fully hydrated) per-plant conductance, the sapwood's over the path from the stem base to the crown:

```math
k_{\text{cond}}=\frac{K_s\,A_{sap}}{H\cdot\text{vessel\_curl}} \qquad(10)
```

```math
K_s=\exp\!\left(2.348-2.455\,\rho+\tfrac{1}{2}\,0.6186\right) \qquad(11)
```

$A_{sap}$ is the sapwood area of the stem allometry, $H$ the plant's height, and `vessel_curl` (1.5, as in ED2) lengthens the path beyond the height. $K_s$ [kg m⁻¹ s⁻¹ MPa⁻¹], the sapwood specific conductivity, is a PFT trait, `[hydraulics].wood_kmax`. The shipped PFTs take it from their wood density $\rho$ [g cm⁻³] through the fit of Xu et al. (2016, New Phytol 212:80) to 72 species of neotropical seasonally dry forests, as ED2 codes it (`plant_hydro_scheme` 2), eq. (11). The fit is computed outside the model, by `scripts/prepare_input/pft_traits`. It is lognormal, so the residual variance 0.6186 turns its median into the mean (×1.36); $\rho$ outside 0.35–0.95 takes the nearest edge. $K_s$ is 6.0 at $\rho$ = 0.35, 3.3 at 0.60 and 1.4 at 0.95. Wood density explains a fifth of the variance ($R^2$ 0.21): one species' $K_s$ can lie a factor of 2 either side. Branches of Panama canopy trees span about 1.6–12, and 2–4 at $\rho$ = 0.6 (Meinzer et al. 2008), and FATES's default is 3.0. For a 35 m BCI canopy tree (72 cm dbh, 0.30 m² of sapwood, 553 m² of leaf) the 0.60 value gives $k_{\text{cond}}$ = 0.019 kg s⁻¹ MPa⁻¹: 1.9 mmol m⁻² s⁻¹ MPa⁻¹ per unit leaf area, and 63 g m⁻² s⁻¹ MPa⁻¹ per unit sapwood, inside the 25–131 measured from soil to branch on Panama canopy trees (Meinzer et al. 2003).

For general $`a\notin\{1,2\}`$ the Kirchhoff integral, eq. (8), is precomputed once into a fixed uniform-grid **lookup table** $G(r)$ over $`0\le r\le 8`$ and read by linear interpolation on the hot path (the closed forms are kept for $`a\in\{1,2\}`$). The table stores the $r$-normalized $G$, so $\psi_{50}$ is a runtime scale. It is built panel by panel, so it rises monotonically, and past $r = 8$ the tail $`\int u^{-a}\,du`$ is added analytically, so the potential keeps approaching its asymptote even for a leaf dried far beyond the vulnerability curve. Each of the table's 512 intervals is one 7-point Gauss–Legendre panel: a single rule over a long range misses the integrand's fall near $u = 1$ (at $r$ = 7800 and $a$ = 4 it gives about 10⁻⁶ instead of 1.11), which would hand the xylem of a desiccated leaf a negative conductance.

**In the flow.** Above zero potential nothing is embolised, so the sapflow takes $\Phi(\psi) = \psi$ there. When the two ends of eq. (9) are within 10⁻⁶ MPa, the conductance and its slopes are taken at their midpoint, where the difference quotient loses its digits. Farther apart, for $a=2$ with both ends at or below zero, the difference is one arctangent, $`\Phi(\psi_{up})-\Phi(\psi_{down})=\psi_{50}\arctan\big((r_{up}-r_{down})/(1+r_{up}\,r_{down})\big)`$, which keeps the digits two close potentials share. The slopes in the solver's Jacobian are those of exactly the $\Phi$ computed (the retained fraction for the closed forms, the slope of the table's linear piece inside the table, the tail's beyond it), so the Jacobian is the derivative of the flow that moves the water. A fully embolised cohort keeps a 10⁻¹² share of $k_{\text{cond}}$, and a plant with no sapwood yet a $k_{\text{cond}}$ of 10⁻¹⁵ kg s⁻¹ MPa⁻¹, so the leaves of every cohort stay joined to its wood and its equations stay solvable.

## 4. Root water uptake

$Q_{\text{root}}$ enters the wood from the soil, along one path per soil layer the cohort roots in (as in ED2). Each path sees its layer's soil potential and **unsaturated** conductivity $`K_{soil}(k)=K(\theta_k)`$, from the layer's retention curve ([soil_biophysics](soil_biophysics.md)).

**Each cohort has its own roots.** It roots to a depth set by its height (ED2's IALLOM 1 allometry, Christoffersen 2013), capped at the soil column's bottom $z_{\text{bot}}$:

```math
D=\min\!\left(b_{1Rd}\,h^{\,b_{2Rd}},\ z_{\text{bot}}\right) \qquad(12)
```

with $`b_{1Rd}`$ = `root_depth_b1` (1.114 m) and $`b_{2Rd}`$ = `root_depth_b2` (0.4223), which give 5 m for a 35 m tree and 1.8 m for a 3 m sapling. Within $D$ its fine roots fall off exponentially, as in ED2's `distrib_root`: the share above depth $d$ is

```math
Y(d)=\frac{1-\beta^{\,\min(d,D)/D}}{1-\beta} \qquad(13)
```

so a layer between $d_{k-1}$ and $d_k$ holds $`Y(d_k)-Y(d_{k-1})`$, the shares sum to one, and a layer below $D$ holds none. Here $\beta$ = `root_beta` (0.1, ED2's default) is the share the profile would leave below $D$ if it went on. The patch's profile weights the root-zone temperature of fine-root respiration and spreads the uptake when no layer supplies any; it is its cohorts' profiles weighted by their fine-root carbon.

**Soil to root, layer by layer: a single root** (Gardner 1960; Cowan 1965). Water flows radially to a root of radius $`r_{\text{root}}`$ from a soil cylinder whose radius is half the distance between roots, so per plant

```math
g_k=\frac{2\pi\,K_{soil}(k)\,L_k}{\ln\!\left(r_{\text{half},k}/r_{\text{root}}\right)},\qquad
L_k=b_{root}\cdot\text{SRL}\cdot Y_k,\qquad
r_{\text{half},k}=\left(\pi\,\frac{\sum_{\text{cohorts}} n_{plant}\,L_k}{\Delta z_k}\right)^{-1/2} \qquad(14)
```

where $L_k$ is the plant's fine-root length in the layer ($b_{root}$ the fine-root carbon, SRL = `specific_root_length` in m kgC⁻¹), $`r_{\text{half},k}`$ comes from the patch's root-length density in the layer, and $`r_{\text{root}}`$ = `fine_root_radius` (the logarithm is floored at ln 2 for very crowded roots). In the code $g_k$ also carries $`n_{plant}`$ and the factor $`\rho_w/`$`grav_head` (`rhizosphere_conductance`), so $`U_k`$ is per m² of ground with the drive in MPa. The conductance grows with root length and falls only by the logarithm as roots crowd, so it does not depend on how the soil is layered or how a stand is split into cohorts. ED2-hydro's per-layer form, $`K\sqrt{\text{RAI}_k}/(\pi\,\Delta z_k)`$ after Katul et al. (2003), depends on both, and connects the thick deep layers one to two orders of magnitude more weakly. Only part of a fine root may absorb water, but a constant absorbing fraction of the root length is exactly a rescaled specific root length: it scales $L_k$ and the root-length density in $`r_{\text{half},k}`$ alike, so `specific_root_length` is the knob to fit. Uptake that grows less than linearly with root length would be a root conductance in series with $g_k$.

Each layer is its own path from the soil to the wood, so the uptake from layer $k$ is

```math
U_k = g_k\,r\!\left(\psi_{soil,k}+\rho g\,z_k-\psi_W\right),\qquad Q_{\text{root}}=\sum_k U_k \qquad(15)
```

with $r$ zero below zero and the identity above a width of 10⁻² MPa, joined smoothly ([soil_plant_water](soil_plant_water.md) eq. (1)). **Hydraulic redistribution is not represented:** a layer drier than the wood gives nothing and takes nothing, so every $U_k\ge 0$, and there is no per-layer efflux or soil re-wetting by roots.

## 5. The solver

The leaves and wood of every cohort are solved with the soil layers and the pond as one network, by implicit Euler with Newton and adaptive sub-steps, in each tissue's symplastic relative water content $R$: [soil_plant_water](soil_plant_water.md). Each cohort's transpiration is its stomata's demand, limited by the water its leaves hold below the closing potential (§4 there). Water is advanced by the flows themselves, so the plant and soil water budgets close to round-off.

## Parameters

These are the PFT file's `[hydraulics]` keys (record `plant_hydraulics_t`). Each key takes one value, for every PFT, or an array with one value per PFT. Every key is required: the model holds no default and derives no trait from another. Each PFT's traits are one record, `pft%hydro(pf)`, and `derive_pft_tables` builds the fast loop's per-PFT table from those records once at load. Each entry carries its **own** Kirchhoff lookup, built from that PFT's `wood_kexp`: sharing one would give every PFT the first one's vulnerability shape while every budget still closes.

| Symbol | Config key | Meaning |
|---|---|---|
| $\pi_{0}$ | `leaf_pi0`, `wood_pi0` | osmotic potential at full turgor [MPa] |
| $\varepsilon$ | `leaf_elastic_mod`, `wood_elastic_mod` | bulk elastic modulus [MPa] |
| $a_f$ | `leaf_apoplast_frac`, `wood_apoplast_frac` | apoplastic fraction [–] |
| $w_{sat}$ | `leaf_water_sat`, `wood_water_sat` | saturated water content [kg H₂O / kgC] |
| $\psi_{50}$ | `wood_psi50` | xylem potential at 50 % loss of conductance [MPa] |
| $a$ | `wood_kexp` | vulnerability-curve shape [–] |
| $K_s$ | `wood_kmax` | sapwood specific conductivity [kg m⁻¹ s⁻¹ MPa⁻¹]; the shipped values from wood density, eq. (11) (`scripts/prepare_input/pft_traits`) |
| — | `vessel_curl` | tortuosity / path-length factor [–] |
| $\beta$ | `root_beta` | the root profile within the rooting depth [–], $0<\beta<1$ (ED2: 0.1) |
| $`b_{1Rd}, b_{2Rd}`$ | `root_depth_b1`, `root_depth_b2` | rooting depth $`b_{1Rd}\,h^{b_{2Rd}}`$ [m], capped at the soil column (ED2 IALLOM 1) |
| SRL | `specific_root_length` | fine-root length per unit fine-root carbon [m kgC⁻¹] |
| $`r_{\text{root}}`$ | `fine_root_radius` | absorbing-root radius [m] |
| $`b_{1SA}, b_{2SA}`$ | `allometry.sapwood_area_b1`, `allometry.sapwood_area_b2` | ED2 sapwood-area allometry; sets the sapwood ring, which with the fine roots is the wood water store, and alone holds the wood thermal store's internal water |

## Code map

| Concept | Routine |
|---|---|
| PV curves / capacitance | `ed_plant_hydr_lib`: one `water_curve_t` record per tissue; `tissue_psi_to_water`, `tissue_capacitance`, `tissue_water_to_psi`, `tissue_psi_tlp`, `tissue_rwc_to_psi`, `tissue_psi_to_rwc` |
| a tissue's water relations in $R$ | `ed_plant_hydr_lib`: `tissue_water_relations`, `tissue_water_to_rwc` |
| vulnerability curve | `ed_plant_hydr_lib`: `conductance_retained_fraction` |
| Kirchhoff potential and mean | `ed_plant_hydraulics`: `kirchhoff_potential` (closed forms, or the PFT's table); the Kirchhoff mean conductance inside `xylem_sapflow`, from `kirchhoff_difference` |
| Kirchhoff table | `kirchhoff_table_t` (`ed_pft_params`), built at load by `build_kirchhoff_table` (`ed_pft_derive`: one 7-point Gauss–Legendre panel per interval, written out in `kirchhoff_panel`); `kirchhoff_tail` past the table (`ed_plant_hydraulics`) |
| xylem and root conductances | `ed_plant_hydr_lib`: `xylem_max_conductance`, `rhizosphere_conductance`, `root_fraction_profile` |
| xylem and root flows | `ed_plant_hydraulics`: `xylem_sapflow`, `root_water_uptake`, `root_water_uptake_profile`, and the leaves' supply limit `supply_limited_transpiration` |
| rooting depth and profile | `ed_plant_hydraulics`: `cohort_root_depth`, `cohort_root_profile`; `crown_gravity_head`, `leaf_node_psi` |
| the solve | `ed_fast_integ_backward_euler`: `water_step` ([soil_plant_water](soil_plant_water.md)) |
| the per-PFT table (`cfg%pft%hydraulics`) | `ed_pft_derive`: `derive_pft_tables` |

## References
- Xu, Medvigy, Powers, Becknell & Guan (2016), *New Phytologist* — X16 hydraulics.
- Bartlett, Scoffoni & Sack (2012); Tyree & Hammel (1972) — pressure–volume theory.
- Gardner (1960) *Soil Science* 89:63–73; Cowan (1965) *J. Appl. Ecol.* 2:221–239 — single-root rhizosphere conductance.
- Christoffersen (2013), PhD thesis, University of Arizona — the rooting-depth allometry (ED2 IALLOM 1).
- ED2 `ED/src/dynamics/plant_hydro.f90`.
