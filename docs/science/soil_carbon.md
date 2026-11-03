# Soil carbon

The decomposable carbon of one patch is a seven-element vector $X$, advanced once per slow step (a day
by default) by a CENTURY-family **carbon matrix ODE**. This page gives the pools, the matrix that moves
carbon between them, the environmental scalar that sets the rates, and how the slow step and the fast
loop share the pools without counting the respiration twice.

The model is a port of ED2's `update_C_and_N_pools` (the same pools, decay rates, and respired and
transferred fractions), written as one linear-algebra object instead of a cascade of hand-written
pool-to-pool transfers. With the decay diagonal $K$, the transfer matrix $A$ and the environmental
scalar $\xi$ as matrices, the exact matrix-exponential step, the semi-analytic steady state (SASU) and
the traceability diagnostics all come from the same operator.

## 1. The carbon matrix ODE

Per patch, with $X$ in $`[\mathrm{kgC\,m^{-2}}]`$, and element by element the mass balance it encodes:

```math
\frac{dX}{dt} = B\,I(t) + A\,\xi(t)\,K\,X(t) \quad [\mathrm{kgC\,m^{-2}\,day^{-1}}],
\qquad
\frac{dX_i}{dt} = u_i + \sum_{j\neq i} a_{ij}\,\xi_j K_j X_j - \xi_i K_i X_i \qquad(1)
```

$`u = B\,I`$ is the litter input (§5); $K$ the diagonal of baseline decay rates, carried in
$`[\mathrm{yr^{-1}}]`$ by the config and divided by `yr_day` at assembly so the kernels run in
$`[\mathrm{day^{-1}}]`$; $\xi$ the diagonal per-pool environmental scalar (§4); $A$ the dimensionless
donor/transfer matrix (§3). The input is litter carbon, never GPP: autotrophic respiration is taken out
upstream, in the plant carbon budget.

## 2. The seven pools and the lignin sub-tracer

The state is `soil_carbon_state_t`, ordered litter → soil organic matter (SOM) → passive. The seven
pools are ED2's CENTURY set. The litter half is a **metabolic/structural × above/below** cross:
chemistry on one axis, position on the other.

| `IP_*` index | pool | name (`SOIL_POOL_NAME`) | `[soil_carbon]` rate, default $`[\mathrm{yr^{-1}}]`$ | gets litter | $\xi$ temperature |
|---|---|---|---|---|---|
| `IP_FAST_GRND` 1 | metabolic litter, above | `fast_grnd` | `k_fast` 11.0 | yes | ground surface |
| `IP_FAST_SOIL` 2 | metabolic litter, below | `fast_soil` | `k_fast` 11.0 | yes | root-weighted soil |
| `IP_STRUCT_GRND` 3 | structural litter + CWD, above | `struct_grnd` | `k_struct` 4.5 | yes | ground surface |
| `IP_STRUCT_SOIL` 4 | structural litter + CWD, below | `struct_soil` | `k_struct` 4.5 | yes | root-weighted soil |
| `IP_MICR` 5 | microbial SOM | `microbial` | `k_micr` 0.0 | no | root-weighted soil |
| `IP_SLOW` 6 | slow / humified SOM | `slow` | `k_slow` 0.2 | no | root-weighted soil |
| `IP_PASSIVE` 7 | passive SOM | `passive` | `k_pass` 0.0 | no | root-weighted soil |

The state holds the pools as `carbon(N_SOIL_POOL)` indexed by `IP_*`, so the kernels read and write the
vector directly.

A pool is **active** when $`K_j>0`$. Under the default `DECOMP_SCHEME_ED2` (ED2 schemes 0–4) the
microbial and passive pools are inert, leaving three active pools: fast, structural and slow.
`DECOMP_SCHEME_CENTURY5` activates all seven, the five-active CENTURY topology. The selector switches
the topology and the temperature form, not the rate constants: those come from `[soil_carbon]` and
default to the scheme 0–4 values. Dead wood has no pool of its own; it enters the structural pools
through `f_labile_stem`, ED2's convention.

The structural pools carry a **lignin sub-tracer**: `lignin(1)` above ground and `lignin(2)` below,
bounded $`0\le L_s\le C_s`$. It is not an eighth pool. Its carbon is already counted inside the
structural pools, so it stays outside the mass vector and keeps its own balance (§6). There is no
nitrogen (§10).

## 3. The transfer matrix, and the respired complement

`decomposition_transfer_matrix` builds $A$, $K$ and the respired fractions together. Each donor loses *all*
of its decayed carbon, so the diagonal is uniform; Rh is exactly what $A$ does *not* transfer:

```math
a_{jj} = -1 \ \ \forall j, \qquad e_{r,j} = 1 - \sum_{i\neq j} a_{ij} = -\sum_i a_{ij},
\qquad R_h = \sum_j e_{r,j}\,\xi_j K_j X_j = -\mathbf{1}^{\mathsf{T}} A\,\xi\,K\,X \qquad(2)
```

with $`a_{ij}\in[0,1]`$ the fraction of $j$'s loss entering pool $i$ and $`e_{r,j}`$ — the **respired
complement** — the fraction leaving the network as CO₂. Keeping $`a_{jj}=-1`$ on inert pools too is
deliberate: $`K_j=0`$ zeroes them anyway, and the column identity $`\sum_i a_{ij}=-e_{r,j}`$ then holds
for every pool, which the tests check.

**Default scheme.** Metabolic litter respires 100% and transfers nothing. Structural litter sends its
complement to slow ($`a_{6,3}=1-e_{r,3}`$, $`a_{6,4}=1-e_{r,4}`$). Slow is the only SOM outlet and
respires 100%. The structural respired fraction is lignin-weighted,
$`e_{r,s}=(1-f_{\ell,s})\,r_{\mathrm{nonlig}}+f_{\ell,s}\,r_{\mathrm{lig}}`$ with
$`f_{\ell,s}=L_s/C_s`$ (0 when $`C_s=0`$) and $`r`$ the keys `er_struct_nonlig` and `er_struct_lig`.
Both default to 0.3, so structural litter sends 70% of its decayed carbon to slow.

**Scheme 5** (the Bolker, Pacala and Parton CENTURY). Fast litter goes to microbial. Structural
non-lignin goes to microbial and lignin to slow. Microbial ⇄ slow and slow ⇄ passive exchange carbon
both ways, controlled by soil texture: sand sets the microbial respired fraction
$`e_{r,\mathrm{micr}}=\texttt{er\_micr\_int}+\texttt{er\_micr\_slp}\cdot x_{\mathrm{sand}}`$, and clay sets
the fractions routed to passive, $`\min(1-e_r,\ \texttt{int}+\texttt{slp}\cdot x_{\mathrm{clay}})`$
(the `fx_micr_pass_*` and `fx_slow_pass_*` keys). The back-transfers sit above the diagonal, so scheme
5's active block is not triangular (§8).

## 4. The environmental scalar, and where it sits in the product

`decomposition_env_scalar` returns one scalar per pool: temperature × moisture and oxygen × (structural
pools only) a lignin brake. The moisture factor reads the relative saturation $`m`$, clipped to
$`[0,1]`$.

```math
f_T = \begin{cases} q_{10}^{\,(T-T_{\mathrm{ref}})/10} & \text{scheme 5, UNCAPPED Q10}\\
\min\!\big(1,\ \exp[a_T(T-T_{\mathrm{sat}})]\big) & \text{schemes 0--4, ED2 capped exp}\end{cases}
\qquad
f_\theta = \begin{cases} \exp[(m-m_{\mathrm{opt}})\,b_{\mathrm{dry}}] & m \le m_{\mathrm{opt}}\\
\exp[(m_{\mathrm{opt}}-m)\,b_{\mathrm{wet}}] & m > m_{\mathrm{opt}}\end{cases}
\qquad(3)
```

```math
m = \frac{\theta-\theta_{\mathrm{dry}}}{\theta_{\mathrm{sat}}-\theta_{\mathrm{dry}}}, \qquad
\xi_j = f_T(T_j)\,f_\theta(\theta)\,\Lambda_j, \qquad
\Lambda_j = \begin{cases}\exp(-e_{\mathrm{lig}}f_{\ell,j}) & j\in\{3,4\}\\
1 & \text{otherwise}\end{cases} \qquad(4)
```

A one-sided exponential on each side of the optimum covers both moisture limbs: substrate diffusion on
the dry side, anoxia on the wet side. The fast Rh kernel uses the same form, so the two timescales
agree at matched inputs. In eq. (3) and eq. (4):
- $`T_j`$ is the ground surface temperature for the above-ground pools (1, 3) and the root-weighted soil
  temperature for the rest (ED2's `A_decomp` / `B_decomp` split);
- $\theta$ is the depth-weighted soil mean, $`\theta_{\mathrm{dry}}`$ and $`\theta_{\mathrm{sat}}`$ the
  top layer's residual and saturated water contents;
- $`f_{\ell,j}=L_j/C_j`$ is the lignin fraction of structural pool $j$ (0 when the pool is empty),
  $`e_{\mathrm{lig}}`$ is `e_lignin`, $`q_{10}`$ is `rh_q10`, $`T_{\mathrm{ref}}`$ is `rh_t_ref`,
  $`a_T`$ is `resp_temp_increase`, $`T_{\mathrm{sat}}`$ is `resp_temp_ref`, $`m_{\mathrm{opt}}`$ is
  `resp_opt_water`, and $`b_{\mathrm{dry}}`$, $`b_{\mathrm{wet}}`$ are `resp_water_below_opt`,
  `resp_water_above_opt`.

Under the default scheme $\xi\le1$; under scheme 5 only $\xi>0$ holds, because its Q10 is uncapped.

**Placement is physics.** The scalar sits inside the transfer, $`A\,(\xi K X)`$, never
$`\xi\,(A K X)`$: each donor decomposes at its own rate $`\xi_j K_j`$, and only then is the loss split
into respired and transferred parts. $\xi$ and $K$ commute with each other but neither commutes with
$A$, so the receiver-scalar ordering silently breaks mass balance and every capacity and residence-time
formula. The tests assert it: the correct form's total loss equals $`R_h`$ of eq. (2); the other's does
not.

## 5. Litter input

Dead plant carbon becomes litter in one place — `necromass_to_litter` (`ed_soil_carbon`), whose
arguments are plain carbon amounts and plain PFT traits, which is what lets both the slow driver and the
demographic operators call it. Each tissue goes where it is:

| Tissue | Labile share | Above / below ground |
|---|---|---|
| leaves; flowers, fruit and seed of reproduction | `f_labile_leaf` | all above |
| fine roots | `f_labile_leaf` | all below |
| wood (stem + coarse roots) | `f_labile_stem` | `aboveground_frac` / the rest |
| storage (non-structural carbon, in sapwood and coarse roots) | `f_labile_leaf` | `aboveground_frac` / the rest |

`aboveground_frac` is the woody split (ED2 `agf_bs`), so it applies to wood and storage only. Only the
structural streams carry lignin: `struct_lignin_frac` of each. The slow driver calls it through
`add_necromass`, which also records the carbon by the tissue it came from. Those records are the litter
outputs (`mean_litter_leaf_py`, `mean_litter_fineroot_py`, `mean_litter_wood_py` for wood and storage,
`mean_litter_repro_py`, each after resorption); with `mean_root_exudate_py` they sum to the carbon
entering pools 1–4.

`build_litter_input` maps the four carbon streams onto $u$ and the two lignin fluxes alongside. It does
no chemistry of its own, because the splits are per-cohort decisions made before summing to the patch.
So only pools 1–4 ever receive litter; microbial, slow and passive are fed by transfers alone.

Carbon arrives by one route: turnover and shedding, background mortality, the seed that did not
establish ([demography](demography.md) §1), and the plants a calendar boundary's cull and treefall kill
all accumulate into `patch%litter_in` through `add_necromass`. The root exudate joins the below-ground
metabolic stream directly. The accumulated litter enters through $u$ at the next soil step, and the
accumulator is then emptied. A boundary's dead carbon therefore reaches the pools with the first soil
step after it.

## 6. The soil step, and why it consumes an accumulated integral

`soil_carbon_step` advances one patch by one slow step. The decay it applies is not an instantaneous
$\xi$. It is the fast loop's **accumulated per-pool integral** $`\xi^{\mathrm{int}}`$ over the step,
with the dimensionless decay fraction $`d`$. The default solver (`DECOMP_STEP_EULER`) is forward Euler
on the frozen operator, as in ED2 ($`a_{ii}=-1`$ removes donor $i$'s own loss, the off-diagonals
transfer):

```math
\xi^{\mathrm{int}}_j = \int_{\mathrm{step}} \xi_j(t)\,dt \ [\mathrm{day}], \qquad
d_j = \xi^{\mathrm{int}}_j K_j \ [-], \qquad
X^{n+1}_i = X^n_i + u_i + \sum_j a_{ij}\,d_j\,X^n_j \qquad(5)
```

For steps where the decay fraction is not small, `DECOMP_STEP_EXPM` takes the exact exponential of that
same frozen operator. An augmented matrix carries the constant input without inverting $A$, which
matters because $`A\,\mathrm{diag}(d)`$ is singular whenever a pool is inert:

```math
\begin{bmatrix} X^{n+1}\\ 1\end{bmatrix}
= \exp\!\left(\begin{bmatrix} A\,\mathrm{diag}(d) & u\\ 0 & 0\end{bmatrix}\right)
\begin{bmatrix} X^{n}\\ 1\end{bmatrix}, \qquad
R_h^{\mathrm{day}} = \sum_i u_i - \Big(\sum_i X^{n+1}_i - \sum_i X^{n}_i\Big) \qquad(6)
```

The respiration $`R_h^{\mathrm{day}}`$ is the actual loss to the atmosphere, computed the same way for
both solvers: inter-pool transfers cancel in the total, so it is what came in minus what stayed. Lignin
decays in lockstep with its host structural carbon, $`L_s\leftarrow L_s\sigma_s+\ell_s`$ with
$`\ell_s`$ the step's lignin input and $`\sigma_s=1-d_s`$ (Euler) or $`\exp(-d_s)`$ (EXPM). So $`f_\ell`$
is unchanged by pure decay and moves only through inputs; the linear form under EXPM would drive $L$
negative at large $`d_s`$.

**Why the integral is the right object.** $\xi$ is a product of nonlinear functions of temperature and
moisture, both of which swing over a day, so by Jensen's inequality
$`\langle f_T(T) f_\theta(\theta)\rangle \neq f_T(\langle T\rangle) f_\theta(\langle\theta\rangle)`$. A
step that recomputed $\xi$ from the mean state would debit the pools by a different amount than the
sub-daily flux delivered to the atmosphere, and ecosystem carbon would drift while every individual
budget still closed. Integrating $\xi$ as it happens removes the gap by construction. The tests run a
diurnal cycle through both paths to show that the mean-state one does not agree.

The step also reports $`\lambda=\max_j d_j`$, the **largest fraction of any pool this step withdraws**.
It says whether freezing the pools across the slow step is sound. It is scale-free: the flux is linear
in the pool, so $\lambda$ does not grow as a pool empties, which is why bare-ground spin-up works. Under
the default scheme $\xi^{\mathrm{int}}\le1$ day, so $\lambda\le$ `k_fast`/`yr_day` for a one-day step
(about 0.03 at the default rate); past 1, forward Euler would drive that pool negative.

## 7. The fast/slow seam

Heterotrophic respiration must reach the canopy air every fast step, but the pools may be written only
once per slow step. Counting the respiration at both rates would double-count it. ED-MORE has both ends
read the same object. At the start of each fast window the fast driver copies the patch's pools into
the fast state, and that **frozen** copy is what every fast step respires, through the same kernel the
slow step's bookkeeping uses:
- each fast step, `heterotrophic_respiration` evaluates eq. (2) on the frozen pools, and the frozen
  build converts $`R_h`$ to $`[\mathrm{\mu mol\,m^{-2}\,s^{-1}}]`$ with `kgcday_2_umols` for the canopy
  air's CO₂;
- the driver integrates $`\xi_j`$ into the per-patch `decomp_xi%xi_days` (the $`\xi^{\mathrm{int}}`$ of
  eq. (5)), and integrates $`R_h`$, still in the pools' kgC m⁻² day⁻¹ so the seam check compares like
  with like, into `decomp_xi%rh_fast`. Both are reset at the start of each fast window.

Within a slow step, `soil_carbon_step` is the only change to the real pools. It debits them with
exactly that $`\xi^{\mathrm{int}}`$, returns the new pools, and `ed_slow_update` writes them. At a
calendar boundary, patch fusion and the treefall gaps blend the pools on area
([demography](demography.md) §4), and the plants the boundary kills reach them through the litter (§5).

The window's total fast Rh therefore equals the soil step's pool debit by construction, and
$`\texttt{rh\_seam\_gap} = R_h^{\mathrm{day}} - \sum_k R_h(t_k)\Delta t_k \approx 0`$ is an
**assertion guard**, never a live correction. It closes to machine precision, and the run reports its
worst value. The two checks differ in what they catch: $\lambda$ (§6) catches the frozen-pool
approximation degrading, and `rh_seam_gap` catches the contract being broken. This is the single
documented fast–slow kernel seam in the model: one kernel called from both timescales, co-located with
the pools it debits.

The boundary's fusion and gaps run before the next fast window opens, so the window always starts from
the patches as they are ([order_of_processes](order_of_processes.md) §5). A nonzero gap on any step is
a broken contract; do not widen the tolerance to absorb it.

**Without the fast loop** (`[fast].fast_biophysics_on = false`, a demography spin-up) there is no soil
climate to integrate and no seam. The soil step takes the fixed scalar `spinup_xi` for every pool,
$`\xi^{\mathrm{int}}_j = \bar\xi\,\Delta t_{\mathrm{slow}}`$ (in days), as growth takes the fixed
`gpp_ref`, and the pools tend to the steady state of §8 for that scalar and the run's litter. Nothing
respires into a canopy air space, so the soil step's own debit is the respiration the slow ledger
declares.

## 8. Steady state (SASU) and traceability diagnostics

Under climatological-mean drivers the equilibrium is a direct linear solve rather than a
multi-millennial integration. The same active-block machinery, evaluated at the current $\xi$ and $u$,
gives three traceability diagnostics (eq. (8)):

```math
0 = \bar u + A\,\bar\xi K X_{ss} \Longrightarrow M X_{ss}[\mathcal{A}] = -\bar u[\mathcal{A}],\ \
M = \big(A\bar\xi K\big)[\mathcal{A},\mathcal{A}],\ \ \mathcal{A}=\{j: K_j>0\} \qquad(7)
```

```math
X_c = -\big(A\,\xi K\big)^{-1}u, \qquad X_p = X_c - X, \qquad
\tau_E = \frac{\sum_i X_{c,i}}{\sum_i u_i}\Big/\texttt{yr\_day} \qquad(8)
```

**The solve runs on the active sub-block only, because the full 7×7 system is singular.** An inert
pool's column of $`A\,\bar\xi K`$ is identically zero, so the inverse does not exist and a blind
back-substitution divides by zero on its diagonal. This is the default scheme, not a corner case. Inert
pools are held at zero. The reduced block is solved by a small dense LU factorization with partial
pivoting, free of allocation and of LAPACK, and valid for scheme 5's non-triangular back-transfers as
well as the default block.

A run can start from eq. (7) with `[soil_carbon].spinup_steady`, in any start but a restart. The start
uses the scalar $`\bar\xi`$ = `spinup_xi` for every pool and a constant litter input
(`spinup_labile_grnd`, `spinup_labile_soil`, `spinup_struct_grnd`, `spinup_struct_soil`
[kgC m⁻² day⁻¹]), on the transfer matrix of empty structural pools (lignin fraction 0). The lignin
tracer starts at zero. Two caveats: those means are prescribed, not fed back from a live canopy, and a
mean-driver equilibrium discards the covariance of periodic forcing.

In eq. (8), $`X_c`$ is the **storage capacity**: the equilibrium the patch is chasing under today's
climate and litter input, a moving attractor rather than a fixed point. $`X_p`$ is the **storage
potential**, the remaining sink (negative means a source). $`\tau_E`$ is the ecosystem **residence
time** in years, stock over throughput. The transient relaxes toward $`X_c`$ at a rate set by
$`A\,\xi K`$, so the three together explain why two runs store different carbon. The model does not
report them; they follow from the pools, the litter input and the scalars.

## 9. Audits, and the difference between conservation and plausibility

Every soil step fills a `soil_carbon_out_t`: the new pools, `rh_seam_gap` (§7), $\lambda$ with the pool
that attained it, and the lignin tracer's residual
$`\max_s|\Delta L_s-(\ell_s-(1-\sigma_s)L_s)|`$, zero by construction for either solver. The carbon
check that can actually fail is the independent one, in the tests: recompute $`R_h`$ from
$`e_r\,\xi^{\mathrm{int}}K X^n`$ and compare it with $`\sum u-\Delta X`$. Comparing eq. (6) with
itself cannot catch dropped litter or a misattributed transfer.

Conservation is not plausibility: a runaway that conserves carbon closes every budget while moving
absurd amounts of it. `soil_carbon_bad_pool` answers what a budget cannot. It returns 0 when every
pool is physically possible, else the index of the first that is not: NaN, negative beyond a round-off
tolerance of $`10^{-9}\,\mathrm{kgC\,m^{-2}}`$, or above a divergence ceiling of
$`10^{4}\,\mathrm{kgC\,m^{-2}}`$ (about 500 times the richest real soil; a temperate forest's holds
20–40). It returns an index, not a flag, so the caller's message can name the pool, and which pool went
bad is the first clue to which operator did it. Round-off negatives must pass, since a daily Euler step
produces them routinely; the ceiling catches divergence, not implausible science.

## 10. What is not here

- **Nitrogen.** The soil is carbon-only: no nitrogen pools, no C:N stoichiometry, no nitrogen
  limitation of decomposition or of plant growth.
- **DAMM, or any other respiration form.** The matrix form (§3, §7) is the only heterotrophic
  respiration. Davidson's dual Arrhenius/Michaelis–Menten scheme (soluble-C and O₂ limitation, a
  unimodal moisture response that emerges rather than being imposed) is not available, nor are separate
  Q10 or capped-exponential respiration kernels: those two temperature responses live inside $\xi$,
  eq. (3).
- **No vertical resolution.** Pools are per-patch scalars, as in ED2; only the environmental scalar
  sees layer-resolved temperature and moisture.
- **No fire**, so nothing consumes the above-ground pools; treefall is the only disturbance.

## Parameters

The `[soil_carbon]` keys; every one is optional, at the default shown in `ed_config_main.toml`.

| Key | Meaning |
|---|---|
| `soil_carbon_on` | master switch (default `true`) |
| `decomp_scheme`, `step_solver` | `"ed2"` (default) or `"century5"`; `"euler"` (default) or `"expm"` |
| `k_fast`, `k_struct`, `k_micr`, `k_slow`, `k_pass` | baseline decay rates $`[\mathrm{yr^{-1}}]`$ |
| `er_fast`, `er_struct_nonlig`, `er_struct_lig`, `er_slow`, `er_pass` | respired fractions, eq. (2) and §3 |
| `e_lignin` | the lignin brake exponent, eq. (4) |
| `er_micr_int`, `er_micr_slp`, `xsand`; `fx_micr_pass_int`, `fx_micr_pass_slp`, `fx_slow_pass_int`, `fx_slow_pass_slp`, `xclay` | scheme-5 texture control: sand sets the microbial respired fraction, clay the routing to passive |
| `rh_q10`, `rh_t_ref`, `resp_temp_increase`, `resp_temp_ref`, `resp_opt_water`, `resp_water_below_opt`, `resp_water_above_opt` | the $\xi$ responses, eq. (3) |
| `spinup_steady`, `spinup_xi`, `spinup_labile_grnd`, `spinup_labile_soil`, `spinup_struct_grnd`, `spinup_struct_soil` | the steady-state start, §8; `spinup_xi` is also the scalar of a run without the fast loop, §7 |

The defaults are the `"ed2"` scheme's. `decomp_scheme = "century5"` changes the topology and the
temperature form but not the rates, so set them with it. ED2's scheme-5 values are `k_fast` 12.0,
`k_struct` 1.5, `k_micr` 6.0, `k_slow` 0.15 and `k_pass` 0.012 yr⁻¹; `er_fast` 0.55,
`er_struct_nonlig` 0.50, `er_struct_lig` 0.20, `er_slow` 0.55 and `er_pass` 0.55; `er_micr_int` 0.60
and `er_micr_slp` 0.17.

The litter-partition traits are per PFT, in the PFT file: `[pft].f_labile_leaf`, `f_labile_stem`,
`struct_lignin_frac`, and `aboveground_frac` (in the allometry block).

**`soil_carbon_on = false` is no soil carbon, not a coarser one.** Litter is discarded at the slow
step and heterotrophic respiration is zero, so ecosystem respiration carries only its autotrophic
part. A prescribed pool respired through an empirical Q10 form would not be a coarser approximation
either: the pools are empty and take no litter, so its CO₂ would be carbon from nothing.

## Code map

| Concept | Routine |
|---|---|
| $\xi$; and $A$, $K$, $`e_r`$ | `ed_respiration`: `decomposition_env_scalar`, `decomposition_transfer_matrix` |
| litter → source term $u$; necromass → destinations | `ed_soil_carbon`: `build_litter_input`, `necromass_to_litter`, `add_necromass` |
| the pool advance (Euler / EXPM) and its audits | `ed_soil_carbon`: `soil_carbon_step`, `soil_carbon_out_t` |
| instantaneous $`R_h`$ (the fast–slow seam kernel) | `ed_respiration`: `heterotrophic_respiration` |
| steady state (SASU) and its cold start | `ed_soil_carbon`: `solve_soil_carbon_steady_state`, `soil_carbon_steady_pools`; `ed_init`: `init_polygon` |
| pool plausibility | `ed_soil_carbon`: `soil_carbon_bad_pool` |
| the pool layout and names | `ed_biogeochem_params`: `N_SOIL_POOL`, `IP_*`, `SOIL_POOL_NAME` |
| the state: pools, lignin, the day's integrals, the litter accumulator | `ed_patch_state_vars`: `soil_carbon_state_t` (`carbon(N_SOIL_POOL)`, `lignin(2)`), `decomp_xi_t`, `litter_input_t` |
| `[soil_carbon]` block and selectors | `ed_biogeochem_params`: `soil_carbon_params_t`, `DECOMP_SCHEME_*`, `DECOMP_STEP_*` |
| the slow step's driver, seam guard and ledger entries | `ed_slow_state`: `advance_soil_carbon`; `ed_slow_update`: `commit_soil_carbon`; `ed_slow_ledger`: `slow_ledger_soil_check` |
| litter accumulation | `ed_slow_state`: `compute_carbon_allocation`, `accumulate_mortality_litter`, `kill_litter`; the kills come from `ed_slow_update`: `terminate_cohorts`, `apply_patch_disturbance` |
| fast-step $`R_h`$ and the $`\xi`$ accumulation | `ed_fast_integ_frozen`: `patch_heterotrophic_respiration`; `ed_fast_state`: `accumulate_fast_to_slow`, `reset_fast_to_slow` |
| output | `mean_soilc_<pool>_py`, `mean_soilc_total_py`, `mean_soilc_total_pa`, `mean_rhet_py`, `mean_litter_*_py`, `mean_root_exudate_py`; the state file's `soilc_<name>_pa`, `soilc_lignin_grnd_pa`, `soilc_lignin_soil_pa` |
| tests | `test/fast/test_respiration.f90` ($\xi$, $A$, $K$, $`e_r`$, $`R_h`$ and the Jensen seam), `test/slow/test_soil_carbon.f90` (the soil step), `test/slow/test_slow_step.f90` (the seam through the driver) |

## References

- ED2 `../ED2/ED/src/`: `dynamics/soil_respiration.f90` (environmental scalars, sub-daily
  accumulation, `update_C_and_N_pools`), `memory/decomp_coms.f90` + `init/ed_params.f90` (pools, decay
  rates, respired/transfer fractions).
- Parton et al. (1987), *SSSAJ* 51:1173 — the CENTURY pool structure. Bolker, Pacala & Parton (1998),
  *Ecol. Appl.* 8:425; Koven et al. (2013), *Biogeosciences* 10:7109 — the five-active topology.
- Luo et al. (2022), *JAMES* — the matrix formulation; Xia et al. (2012), *GMD* 5:1259 — semi-analytic
  spin-up; Luo et al. (2017), *Biogeosciences* 14:145 — capacity and potential; Davidson et al. (2012),
  *Glob. Change Biol.* 18:371 — DAMM.
