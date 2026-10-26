# PFT traits from published fits

ED-MORE reads every PFT value from the PFT file. It holds no rule that derives one trait from another,
because such rules are fits to particular data sets, not laws. [`pft_traits.py`](pft_traits.py)
computes the traits that four published relations give, so a PFT file can be written from them:

| PFT-file key | From | Fit |
|---|---|---|
| `[hydraulics].wood_kmax` [kg m⁻¹ s⁻¹ MPa⁻¹] | wood density ρ | Xu et al. (2016, New Phytol 212:80): ln K<sub>s</sub> = 2.348 − 2.455 ρ, lognormal over 72 neotropical dry-forest species; the mean exp(2.348 − 2.455 ρ + 0.6186/2), with ρ outside 0.35–0.95 g cm⁻³ taking the nearest edge ([plant_hydraulics](../../../docs/science/plant_hydraulics.md) eq. (11)) |
| `[pft].mort_gamma`, `mort_alpha`, `mort_beta` | wood density ρ | Camac et al. (2018, PNAS), fitted at Barro Colorado Island: each parameter is param₀ (ρ / 0.6)<sup>exp</sup> |
| `[pft].kplastic_vcmax`, `kplastic_rd`, `kplastic_sla` | Vcmax25 | Lloyd et al. (2010) as ED2's `trait_plasticity_scheme = 2`: k<sub>Vcmax</sub> = −exp(−2.788 + 0.01439 Vcmax25), Rd follows Vcmax, SLA rises by 1.18/1.10 of it ([plant_traits](../../../docs/science/plant_traits.md) eq. (3)) |
| `[pft].kplastic_leaf_turnover` | top-of-canopy leaf turnover rate | minus ED2's lifespan slope, 0.062 ln(12 / rate) − 0.2126, capped at 0 so a slow-turnover PFT keeps one rate through the canopy ([plant_traits](../../../docs/science/plant_traits.md) eq. (4)) |

```bash
python pft_traits.py --pft-file ../../../ed_config_pft.toml
python pft_traits.py --wood-density 0.40 0.60 0.85 --vcmax25 60 45 40 --leaf-turnover 1.0 0.5 0.33
```

With `--pft-file`, the inputs are the file's `[pft]` `wood_density`, `vcmax25` and
`leaf_turnover_rate_toc`; a value given on the command line takes precedence, and every input needs one
value per PFT. The output is the TOML lines for the PFT file's `[pft]` and `[hydraulics]` blocks, one
value per PFT, each printed to the digits that read back as the same number.

The shipped PFT files carry these fits' values to about 10⁻¹⁵. Example 05's PFT files carry their own
plasticity slopes (Ma et al. 2025, from Barro Colorado's crane), and the Barro Colorado file its own
mortality hazard, fitted to the censuses.

A fit's coefficients are the constants at the top of the script. A site with its own fit changes them
there, or writes its PFT file's values directly: ED-MORE reads the values, not the fit.
