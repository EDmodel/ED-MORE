# A site's soil horizons

A site's `[soil]` may give each layer its own texture: a list per key, one value per layer
([configuration](../../../docs/configuration.md), "The soil layers"). [`soil_horizons.py`](soil_horizons.py)
computes those values from a soil pit's horizons by the model's own formulas, so that a site's script
prints what ED-MORE would derive from the same texture:

| Function | Gives | From |
|---|---|---|
| `horizon_hydraulics` | `theta_sat`, `theta_res`, `curve_par_a`, `curve_par_n`, `ksat` | sand, silt, clay, organic carbon, CEC, pH and bulk density: Hodnett & Tomasella (2002) and Ottoni et al. (2019), the coarse fragments after Saxton & Rawls (2006), as `soil_texture_to_hydraulics` (`ed_soil_hydr_lib.f90`) |
| `texture_thermal` | `solid_conductivity`, `dry_conductivity`, `dry_heat_capacity` | sand, clay and porosity, as `soil_texture_to_thermal` (`ed_therm_lib.f90`): the solids conduct by their quartz, taken as the sand (Johansen 1975; Peters-Lidard et al. 1998) |
| `organic_fraction`, `organic_thermal` | the three thermal keys of an organic horizon | the mineral's values blended with peat's by the organic fraction, organic matter per volume over 130 kg m⁻³, as CLM5 does (Lawrence & Slater 2008) |
| `read_megapit` | each horizon's depths, texture, chemistry and bulk density | NEON's soil megapit (DP1.00096.001) |

The sites' scripts import it: example 04's [`harvard_traits.py`](../../../examples/example04_column_biophysics/sites/harvard/harvard_traits.py)
and [`toolik_traits.py`](../../../examples/example04_column_biophysics/sites/toolik/toolik_traits.py).
The tests (`tests/`) check the formulas against the model's and the blend's two ends.
