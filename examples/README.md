# Examples

Six runnable examples of ED-MORE: one module on its own (01–03), the coupled column at a flux tower (04),
and the whole model at a site and over a region (05–06). Each folder holds the config, the scripts, the
figures and a README with the commands that reproduce it, run from the repository root. Data the
repository does not hold (the processed ERA5-Land archive, the BCI census tables, the soil maps) is
named on the command line or by an environment variable such as `ERA5LAND_ARCHIVE`; each README says
which.

| Example | What it shows |
|---|---|
| [`example01_leaf_gas_exchange/`](example01_leaf_gas_exchange/) | The leaf gas-exchange module on its own, driven from Python: C3 photosynthesis, Medlyn stomata and the coupled Cᵢ solve, against the A–Cᵢ curves and leaf-temperature responses of four tropical tree species (Slot & Winter 2017). |
| [`example02_canopy_phenology/`](example02_canopy_phenology/) | The leaf-phenology module on its own, driven from Python: one kernel, four forests, four leaf habits, each site's light from ERA5-Land. Deciduous broadleaf (Harvard Forest), evergreen pine (Hyytiälä), a tropical tree that renews its leaves in the dry season (Barro Colorado Island) and a drought-deciduous dry forest (Palo Verde). The sites differ only in parameter values. |
| [`example03_demography/`](example03_demography/) | The demography module on its own, driven from Python: growth, mortality (Camac et al. 2018) and recruitment laws fitted to the censuses of the Barro Colorado Island 50-ha plot, run through the cohort-and-patch engine from the 1985 census and from near-bare ground. |
| [`example04_column_biophysics/`](example04_column_biophysics/) | The coupled column (canopy radiation, leaf gas exchange, energy balances, plant hydraulics, soil water) every 15 minutes at the Barro Colorado Island flux tower, Panama. It starts from the plot's 2010 census, calibrates the fast parameters against the tower, compares the default and calibrated runs with the tower's fluxes, and shows ten days at half-hourly output. |
| [`example05_forest_regeneration/`](example05_forest_regeneration/) | The whole model at once: the column of example 04 with the demography of example 03, regrowing forests from bare ground on repeated ERA5-Land weather and the CO₂ of its history. Barro Colorado Island from 1600 to 2020 (three successional PFTs, checked against the census and the carbon budget), and Palo Verde's dry forest for 150 years (four leaf habits after Xu et al. 2016, set against the plots' 2009 inventory). |
| [`example06_regional_run/`](example06_regional_run/) | The whole model over a region: example 05's four dry-forest leaf habits regrow from bare ground in every cell of a Central American box, run as ED-MORE region runs with each cell's own ERA5-Land weather, soil and height allometry. The leaf area's seasonality is set against MODIS and the biomass against Neotropical chronosequences (Poorter et al. 2016). |
