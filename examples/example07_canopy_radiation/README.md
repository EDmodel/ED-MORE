# Example 07 — Canopy radiation against vertical PAR sensors (Harvard Forest) — PROTOTYPE

ED-MORE's two-stream canopy radiation ([`docs/science/canopy_radiation_transfer.md`](../../docs/science/canopy_radiation_transfer.md)),
driven on its own from Python, against the light measured inside a forest. The site is the NEON tower at
Harvard Forest (Massachusetts), with PAR sensors at six heights and three on the forest floor, inside the
35-ha ForestGEO plot. One year, 2019, every daytime half hour.

This is a prototype. The pipeline runs end to end, but the fetch script, the choices below and the README
are not final.

## Run

```bash
cd examples/example07_canopy_radiation
python run_example.py --data <data> --hgt-max 28    # about 30 s
python plot_example.py                               # canopy_radiation.png in output/
```

It needs a built `libedmore.so` (`build-dev/`, or `EDMORE_LIB`), numpy, pandas and matplotlib.

## The data

`<data>` holds what the fetch scripts download:
- **NEON HARV:** a NEON API token is needed since 2026; the data are CC BY 4.0.
  - PAR at the tower levels, 30 min (DP1.00024.001): 0.36, 5.13, 16.13, 23.63, 29.37 and 39.09 m in 2019.
  - The forest-floor line-quantum sensors of soil plots 1, 3 and 5 (DP1.00066.001), 41–213 m from the tower.
- **Harvard Forest data archive (CC0):**
  - the EMS tower's direct and diffuse PAR at 29 m (HF004), 134 m from the NEON tower;
  - the 2019 census of the ForestGEO plot (HF253);
  - the phenology of tagged trees (HF003).

## The model

- **Drivers.** The measured PAR above the canopy (39 m), split into beam and diffuse by the EMS tower's
  direct/diffuse sensor, with the model's own sun.
  - NEON's own direct/diffuse sensor was shaded by the lightning rod until 2025-06.
- **Stand.** Around each sensor, every live stem within 40 m of the 2019 census becomes a cohort in its 20 m
  quadrat's patch. Two PFTs, deciduous and evergreen (the conifers), with the default optics.
  - Height, leaf area and wood area come from ED-MORE's default allometry, with the height capped at 28 m. The
    deciduous trees NEON measured reach 28.5 m, and lidar puts the canopy top at 24–26 m around the tower.
- **The deciduous leaves through the year.** Each species' leaf area follows the ground observations of
  tagged trees (HF003):
  - in spring, the share of buds broken times the leaves' size;
  - in autumn, the share of leaves not yet fallen.
- **What a sensor reads.** A sensor at height z reads the downwelling PAR at the top of the tallest cohort
  layer at or below z, or the ground's when none is. The value is averaged over the patches by area: the
  horizontal mean of a 1-D layered canopy.
  - `edmore.canopy.Canopy.radiation` returns this light at every cohort's top as `cohort_incid_vis`.

## Results (2019)

![The modelled PAR at NEON's sensors against what they measured](canopy_radiation.png)

PAR at a sensor over the PAR above the canopy, with the sun above about 20°, observed / modelled:

| | April (leaf-off) | July–August | October |
|---|---|---|---|
| 29 m | 1.02 / 1.00 | 1.01 / 1.00 | 1.05 / 1.00 |
| 24 m | 0.95 / 0.78 | 0.96 / 0.23 | 0.90 / 0.25 |
| 16 m | 0.71 / 0.60 | 0.16 / 0.06 | 0.19 / 0.10 |
| 5 m | 0.47 / 0.50 | 0.034 / 0.024 | 0.09 / 0.06 |
| 0.4 m | 0.39 / 0.50 | 0.030 / 0.023 | 0.08 / 0.06 |
| floor, 3 sensors | 0.27–0.33 / 0.22–0.47 | 0.010–0.025 / 0.009–0.019 | 0.05–0.06 / 0.04–0.05 |

What the prototype shows:
- **The understory and the floor follow the seasons:** leaf-out in late May to June, leaf fall in October.
  In summer they are within a factor of ~1.5, darker than measured. The allometry's leaf area (LAI
  6.3–7.7 around the sensors) is above the LAI-2000 plant area at the EMS tower (4.1–4.6). The next step
  is the leaf area from litterfall × SLA (HF069).
- **The patch size matters as much as the optics.** With one patch per 10 m cell, the summer floor is 2–4×
  brighter than measured: crowns spread over more than 10 m, so cells with few stems become gaps.
  With one uniform patch it is slightly darker than measured. 20 m, the census's quadrat, matches best.
- **Inside the crowns (16 and 24 m) the model is far too dark.** A cohort's leaves sit in a thin layer at
  the tree's top. A sensor at 16 m therefore has all the leaves of a 17–28 m tree above it, where the
  real crown reaches below the sensor. The 24 m sensor also sits by the tower's own gap.
- **The model ranks the skies the wrong way round.** Measured, overcast light reaches 16 m and the floor
  as well as or better than clear light: in October 2.3× better at 16 m and 2× at the floor. The model
  passes the beam better than isotropic diffuse light.
  - The candidates are wood with the leaves' angle distribution, an isotropic sky, and no crown-scale
    gaps. These are the open questions of sunlit and shaded leaves and of crown-scale darkening.
- **MODIS LAI** as the leaf-out driver (`--phenology modis`) puts leaves out in April and takes them
  off in mid-October, while the oaks keep them to November. The ground observations fix both.

## Files

| File | What it is |
|---|---|
| [`run_example.py`](run_example.py) | the run: drivers, stands, the two-stream, the comparison table |
| [`harvard_data.py`](harvard_data.py) | the readers: NEON sensors, the EMS split, the phenology, the census stand |
| [`plot_example.py`](plot_example.py) | the figure |
