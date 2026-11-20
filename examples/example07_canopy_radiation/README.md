# Example 07 — Canopy radiation against vertical PAR sensors (Harvard Forest)

The canopy radiation of ED-MORE on its own: the multi-layer two-stream that every patch solves each step,
driven from Python through [`edmore.canopy`](../../python/edmore/canopy/__init__.py) and compared with
the light measured inside a forest. The equations are in
[`docs/science/canopy_radiation_transfer.md`](../../docs/science/canopy_radiation_transfer.md).

The site is NEON's tower at Harvard Forest (Massachusetts), inside the ForestGEO 35-ha plot. It has
up-facing PAR sensors at 0.4, 5, 16, 24, 29 and 39 m, and three line sensors on the forest floor. Two
windows of 2019 isolate the canopy's two states:
- **leaf-off, 1–25 April:** after the snow, before bud break; only the conifers carry leaves;
- **leaf-on, July.**

The same windows of 2021 test the result on a year the calibration did not see.

## Run

With a built `libedmore.so` (`build-dev/`, as in [example 01](../example01_leaf_gas_exchange/README.md)):

```bash
cd examples/example07_canopy_radiation
python fetch_data.py <data>                       # once, ~60 MB
python run_example.py --data <data>               # ~10 s: output/model_vs_sensors.csv and the table
python run_example.py --data <data> --year 2021   # the held-out year, in output/2021/
python plot_example.py                            # output/canopy_radiation.png
python calibrate_optics.py --data <data>          # optional, ~9 min: rewrites calibrated_optics.json
python fit_allometry.py --data <data>             # optional: derives HARVARD_PFT in run_example.py
```

`fetch_data.py` needs a NEON API token (free with a NEON account) in `$NEON_API_TOKEN`, or in a file
given by `--token-file`. It downloads:
- **NEON HARV** (CC BY 4.0): PAR at the tower levels and on the floor (DP1.00024.001, DP1.00066.001),
  vegetation structure (DP1.10098.001) and foliar traits (DP1.10026.001);
- **the Harvard Forest data archive** (CC0): the EMS tower's direct and diffuse PAR (HF004), the
  ForestGEO census (HF253), and litterfall and LAI-2000 (HF069);
- **MODIS LAI** (MCD15A3H), from the ORNL DAAC's subsets.

## Setup

| Input | Value |
|---|---|
| Light above the canopy | NEON's 39 m PAR, split into beam and diffuse by the EMS tower's direct/diffuse sensor 134 m away; the model's own sun |
| Stand | every live tree of the 2019 census within 40 m of a sensor, one cohort per tree (its stems lumped at equal basal area), one patch per 20 m quadrat; deciduous and evergreen (conifer) PFTs |
| Height | ED-MORE's log-quadratic curve per PFT, fitted to 3,400 trees NEON measured at the site, capped at 28 m (the lidar canopy top) |
| Leaf area | deciduous `lai_b1` scaled to a July LAI of 3.95, the median of litterfall × leaf mass per area (4.5), MODIS (3.8) and LAI-2000 (3.95) |
| Wood angle | the default 75° ± 20°: the woody surface's inclination in unpublished TLS data |
| Optics | leaf clumping 0.93, mean leaf angle 34°, wood clumping 1.0, calibrated to the sensors below the crowns in 2019 ([`calibrated_optics.json`](calibrated_optics.json)); every other key the default |
| A sensor's light | the downwelling PAR at the top of the tallest cohort layer at or below the sensor (`cohort_incid_vis`), the patches' area mean |

The 20 m patches already hold the gaps between crowns, so the clumping within a layer comes out near 1.
The sensors at 16 and 24 m are left out of the calibration: they sit in the crowns, where the layers'
geometry decides (Results).

## Results

![ED-MORE's two-stream against NEON's PAR sensors at Harvard Forest, 2019](canopy_radiation.png)

PAR over the PAR above the canopy, sun above ~20°: observed / calibrated / default optics, 2019.

| | April, leaf-off | July, leaf-on |
|---|---|---|
| 24 m | 0.96 / 0.84 / 0.88 | 0.96 / 0.60 / 0.66 |
| 16 m | 0.72 / 0.57 / 0.70 | 0.16 / 0.07 / 0.12 |
| 5 m | 0.47 / 0.45 / 0.59 | 0.033 / 0.030 / 0.063 |
| 0.4 m | 0.39 / 0.44 / 0.58 | 0.025 / 0.027 / 0.058 |
| floor | 0.26–0.32 / 0.22–0.43 / 0.33–0.57 | 0.009–0.017 / 0.009–0.012 / 0.026–0.033 (two sensors) |

- **(a, b) Below the crowns** the calibrated model holds the light in both seasons; the default optics let
  through about twice as much in July. **In 2021**, which the calibration did not see, the fit's cost
  (mean squared log of model over observed) is 0.20, against 0.16 in 2019 and 0.33 with the default
  optics.
- **In the crowns (16, 24 m) the model is too dark.** ED-MORE puts each cohort's leaves and wood in a thin
  layer at the tree's top, with no crown depth, so a sensor at 16 m has all of a 17–25 m tree's foliage
  above it. The 24 m sensor also sits by the tower's own gap.
- **(d, e) Under a clear sky** the light falls with the sun's height, more in the data than in the model.
  Leaf-off, from cos z 0.85 to 0.35, the light below the crowns falls to ×0.24–0.56 of its value
  observed, and to ×0.48–0.59 modelled.
  - The near-vertical wood accounts for part of the drop: with the wood at the leaves' angle, the model
    falls only to ×0.66–0.77.
  - The rest is the layers' horizontal uniformity: a low sun crosses whole crowns and trunks.
- **(c) Under overcast skies** the light below the crowns is as high as under a clear sky, and the
  model's is lower: in July, 0.020 against 0.030 at 5 m. The two-stream's diffuse light comes from an
  isotropic sky, while an overcast sky is brightest overhead.
- **(f) A sensor is a point; the model a horizontal mean.** On clear days the light at 5 m swings between
  0.1 and 0.9 as the sun's path crosses gaps and trunks above the sensor. The pattern repeats from day to
  day (r = 0.97); under overcast skies the sensor varies by 4 %. So the comparison is made over each
  window's sums, not half hour by half hour. The patches' 10–90 % range holds 89–94 % of the half hours
  at 5 and 0.4 m, and 55–81 % on the floor.

## Files

- [`fetch_data.py`](fetch_data.py) — the downloads, in the layout the readers expect.
- [`run_example.py`](run_example.py) — the drivers, the stands, the two-stream in both optics, the table.
- [`harvard_data.py`](harvard_data.py) — the readers: NEON sensors, the EMS split, the census stand.
- [`calibrate_optics.py`](calibrate_optics.py) and [`calibrated_optics.json`](calibrated_optics.json) — the optics' fit and its result.
- [`fit_allometry.py`](fit_allometry.py) — the Harvard height and leaf-area allometry.
- [`plot_example.py`](plot_example.py) — the figure.
