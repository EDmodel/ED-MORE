# Example 07 — Canopy radiation against vertical PAR sensors (Harvard Forest) — PROTOTYPE

ED-MORE's two-stream canopy radiation ([`docs/science/canopy_radiation_transfer.md`](../../docs/science/canopy_radiation_transfer.md)),
driven on its own from Python, against the light measured inside a forest. The site is NEON's tower at
Harvard Forest (Massachusetts), with up-facing PAR sensors at 0.36, 5.1, 16.1, 23.6, 29.4 and 39.1 m and three
line sensors on the forest floor, all inside the ForestGEO 35-ha plot.

Two windows of 2019 isolate the canopy's two states:
- **1–25 April, leaf-off:** after the snow (gone on 30 March), before bud break. Only the conifers carry
  leaves, and the bare wood intercepts the rest.
- **July, leaf-on:** full leaves.

March would be leafless too, but in 2019 snow buried the floor sensors and made the ground far brighter than
the model's bare soil. This is a prototype: the pipeline runs end to end, but the fetch script and the README
are not final.

## Run

```bash
cd examples/example07_canopy_radiation
python run_example.py --data <data>      # about 10 s; writes output/
python plot_example.py                   # output/canopy_radiation.png
python fit_allometry.py --data <data>    # how the Harvard allometry below was derived (optional)
```

It needs a built `libedmore.so` (`build-dev/`, or `EDMORE_LIB`), numpy, pandas and matplotlib. `<data>`
holds what the fetch scripts download:
- **NEON HARV:** a NEON API token is needed since 2026; the data are CC BY 4.0.
  - tower PAR (DP1.00024.001) and floor line-quantum PAR (DP1.00066.001), 30 min;
  - vegetation structure (DP1.10098.001), for tree heights;
  - foliar traits (DP1.10026.001), for leaf mass per area.
- **Harvard Forest data archive (CC0):**
  - the EMS tower's direct and diffuse PAR (HF004), 134 m from the NEON tower;
  - the 2019 ForestGEO census (HF253);
  - litterfall and LAI-2000 (HF069).
- **MODIS LAI** (MCD15A3H, ORNL DAAC subset).

## The model

- **Drivers.** The PAR measured above the canopy (39 m), split into beam and diffuse by the EMS tower's
  direct/diffuse sensor, with the model's own sun. NEON's own direct/diffuse sensor was shaded by the
  lightning rod until 2025-06.
- **Stand.** Around each sensor, every live stem of the 2019 census within 40 m becomes a cohort in its 20 m
  quadrat's patch. Two PFTs: deciduous and evergreen (the conifers). In April the deciduous trees have no
  leaves; in July all their leaves.
- **Harvard allometry** (`fit_allometry.py`; set on the repository's PFT file in `run_example.py`):
  - **height:** ED-MORE's log-quadratic curve per PFT, fitted to 3,400 trees NEON measured at HARV, each
    DBH class weighted equally, capped at 28 m. The fit is within ~1 m of the measured mean in every DBH
    class. The tallest tree per patch is 25 m, against 24–27 m in the lidar.
    - The default pan-tropical curve put 60–80 cm oaks at 32–37 m and the conifers 3–4 m too tall.
  - **deciduous leaf area:** scaled (`lai_b1` 0.234 → 0.171) to a July LAI of 3.95 around the towers, the
    median of three sources:
    - leaf litterfall over leaf mass per area: 4.5 (2.95 at sun-leaf LMA);
    - MODIS July minus leaf-off: 3.8;
    - LAI-2000 over the default clumping: 3.95.
  - **evergreen and wood area:** the defaults (conifer LAI 0.55–1.65 around the sensors, WAI 0.57–0.69).
- **Optics: two settings.** Every key is the default except the clumping:
  - **foliage at random** within each cohort layer (leaf and wood clumping 1.0), since the 20 m patches
    already carry the gaps between crowns;
  - **the default clumping** (leaf 0.8, wood 0.5).
- **What a sensor reads.** A sensor at height z reads the downwelling PAR at the top of the tallest cohort
  layer at or below z, or the ground's when none is. The value is averaged over the patches by area: the
  horizontal mean of a 1-D layered canopy. `edmore.canopy.Canopy.radiation` returns this light at every
  cohort's top as `cohort_incid_vis`.

## Results (2019)

![The modelled PAR at NEON's sensors against what they measured](canopy_radiation.png)

PAR at a sensor over the PAR above the canopy, with the sun above ~20°: observed / foliage at random /
default clumping:

| | April, leaf-off | July, leaf-on |
|---|---|---|
| 29 m | 1.02 / 1.00 / 1.00 | 1.01 / 1.00 / 1.00 |
| 24 m | 0.96 / 0.83 / 0.88 | 0.96 / 0.60 / 0.66 |
| 16 m | 0.72 / 0.55 / 0.69 | 0.16 / 0.07 / 0.12 |
| 5 m | 0.47 / 0.42 / 0.57 | 0.033 / 0.030 / 0.063 |
| 0.4 m | 0.39 / 0.42 / 0.57 | 0.025 / 0.027 / 0.058 |
| floor (3 sensors) | 0.26–0.32 / 0.20–0.40 / 0.32–0.56 | 0.009–0.017 / 0.010–0.012 / 0.026–0.033 |

- **Below the crowns the two-stream reproduces the light in both seasons,** once the foliage is placed at
  random within each layer.
  - The default clumping lets twice the light through, because the patches already hold the crown-scale gaps.
  - The model's effective plant area with random foliage (~5.1 at the tower) matches LAI-2000's 4.8.
- **The floor follows the sun under clear skies** in July (panel e), with the model's slope.
- **Leaf-off, low sun (panel d):** the observed floor light falls ~3.5× from a high to a low sun; the
  model's falls 1.7×. In the model the wood takes the leaves' angle distribution, but trunks and branches
  stand nearer vertical and block a low sun far more.
- **Overcast skies (panel c):** the measured light below the crowns is as high as under a clear sky (July
  floor: 0.014 and 0.014). The model's is half: the two-stream's diffuse light comes from an isotropic sky,
  while an overcast sky is brightest overhead, and zenith light penetrates furthest.
- **Inside the crowns (16 and 24 m) the model is too dark.** This is an assumption of ED-MORE: each
  cohort's leaves and wood sit in a thin layer at the tree's top, with no crown depth, so a sensor at 16 m
  has all of a 17–25 m tree's foliage above it. The 24 m sensor also sits beside the tower's own gap.
- **A point sensor and a horizontal mean** (panel f): one sensor swings between shade and sunfleck from
  half hour to half hour, so the comparison is made over a window's sums, not hour by hour.

Not in this example: the near-infrared and thermal bands (the same solver), which a full-band comparison
may add later.

## Files

| File | What it is |
|---|---|
| [`run_example.py`](run_example.py) | the run: drivers, stands, the two-stream in both optics, the comparison table |
| [`harvard_data.py`](harvard_data.py) | the readers: NEON sensors, the EMS split, the census stand with ED-MORE's allometry |
| [`fit_allometry.py`](fit_allometry.py) | the Harvard height and leaf-area allometry, from NEON, HF069 and MODIS |
| [`plot_example.py`](plot_example.py) | the figure |
