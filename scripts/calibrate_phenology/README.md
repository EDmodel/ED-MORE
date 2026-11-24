# calibrate_phenology — a stand's leaf phenology from observations, with each cohort's own cues

Fits the phenology traits of a stand's PFTs ([`docs/science/plant_phenology.md`](../../docs/science/plant_phenology.md))
to observations of its canopy. These include satellite LAI, tagged trees, litter baskets and bud-break dates.
Every cohort keeps the cues it reads in the coupled model:
- the temperature of the canopy air around it;
- the hours of light at its top.

A cohort in a closed patch sits in cooler canopy air on a sunny day than one in a gap, and an
understorey cohort is darker than the trees above it. So each keeps its own calendar, and the fit sees it do so.

## The method

The cues depend on the canopy, and the canopy depends on the phenology. A fit whose every trial was a coupled run
would cost a full simulation each time. Instead, the tool alternates between two steps.

1. **A coupled run** of the site with the current values. The stand is frozen
   (`[demography].demography_on = false`), so every cohort keeps its size and its place; the leaves still follow
   the season. Each cohort's daily cues come from its output (`mean_pheno_temp_co`, `mean_pheno_light_hours_co`,
   `mean_dmax_psi_leaf_co`), with its leaf turnover rate and its leaf area.
2. **The closure check.** The offline stand runs every cohort as a column through the model's own kernel and leaf rule
   (`edmore.plant.pheno.run`). Run with the coupled run's own values, it must give that run's leaf cover to within
   `[fit].closure`, every day and PFT, or the tool stops.
3. **The offline fit.** Differential evolution searches the fitted keys' ranges. Each trial runs the offline stand with its values on the coupled run's
   cues and scores it on the targets. A trial costs a fraction of a second where a coupled run costs minutes.
4. **Again**: a coupled run with the fitted values gives new cues, and the next fit starts from those values. The
   iteration stops when no key moves more than `[fit].settle` of its range between two fits (the canopy the cues
   came from no longer changes), or after `[fit].iterations` fits. A cue the canopy feeds back on can keep the
   fits from settling: a grass's leaves dry the soil whose potential is its WATER cue, so each fit, made on
   the last run's cues, overshoots on its own. The tool then keeps the values whose own coupled run scored
   best (each run's score of the values it was made with, `coupled_score` in `iterations.json`), not the
   last fit's, which no run has tried.

The tool writes the site's PFT file with the fitted values, every other line and comment kept:
`<work>/pft_parameters_fitted.toml`. Each fit writes `fit.json` (the values, the scores before and after, the
closure) and `targets.csv` (each observation beside the model) in its run's directory.

## Run

```
calibrate_phenology.py iterate --config phenology.toml --work runs --workers 40   # the method above
calibrate_phenology.py run     --config phenology.toml --work runs                # one coupled run and its closure
calibrate_phenology.py fit     --config phenology.toml --work runs --workers 40   # one fit on that run
```

`--values fit.json` makes `run` and `fit` start from a fit's values instead of the PFT file's. A coupled run is
reused when its directory holds one made with the same values. `--workers` counts processes; each loads the
cues once (about 30 MB for 1,000 cohorts and five years). On Slurm, ask for the whole node (`--mem=0`).

## The calibration file

Every setting is listed in [`phenology_reference.toml`](phenology_reference.toml) with its default and the reason
for it; that file is also the schema, so a setting it does not list stops the tool. `[run].config` names the
site's main file and `[run].score_from` the first day scored. The year before it spins up the cue memory and the
canopy from the census, which starts in full leaf. A site whose run starts from a spin-up's state file
(`[init].init_mode = 2`, as Toolik's) needs no such year: each cohort starts from the file's phenology memory, and
the offline stand starts its column from the same memory and the run's first-day leaf cover.

- **The fitted keys** are named `key_pftN` (`flush_degree_days_pft1`). The PFT file's cue masks say which cues a
  PFT has; the fit sets their values. `leaf_turnover_rate` is the PFT file's `[pft].leaf_turnover_rate_toc`. A
  cohort's own rate is that times its light gradient, so offline each cohort's coupled rate is scaled by the
  trial's top-of-canopy rate.
- **`par_min` is not fitted.** The coupled run counts each cohort's hours of light at the PFT file's `par_min`, and the
  hour thresholds (`flush_light_hours`, `shed_light_hours`) absorb a different choice. At a low `par_min` (a few
  µmol m⁻² s⁻¹) the hours are the photoperiod the cohort senses at its height.
- **The WATER cue** sums each cohort's predawn leaf water potential above and below the PFT's turgor-loss point,
  which the offline stand derives from the leaf's pressure-volume curve as the model does
  (`[hydraulics].leaf_pi0` and `leaf_elastic_mod`). Its sums and sharpnesses are fitted; the curve is not.

## Targets

A `[[target]]` compares a series of observations with the same quantity of the offline stand. It sees the
cohorts of its `pfts` from `min_height` to `max_height`, each counted by its full-canopy leaf area. Its score is
the mean squared error in units of `sigma`, times `weight`; the fit minimises the sum.

| kind | observations (file columns) | the model's side |
|---|---|---|
| `fraction` | how full the canopy is, 0 to 1 (`date, value`) | leaf / full leaf, or (`scale = "year_minmax"`) the leaf area scaled between each year's least and most, as MODIS LAI is; averaged over `window_days` |
| `fallen` | the share of the leaves fallen, as on tagged trees (`date, value`) | the share of the year's fullest canopy gone, from `from_doy` |
| `litter_fraction` | the share of each leaf year's litter caught by each date (`date, value`) | the cumulative leaf loss from `year_start` |
| `annual_loss` | canopies lost a year (`value`) | the leaf loss of each scored year over the full canopy |
| `date` | the day of the year of an event, as 50 % bud break (`year, doy`) | the day the canopy, scaled between the year's least and most, crosses `level` |

The height range ties a target to the trees it was observed on. Tagged trees and litter baskets see the canopy
trees; a sensor below the canopy (LAI-2000) or a satellite sees every layer. An evergreen's season shows in its
needle litter and in its spring flush rising from its floor (`min_leaf_cover`). Its annual loss is the inverse of
its needles' lifespan.

## The closure check

Offline, carbon never limits the flush. In the coupled run a cohort short of storage refills its leaves more
slowly. A closure that fails means the offline stand is not the coupled run. The cause is a carbon-limited flush,
a cohort that changed (a stand that was not frozen), or cues misaligned in time. The tool reports each PFT's
largest daily difference and its day, and stops above `[fit].closure`.

## Worked example: Harvard Forest

[Example 02, part 2](../../examples/example02_canopy_phenology/README.md#part-2--fitting-a-stand-harvard-forest)
([`harvard_stand/phenology.toml`](../../examples/example02_canopy_phenology/harvard_stand/phenology.toml))
fits a cold-deciduous broadleaf and an evergreen conifer on example 04's Harvard stand, 2019–2024, to
HF003's tagged trees, HF069's LAI-2000 and litter baskets, MODIS, and the conifers' climatologies.
Four lessons from it hold at any site:

- **Scaled targets need an anchor.** A target scaled within each year (MODIS, a crossing day) says nothing
  about how full the canopy gets. With only such targets a fit can lower the whole cycle, leaving a canopy
  that never fills while every target matches. Add a `fraction` target with `scale = "none"` where the
  canopy is known to be full: the census gives each cohort its full canopy.
- **Read what an observation is.** HF003's `l75` is a code that tops out at 80, not a share of the
  canopy. Use the day it reaches half its plateau (a `date` target), not `l75 / 100`.
- **An evergreen's winter dip needs a constraint.** Its annual loss can come from an autumn dip and a
  spring refill, or from an exchange: falling while the canopy re-flushes. Timing targets cannot tell the
  two apart. Set `min_leaf_cover` from the needles' lifespan (two thirds stay when they live three years),
  and give the flush a sharp day-length gate.
- **The light cue reads the light where the cohort is.** Under a closed canopy a cohort counts fewer
  hours above `par_min`. If its flush comes after the canopy closes (a conifer in June), a `par_min` set
  for open-grown trees can keep it from ever refilling. A photoperiod is sensed at very low light.

## Tests

`tests/` (CTest `calibrate_phenology`) checks the settings, the fitted PFT file, each target kind on series made
by hand, and the daily output read into cues. An offline fit on synthetic cohorts must find the values its
observations were made with, using two worker processes.
