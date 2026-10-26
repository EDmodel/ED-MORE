# MEDS documentation

**Coming from ED2?** Start with **[`ed2_comparison.md`](ed2_comparison.md)**: what carried over,
what changed and what is missing, keyed to the `ED2IN` options you already know.

**Setting up a run?** Start with [`configuration.md`](configuration.md): the files a run needs, how it
starts, and what each block of the config controls.

**Want to change the code?** Start with [`../src/README.md`](../src/README.md): the source layout,
the placement rules and the library graph.

## What is here

- **[`science/`](science/)**: the model's **equations and their rationale**, one page per process,
  versioned with the code that implements them.
- **[`configuration.md`](configuration.md)**: the two TOML files, the parameter philosophy, the three
  ways a run starts, output selection, regional runs and several sites per polygon.
- **[`building.md`](building.md)**: compilers, dependencies, the four build kinds, the tests and the
  parallel builds.
- **[`ed2_comparison.md`](ed2_comparison.md)**: MEDS against ED2, process by process.
- **[`ROADMAP.md`](ROADMAP.md)**: everything MEDS has deliberately deferred, each item with its issue.
- **[`dev_plans/`](dev_plans/)**: dated **design and implementation records**, working documents rather
  than reference; [`dev_plans/README.md`](dev_plans/README.md) says which are still live, and
  [`dev_plans/archive/`](dev_plans/archive/) holds the ones whose work is done.

The changelog is at the repository root: [`../CHANGELOG.md`](../CHANGELOG.md).

## The science pages

**The whole model**

| Page | Covers |
|---|---|
| [`order_of_processes.md`](science/order_of_processes.md) | Which process runs when: a run, one slow step, the fast loop's sub-steps, the slow dynamics, the calendar boundary's restructuring, and output. |

**Within the fast loop (every `dt_fast`, 15 minutes in the shipped configs)**

| Page | Covers |
|---|---|
| [`coupled_biophysics.md`](science/coupled_biophysics.md) | How a patch's stores (canopy air, leaves, wood, ground, snow, soil) are coupled and advanced together, and what each conserves. |
| [`numerical_scheme.md`](science/numerical_scheme.md) | How the fast loop is advanced in time: what is held fixed over a step, the two integrators, how they differ, and how to choose them and `dt_fast`. |
| [`canopy_radiation_transfer.md`](science/canopy_radiation_transfer.md) | The ED2 two-stream solver, leaf-angle distributions, and the multi-band solution. |
| [`canopy_aerodynamics.md`](science/canopy_aerodynamics.md) | The Monin–Obukhov surface layer, the leaf and wood boundary layers, and the wind inside the canopy. |
| [`canopy_air_space_biophysics.md`](science/canopy_air_space_biophysics.md) | The canopy air's three prognostic variables: enthalpy, humidity and CO₂. |
| [`leaf_gas_exchange.md`](science/leaf_gas_exchange.md) | Photosynthetic demand, three stomatal-conductance models, water stress in two limbs, and the coupled Cᵢ solver. |
| [`plant_hydraulics.md`](science/plant_hydraulics.md) | Plant water: pressure–volume curves, xylem vulnerability and the Kirchhoff conductance, root uptake by layer. |
| [`soil_plant_water.md`](science/soil_plant_water.md) | The pond, the soil layers and every cohort's wood and leaves solved as one water network. |
| [`vegetation_energy_dynamics.md`](science/vegetation_energy_dynamics.md) | Leaf and wood temperature as an exact relaxation toward each tissue's energy balance. |
| [`soil_biophysics.md`](science/soil_biophysics.md) | Soil water (mixed-form Richards) and soil heat as conserved quantities, frozen soil, and the soil parameters from texture and maps. |
| [`snow_biophysics.md`](science/snow_biophysics.md) | The snow and temporary surface-water store. |
| [`forcing.md`](science/forcing.md) | Meteorological forcing: the file format, interpolation, shortwave disaggregation and partition, the corrections to the site and the canopy-air top, calendar recycling, and CO₂. |

**Once per slow step (`dt_slow`, a day by default)**

| Page | Covers |
|---|---|
| [`plant_phenology.md`](science/plant_phenology.md) | Leaf phenology: temperature, hours-of-light and water cues → flush and senescence tendencies → leaf growth and loss. |
| [`plant_traits.md`](science/plant_traits.md) | Light-driven plasticity of leaf traits and its carbon consequences. |
| [`plant_carbon_allocation.md`](science/plant_carbon_allocation.md) | The PARTEH-H1 carbon budget: net carbon → priority ladder → growth of each pool. |
| [`plant_respiration.md`](science/plant_respiration.md) | The autotrophic budget: leaf, stem and fine-root maintenance, and growth respiration on realized growth. |
| [`soil_carbon.md`](science/soil_carbon.md) | The CENTURY decomposition network as a carbon matrix ODE, and the seam between fast and slow respiration. |
| [`demography.md`](science/demography.md) | Growth, mortality, recruitment and treefall disturbance as rates, and the operators that apply them: cohort and patch fusion, fission, the cull and termination. |
| [`conservation.md`](science/conservation.md) | The slow-step ledger: the site's carbon, water and energy across each slow-step operator, what each declares, and plausibility. |

**Output**

| Page | Covers |
|---|---|
| [`diagnostics.md`](science/diagnostics.md) | What is written, on which axes and tiers, how it is aggregated across the hierarchy and across time, and how to choose variables. |

For the API of a routine, read its doc comment in the Fortran source.

## Writing a science page

A science page states equations and their rationale in the present tense. What changed and when
goes in [`../CHANGELOG.md`](../CHANGELOG.md); what is deferred goes in [`ROADMAP.md`](ROADMAP.md).

- **Shape.** Open with what the page covers and what the reader gets from it, then the numbered
  sections, then `## Parameters`, `## Code map` and `## References`, in that order, as the page needs
  them.
- **Equation numbers are cited from the code.** Every display equation ends in `\qquad(N)`, numbered
  in reading order. Source comments cite `docs/science/<page>.md eq. (N)` and `§N`, so inserting or
  removing an equation, or renumbering a section, means updating every citation in the same change.
- **Math.** GitHub's Markdown sanitizer runs *inside* `$…$` before MathJax, so display equations go in
  ` ```math ` fenced blocks, never `$$…$$`, and inline math with backslash spacing or paired brace
  subscripts needs the dollar-backtick form. The full rules are in
  [`.claude/rules/docs-math.md`](../.claude/rules/docs-math.md).

There is **no local preview**. Validate structurally (no `$$` left, fences balanced, every
conflicting construct protected) and look at the rendered file on a branch before merging.
