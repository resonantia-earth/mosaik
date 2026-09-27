# mosaik

<!-- badges: start -->
<!-- badges: end -->

**mosaik** is categorical raster algebra for landscape analysis.

## Core idea

Landscape-metric software has largely been a matter of fixed metric lists: you
get the metrics the author implemented, and a metric that is not on the list is
a feature request. mosaik takes the other route. It exposes the **measurement
primitives and the means of composing them**, so that a new metric is an
expression rather than a new function.

Nearly every published landscape metric decomposes into five primitives —
adjacency, area, perimeter, distance, dissimilarity — measured at patch, class,
or landscape scale. `msr()` combines them from an equation written in
`metric.scale` notation.

Four things distinguish it:

- **Categorical-first.** Cell values *are* group IDs, with `@categories` and
  `@patches` as first-class slots. terra is numeric-first and treats categories
  as levels bolted on.
- **Native provenance.** `@provenance` is a slot that survives every operation,
  recorded in the vocabulary of the W3C PROV ontology.
- **Composability.** Mosaic in, mosaic out, so everything chains with the pipe.
  A recorded sequence of operations is a recipe, replayable on any raster.
- **Five primitives**, out of which the rest is composed.

The `mspa` vignette is the demonstration: a published segmentation algorithm
reproduced at full fidelity out of the primitives alone.

## Architecture

| Prefix | Purpose | Example |
|--------|---------|---------|
| `syn_*` | Synthesise an abstract field | `syn_noise()`, `syn_texture()`, `syn_gradient()` |
| `mdf_*` | Modify layers | `mdf_erode()`, `mdf_blend()`, `mdf_mask()` |
| `msr_*` | Measure primitives | `msr_area()`, `msr_perimeter()`, `msr_adjacency()` |
| `msr()` | Compose a metric from an equation | `msr(m, equation = "perimeter.class / area.class", label = "edge_density")` |
| `msk_*` | Utilities and accessors | `msk_vis()`, `msk_terra()`, `msk_extent()` |

For **generating whole synthetic landscapes** — climate, soil, vegetation, land
use, tenure — see the companion package
[mundus](https://github.com/resonantia-earth/mundus), which builds on this class.

## Installation

``` r
# install.packages("remotes")
remotes::install_github("resonantia-earth/mosaik")
```

## Quick start

``` r
library(mosaik)

# a synthetic field to analyse
m <- mosaik(extent = c(0, 1000, 0, 1000), res = 10) |>
  syn_noise(type = "perlin", name = "field") |>
  mdf_binarise(thresh = 0.5, layer = "field") |>
  mdf_componentise(layer = "field")

# measure the primitives
m <- m |>
  msr_area(scale = "class") |>
  msr_perimeter(scale = "class")

# compose a metric from them
m <- msr(m, equation = "perimeter.class / area.class", label = "edge_density")

msk_categories(m)
msk_vis(m)
```

## Interactive metrics explorer

An interactive d3-based dependency graph with a built-in decision guide is
hosted at
[mosaik.resonantia.earth/metrics](https://mosaik.resonantia.earth/metrics/). It
shows how the ~50 supported landscape metrics decompose into the measurement
primitives, with a decision tree to find the right metric for your analysis and
a "Show code" view giving the exact `msr()` call for each.

## License

GPL (>= 3)
