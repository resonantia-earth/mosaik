# mosaik: reference for language models

This file lets a language model write correct mosaik code without reading the
source. It lists every exported function with its real signature, where results
are stored, and the mistakes that recur. If it disagrees with the help pages,
the help pages are right.

Read "Common mistakes" first.

## Common mistakes

1. **There is no `derive()`.** A metric built from primitives is computed with
   `msr(obj = NULL, equation, label, layer)`, where `equation` is a character string in
   `metric.scale` notation:
   `msr(m, "perimeter.class / area.class", "edge_density", layer = "cover")`.
2. **Measure the primitives before `msr()`.** The equation reads results that
   `msr_*` already stored. `area.landscape` exists only after
   `msr_area(scale = "landscape")`.
3. **The name in an equation is the stored column name.** `msr_adjacency(type =
   "like")` stores `likeAdj`, so the equation uses `likeAdj.class`, not
   `adjacency.class`. See the storage table below.
4. **`msr_number()` names its scale by what is counted.** `scale = "patch"`
   counts patches per class and stores `number` in the categories (use
   `number.class`). `scale = "class"` counts classes and stores `number` in the
   global slot (use `number.landscape`).
5. **`msk_vis()` takes layer specifications, not a `layer` argument.** Write
   `msk_vis(m, .layer("cover"))`, not `msk_vis(m, layer = "cover")`.
6. **`msr_*` results are not returned as tables.** Every function returns the
   mosaik with results attached. Read them with `msk_categories()`,
   `msk_patches()` or `msk_global()`.
7. **A recipe is built by calling functions without an object.** Called
   without one, an `mdf_*` or `msr_*` function (and `msr()`) records itself
   instead of running. Chain such calls, then apply the recipe with `mdf()`.
   Arguments are stored by value when recorded, except the expressions of
   `mdf_filter(expr)` and `mdf_loop(until)`, which refer to layers of the map
   the recipe is later applied to.
8. **Functions that do not exist:** `derive`, `as_terra`, and anything with the
   prefixes `mk_`, `gnrt_`, `mg_`, `sim_`, `make_`. Terrain, climate, soil,
   vegetation, land-use and simulation functions belong to the separate package
   mundus, not to mosaik. `msr_distance` does not exist either: distance is
   `msr_cost()` without a cost surface.
9. **C++ functions are internal.** Do not call `morphCpp`, `distanceCpp` and so
   on; use the R functions.
10. **Number the patches before measuring them.** Every patch-level measure
    (`msr_area/msr_perimeter/msr_cost/msr_adjacency(scale = "patch")`,
    `msr_number(scale = "patch")`, and `msr()` with a `.patch` or `.cell`
    variable) measures the patches `mdf_componentise(layer = ...)` numbered
    on that layer, and stops if there are none. The measures never find
    patches themselves, because the connectivity (4 or 8) must be stated.

## The object

One S4 class, `mosaik`, with these slots:

| slot | holds |
|---|---|
| `extent` | `c(xmin, xmax, ymin, ymax)` |
| `dims` | `c(ncols, nrows)` |
| `layers` | named list of flat row-major vectors, compressed with `rle()` when smaller |
| `categories` | per layer: `gid` (class IDs) and `val` (labels) for a categorical layer, plus class-level results |
| `patches` | per layer, written by `mdf_componentise`: `ids` (the layer holding the patch numbers), `connectivity`, `class` and `patch` per patch; then one field per patch-level metric |
| `global` | landscape-level results |
| `crs` | a CRS string, or `NA` |
| `provenance` | one entry per operation, in order |

Cell values of a categorical layer are the class IDs themselves. Read layer
values with `msk_pull()`, never from `obj@layers` directly.

Every function takes a mosaik and returns a mosaik, so everything chains with
`|>`. Most functions take `layer` (default: the first layer) and, for `mdf_*`,
`add` (`NULL` overwrites `layer`, a string writes a new layer of that name).

## Create, read, write

- `mosaik(extent = NULL, res = NULL, crs = NA_character_, vals = NULL, rast = NULL, group = FALSE)`:
  an empty grid from `extent` and `res`; with `vals`, a matrix or named list of
  layers; or from a terra `SpatRaster` via `rast`.
- `msk_terra(obj)`: to a `SpatRaster`.
- `msk_pull(obj, layer = NULL)`: a layer's values as a vector.
- `msk_spaghettify(x)`: sf or terra vectors as a table of vertices `x, y, id,
  part`, the input of `msk_rasterise()`.
- `msk_rasterise(obj, geom, type = "polygon", name = "values")`: write points,
  lines or polygons, given as a data frame `x, y, id` (and optionally `part`),
  into a new layer as their `id`; `NA` elsewhere. Rings combine by even-odd, so
  holes are just further parts. Works on an empty grid.
- `msk_add(obj, from = NULL, ..., rename = NULL)`: copy layers from another
  mosaik on the same grid (extent, dims and CRS are checked), or add vectors of
  cell values as `name = values`; category tables come along.
- `msk_select(obj, ...)`, `msk_remove(obj, ...)`: keep or drop layers, named
  unquoted.
- Accessors: `msk_extent`, `msk_dims`, `msk_res`, `msk_ncells`, `msk_crs`,
  `msk_names`, `msk_categories`, `msk_patches`, `msk_global`, `msk_provenance`.

## drw_*: draw a field

All take an existing mosaik (for its grid) and add one layer named `name`;
the random ones accept `seed`.

- `drw_noise(obj, type = "white", frequency = 4, name = "values", seed = NULL)`:
  `type` is `"white"`, `"perlin"` or `"simplex"`.
- `drw_texture(obj, type = "diamondSquare", base = "perlin", hurst = 0.7, octaves = 6L, lacunarity = 2, frequency = 4, startDev = 1, name = "values", seed = NULL)`:
  `type` is `"diamondSquare"`, `"fbm"`, `"billow"` or `"ridged"`; `base` is
  `"perlin"` or `"simplex"`.
- `drw_gradient(obj, type = "planar", angle = 0, position = c(0.5, 0.5), size = 0.3, origin = NULL, invert = FALSE, name = "values")`:
  `type` is `"planar"`, `"point"`, `"line"`, `"circle"` or `"square"`;
  `origin` names a binary layer of `obj` to measure the gradient from.
- `drw_pattern(obj, type = "checkerboard", frequency = 10, n = NULL, angle = 0, name = "values")`:
  `type` is `"checkerboard"`, `"stripes"`, `"rings"`, `"waves"`, `"grid"` or
  `"hexagonal"`.
- `drw_cluster(obj, type = "percolation", p = 0.5, n = 3L, name = "values", seed = NULL)`:
  `type` is `"percolation"` or `"randomCluster"`.
- `drw_tessellation(obj, type = "voronoi", n = 20L, interaction = 0.1, size = c(5L, 20L), name = "values", seed = NULL)`:
  `type` is `"voronoi"`, `"rectangle"` or `"gibbs"`.

## mdf_*: modify layers

Generic operations; none knows what a layer represents. All have
`obj = NULL` first: called without an object, they record into a recipe.

Values:

- `mdf_binarise(obj, thresh = NULL, match = NULL, layer, add)`: 1 where the value
  is above `thresh` or in `match`, else 0.
- `mdf_categorise(obj, breaks = NULL, n = NULL, shares = NULL, layer, add)`: bin
  into classes, by break points, into `n` classes of equal width, or so that
  the classes cover the given `shares` of the cells (lowest values first;
  must sum to 1; an error if a cut falls among cells of one value). Exactly
  one of the three.
- `mdf_replace(obj, old, new, layer, add)`, `mdf_range(obj, lower, upper, background = NA, layer, add)`,
  `mdf_scale(obj, range, layer, add)`, `mdf_offset(obj, fun = "+", value = 1, layer, add)`,
  `mdf_perturb(obj, sd = 1, layer, add)`.
- `mdf_permute(obj, type = "invert", by = NULL, layer, add)`: `type` is
  `"invert"`, `"revert"`, `"descending"`, `"ascending"` or `"cycle"`.
- `mdf_filter(obj, expr, value = FALSE, layer, add)`: where an expression over
  layer names holds, e.g. `seed == 1 & forest == 1`: a 1/0 mask, or with
  `value = TRUE` the values of `layer` there and `NA` elsewhere.

Shape and morphology:

- `mdf_morph(obj, struct, blend, merge, rotate = FALSE, strict = FALSE, background = NA, layer, add)`:
  the general kernel operation. `blend`: `identity`, `equal`, `lower`,
  `greater`, `plus`, `minus`, `product`. `merge`: `min`, `max`, `all`, `any`,
  `!all`, `!any`, `sum`, `mean`, `median`, `sd`, `cv`, `sumNa`. A 0 in the
  pattern is outside the neighbourhood, except for `blend = "equal"`, where it
  means "must be 0".
- `mdf_dilate`, `mdf_erode`, `mdf_interpolate` (smoothing),
  `mdf_match(obj, struct, rotate = TRUE, ...)` (hit-or-miss): convenience forms
  of `mdf_morph`, all `(obj, struct = NULL, layer, add)`.
- `msk_struct(type = "disc", width = 3, height = 3, rotate = FALSE, background = NA, custom = NULL)`:
  a kernel; `type` is `"disc"`, `"box"`, `"diamond"` or `"cross"`.
- `mdf_componentise(obj, connectivity = 4L, background = NA, layer, add)`:
  number the patches, the connected cells of equal value; 0 and `NA` form none.
  Works on binary and categorical layers, and records the patches with `layer`
  for the patch-level measures.
- `mdf_fill(obj, connectivity = 4L, layer, add)`: set the holes of a binary
  layer to 1.
- `mdf_skeletonise(obj, background = NA, anchor = NULL, method = "zhangSuen", layer, add)`:
  `method` is `"zhangSuen"` or `"homotopic"`.
- `mdf_centroid(obj, background = NA, layer, add)`, `mdf_tesselate(obj, layer, add)`.
- `mdf_distance(obj, source = "foreground", method = "euclidean", layer, add)`:
  distance map; `method` is `"euclidean"`, `"manhattan"` or `"chessboard"`.
  Distance from points or lines: `msk_rasterise()` them first.

Combining and zones:

- `mdf_blend(obj, layers = NULL, fun = "+", weights = NULL, add)`: combine layers
  of one object cell by cell; `fun` is an arithmetic operator, one of the
  summaries of `mdf_summarise()`, or a function. Bring layers in with
  `msk_add()` first.
- `mdf_summarise(obj, by = NULL, fun = "sum", neighbours = FALSE, connectivity = 8L, background = NA, layer, add)`:
  summarise a layer within (or around) each zone of `by`; `fun` is a function or
  `"max"`, `"min"`, `"sum"`, `"mean"`, `"median"`, `"any"`, `"all"`, `"n"`,
  `"n_distinct"`, `"unique"` (the one value, `NA` if several).

Grid:

- `mdf_crop(obj, extent)`, `mdf_pad(obj, width = 1L, sides, value = NA, layer)`,
  `mdf_resize(obj, factor, layer)`, `mdf_rotate(obj, angle = 90L, layer)`,
  `mdf_transpose(obj, layer)`.

Recipes:

- `mdf(obj, recipe)`: apply a recipe to a real mosaik.
- `mdf_loop(obj, recipe, times = Inf, until = NULL, stable = FALSE, layer)`:
  repeat a recipe a number of times, until a condition holds, or until the layer
  stops changing (`stable = TRUE`).

```r
core <- mdf_binarise(match = 47, layer = "cover", add = "forest") |>
  mdf_erode(layer = "forest", add = "core")
result <- mdf(landscape, core)
```

## msr_*: measure the primitives

Every result is stored per layer: `@patches[[layer]]`, `@categories[[layer]]`,
`@global[[layer]]`. Read them with `msk_patches(m, layer)`,
`msk_categories(m, layer)`, `msk_global(m, layer)`. Rewriting a layer drops
its patch and landscape results; rewriting the layer that holds the patch
numbers drops the patches too. Patch-level calls need
`mdf_componentise(layer = ...)` first.

| call | stored in | name in an equation |
|---|---|---|
| `msr_area(scale = "patch")` | `@patches[[layer]]$area` | `area.patch` |
| `msr_area(scale = "class")` | `@categories[[layer]]$area` | `area.class` |
| `msr_area(scale = "landscape")` | `@global[[layer]]$area` | `area.landscape` |
| `msr_perimeter(...)` | as `msr_area`, column `perimeter` | `perimeter.patch/.class/.landscape` |
| `msr_number(scale = "patch")` | `@categories[[layer]]$number` (patches per class) | `number.class` |
| `msr_number(scale = "class")` | `@global[[layer]]$number` (number of classes) | `number.landscape` |
| `msr_adjacency(type = "like")` | `@categories[[layer]]$likeAdj` | `likeAdj.class` |
| `msr_adjacency(type = "pairedSum")` | `@categories[[layer]]$pairedSum` | `pairedSum.class` |
| `msr_adjacency(type = "paired")` | `@global[[layer]]$adjacency` (a matrix) | not usable in equations |
| `msr_dissimilarity(contrast, scale = "class")` | `@categories[[layer]]$dissimilarity` | `dissimilarity.class` |
| `msr_dissimilarity(contrast, scale = "landscape")` | `@global[[layer]]$dissimilarity` | `dissimilarity.landscape` |
| `msr_cost(scale = "patch")` | `@patches[[layer]]$distance` (per class, patch-to-patch matrix) | `distance.patch` |
| `msr_cost(scale = "cell")` | internal layer `_distance_<layer>` | `distance.cell` |

Signatures:

- `msr_area(obj = NULL, scale = "patch", unit = "cells", layer = NULL)`, `unit` is
  `"cells"` or `"map"`; `msr_perimeter` the same.
- `msr_number(obj = NULL, scale = "class", layer = NULL)`, `scale` is `"class"` or
  `"patch"`.
- `msr_adjacency(obj = NULL, scale = "class", type = "like", count = "double", connect = 4, layer = NULL)`.
- `msr_dissimilarity(obj = NULL, contrast, scale = "class", layer = NULL)`, `contrast`
  a symmetric matrix with class IDs as row and column names.
- `msr_cost(obj = NULL, scale = "patch", cost = NULL, routing = "cheapest", accumulate = "sum", layer = NULL)`:
  without `cost`, the cost is distance in metres; with a layer of per-cell
  traversal costs, any other cost. `routing` is `"cheapest"` or `"straight"`;
  `accumulate` is `"sum"`, `"max"`, `"min"`, `"product"` or `"mean"`.

## msr(): compose a metric

`msr(obj = NULL, equation, label, layer = NULL)`. The result goes where its length says:
one value per class to the categories, one per patch to the patches, a single
value to the global slot, always under `layer`. Scales in names: `class`,
`patch`, `landscape`, `cell`. A variable may end in `_<layer>` to read another
layer (`area.class_core / area.class_forest`); without it, it reads `layer`.
Class values of two layers combine only if both have the same classes, patch
values only if both have the same patches; otherwise `msr()` stops (relate
patches of different layers with `mdf_summarise()`).
With `distance.patch`, the equation is evaluated once per patch on that patch's
row of the distance matrix (self-distance is `Inf`), so `"min(distance.patch)"`
is the nearest-neighbour distance. Any `.patch` or `.cell` variable needs the
patches numbered by `mdf_componentise()` on `layer`:

```r
m <- landscape |>
  mdf_componentise(connectivity = 8L, layer = "cover", add = "patch") |>
  msr_cost(scale = "patch", layer = "cover") |>
  msr("min(distance.patch)", "enn", layer = "cover")
```

Costs are named after what they measure, not after the primitive: `distance`
without a cost surface, otherwise the name of the cost layer
(`friction.patch`, `friction.cell`).

```r
m <- landscape |>
  msr_area(scale = "class", layer = "cover") |>
  msr_area(scale = "landscape", layer = "cover") |>
  msr(equation = "area.class / area.landscape * 100", label = "pland",
      layer = "cover")
msk_categories(m, layer = "cover")$pland

m <- msr(m, layer = "cover", label = "shannon",
         equation = "-sum(area.class / area.landscape * log(area.class / area.landscape))")
msk_global(m, layer = "cover")$shannon
```

## Plotting

- `msk_vis(obj, ..., title = NULL, shared_scale = FALSE, window = NULL, theme = NULL, trace = FALSE)`:
  with no `...`, every layer in its own panel.
- `.layer(layer, panel = NULL, colours = NULL, limits = NULL, bins = NULL, hillshade = NULL, legend = TRUE, title = NULL)`:
  how one layer is drawn. Layers with the same `panel` are drawn into one
  panel, bottom to top. `colours` is a palette name, colour stops, or a named
  vector of colours per category label. `hillshade = list(layer = "dem",
  exaggeration = 5)` shades by relief.
- `.theme(from = NULL, ...)`: graphical settings for the whole plot.

```r
msk_vis(landscape, .layer("cover", colours = "viridis"))
```

## Data

- `landscape`: a 60 by 56 grid with layers `cover` (categorical, labelled:
  1 river, 11 arable land, 21 intensive grassland, 24 extensive grassland,
  27 fallow, shrub and clear-cuts, 31 settlement, 35 road, 41 wetland,
  44 orchard, 47 forest and hedgerows) and `canopy` (canopy height in
  metres). Built by `data-raw/landscape.R`; its structures (glades, a
  clear-cut, corridors, a road beside the river, hedgerows) are placed so the
  vignettes find something known.
- `mspa`: the test pattern from the original MSPA paper, used in
  `vignette("mspa")`.

## Further reading

`vignette("mosaik")` (getting started), `vignette("techniques")` (how the
operations combine), `vignette("mspa")` (a published algorithm rebuilt from the
primitives).
