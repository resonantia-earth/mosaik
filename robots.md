# mosaik: reference for language models

This file lets a language model write correct mosaik code without reading the
source. It lists every exported function with its real signature, where results
are stored, and the mistakes that recur. If it disagrees with the help pages,
the help pages are right.

Read "Common mistakes" first.

## Common mistakes

1. **There is no `derive()`.** A metric built from primitives is computed with
   `msr(obj = NULL, equation, label, layer)`, where `equation` is a character string in
   `name.focus` notation:
   `msr(m, "perimeter.self / area.self", "para", layer = "cover")`.
2. **Measure the primitives before `msr()`.** The equation reads results that
   `msr_*` already stored. There is no landscape scale: the area of the layer
   is `sum(area.all)`.
3. **There are no scales.** No `msr_*` takes `scale`, no variable ends in
   `.class`, `.patch`, `.landscape` or `.cell`, and a label has no level. The
   focus words `self`, `others`, `all` decide everything (see `msr()` below).
4. **Patches are a layer, not a scale.** To measure patches, measure the layer
   `mdf_componentise()` wrote: `msr_area(m, layer = "patch")`. There is no
   `msr_number()` and no `msk_patches()`; the number of patches is
   `length(area.all)` on that layer.
5. **`msk_vis()` takes layer specifications, not a `layer` argument.** Write
   `msk_vis(m, .layer("cover"))`, not `msk_vis(m, layer = "cover")`.
6. **`msr_*` results are not returned as tables.** Every function returns the
   mosaik with results attached. Read them with `msk_table(m, layer)`: a
   list with a print method, so `msk_table(m, "cover")$area`.
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
   `msr_distance()` without a cost surface.
9. **C++ functions are internal.** Do not call `morphCpp`, `distanceCpp` and so
   on; use the R functions.
10. **Number the patches before measuring them.** The measures never find
    patches themselves, because the connectivity (4 or 8) must be stated:
    `mdf_componentise(connectivity = 8L, layer = "forest", add = "patch")`,
    then measure `layer = "patch"`.
11. **Reserved names.** No layer and no label may be called `area`,
    `perimeter`, `adjacency`, `distance`, `dissimilarity`, `gid`, `complete`,
    `val`, `colour`, `x` or `y`. Labels, and layers read by their name in an
    equation (`canopy.self`), must not contain `.` or `_`: `msr()` reads the
    `.` as the start of the focus and the `_` as the start of the layer.

## The object

One S4 class, `mosaik`, with these slots:

| slot | holds |
|---|---|
| `extent` | `c(xmin, xmax, ymin, ymax)` |
| `dims` | `c(ncols, nrows)` |
| `layers` | named list of flat row-major vectors, compressed with `rle()` when smaller |
| `categories` | per layer: `gid` (class IDs) and `val` (labels) for a categorical layer, plus the values measured per class |
| `global` | per layer: values for the whole layer, from `msr()` equations with only `.all` |
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
  `msk_names`, `msk_provenance`; what is measured on a layer: `msk_table(obj, layer)`, a list (per class and overall values by name) with a print method.

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
  Works on binary and categorical layers. Each class of the written layer is a
  patch; cells outside every patch are `background` (keep `NA`, or the
  background is measured as one more class).
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
core <- mdf_filter(expr = cover == 47, add = "forest") |>
  mdf_erode(layer = "forest", add = "core")
result <- mdf(landscape, core)
```

## msr_*: measure the primitives

Every primitive is measured for each class of `layer` and stored in its class
table, `@categories[[layer]]`; read it with `msk_table(m, layer)`. On a
layer of patch numbers (from `mdf_componentise()`), each class is a patch, so
the same call measures the patches. Rewriting a layer drops what was measured
on it. Values describe the part of a class on the map, also for a class the
map border cuts; `complete` tells which classes lie wholly on the map.

| call | stored in `@categories[[layer]]` | name in an equation |
|---|---|---|
| `msr_area()` | `$area` | `area` |
| `msr_perimeter()` | `$perimeter`; edges to `NA` cells count, edges along the map border do not | `perimeter` |
| `msr_adjacency()` | `$adjacency` (bordering cell pairs, both sides counted), class by class matrix | `adjacency` |
| `msr_dissimilarity(contrast)` | `$dissimilarity` | `dissimilarity` |
| `msr_distance()` | `$distance`, or under `name`; class by class matrix, `Inf` on the diagonal | `distance` or the `name` |
| any layer, e.g. `canopy` | the layer itself | `canopy` (cells) |

Signatures:

- `msr_area(obj = NULL, unit = "cells", layer = NULL)`, `unit` is `"cells"` or
  `"map"`; `msr_perimeter` the same.
- `msr_adjacency(obj = NULL, connect = 4, layer = NULL)`.
- `msr_dissimilarity(obj = NULL, contrast, layer = NULL)`, `contrast` a
  symmetric matrix with class IDs as row and column names.
- `msr_distance(obj = NULL, cost = NULL, routing = "cheapest", name = "distance", layer = NULL)`:
  between the nearest cells of every two classes. Without `cost`, every cell
  counts the same; with a layer of per-cell costs, the distance is the sum of
  the costs of the cells crossed (`NA` = impassable, unreachable = `NA`).
  `routing` is `"cheapest"` or `"straight"`. `name` must not be a layer.

## msr(): compose a metric

`msr(obj = NULL, equation, label, layer = NULL)`. The equation is evaluated
with one class of `layer` in focus at a time. A variable is
`<name>.<focus>[_<layer>]`:

| focus | reads | areas A 10, B 20, C 5, A in focus |
|---|---|---|
| `.self` | the class in focus | 10 |
| `.others` | every other class | 20, 5 |
| `.all` | all classes | 10, 20, 5 |

- With `.self` or `.others` anywhere, the equation runs once per class and
  the result goes to `@categories[[layer]]`, also when it contains `.all`
  (`area.self / sum(area.all)`: `sum(area.all)` is the same in every run).
  With only `.all`, it runs once and goes to `@global[[layer]]`. Each run
  gives the value of one class, or of the whole layer, so it must come to one
  number: reduce `.others` and `.all` with `sum()`, `max()` and the like;
  `area.all * 2` is an error. The label is a plain name (`"pland"`, not
  `"pland.class"`); it must not be a layer, a reserved name or already stored.
- A name that is a layer reads cells: `canopy.self` = the canopy values in the
  cells of the class in focus. Otherwise the name is a stored value. A value
  stored for the whole layer is read with `.all`.
- Needs no measuring: `gid`, `complete` (no cell on the map border), `x`, `y`
  (cell centres in map units).
- A class by class matrix is read along the row of the class in focus:
  `distance.self` = the diagonal entry (`Inf`), `distance.others` = the rest
  of the row, so `"min(distance.others)"` is the nearest neighbour;
  `distance.all` = the whole row. With only `.all`, the whole matrix.
- `_<layer>` reads another layer: with `.all` its classes
  (`sum(area.all_cover)` = area of the map), with `.self`/`.others` its class
  values through the cells of the classes in focus: after
  `msr_area(layer = "cover")`, `area.self_cover` on the patch layer = for each
  cell of the patch, the area of its cover class.
- `pi` and other base constants may appear.

```r
m <- landscape |>
  mdf_filter(cover == 47, add = "forest") |>
  mdf_componentise(connectivity = 8L, layer = "forest", add = "patch") |>
  msr_area(layer = "patch") |>
  msr_distance(layer = "patch") |>
  msr("min(distance.others)", "enn", layer = "patch") |>
  msr("sum(area.others / distance.others^2 * (distance.others <= 10))", "prox",
      layer = "patch") |>
  msr("max(area.all[complete.all])", "largest", layer = "patch")
msk_table(m, layer = "patch")$enn

m <- landscape |>
  msr_area(layer = "cover") |>
  msr(equation = "area.self / sum(area.all) * 100", label = "pland",
      layer = "cover") |>
  msr(equation = "-sum(area.all / sum(area.all) * log(area.all / sum(area.all)))",
      label = "shannon", layer = "cover")
msk_table(m, layer = "cover")$shannon
```

A distance for every cell is a layer: `mdf_distance()`, with `cost` for a cost
layer.

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
