# mosaik: reference for language models

This file lets a language model write correct mosaik code without reading the
source. It lists every exported function with its real signature, where results
are stored, and the mistakes that recur. If it disagrees with the help pages,
the help pages are right.

Read "Common mistakes" first.

## Common mistakes

1. **There is no `derive()`.** A metric built from primitives is computed with
   `msr(obj, equation, label, layer)`, where `equation` is a character string in
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
7. **Recipes hold `mdf_*` steps only.** Called without an object, an `mdf_*`
   function records itself into a recipe. `msr_*` functions always run
   immediately and cannot be recorded. Apply the recipe with `mdf()`, then
   measure.
8. **Functions that do not exist:** `derive`, `as_terra`, and anything with the
   prefixes `mk_`, `gnrt_`, `mg_`, `sim_`, `make_`. Terrain, climate, soil,
   vegetation, land-use and simulation functions belong to the separate package
   mundus, not to mosaik. `msr_distance` does not exist either: distance is
   `msr_cost()` without a cost surface.
9. **C++ functions are internal.** Do not call `morphCpp`, `distanceCpp` and so
   on; use the R functions.

## The object

One S4 class, `mosaik`, with these slots:

| slot | holds |
|---|---|
| `extent` | `c(xmin, xmax, ymin, ymax)` |
| `dims` | `c(ncols, nrows)` |
| `layers` | named list of flat row-major vectors, compressed with `rle()` when smaller |
| `categories` | per layer: `gid` (class IDs) and `val` (labels) for a categorical layer, plus class-level results |
| `patches` | patch-level results: `class`, `patch`, then one column per metric |
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
- `msk_add(obj, from, ..., rename = NULL)`: copy layers from another mosaik on
  the same grid (extent, dims and CRS are checked); category tables come along.
- `msk_select(obj, ...)`, `msk_remove(obj, ...)`: keep or drop layers, named
  unquoted.
- `msk_set(obj, layer, values, prov = NULL, keep = TRUE, gid = NULL, val = NULL)`:
  write raw values into a layer. Meant for other packages; users combine
  objects with `msk_add()`. `gid` and `val` make the layer categorical.
- `msk_prov(fn, args, derivedFrom = NULL, activity = list(), step = FALSE)`:
  build a provenance entry for `msk_set()`.
- Accessors: `msk_extent`, `msk_dims`, `msk_res`, `msk_ncells`, `msk_crs`,
  `msk_names`, `msk_categories`, `msk_patches`, `msk_global`, `msk_provenance`.

## syn_*: synthesise a field

All take an existing mosaik (for its grid), add one layer named `name`, and
accept `seed`.

- `syn_noise(obj, type = "white", frequency = 4, name = "values", seed = NULL)`:
  `type` is `"white"`, `"perlin"` or `"simplex"`.
- `syn_texture(obj, type = "diamondSquare", base = "perlin", hurst = 0.7, octaves = 6L, lacunarity = 2, frequency = 4, startDev = 1, name = "values", seed = NULL)`:
  `type` is `"diamondSquare"`, `"fbm"`, `"billow"` or `"ridged"`; `base` is
  `"perlin"` or `"simplex"`.
- `syn_gradient(obj, origin = NULL, type = "planar", angle = 0, position = c(0.5, 0.5), size = 0.3, invert = FALSE, name = "values", seed = NULL)`:
  `type` is `"planar"`, `"point"`, `"line"`, `"circle"`, `"rectangle"`,
  `"square"`, `"polygon"`, `"ellipse"`, `"triangle"` or `"hexagon"`.
- `syn_pattern(obj, type = "checkerboard", frequency = 10, n = NULL, angle = 0, name = "values", seed = NULL)`:
  `type` is `"checkerboard"`, `"stripes"`, `"rings"`, `"waves"`, `"grid"` or
  `"hexagonal"`.
- `syn_cluster(obj, type = "percolation", p = 0.5, n = 3L, name = "values", seed = NULL)`:
  `type` is `"percolation"` or `"randomCluster"`.
- `syn_tessellation(obj, type = "voronoi", n = 20L, interaction = 0.1, name = "values", seed = NULL)`:
  `type` is `"voronoi"`, `"rectangle"` or `"gibbs"`.

## mdf_*: modify layers

Generic operations; none knows what a layer represents. All have
`obj = NULL` first: called without an object, they record into a recipe.

Values:

- `mdf_binarise(obj, thresh = NULL, match = NULL, layer, add)`: 1 where the value
  is above `thresh` or in `match`, else 0.
- `mdf_categorise(obj, breaks = NULL, n = NULL, layer, add)`: bin into classes.
- `mdf_replace(obj, old, new, layer, add)`, `mdf_range(obj, lower, upper, background = NA, layer, add)`,
  `mdf_scale(obj, range, layer, add)`, `mdf_offset(obj, fun = "+", value = 1, layer, add)`,
  `mdf_perturb(obj, sd = 1, layer, add)`.
- `mdf_permute(obj, type = "invert", by = NULL, layer, add)`: `type` is
  `"invert"`, `"revert"`, `"descending"`, `"ascending"` or `"cycle"`.
- `mdf_filter(obj, expr, value = FALSE, background = NA, layer, add)`: keep cells
  where an expression over layer names holds, e.g. `seed == 1 & forest == 1`.

Shape and morphology:

- `mdf_morph(obj, struct, blend, merge, rotate = TRUE, strict = TRUE, background = NA, layer, add)`:
  the general kernel operation. `blend`: `identity`, `equal`, `lower`,
  `greater`, `plus`, `minus`, `product`. `merge`: `min`, `max`, `all`, `any`,
  `sum`, `mean`, `median`, `sd`, `cv`, `one`, `zero`, `na`.
- `mdf_dilate`, `mdf_erode`, `mdf_interpolate` (smoothing),
  `mdf_match(obj, struct, rotate = TRUE, ...)` (hit-or-miss): convenience forms
  of `mdf_morph`, all `(obj, struct = NULL, layer, add)`.
- `msk_struct(type = "disc", width = 3, height = 3, rotate = FALSE, background = NA, custom = NULL)`:
  a kernel; `type` is `"disc"`, `"box"`, `"diamond"` or `"cross"`.
- `mdf_componentise(obj, connectivity = 4L, background = NA, layer, add)`: label
  connected patches.
- `mdf_fill(obj, background = 0, value = NULL, connectivity = 4L, layer, add)`:
  fill holes enclosed by a patch.
- `mdf_skeletonise(obj, background = NA, anchor = NULL, method = "zhangSuen", layer, add)`:
  `method` is `"zhangSuen"` or `"homotopic"`.
- `mdf_centroid(obj, background = NA, layer, add)`, `mdf_tesselate(obj, layer, add)`.
- `mdf_distance(obj, source = "foreground", coords = NULL, snap = TRUE, method = "euclidean", layer, add)`:
  distance map; `method` is `"euclidean"`, `"manhattan"` or `"chessboard"`.

Combining and zones:

- `mdf_blend(obj, layers = NULL, fun = "+", weights = NULL, add)`: combine layers
  of one object cell by cell. Bring layers in with `msk_add()` first.
- `mdf_mask(obj, by, background = NA, layer, add)`: keep cells where layer `by`
  is non-zero.
- `mdf_zonal(obj, by = NULL, fun = "sum", neighbours = FALSE, connectivity = 8L, background = NA, layer, add)`:
  summarise a layer within (or around) each zone of `by`; `fun` is a function or
  `"max"`, `"min"`, `"sum"`, `"mean"`, `"median"`, `"any"`, `"all"`, `"n"`,
  `"n_distinct"`.
- `mdf_zonify(obj, geom, value, background = NA, layer, add)`: burn polygons,
  given as closed two-column coordinate matrices, into a layer.
- `mdf_layerise(obj, by = NULL, flatten = FALSE, background = NA, layer)`: one
  layer per value.

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

| call | stored in | name in an equation |
|---|---|---|
| `msr_area(scale = "patch")` | `@patches$area` | `area.patch` |
| `msr_area(scale = "class")` | `@categories[[layer]]$area` | `area.class` |
| `msr_area(scale = "landscape")` | `@global$area` | `area.landscape` |
| `msr_perimeter(...)` | as `msr_area`, column `perimeter` | `perimeter.patch/.class/.landscape` |
| `msr_number(scale = "patch")` | `@categories[[layer]]$number` (patches per class) | `number.class` |
| `msr_number(scale = "class")` | `@global$number` (number of classes) | `number.landscape` |
| `msr_adjacency(type = "like")` | `@categories[[layer]]$likeAdj` | `likeAdj.class` |
| `msr_adjacency(type = "pairedSum")` | `@categories[[layer]]$pairedSum` | `pairedSum.class` |
| `msr_adjacency(type = "paired")` | `@global$adjacency` (a matrix) | not usable in equations |
| `msr_dissimilarity(contrast, scale = "class")` | `@categories[[layer]]$dissimilarity` | `dissimilarity.class` |
| `msr_dissimilarity(contrast, scale = "landscape")` | `@global$dissimilarity` | `dissimilarity.landscape` |
| `msr_cost(scale = "patch")` | `@patches$distance` (per class, patch-to-patch matrix) | `distance.patch` |
| `msr_cost(scale = "cell")` | internal layer `_distance` | `distance.cell` |

Signatures:

- `msr_area(obj, scale = "patch", unit = "cells", layer = NULL)`, `unit` is
  `"cells"` or `"map"`; `msr_perimeter` the same.
- `msr_number(obj, scale = "class", layer = NULL)`, `scale` is `"class"` or
  `"patch"`.
- `msr_adjacency(obj, scale = "class", type = "like", count = "double", connect = 4, layer = NULL)`.
- `msr_dissimilarity(obj, contrast, scale = "class", layer = NULL)`, `contrast`
  a symmetric matrix with class IDs as row and column names.
- `msr_cost(obj, scale = "patch", cost = NULL, routing = "cheapest", accumulate = "sum", layer = NULL)`:
  without `cost`, the cost is distance in metres; with a layer of per-cell
  traversal costs, any other cost. `routing` is `"cheapest"` or `"straight"`;
  `accumulate` is `"sum"`, `"max"`, `"min"`, `"product"` or `"mean"`.

## msr(): compose a metric

`msr(obj, equation, label, layer = NULL)`. The result goes where its length says:
one value per class to the categories, one per patch to the patches, a single
value to the global slot. Scales in names: `class`, `patch`, `landscape`, `cell`.
With `distance.patch`, the equation is evaluated once per patch on that patch's
row of the distance matrix (self-distance is `Inf`), so `"min(distance.patch)"`
is the nearest-neighbour distance. Any `.patch` variable needs the patch table,
which only `msr_area(scale = "patch")` or `msr_perimeter(scale = "patch")`
create; `msr_cost()` alone does not:

```r
m <- landscape |>
  msr_area(scale = "patch", layer = "cover") |>
  msr_cost(scale = "patch", layer = "cover") |>
  msr("min(distance.patch)", "enn", layer = "cover")
```

```r
m <- landscape |>
  msr_area(scale = "class", layer = "cover") |>
  msr_area(scale = "landscape", layer = "cover") |>
  msr(equation = "area.class / area.landscape * 100", label = "pland",
      layer = "cover")
msk_categories(m)$cover$pland

m <- msr(m, layer = "cover", label = "shannon",
         equation = "-sum(area.class / area.landscape * log(area.class / area.landscape))")
msk_global(m)$shannon
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

- `landscape`: a 60 by 56 grid with layers `cover` (categorical, classes 1, 11,
  21, 24, 27, 31, 41, 44, 47) and `intensity` (continuous).
- `mspa`: the test pattern from the original MSPA paper, used in
  `vignette("mspa")`.

## Further reading

`vignette("mosaik")` (getting started), `vignette("recipes")` (pipelines),
`vignette("measurement")` (the argument behind the design), `vignette("mspa")`
(a published algorithm rebuilt from the primitives).
