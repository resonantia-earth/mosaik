# mosaik — reference for LLMs

This document is the single source of truth for the mosaik package. It
exists so an LLM can write correct mosaik code from this file alone,
without opening source. When in doubt, this document overrides any
other description. Update this file whenever a function's signature,
behaviour, or the package structure changes.

If you are an LLM reading this: **read the "Common mistakes" section
first**. Most failures with mosaik come from a small set of recurring
errors that the section names directly.

---

## Common mistakes

These are the mistakes that recur most often. Read them first.

### 1. `derive()` takes a character string, not a formula.

`derive()` does not accept a formula. The `equation` argument is a
**character string** that uses `metric.scale` dot notation. The
`label` is a separate argument that names where the result is stored.

```r
# WRONG — formula not accepted
derive(m, edge_density ~ perimeter / area)

# CORRECT
derive(m, equation = "perimeter.class / area.class",
       label = "edge_density")
```

Allowed scale suffixes: `class`, `patch`, `landscape`, `cell`.

### 2. `derive()` requires `msr_*` to have been called first.

The equation reads pre-computed primitives from `@categories[[layer]]`,
`@patches`, or `@global`. There is no fall-through to raw layer
values. Run the matching `msr_*(scale = ...)` first.

```r
# WRONG — area not computed yet
m <- derive(m, "area.class / area.landscape", "pland")

# CORRECT
m <- msr_area(m, scale = "class")
m <- msr_area(m, scale = "landscape")
m <- derive(m, "area.class / area.landscape", "pland")
```

### 3. The metric names in `derive()` are the **column** names, not the function names.

When you call `msr_adjacency(type = "like")`, it stores the result
under `@categories[[layer]]$likeAdj`, **not** `$adjacency`. Use
`likeAdj.class`, not `adjacency.class`.

| `msr_*` call | Stored as | In `derive()` use |
|---|---|---|
| `msr_area(scale="class")` | `$area` | `area.class` |
| `msr_perimeter(scale="patch")` | `$perimeter` | `perimeter.patch` |
| `msr_number(scale="patch")` | `$number` | `number.class` |
| `msr_adjacency(type="like")` | `$likeAdj` | `likeAdj.class` |
| `msr_adjacency(type="paired")` | `@global$adjacency` (matrix) | not directly usable |
| `msr_adjacency(type="pairedSum")` | `$pairedSum` | `pairedSum.class` |
| `msr_distance(scale="patch")` | `$distance` (list of matrices) | `distance.patch` (row-wise mode) |
| `msr_distance(scale="cell")` | `_distance` internal layer | `distance.cell` (per-patch mode) |
| `msr_dissimilarity` | `$dissimilarity` | `dissimilarity.class` or `.landscape` |

### 4. `mg_propensity` is built before `mg_run` and passed in.

It is its own object, not an argument bundle.

```r
# WRONG — propensity is not a list of inline parameters
mg_run(m, types, targets, propensity = list(forest = "score"))

# CORRECT
prop <- mg_propensity(forest = "forest_score", crop = "crop_score")
mg_run(m, types, targets, propensity = prop)
```

### 5. Targeting a derived metric requires `derive()` first.

`mg_target(name = "edge_density")` only resolves at run time if a
prior `derive(label = "edge_density")` is recorded in
`@provenance`. The engine looks up the equation and re-evaluates it
each iteration.

### 6. Inventing functions that don't exist.

The following are NOT in the package: `mk_combine`, `mk_from_terra`,
`mk_from_vector`, `mk_report`, `mk_save`, `mk_load`, `card`, `punch`,
`sim_random`. Do not write code that calls them.

To ingest a SpatRaster, use `mosaik(rast = r)`.
To emit a SpatRaster, use `as_terra(m)`.
There is no built-in stack-two-mosaiks utility.

### 7. `mdf_zonify` takes raw matrices, not sf objects.

`geom` is a closed two-column matrix or a list of such matrices. The
first row must equal the last row. sf geometries are not accepted.

### 8. `gnrt_*` functions create a new mosaik; `sim_*` modify one.

The `gnrt_*` family takes `extent` and `dims` and returns a fresh
mosaik. The `sim_*` family takes an existing mosaik and adds or
overwrites layers. Do not chain `gnrt_gradient(m, ...)` — pass
`extent` and `dims`, not a mosaik.

### 9. Predecessor packages are gone.

Do not document or call `loom`, `weave_*`, `wv_*`, `fabric()`,
`geomio`, `geometr`. Use mosaik's verbs.

### 10. `landscapemetrics`-style return shape is wrong.

`msr_*` does not return a tibble of metrics. It returns a mosaik
with metrics attached at `@global`, `@categories`, or `@patches`.
Read the metric out of those slots after the call.

---

## Object model

mosaik is one S4 class, `mosaik`, with eight slots:

| Slot | Type | Holds |
|---|---|---|
| `extent` | `numeric(4)` | `c(xmin, xmax, ymin, ymax)` |
| `dims` | `integer(2)` | `c(ncols, nrows)` |
| `layers` | named `list` | flat row-major vectors, RLE-compressed when smaller |
| `categories` | named `list` | per-layer list keyed by layer name; each entry has at least `gid` (the class IDs present in the layer), optionally `val` (labels), `role` (semantic tag), and class-level metrics added by `msr_*(scale = "class")` |
| `patches` | `list` | patch-level metrics added by `msr_*(scale = "patch")`. After any patch-level `msr_*`, contains `class`, `patch`, and the metric columns |
| `global` | `list` | landscape-level metrics added by `msr_*(scale = "landscape")` |
| `crs` | `character(1)` | proj4 string or `NA_character_` |
| `provenance` | `list` | one entry per operation; each is a named list produced by `.make_prov(fn_name, args)` |

A second S4 class, `struct`, holds a structuring element (kernel) for
morphological `mdf_*` operators. Slots: `pattern` (matrix), `rotate`
(logical), `background` (numeric).

A third S4 class, `mkTheme`, holds visual parameters for `mk_vis`.
See `mk_theme()`.

### Layer values

Layer values are stored as flat row-major vectors. `mk_pull(obj, "x")`
returns the vector; `inverse.rle()` is applied transparently when the
storage is RLE-compressed. Use `mk_pull` for read access; do not
touch `obj@layers[[name]]` directly.

### Categories and the `role` tag

`@categories` is keyed by layer name. The minimal entry is a list
with `gid` (integer class IDs). Continuous layers may have only
`role` (a semantic tag like `"surface"`, `"precipitation"`,
`"temperature"`, `"water"`, `"vegetation"`, `"settlement"`,
`"propensity"`). Several `sim_*` functions auto-detect their inputs
by scanning categories for these roles via the internal
`.find_role(obj, role_name)` lookup.

### Patches roster

Once any `msr_*(scale = "patch")` has been called, an internal
layer named `_patches` is created carrying globally unique patch
IDs, and `@patches` carries at least `class` and `patch` columns.
The roster is reused across patch-scale `msr_*` calls; you do not
need to call `mdf_componentise` first.

### Provenance

Each operation appends a named list to `@provenance`. `derive()`
records the equation and label so `mg_target()` / `mg_run()` can
look it up later. There is currently no provenance rendering or
diffing function.

---

## Minimum working example

```r
library(mosaik)

# 1. start from the bundled example landscape
m <- landscape

# 2. transform: keep only forest cells (class 47)
forest <- mdf_binarise(m, match = 47, layer = "cover")

# 3. measure: number of forest patches and their areas
forest <- forest |>
  mdf_componentise() |>
  msr_area(scale = "patch") |>
  msr_perimeter(scale = "patch")

# 4. derive: shape index per patch
forest <- derive(forest,
                 equation = "perimeter.patch / sqrt(area.patch)",
                 label = "shape_index")

forest@patches$shape_index

# 5. simulate change: convert some forest to non-forest
types <- mg_landtype(name = "forest", resistance = 0.3,
                     produces = list(timber = 5),
                     allow = list(other = list(strategy =
                                  mg_strategy("frontier")))) |>
         mg_landtype(name = "other", resistance = 0.1,
                     produces = list(timber = 0),
                     allow = list())

targets <- mg_target(name = "timber", target = 200, match = "minimum")

result <- mg_run(forest, types = types, targets = targets,
                 layer = "cover")

# 6. visualise side by side
mk_vis(before = forest, after = result, layer = "cover")
```

---

## Constructors and I/O

### `mosaik(extent, res, crs, vals, rast, group)`

User-facing constructor. Two paths:

- **From scratch.** Provide `extent` (numeric(4): xmin, xmax, ymin,
  ymax) and `res` (numeric(1) or numeric(2)). `vals` is `NULL`,
  a matrix, or a named list of matrices/vectors. `crs` is a proj4
  string or `NA_character_` (default).
- **From a SpatRaster.** Provide `rast = your_spatraster`. `extent`,
  `res`, `crs`, `vals` are ignored. `group = TRUE` stores unique
  values as categories when the raster has none. Requires the
  `terra` package.

```r
# empty grid
m <- mosaik(extent = c(0, 100, 0, 100), res = 1)

# multi-layer
m <- mosaik(extent = c(0, 10, 0, 10), res = 1,
            vals = list(cover = sample(1:5, 100, replace = TRUE),
                        ndvi  = runif(100)))

# from terra
m <- mosaik(rast = my_spatraster, group = TRUE)
```

### `new_mosaik(extent, dims, layers, categories, patches, global, crs, provenance)`

Low-level constructor. Validates inputs and auto-RLE-compresses
layers. Used internally; user code should prefer `mosaik()`.

### `as_terra(obj)`

Emits a `SpatRaster` with one layer per mosaik layer. CRS and
categories are transferred. Provenance is dropped (with a warning).

---

## sim_* — generate

Each `sim_*` takes an existing mosaik and adds or overwrites layers.
Most call `.find_role(obj, ...)` to locate inputs and tag their
output with a semantic role.

### `make_climate(obj, latitude, continentality, wind_direction, wind_strength, base_precip, base_temp, temp_offset, precip_regime, temp_monthly, precip_monthly, layer, add, seed)`

Generate temperature and precipitation surfaces from a
`"surface"`-tagged terrain plus geographic context.

- `latitude` -90 to 90, default 45 — sets base temperature and the
  seasonal amplitude; the sign selects the hemisphere
- `continentality` 0 (fully maritime) to 1 (continental interior),
  default 0.5 — sets the annual temperature amplitude and **nothing
  else**. Not distance to coast: a cold-current coast is continental in
  behaviour.
- `wind_direction` degrees 0–360, default 270 (westerlies)
- `wind_strength` 0–1, default 0.5 — orographic rain shadow intensity
- `base_precip` mm/year, default 800 — a stated boundary condition, so
  it is the realised landscape mean and nothing rescales it
- `base_temp` °C at sea level; if NULL, derived from latitude
- `temp_offset` °C added to every month, default 0
- `precip_regime` one of `"uniform"`, `"summer_max"`, `"winter_max"`,
  `"dry_summer"`, `"monsoon"`, `"bimodal"`, `"polar"`
- `temp_monthly` 12 absolute °C, overrides the computed cycle
- `precip_monthly` 12 proportions summing to 1, overrides the regime

Adds layers `temperature` (role `"temperature"`, °C) and
`precipitation` (role `"precipitation"`, mm/year). Temperature uses the
environmental lapse rate (-6.5 °C/1000 m) plus a latitude gradient;
precipitation uses the linear-orographic upslope-wind model, normalised
so the landscape mean equals `base_precip`.

Seasonality lives in `@categories`, not in layers:
`@categories$temperature$monthly` is twelve absolute °C (a cosine, since
insolation is sinusoidal over the year) and
`@categories$precipitation$monthly` is twelve proportions of the annual
total (a von Mises per regime). Landscape-level because circulation sets
rain timing at scales far larger than a mosaik. `make_vegetation` reads
both to compute a monthly water balance.

### `sim_flow(obj, precipitation, infiltration, sea_level, method, layer, add, diagnostics, seed)`

Generate a drainage network from a `"surface"`-tagged terrain layer
using depression-aware multi-flow routing (Holmgren 1994; Qin 2007),
priority-flood depression filling (Barnes 2014), Manning's equation
for water depth, and Leopold–Maddock hydraulic geometry for channel
visibility.

- `precipitation` numeric default 1.0, or a layer name (auto-detects
  a `"precipitation"`-tagged layer if default is used)
- `infiltration` numeric in 0–1 default 0.5, or a layer name
- `method` `"depression"` (default, routes on original terrain with
  overflow-aware filling) or `"filled"` (routes on filled surface for
  continuous streams)
- `diagnostics` if TRUE, also adds layers `accumulation`, `volume`,
  `depth`, `fill_ratio`, `sinks`, `pour_points`

Output layer tagged with `role = "water"`. Cells classified as
`"water"` (river or lake) where the channel/lake fills a significant
fraction of the cell.

### `make_vegetation(obj, targets, layer, add, seed)`

Compute five continuous structural surfaces describing vegetation:
`height` (m), `cover` (0–1), `woodiness` (0–1), `evergreenness` (0–1),
`LAI` (0–~10). One computation parameterised by a disturbance regime.
With no pressure present it yields the potential natural vegetation
(climate-driven equilibrium); with a pressure signal present it blends
each field toward a use's target, `field = (1 - s)·equilibrium +
s·target`. Auto-detects layers tagged `"temperature"`,
`"precipitation"`, `"temp_seasonality"`, `"precip_seasonality"`,
`"surface"`, and optionally `"water"`. With no climate at all, falls
back to spatially autocorrelated surfaces per field (neutral mode).

The pressure signal is read from two role-tagged layers when present:
`"pressure_type"` (categorical, value = a `targets` key) and
`"pressure_strength"` (continuous 0–1). Pressure functions
(`make_production`, `make_urban`, `make_disturbance`) write these.

- `targets` data.frame mapping pressure `type` to a target structural
  vector (columns `height`, `cover`, `woodiness`, `LAI`; optional
  `evergreenness`). NULL = built-in table (arable, pasture, forestry,
  urban, clearcut).

Classification into biomes/habitats is a separate downstream step
(`mdf_classify`). Each output layer is tagged with `role` equal to its
field name (`"height"`, `"cover"`, …). Equilibrium response curves are
literature-shaped but uncalibrated (pipeline-level calibration pass).

### `sim_urban(obj, cover, productivity, delta_cost, speed, seeds, cell_size, add, seed)`

Grow settlements via iterative cheapest-cell allocation on a cost
surface combining terrain resistance (1 - propensity), construction
cost, and agglomeration attraction. Grounded in bid-rent theory
(Alonso 1964), Marchetti's constant (30-min travel budget), and urban
scaling laws (Bettencourt 2007).

- `cover` required, 0–1 — target fraction of viable cells.
- `productivity` 0.2–1.2 controls kernel amplitude (high = compact
  cities).
- `delta_cost` scalar or layer name — densification cost.
- `speed` km/h sets kernel radius via Marchetti's constant:
  walking = 5, cart = 15, car = 40.
- `seeds` initial nucleation count (only on the first call; on
  subsequent calls the function continues from the existing layer).

Calls can be chained: walking → cart → car for staged historical
development.

Adds three layers: `<add>` (settlement count, role `"settlement"`),
`<add>_order` (allocation tick, role `"allocation_order"`),
`<add>_nucleus` (nucleus ID, role `"nucleus"`).

### `sim_production(obj, class, types, weights, gradient, size, size_sd, geometry, fill, layer, add, seed)`

Fill available space with managed parcels (cropland, pasture, or
forestry) via a Voronoi/Delaunay merge process (Gaucherel 2008).
Respects barriers from layers tagged `"settlement"`, `"network"`,
`"water"`, `"place"`. Assigns sub-types across parcels by ranking
parcels along a spatial `gradient` and partitioning the rank space
according to per-type `weights`.

- `class` one of `"cropland"`, `"pasture"`, `"forestry"`; obligatory.
  Drives default `size`, default `geometry`, default `types`, and the
  `role = "production"` tag.
- `types` character vector of sub-type labels. Order is load-bearing:
  `types[1]` lands on the parcels ranked lowest by `gradient` (closest
  to settlements/roads under Thünen); the last entry lands on the
  highest-ranked parcels. `NULL` → class default vocabulary.
- `weights` numeric vector, one per `types` entry. Relative shares of
  the rank space. `NULL` → equal shares.
- `gradient` one of `"thunen"` (default), `"uniform"`, or a numeric
  vector of length `ncells(obj)`. Thünen degrades to uniform when
  neither `"settlement"` nor `"network"` is present.
- `size` numeric ≥ 1 — mean parcel size in cells. `NULL` → class
  default: cropland 50, pasture 100, forestry 200.
- `size_sd` numeric ≥ 0, default 0.5 — coefficient of variation.
- `geometry` one of `"regular"`, `"irregular"`, `"elongated"`. `NULL`
  → class default: cropland and forestry use `"regular"`; pasture uses
  `"irregular"`.
- `fill` logical, default TRUE — fill all available space.
- `add` `NULL` overwrites `layer`, otherwise writes to that new layer.

Output is tagged with `role = "production"`. Default vocabularies:
- cropland: `intensive_crop`, `extensive_crop`, `pasture`, `fallow`, `orchard`
- pasture: `intensive_pasture`, `extensive_pasture`, `rough_grazing`, `meadow`, `wood_pasture`
- forestry: `coppice`, `plantation`, `mixed_stand`, `selective_cut`, `conservation`

Mixed systems (agroforestry, silvopasture, Hutewald) are not separate
classes — pass labels like `"agroforestry_alley"` via `types`.

---

## gnrt_* — pattern primitives

Each `gnrt_*` returns a fresh `mosaik` from `extent` and `dims`. All
have `name` (default `"values"`), optional `role`, and `seed`.

### `gnrt_gradient(extent, dims, origin, type, angle, position, size, invert, name, role, seed)`

Distance-from-origin scaled to [0, 1].

- `origin` optional binary mosaik whose foreground cells are the
  origin. Overrides `type` if given.
- `type` one of `"planar"` (default), `"point"`, `"line"`,
  `"circle"`, `"rectangle"`, `"square"`, `"polygon"`, `"ellipse"`,
  `"triangle"`, `"hexagon"`
- `angle` rotation in degrees (planar gradient direction or line angle)
- `position` `c(x, y)` in 0–1 grid-relative coordinates, default centre
- `size` shape size as fraction of shorter grid dimension, default 0.3
- `invert` TRUE = high values near origin

Deterministic; `seed` is accepted but ignored.

### `gnrt_noise(extent, dims, type, frequency, name, role, seed)`

- `type` one of `"white"` (uniform random), `"perlin"`, `"simplex"`
- `frequency` ≥ 0.01, default 4 — higher = finer detail (Perlin/Simplex)

### `gnrt_pattern(extent, dims, type, frequency, n, angle, name, role, seed)`

Deterministic geometric patterns.

- `type` one of `"checkerboard"`, `"stripes"`, `"rings"`, `"waves"`,
  `"grid"`, `"hexagonal"`
- `frequency` block size (checkerboard, grid), stripe width (stripes),
  number of rings (rings), wave cycles (waves), tile size (hexagonal)
- `n` number of distinct classes for `"checkerboard"` (default 2) and
  `"hexagonal"` (default 3). `n = 0` gives each tile a unique sequential
  ID — useful for zoning. Ignored for other types.
- `angle` rotation in degrees (stripes, waves only)

### `gnrt_cluster(extent, dims, type, p, n, name, role, seed)`

Stochastic spatial growth.

- `type` `"percolation"` (binary random map at probability `p`) or
  `"randomCluster"` (Saura & Martínez-Millán 2000: percolate, label
  components, assign random class IDs, fill remaining background by
  iterative neighbour propagation)
- `p` percolation probability, default 0.5
- `n` number of classes for `"randomCluster"`, default 3 (ignored for
  percolation)

### `gnrt_tessellation(extent, dims, type, n, interaction, name, role, seed)`

- `type` `"voronoi"` (random seeds), `"rectangle"` (random rectangular
  blocks), `"gibbs"` (Strauss-process inhibition before Voronoi)
- `n` number of regions/seeds, default 20
- `interaction` for Gibbs only: inhibition radius as fraction of the
  grid diagonal, default 0.1

### `gnrt_texture(extent, dims, type, base, hurst, octaves, lacunarity, frequency, startDev, name, role, seed)`

Multi-scale spatially autocorrelated continuous surfaces.

- `type` `"diamondSquare"` (default), `"fbm"`, `"billow"` (abs-value
  fBm), `"ridged"` (inverted ridged fBm)
- `base` octave noise function: `"perlin"` (default) or `"simplex"`.
  Ignored for diamondSquare.
- `hurst` 0–1 roughness exponent, default 0.7 (0 = very rough,
  1 = very smooth)
- `octaves` 1–16, default 6 (fBm/billow/ridged only)
- `lacunarity` ≥ 1, default 2.0 (frequency multiplier per octave)
- `frequency` base frequency, default 4
- `startDev` initial standard deviation, default 1 (diamondSquare only)

Output scaled to [0, 1].

---

## mdf_* — transform

Each `mdf_*` takes a mosaik, modifies one or more layers, records
the operation, and returns a mosaik. The `add` argument is universal:
`NULL` overwrites the input layer, a string writes to a new layer.
The `layer` argument defaults to the first layer.

### Geometric

#### `mdf_resize(obj, factor, layer)`
Up- or down-scale via nearest-neighbour. `factor` integer > 1 (up) or
inverse integer 1/n (down). Other values are coerced to the nearest
allowed value.

#### `mdf_rotate(obj, angle, layer)`
Clockwise rotation by 90, 180, or 270 degrees. Adjusts extent and
dims. `layer = NULL` rotates all layers. R-loop implementation;
slow on very large grids.

#### `mdf_transpose(obj, layer)`
Mirror along the main diagonal. Swaps dims, adjusts extent,
preserves categories. R-loop implementation.

#### `mdf_crop(obj, extent)`
Subset to a smaller `extent = c(xmin, xmax, ymin, ymax)`. Snaps
outward to cell boundaries. Preserves categories, patches, global
(though the latter two may reference values no longer present in
the cropped region).

#### `mdf_pad(obj, width, sides, value, layer)`
Add (positive `width`) or remove (negative `width`) cells on
specified sides. `sides` is a subset of `c("left", "right", "top",
"bottom")`, default all. Added cells filled with `value` (default NA).
`layer = NULL` pads all layers.

#### `mdf_offset(obj, fun, value, layer, add)`
Apply an arithmetic operator to all cell values. `fun` one of
`"+"`, `"-"`, `"*"`, `"/"`, `"%%"`, `"%/%"`, `"^"`, `"**"`. `value`
must be integerish (despite the name).

#### `mdf_centroid(obj, background, layer, add)`
Replace each foreground patch's cells with `background` except for
the patch's centroid cell, which carries the patch ID. Patch IDs
are taken from the layer values; binarise + componentise first if
needed.

### Value-level

#### `mdf_binarise(obj, thresh, match, layer, add)`
Set cells to 0/1. With `thresh`, cells with value > `thresh` become 1.
With `match`, cells whose value is in `match` become 1. Use exactly
one of the two.

#### `mdf_categorise(obj, breaks, n, layer, add)`
Bin values into integer categories. Provide `breaks` (numeric vector)
or `n` (number of equal-width bins). Attaches a categories entry
with `gid` and human-readable `val` labels like `"0-25"`.

#### `mdf_replace(obj, old, new, layer, add)`
Substitute values. `old` and `new` are vectors; if lengths differ,
`new` is recycled. `old = NA` replaces NA cells with `new`.

#### `mdf_range(obj, lower, upper, background, layer, add)`
Set cells outside `[lower, upper]` to `background` (default NA).
Either bound can be NULL.

#### `mdf_scale(obj, range, layer, add)`
Linearly rescale values to `range = c(min, max)`.

#### `mdf_perturb(obj, sd, layer, add)`
Add Gaussian noise N(0, `sd`) to cell values. Preserves NA.

#### `mdf_permute(obj, type, by, layer, add)`
Reassign distinct values among themselves.

- `type` `"invert"` (max - v), `"revert"` (reverse sorted order),
  `"descending"`, `"ascending"`, `"cycle"` (shift by `by`, wrapping)

#### `mdf_distance(obj, source, coords, snap, method, layer, add)`
Distance transform.

- `source` `"foreground"` (default, distance to nearest 1-cell;
  binary input required), `"background"` (distance to nearest 0-cell;
  binary input required), or a layer name (cells are matched by
  patch ID against target points in the named layer)
- `coords` two-column matrix of source coordinates in CRS space.
  When provided, `source` is ignored and the input does not need to
  be binary.
- `snap` only used with `coords`. TRUE (default) snaps to nearest
  centroid and uses the fast Meijster algorithm. FALSE computes
  exact Euclidean distances from the original coordinates.
- `method` `"euclidean"`, `"manhattan"`, `"chessboard"`. Only used
  with the Meijster algorithm.

Returned distances are in cell units (Meijster) or CRS units
(`snap = FALSE`).

### Morphological

#### `mdf_morph(obj, struct, blend, merge, rotate, strict, background, layer, add)`
The core primitive. Iterates a structuring element over the grid;
at each position, blends covered cell values with the kernel pairwise
(`blend`), then merges them into one scalar (`merge`).

- `struct` a `struct` object (see `mk_struct`); default 3×3 disc
- `blend` one of `"identity"`, `"equal"`, `"lower"`, `"greater"`,
  `"plus"`, `"minus"`, `"product"`
- `merge` one of `"min"`, `"max"`, `"all"`, `"any"`, `"!all"`,
  `"!any"`, `"sum"`, `"mean"`, `"median"`, `"sd"`, `"cv"`, `"sumNa"`
- `rotate` apply all rotations of the kernel, default TRUE
- `strict` strict matching, default TRUE

Most other morphological functions are convenience wrappers around
this.

#### `mdf_dilate(obj, struct, layer, add)`
Default kernel = 3×3 disc. Increases values matching the kernel.

#### `mdf_erode(obj, struct, layer, add)`
Default kernel = 3×3 disc. Decreases values matching the kernel.

#### `mdf_skeletonise(obj, background, layer, add)`
Topology-preserving morphological skeleton. Input must be binary.

#### `mdf_componentise(obj, connectivity, background, layer, add)`
Connected-component labelling on a binary input. `connectivity` is
4 (rook, default) or 8 (queen).

#### `mdf_match(obj, struct, rotate, background, layer, add)`
Hit-or-miss transform. Convenience wrapper for `mdf_morph` with
`blend = "equal"`, `merge = "sumNa"`.

#### `mdf_interpolate(obj, struct, layer, add)`
Spatial smoothing (weighted mean filter). Convenience wrapper for
`mdf_morph` with `blend = "product"`, `merge = "mean"`,
`rotate = FALSE`, `strict = FALSE`. Default kernel = 3×3 disc.

### Structural

#### `mdf_layerise(obj, by, flatten, background, layer)`
Segregate distinct values into one layer per value. `flatten = TRUE`
sets foreground to 1; FALSE keeps original values. `by` optional
mosaik whose values determine the segregation; if NULL, uses `obj`'s
own values. Output layers are RLE-compressed.

#### `mdf_mask(obj, by, background, layer, add)`
Set cells where `by` (a binary single-layer mosaik of the same
dimensions) is 0 to `background` (default NA).

#### `mdf_blend(obj, layers, fun, weights, add)`
Combine two or more layers of one mosaik into a single layer.

- `layers` the layers to blend, at least two; NULL blends all layers
- `fun` arithmetic operator string (vectorised pairwise reduction) or
  a custom function of one value per layer (uses `reduceCpp`)
- `weights` per-layer multipliers applied before blending
- `add` NULL overwrites `layers[1]`, otherwise names a new layer

Blending is within one object: bring a layer in from another mosaik
with `msk_add()` first, then blend. Replaces the former `mdf_reduce()`.

#### `mdf_tesselate(obj, layer, add)`
Voronoi tessellation: assign every cell to the nearest patch
centroid in the input. Requires ≥ 2 distinct patch values.

#### `mdf_zonify(obj, geom, value, background, layer, add)`
Burn polygon zones into a layer.

- `geom` two-column matrix `cbind(x, y)` with first row = last row,
  or a list of such matrices
- `value` integerish; recycled or matched per polygon
- `background` value for cells outside all polygons (default NA);
  set NULL to leave existing values untouched

Cell-centroid-in-polygon test. Later polygons overwrite earlier ones.

---

## msr_* — measure

Each `msr_*` takes a mosaik and writes results to one of `@global`,
`@categories[[layer]]`, or `@patches`, depending on `scale`. The
input mosaik is returned (with metrics attached). All accept
`layer = NULL` (defaults to first).

The internal `.ensure_patch_layer(obj, layer)` is called by
patch-scale `msr_*` and creates a hidden `_patches` layer plus
`@patches$class` and `@patches$patch` columns. Subsequent
patch-scale calls reuse it.

### `msr_area(obj, scale, unit, layer)`
- `scale` `"patch"` (default), `"class"`, or `"landscape"`
- `unit` `"cells"` (default) or `"map"` (multiplied by cell area)

Stores: `@patches$area`, `@categories[[layer]]$area`, or
`@global$area`.

### `msr_perimeter(obj, scale, unit, layer)`
Same scales and units as `msr_area`. Stores `$perimeter`.

### `msr_number(obj, scale, layer)`
- `scale` `"class"` (default; counts distinct classes, stored as
  `@global$number`) or `"patch"` (counts patches per class, stored as
  `@categories[[layer]]$number`)

### `msr_adjacency(obj, type, count, layer)`
Cell-adjacency matrix.

- `type` `"like"` (default; diagonal of adjacency matrix, stored as
  `@categories[[layer]]$likeAdj`), `"paired"` (full matrix, stored as
  `@global$adjacency`), `"pairedSum"` (row sums, stored as
  `@categories[[layer]]$pairedSum`)
- `count` `"single"` (right and bottom neighbours only) or
  `"double"` (default; also left and top)

### `msr_distance(obj, scale, source, layer)`
- `scale = "patch"` (default): pairwise nearest-edge Euclidean
  distance between all patches of each class. Stored as
  `@patches$distance` — a named list of symmetric matrices per
  class, with `Inf` on the diagonal.
- `scale = "cell"`: per-cell distance stored in an internal layer
  named `_distance`. Sub-mode by `source`:
  - `"centroid"` (default): distance from each cell to its own
    patch's centroid. Computes a centroid layer if not present.
  - `"background"`: distance from each foreground cell to the
    nearest background. Binary input required.
  - layer name: passed to `mdf_distance(source = ...)`.

In `derive()`, `distance.patch` triggers per-patch row-wise
evaluation; `distance.cell` triggers per-patch cell-grouped
evaluation.

### `msr_dissimilarity(obj, contrast, scale, layer)`
Contrast-weighted edge dissimilarity.

- `contrast` symmetric numeric matrix with row and column names
  matching class values (as character). Diagonal should be 0.
- `scale` `"class"` (default; per-class row sum, stored as
  `@categories[[layer]]$dissimilarity`) or `"landscape"` (total /
  2, stored as `@global$dissimilarity`)

### `msr_values(obj, ...)`
Reads layer values into the attribute tables. Used as a primitive
in extensions of `derive()`. See source for argument details.

---

## mg_* — simulate land-use change

CLUMondo-style allocation: each cell is scored for every reachable
type based on propensity, conversion resistance, and demand
elasticity; the cell becomes whichever scores highest. The engine
iterates until supply meets demand within tolerance.

### `mg_landtype(obj, name, resistance, produces, allow)`
Define one land type. Pipe to add more.

- `name` unquoted (or quoted) type name
- `resistance` 0–1, conversion inertia (0 = easily converted,
  1 = permanent), default 0
- `produces` named list of goods → per-cell production rate, e.g.
  `list(timber = 5, carbon = 10)`. Default `list()`
- `allow` named list of permitted transitions. Each element is keyed
  by target type and contains:
  - `strategy` an `mg_strategy` object (default
    `mg_strategy("frontier")`)
  - `min` minimum timesteps before this transition (default 0)
  - `max` maximum timesteps before forced (default `Inf`)
  - `auto` timestep at which to auto-convert (default NULL)

Use `character(0)` or `list()` for immutable types.

Returns an `mg_landtypes` collection (S3 list).

### `mg_strategy(type, reach, expand, min_distance, density, exposure, attractors, decay, threshold)`
Selection strategy for a transition.

- `type` `"frontier"` (only cells within `reach` of existing
  target-type patches; `expand` = TRUE auto-grows the frontier when
  demand isn't met), `"suitability"` (global score ranking),
  `"nucleation"` (new patches at `min_distance` from existing
  same-type, `density` fraction returned as candidates),
  `"retreat"` (cells with ≥ `exposure` non-same-type 8-neighbours),
  `"gravity"` (interaction potential between attractor features in
  layer `attractors`, with `decay` exponent and `threshold` minimum)

Returns an `mg_strategy` object (S3 list).

### `mg_propensity(...)`
Declare which mosaik layer(s) encode per-type suitability. Two
shapes:

- One unnamed argument: a single shared layer for all types
- Multiple named arguments: per-type layer names. Types not listed
  get uniform propensity = 1.

```r
mg_propensity("overall_score")                          # shared
mg_propensity(forest = "forest_score", crop = "crop_score")  # per-type
```

Returns an `mg_propensity` object (S3 list).

### `mg_target(obj, name, target, match, class)`
Add one target to an `mg_targets` collection.

- `name` for goods targets, the goods name; for metric targets, the
  metric name
- `target` single number (global) or data frame with columns
  `region` (integer region IDs) and `value` (per-region targets)
- `match` `"exact"`, `"minimum"`, or `"maximum"`
- `class` NULL = goods target; otherwise the class name (or
  `"landscape"`) for a metric target. Metric names can be primitives
  (`area`, `number`, `perimeter`, `adjacency`, `distance`,
  `dissimilarity`) or labels created by a prior `derive()` call.

Returns an `mg_targets` collection (S3 list with `goods` and
`metrics`).

### `mg_run(obj, types, targets, n, record, propensity, regions, noise, contrast, layer, age_layer, max_iter, tolerance, speed, dampen)`
Run the engine.

- `obj` mosaik
- `types` `mg_landtypes` from `mg_landtype()`
- `targets` `mg_targets` from `mg_target()` **or** a function
  `function(t, obj)` returning `mg_targets` per timestep (dynamic)
- `n` integer ≥ 1, default 1 — number of timesteps
- `record` if TRUE and n > 1, returns a list of mosaiks (one per
  timestep). Otherwise returns the final mosaik only.
- `propensity` an `mg_propensity` object, a single layer name
  (character), or NULL (uniform)
- `regions` name of an integer layer defining accounting regions
  for per-region targets. NULL = all targets are global.
- `noise` standard deviation of per-cell score noise. 0 =
  deterministic, default 0.
- `contrast` optional contrast matrix for dissimilarity metrics
  (same format as `msr_dissimilarity`).
- `layer` the cover layer to allocate over (default first layer)
- `age_layer` name of the age tracking layer, default `"_age"`,
  created automatically
- `max_iter`, `tolerance`, `speed`, `dampen` convergence parameters

Returns a single mosaik (or a list when `n > 1` and
`record = TRUE`). Provenance entry records `iter`, `converged`,
goods names, and metric names.

---

## derive()

`derive(obj, equation, label, layer)`.

The equation is a character string using `metric.scale` dot
notation. Allowed scales: `class`, `patch`, `landscape`, `cell`.

```r
m <- msr_area(m, scale = "class")
m <- msr_perimeter(m, scale = "class")
m <- derive(m, equation = "perimeter.class / area.class",
            label = "edge_density")
```

Storage level is determined by result length: equal to
`length(@categories[[layer]]$gid)` → class, equal to
`nrow(@patches)` → patch, scalar → landscape.

Two specialised modes:

- **Matrix-valued patch primitives.** If any variable is
  `distance.patch` (a per-class list of n×n matrices), evaluation is
  row-wise per patch: that variable resolves to the patch's row,
  other `.patch` variables resolve to the full vector for the
  patch's class, `.class` and `.landscape` pass through. The Inf
  diagonal naturally excludes self-terms (e.g. `1/Inf = 0`). Must
  reduce to a scalar per patch.
- **Cell scale.** A `metric.cell` variable reads the internal layer
  `_<metric>`, groups cell values by patch ID, and evaluates the
  equation per patch. Must reduce to a scalar per patch.

Provenance records the full equation under
`provenance[[k]]$derive$args` so `mg_target()` / `mg_run()` can
look it up by `label`.

---

## mk_* — utilities

### Accessors (`R/004_accessors.R`)
- `mk_extent(x)` — `numeric(4)`
- `mk_dims(x)` — `c(ncols, nrows)`
- `mk_res(x)` — `c(res_x, res_y)`
- `mk_ncells(x)` — total cells
- `mk_crs(x)` — proj4 string or NA
- `mk_provenance(x)` — the provenance list
- `mk_names(x)` — layer names

### Tabular utilities (`R/005_utils.R`)
- `mk_select(obj, ...)` — keep specified layers (unquoted or
  character names; e.g. `mk_select(m, cover, intensity)`).
  Categories for non-selected layers are dropped.
- `mk_filter(obj, expr, layer)` — mask cells where `expr` evaluates
  to FALSE to NA. The expression is evaluated with all layer
  values available as variables (e.g. `mk_filter(m, cover == 47)`).
  Sets cells in `layer` to NA where the predicate fails.
- `mk_pull(obj, layer)` — extract a layer's values as a
  numeric/integer vector. RLE storage is decompressed transparently.

### Visualisation
- `msk_vis(..., layer, overlay, hillshade, window, theme,
  trace, new)` — plot one or more mosaiks. Multiple mosaiks become
  multiple panels (named via the `...` argument names if provided).
  `layer = NULL` shows all layers as separate panels per mosaik;
  pass a layer name to show only that one. `overlay = list(layer =,
  colours =)` draws a layer on top (rivers over terrain).
  `hillshade = list(layer =, azimuth =, altitude =, exaggeration =,
  intensity =)` shades the panel by an elevation layer's relief and
  is applied before the overlay. `trace = TRUE` prints provenance.
- `mk_theme(from, title, box, xAxis, yAxis, grid, legend, scale,
  parameters)` — build or modify an `mkTheme` for `mk_vis`.

### Structuring elements
- `mk_struct(type, width, height, rotate, background, custom)` —
  create a `struct` for morphological operators. Default `type =
  "disc"`, `width = 3`, `height = 3`. Pass a matrix to `custom` to
  override the built-in patterns.

---

## Recipes

Direct mappings from intent to call.

### Compute landscape Shannon diversity
```r
m <- msr_area(landscape, scale = "class")
m <- msr_area(m, scale = "landscape")
m <- derive(m,
            equation = "-sum((area.class / area.landscape) *
                            log(area.class / area.landscape))",
            label = "shannon")
m@global$shannon  # stored at landscape level since result is scalar
```

### Compute proportion of landscape per class (PLAND)
```r
m <- msr_area(landscape, scale = "class")
m <- msr_area(m, scale = "landscape")
m <- derive(m, equation = "area.class / area.landscape * 100",
            label = "pland")
m@categories$cover$pland
```

### Compute mean patch size per class
```r
m <- msr_area(landscape, scale = "class")
m <- msr_number(m, scale = "patch")
m <- derive(m, equation = "area.class / number.class",
            label = "mean_patch_size")
```

### Compute edge density per class
```r
m <- msr_perimeter(landscape, scale = "class")
m <- msr_area(m, scale = "landscape")
m <- derive(m, equation = "perimeter.class / area.landscape",
            label = "edge_density")
```

### Compute shape index per patch
```r
m <- msr_area(landscape, scale = "patch")
m <- msr_perimeter(m, scale = "patch")
m <- derive(m, equation = "perimeter.patch / sqrt(area.patch)",
            label = "shape_index")
```

### Find largest patch per class
```r
m <- msr_area(landscape, scale = "patch")
sapply(split(m@patches$area, m@patches$class), max)
```

### Compute nearest-neighbour distance per patch
```r
m <- msr_distance(landscape, scale = "patch")
m <- derive(m, equation = "min(distance.patch)", label = "enn")
m@patches$enn
```

### Rasterise a polygon into a layer
```r
poly <- matrix(c(10,10, 50,10, 50,40, 10,40, 10,10),
               ncol = 2, byrow = TRUE)
m <- mdf_zonify(landscape, geom = poly, value = 99L,
                layer = "cover", add = "zone")
```

### Build a generated landscape from terrain to vegetation
```r
m <- gnrt_texture(extent = c(0, 100, 0, 100), dims = c(100L, 100L),
                  type = "diamondSquare",
                  name = "terrain", role = "surface", seed = 42)
m <- sim_climate(m, latitude = 50)
m <- sim_flow(m, add = "water")
m <- sim_vegetation(m, add = "vegetation")
mk_vis(m)
```

### Run a simple multi-objective land-use scenario
```r
types <- mg_landtype(name = "forest", resistance = 0.5,
                     produces = list(timber = 5, carbon = 10),
                     allow = list(crop = list(strategy =
                                  mg_strategy("frontier")))) |>
         mg_landtype(name = "crop", resistance = 0.1,
                     produces = list(food = 8),
                     allow = list(forest = list(auto = 25)))

targets <- mg_target(name = "food",   target = 1000, match = "minimum") |>
           mg_target(name = "carbon", target = 5000, match = "minimum")

prop <- mg_propensity(forest = "forest_suit", crop = "crop_suit")

result <- mg_run(landscape, types = types, targets = targets,
                 propensity = prop, n = 5, record = TRUE)
```

### Smooth a continuous layer
```r
m <- mdf_interpolate(landscape, layer = "intensity")  # 3x3 mean
m <- mdf_interpolate(landscape, layer = "intensity",
                     struct = mk_struct("square", width = 5))
```

### Detect edges of a class
```r
forest <- mdf_binarise(landscape, match = 47, layer = "cover")
edges <- mdf_morph(forest, struct = mk_struct("disc", width = 3),
                   blend = "equal", merge = "all", strict = TRUE,
                   background = 0)
```

### Compute distance from each cell to nearest road
```r
roads <- mdf_binarise(landscape, match = c(11, 21), layer = "cover")
dist  <- mdf_distance(roads, source = "foreground",
                      method = "euclidean")
```

### Visualise two mosaiks side by side
```r
mk_vis(before = original, after = result, layer = "cover")
```

---

## Conventions

- **Pipe-based.** Every verb takes a mosaik in and returns a mosaik
  out. Metrics attach to the returned object; they are not
  separate return values.
- **`layer` parameter.** Most functions accept `layer = NULL` and
  default to the first layer name. Set explicitly to operate on a
  named layer in a multi-layer mosaik.
- **`add` parameter.** Most `mdf_*` and `sim_*` functions accept
  `add = NULL` (overwrite the input layer) or `add = "newname"`
  (write to a new layer).
- **Provenance accumulates.** Every operation appends to
  `@provenance` via the internal `.make_prov(fn_name, args)`.
- **RLE compression.** Layers are auto-compressed with `rle()` when
  the compressed form is smaller. `mk_pull()` decompresses
  transparently. Do not access `obj@layers[[name]]` directly; use
  `mk_pull`.
- **S4 with formal validity.** `setValidity("mosaik", ...)` checks
  extent, dims, layer lengths, and slot types. Constructors run
  this via `new()`.
- **Role tags.** Continuous layers can carry a semantic `role`
  (e.g. `"surface"`, `"precipitation"`, `"temperature"`, `"water"`,
  `"vegetation"`, `"settlement"`, `"propensity"`). `sim_*` functions
  auto-detect inputs by scanning for these roles.
- **C++ backend in `src/`.** Internal, called via `.Call()` through
  R wrappers. Includes `morphCpp`, `componentsCpp`, `binariseCpp`,
  `distanceCpp`, `skeletoniseCpp`, `rescaleGridCpp`, `reduceCpp`,
  `tesselateCpp`, `pointInPolyCpp`, `countCellEdgesCpp`,
  `countCellAdjacenciesCpp`, `countCellValuesCpp`, `perlinCpp`,
  `simplexCpp`, `diamondSquareCpp`, `percolationCpp`,
  `mergeRegionsCpp`, `flowSurfaceCpp`, `fillDepressionsCpp`,
  `settleGrowCpp`, `distanceFromPointsCpp`, `isBinaryCpp`. All take
  flat vectors (not matrices) for layer data.
- **Dependencies.** `Imports`: checkmate, crayon, digest, grDevices,
  grid, methods, Rcpp, utils. `Suggests`: knitr, rmarkdown, sf,
  terra, testthat. Tidyverse, mmand, dplyr, rlang are deliberately
  removed.

## Internal vs. external boundary

- **Public API** (export, document, recommend): the functions named
  in this document.
- **Bundled data**: `landscape` (a 60 × 56 example mosaik with
  `cover` and `intensity` layers; 9 cover classes: 1, 11, 21, 24,
  27, 31, 41, 44, 47).
- **Internal helpers** (do not call): functions with leading dot —
  `.find_role`, `.make_prov`, `.replace_layer`,
  `.replace_layer_categorical`, `.resolve_add`, `.next_type_values`,
  `.ensure_patch_layer`, `.find_derive_provenance`, and the
  `.derive_*` helpers.
- **C++ entry points** (use the R wrapper, not `.Call` directly).
- **`R/deprecated.R`** holds `mk_from_terra` and `mk_as_terra` shims
  that should not be referenced. Use `mosaik(rast = ...)` and
  `as_terra()` instead.
- **`R/RcppExports.R`** is auto-generated.
- **`R/001_show.R`** defines the `show` method for `mosaik`.

## Bundled example dataset

`landscape` — a 60 × 56 example mosaik with extent
`c(0, 60, 0, 56)`, resolution 1, two layers: `cover` (integer
land-cover class with 9 classes: 1, 11, 21, 24, 27, 31, 41, 44, 47)
and `intensity` (continuous values in 0–100). Used in most examples
in this document. Loaded automatically with `library(mosaik)`.
