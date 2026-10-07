# Groups and focus: the next shape of msr()

Discussed and built 2026-10-06 on the branch `groups`; see "As built" at the
end for what the build decided. Read this before touching `msr()`, the `msr_*`
primitives or the metric app.

## Why

Two problems surfaced on 2026-10-06:

1. PROX (and SIMI, and "the area of the nearest patch") need the values of
   the *other* patches next to a patch's own distances. Under the group rule
   (the label's scale is the group, every variable holds the group's own
   values) no variable reaches them. PROX has no code in the app since.
2. Patches are not paired consistently across primitives. `msr_area(scale =
   "patch")` stores one list over all patches of the layer, `msr_distance`
   one table per class. On a map with several classes the two do not line
   up. The habitat form (one class) hides this.

Classes and cells are intrinsic and need no assumptions. Patches are
arbitrary groups, so they should not be built into the measuring.

## Decided direction

The user's name for the whole mechanism: a **spatio-thematic group_by**. A
grouping can be spatial (patches, zones, tessellation cells) or thematic
(land cover classes); `msr()` treats both as `group_by()` followed by an
equation. The focus (`self`, `others`) corresponds to dplyr's
`cur_group_rows()` (dplyr >= 1.1), which marks the current group within the
whole table; in mosaik `.others` comes ready as a value instead of being
reached by hand.

**Groups are the classes of a group layer.** A grouping is a layer of
numbers: every cell carries the number of its group. `mdf_componentise`
makes one kind (connected areas), `msk_rasterise` (zones), `mdf_tesselate`
or the user's own layer make others. Several groupings live side by side
under their own names. Their attribute table is the class table of that
layer (`@categories[["patch"]]`). Groups within classes is something a user
builds, not a default.

Consequences:

- The patch scale goes: `scale = "patch"` of every `msr_*`, the `@patches`
  slot, `msk_patches()`, the `.patch` handling in `msr()`.
- The landscape scale of the primitives goes too: every landscape value they
  measure is a sum or count over the classes (`area.landscape` =
  `sum(area.class)`, `perimeter.landscape` = `sum(perimeter.class)`,
  `number.landscape` = `length(gid.class)`). The primitives always measure
  per class of a layer and lose their `scale` argument.
- `msr_distance` measures between the classes of a layer: one table, in the
  order of the class table, so it lines up with the other class values. On
  a group layer these are distances between groups (also across what used to
  be classes, which SIMI needs); on `cover`, between land cover classes.
- The metric app changes throughout: its middle column (each primitive at
  patch, class and landscape scale) has nothing left to show. It could show
  the focus instead (`area.self`, `area.others`, `area.all`), which tells
  which metrics look only at the class in focus and which at the others.
  All metric code moves to the new notation; the decision tree stays.

**Variables name a value and a focus, not a scale.** Notation
`<name>.<focus>_<layer>`; the scale moves out of the variable. Within one
evaluation (one class in focus):

| focus     | meaning                       | A in focus, areas A 10, B 20, C 5 |
|-----------|-------------------------------|-----------------------------------|
| `.self`   | the class in focus            | 10                                |
| `.others` | every other class of the layer| 20, 5                             |
| `.all`    | all classes                   | 10, 20, 5                         |

The label needs no level (user, 2026-10-06): the focus words decide it. An
equation with `.self` or `.others` runs once per class and is stored in the
class table; an equation with only `.all` runs once and is stored as one
value for the layer. There is no third case (only `.all` evaluated per class
would give every class the same value). The label is then just a name
(`prox`, not `prox.class`), and the old ambiguity of storing by the number
of values is gone.

| metric | equation | label | result |
|---|---|---|---|
| PARA | `perimeter.self / area.self` | `para` | one per class |
| PLAND | `area.self / sum(area.all) * 100` | `pland` | one per class |
| ENN | `min(distance.others)` | `enn` | one per class |
| PROX | `sum(area.others / distance.others^2)` | `prox` | one per class |
| LPI | `max(area.all) / sum(area.all) * 100` | `lpi` | one for the layer |
| SHDI | `-sum(area.all / sum(area.all) * log(area.all / sum(area.all)))` | `shdi` | one for the layer |
| mean canopy per patch | `mean(canopy.self)` | `canopy` | one per class |

No `[]`, no scale in the variable, every variable means the same in every
equation.

**A name is either a layer or a stored value.** `canopy.self` reads cells
(the canopy values in the cells of the class in focus) because `canopy` is
a layer; `area.self` reads a stored value because `area` was measured. This
is the strong point of the design, and it is not new: `msr_distance(cost =
"friction")` already lets a layer's name become the name of the measured
value (`friction.patch`). Any layer enters the interface, from a
`drw_` field to a table tied to coordinates (`msk_rasterise` puts records on
cells, `mdf_replace` turns a column into a layer). Tables tied to classes
(a cost per land cover class) are columns of the class table and read the
same way (`cost.self`).

Names (user, 2026-10-06): the names of the primitives are reserved, no layer
may take them, and layer names contain no `.` or `_`; documenting this is
enough. A name in an equation that matches a layer always means that layer;
otherwise it means a stored value. No further rule is needed. Weighted
distances (`msr_distance` with a cost layer) are stored under a name the
user chooses, not automatically under the cost layer's name as today.

Reserved names (no layer may take them):

| name | what it holds |
|---|---|
| `area`, `perimeter`, `adjacency`, `distance`, `dissimilarity` | the primitives (`number` only if `msr_number` survives; without scales, the number of classes is `length(gid.all)`) |
| `gid` | the code of each class |
| `complete` | TRUE if a class has no cell on the map border |
| `val`, `colour` | label and colour of each class (class table) |
| `x`, `y` | the coordinates of every cell's centre, in map units |

Names R knows as constants (`pi`, `T`, `F`, `LETTERS`, ...) are read as
constants, so a layer with such a name could not be reached either.

`x` and `y` are cell values provided without a layer. The radius of
gyration of each patch, without adding coordinate layers by hand:
`mean(sqrt((x.self - mean(x.self))^2 + (y.self - mean(y.self))^2))`, on the
group layer `patch`, one value per patch.

## Prototype test, 2026-10-06

`inst/design/archive/groups_prototype/groups_prototype.R` (archived; it ran
against the old API and no longer runs) implemented the design outside the
package: class values per layer (area, perimeter,
adjacency table, distance table between classes, border flag), an
evaluator for `<name>.<self|others|all>[_layer]` with the label's level as
focus, and all 42 metrics of the app rewritten in the new notation
(`groups_prototype_today.R` holds today's code for comparison). It assumes
the open border rule as "leave cut groups out of the population".

Result on `landscape` (forest patches as the group layer `patch`):

- **37 metrics give exactly today's values**: AREA PERIM GYRATE CA PLAND LPI
  TE ED PARA SHAPE FRAC PAFRAC LSI CORE NCORE CAI TCA CPLAND CWED TECI PLADJ
  AI CONTAG IJI DIVISION SPLIT MESH COHESION PR PRD RPR SHDI SIDI MSIDI SHEI
  SIEI MSIEI. PROX is back and equals the old values:
  `sum(area.others / distance.others^2 * (distance.others <= 10))`.
- **4 differ, all for one reason, the border rule**: leaving the 3 cut
  patches out of the population makes NP 14 instead of 17 (and PD), CONNECT
  28.6 instead of 21.7, and ENN of patch 8 3.61 instead of 2.24 (its nearest
  neighbour is a cut patch). So "leave cut groups out" is wrong. A cut group
  must stay in the population (it exists, it can be a neighbour, it counts);
  only its own result is unknown. That is today's rule (own row NA, the way
  to it measured), and it points the open decision below in that direction.
- **The name clash happened at once**: the CORE metric was first labelled
  `core.class`, and then `core.self` in CAI read that stored value instead of
  the layer `core` (CAI 1.18 instead of 53.39). The prototype looked up
  stored values before layers; with layers first (see Names above) the
  clash does not arise.
- **A group layer's outside cells must count as "outside", not as no data**:
  the prototype sets them to 0 before measuring; with NA, the edges along
  the outside would not be counted and every perimeter would be too short.

## Classes cut by the map border (decided 2026-10-06)

Every class of every layer gets a column `complete` in its class table: TRUE
if none of its cells lies on the map border, FALSE otherwise. It is measured
like any value, for patches and land cover classes alike, so no function has
to know whether a layer's classes are groups. Primitives store what they
measure, no NA, and `msr()` leaves nothing out on its own. The user decides
per equation, e.g. the largest complete patch:
`max(area.all[complete.all])`. Cut groups stay in the population as
neighbours and in counts (the prototype showed leaving them out is wrong).
`complete` is a reserved name like the primitives.

The NA-for-cut-patches built on 2026-10-05/06 (in `msr_area`,
`msr_perimeter`, `msr_adjacency`, `msr_distance`, `msr()`) goes when this is
built.

## Rejected on the way (do not propose again)

- `.partner` notation (partner is not a scale), a `partner()` function.
- A `where =` argument (filters every term of the equation at once;
  filtering belongs in the layer).
- A plural scale (`area.patches`) as population mark.
- `self` as a TRUE/FALSE value multiplied into a term (`area * !self` sets
  the focus to 0 instead of removing it: wrong for mean, min, max), and
  `area.patch[!self]` (superseded by `.others`).
- `area.patch` meaning "the class's other patches" next to `distance.patch`
  (the old code; a meaning that depends on context).

## As built (2026-10-06)

The new `msr()` reproduces the prototype's 42 metrics exactly, with cut
patches kept in. Decisions taken during the build (user-approved):

- **Primitives.** `msr_area`, `msr_perimeter`, `msr_adjacency`,
  `msr_distance`, `msr_dissimilarity` measure per class of `layer` only and
  store in `@categories[[layer]]` (helper `.store_class`, in the order of the
  class table). `msr_number` is deleted (`length(area.all)`).
- **`msr_adjacency`** lost `type` as well as `scale`: it always stores the full
  `adjacency` matrix (both sides counted). `regions` (in how many separate
  places two patches touch; built for MSPA's loops, which was then solved
  without it) and `patchAdjacencyCpp` were REMOVED 2026-10-07: nothing used
  them. They are in the git history if a recipe needs loops (network
  redundancy). `like` = `adjacency.self`;
  the per-class total is `adjacency.self + sum(adjacency.others)`, NOT
  `sum(adjacency.all)`, because an equation with only `.all` runs once for the
  layer and then reads the whole matrix.
- **`msr_perimeter`** counts the edge to an `NA` cell (an NA cell is not of
  the class); edges along the map border are still not counted. This is what
  gives a patch layer (outside = NA) its full outlines without any function
  knowing the layer is a grouping.
- **`msr_distance`** measures between the nearest cells of every two classes,
  one matrix, `Inf` on the diagonal; `name =` (default `distance`) is the
  user-chosen name, refused if it is a layer name.
- **`complete`, `gid`, `x`, `y`** are worked out by `msr()` from the layer,
  not stored by a primitive: nothing has to run first.
- **Other layers.** With `.all`, `_<layer>` reads that layer's classes
  (`sum(area.all_cover)`); with `.self`/`.others` it reads that layer's class
  values through the cells of the classes in focus (`mean(cost.self_cover)`).
  The old rule (class values of two layers combine only with identical
  classes) is gone. A layer name takes no `_layer` suffix; a matrix of another
  layer is read only as a whole.
- **No patch record.** `mdf_componentise` writes only the layer of patch
  numbers; the `@patches` slot, `msk_patches`, `.patches_of` and
  `.unlink_patches` are gone. The source class of a patch is read through its
  cells (`msr(m, "cover.self[1]", "source", layer = "patch")`), and a patch can
  hold several classes. With `background` other than `NA`, the background is
  one more class.
- **Rewriting a layer drops its class table.** `.update_mosaik` used to keep
  every non-`gid` field on a rewrite, so `mdf_erode` left the old `area` and
  `msr()` read it as current. Now only the fields of a layer without classes
  (attached by mundus) survive.
- Labels: no `.` or `_`, not a layer, not reserved, not already stored.
  Reserved: `area perimeter adjacency distance dissimilarity gid
  complete val colour x y`.
- Known risk, untested at size: `msr_distance(cost =, routing = "straight")`
  compares every cell of one class with every cell of the other
  (`.straight_path`), which may exhaust memory on large land cover classes.

## Decided 2026-10-07 (user)

- **Map-border edges are counted nowhere.** FRAGSTATS counts them in PERIM,
  LSI and PLADJ but not by default in TE; mosaik does not copy that
  inconsistency. Edges to NA cells count (see As built), edges along the map
  border never; `complete` flags the classes the border cuts.
- **A rewrite drops the class table** (an operator writing new values into an
  existing layer, e.g. `mdf_erode(layer = "forest")` without `add`); the
  geometry operators (crop, pad, resize) keep the values and warn, so a user
  can compare before and after on purpose.
- **Classes without core** get 0 in TCA/CPLAND now (every class is
  evaluated); the old difference to FRAGSTATS is gone.
- **ENN of a patch without neighbour** is the equation's business, not the
  package's: `min()` of nothing is `Inf` in R. The metric app writes
  `if (length(distance.others)) min(distance.others) else NA`, because there
  is no nearest neighbour, not an infinitely far one.
