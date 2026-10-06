# Groups and focus: the next shape of msr()

Discussed 2026-10-06, nothing built. Work on this happens on the branch
`groups`. Read this before touching `msr()`, the `msr_*` primitives or the
metric app.

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
otherwise it means a stored value. No further rule is needed.

## Prototype test, 2026-10-06

`inst/design/groups_prototype.R` (run from the package root) implements the
design outside the package: class values per layer (area, perimeter,
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
