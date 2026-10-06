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

The label keeps its level and with it the focus and the storage: `.class`
evaluates once per class, `.landscape` once for the whole layer (only
`.all` makes sense there).

| metric | equation | label |
|---|---|---|
| PARA | `perimeter.self / area.self` | `para.class` |
| PLAND | `area.self / sum(area.all) * 100` | `pland.class` |
| ENN | `min(distance.others)` | `enn.class` |
| PROX | `sum(area.others / distance.others^2)` | `prox.class` |
| LPI | `max(area.all) / sum(area.all) * 100` | `lpi.landscape` |
| SHDI | `-sum(area.all / sum(area.all) * log(area.all / sum(area.all)))` | `shdi.landscape` |
| mean canopy per patch | `mean(canopy.self)` | `canopy.class` |

No `[]`, no scale in the variable, every variable means the same in every
equation.

**A name is either a layer or a stored value.** `canopy.self` reads cells
(the canopy values in the cells of the class in focus) because `canopy` is
a layer; `area.self` reads a stored value because `area` was measured. The
user sees this as the strong point: any layer enters the interface, from a
`drw_` field to a table tied to coordinates (`msk_rasterise` puts records on
cells, `mdf_replace` turns a column into a layer). Tables tied to classes
(a cost per land cover class) are columns of the class table and read the
same way (`cost.self`). Layer names and stored value names must not clash.

## Open

- **Groups cut by the map border.** Primitives store what they measure, no
  NA. `msr()` leaves out the classes that touch the map border when it
  evaluates an equation. On `cover` almost every class touches the border,
  so `msr()` must know whether a layer's classes are groups or land cover
  classes (marked by `mdf_componentise`, or said by the user). Not decided.
- The NA-for-cut-patches built on 2026-10-05/06 (in `msr_area`,
  `msr_perimeter`, `msr_adjacency`, `msr_distance`, `msr()`) is replaced by
  the rule above once it is decided.

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
