# The metric app's codes. HAB is the habitat head every patch metric starts with.
HAB <- '# habitat: the classes of the map the species lives in
mosaik(...) |>
  mdf_filter(cover %in% habitat, add = "habitat") |>
  mdf_componentise(connectivity = 8L, layer = "habitat", add = "patch") |>'

CORE <- '  mdf_distance(source = "background", method = "manhattan", layer = "habitat",
               add = "depth") |>
  mdf_filter(depth >= 2, add = "core") |>'

codes <- list(
AREA = paste0(HAB, '
  msr_area(scale = "patch", unit = "map", layer = "habitat") |>
  msr(
    equation = "area.patch / 10000",
    label = "area_ha.patch", layer = "habitat"
  )'),
PERIM = paste0(HAB, '
  msr_perimeter(scale = "patch", unit = "map", layer = "habitat") |>
  msr(
    equation = "perimeter.patch",
    label = "perim.patch", layer = "habitat"
  )'),
GYRATE = paste0('m <- mosaik(...)
dims <- msk_dims(m)
# the x and y coordinate of every cell, as the layers x and y
x <- rep(seq_len(dims[1]), dims[2]) * msk_res(m)[1]
y <- rep(rev(seq_len(dims[2])), each = dims[1]) * msk_res(m)[2]
', sub("mosaik(...)", "m", HAB, fixed = TRUE), '
  msk_add(x = x, y = y) |>
  msr(
    equation = "mean(sqrt((x.cell - mean(x.cell))^2 +
                          (y.cell - mean(y.cell))^2))",
    label = "gyrate.patch", layer = "habitat"
  )'),
CA = 'mosaik(...) |>
  msr_area(scale = "class", unit = "map", layer = "cover") |>
  msr(
    equation = "area.class / 10000",
    label = "ca.class", layer = "cover"
  )',
PLAND = 'mosaik(...) |>
  msr_area(scale = "class", layer = "cover") |>
  msr_area(scale = "landscape", layer = "cover") |>
  msr(
    equation = "area.class / area.landscape * 100",
    label = "pland.class", layer = "cover"
  )',
LPI = paste0(HAB, '
  msr_area(scale = "patch", layer = "habitat") |>
  msr_area(scale = "landscape", layer = "habitat") |>
  msr(
    equation = "max(area.patch, na.rm = TRUE) / area.landscape * 100",
    label = "lpi.landscape", layer = "habitat"
  )'),
TE = 'mosaik(...) |>
  msr_perimeter(scale = "class", unit = "map", layer = "cover") |>
  msr(
    equation = "perimeter.class",
    label = "te.class", layer = "cover"
  )',
ED = 'mosaik(...) |>
  msr_perimeter(scale = "class", unit = "map", layer = "cover") |>
  msr_area(scale = "landscape", unit = "map", layer = "cover") |>
  msr(
    equation = "perimeter.class / area.landscape * 10000",
    label = "ed.class", layer = "cover"
  )',
PARA = paste0(HAB, '
  msr_perimeter(scale = "patch", unit = "map", layer = "habitat") |>
  msr_area(scale = "patch", unit = "map", layer = "habitat") |>
  msr(
    equation = "perimeter.patch / area.patch",
    label = "para.patch", layer = "habitat"
  )'),
SHAPE = paste0(HAB, '
  msr_perimeter(scale = "patch", layer = "habitat") |>
  msr_area(scale = "patch", layer = "habitat") |>
  # the smallest perimeter a patch of this many cells can have
  msr(
    equation = "4 * floor(sqrt(area.patch)) +
      2 * (area.patch > floor(sqrt(area.patch))^2) +
      2 * (area.patch - floor(sqrt(area.patch))^2 > floor(sqrt(area.patch)))",
    label = "minperim.patch", layer = "habitat"
  ) |>
  msr(
    equation = "perimeter.patch / minperim.patch",
    label = "shape.patch", layer = "habitat"
  )'),
FRAC = paste0(HAB, '
  msr_perimeter(scale = "patch", unit = "map", layer = "habitat") |>
  msr_area(scale = "patch", unit = "map", layer = "habitat") |>
  msr(
    # a single cell has no shape to measure; FRAGSTATS sets it to 1
    equation = "ifelse(area.patch == 1, 1,
                       2 * log(0.25 * perimeter.patch) / log(area.patch))",
    label = "frac.patch", layer = "habitat"
  )'),
PAFRAC = paste0(HAB, '
  msr_perimeter(scale = "patch", unit = "map", layer = "habitat") |>
  msr_area(scale = "patch", unit = "map", layer = "habitat") |>
  msr(
    equation = "2 / (cov(log(area.patch), log(perimeter.patch),
                         use = \\"complete.obs\\") /
                     var(log(perimeter.patch), na.rm = TRUE))",
    label = "pafrac.landscape", layer = "habitat"
  )'),
LSI = 'mosaik(...) |>
  msr_perimeter(scale = "class", layer = "cover") |>
  msr_area(scale = "class", layer = "cover") |>
  # the smallest perimeter a class of this many cells can have
  msr(
    equation = "4 * floor(sqrt(area.class)) +
      2 * (area.class > floor(sqrt(area.class))^2) +
      2 * (area.class - floor(sqrt(area.class))^2 > floor(sqrt(area.class)))",
    label = "minperim.class", layer = "cover"
  ) |>
  msr(
    equation = "perimeter.class / minperim.class",
    label = "lsi.class", layer = "cover"
  )',
CORE = paste0(HAB, '\n', CORE, '
  msr(
    equation = "mean(core.cell)",
    label = "coreshare.patch", layer = "habitat"
  ) |>
  msr_area(scale = "patch", unit = "map", layer = "habitat") |>
  msr(
    equation = "coreshare.patch * area.patch / 10000",
    label = "core_ha.patch", layer = "habitat"
  )'),
NCORE = paste0(HAB, '\n', CORE, '
  mdf_componentise(connectivity = 8L, layer = "core",
                   add = "coreid") |>
  msr(
    equation = "length(unique(coreid.cell[!is.na(coreid.cell)]))",
    label = "ncore.patch", layer = "habitat"
  )'),
CAI = paste0(HAB, '\n', CORE, '
  msr(
    equation = "mean(core.cell) * 100",
    label = "cai.patch", layer = "habitat"
  )'),
TCA = paste0(HAB, '\n', CORE, '
  msr_area(scale = "class", unit = "map", layer = "core") |>
  msr(
    equation = "area.class / 10000",
    label = "tca.class", layer = "core"
  )'),
CPLAND = paste0(HAB, '\n', CORE, '
  msr_area(scale = "class", layer = "core") |>
  msr_area(scale = "landscape", layer = "habitat") |>
  msr(
    equation = "area.class / area.landscape_habitat * 100",
    label = "cpland.class", layer = "core"
  )'),
CWED = '# my_contrast: how different each pair of classes is, from 0 to 1
mosaik(...) |>
  msr_dissimilarity(contrast = my_contrast, scale = "class", layer = "cover") |>
  msr_area(scale = "landscape", unit = "map", layer = "cover") |>
  msr(
    equation = "dissimilarity.class / area.landscape * 10000",
    label = "cwed.class", layer = "cover"
  )',
TECI = '# my_contrast: how different each pair of classes is, from 0 to 1
mosaik(...) |>
  msr_dissimilarity(contrast = my_contrast, scale = "class", layer = "cover") |>
  msr_perimeter(scale = "class", unit = "map", layer = "cover") |>
  msr(
    equation = "dissimilarity.class / perimeter.class * 100",
    label = "teci.class", layer = "cover"
  )',
PLADJ = 'mosaik(...) |>
  msr_adjacency(type = "paired", layer = "cover") |>
  msr(
    equation = "sum(diag(adjacency.class)) / sum(adjacency.class) * 100",
    label = "pladj.landscape", layer = "cover"
  )',
AI = 'mosaik(...) |>
  msr_adjacency(type = "like", layer = "cover") |>
  msr_area(scale = "class", layer = "cover") |>
  # the most like adjacencies a class of this many cells can have
  msr(
    equation = "2 * floor(sqrt(area.class)) * (floor(sqrt(area.class)) - 1) +
      (area.class > floor(sqrt(area.class))^2) *
        (2 * (area.class - floor(sqrt(area.class))^2) - 1) -
      (area.class - floor(sqrt(area.class))^2 > floor(sqrt(area.class)))",
    label = "maxadj.class", layer = "cover"
  ) |>
  msr(
    equation = "likeAdj.class / 2 / maxadj.class * 100",
    label = "ai.class", layer = "cover"
  )',
CONTAG = 'mosaik(...) |>
  msr_adjacency(type = "paired", layer = "cover") |>
  msr_area(scale = "class", layer = "cover") |>
  msr_area(scale = "landscape", layer = "cover") |>
  msr_number(scale = "landscape", layer = "cover") |>
  msr(
    equation = "(1 + sum(
        area.class / area.landscape * adjacency.class / rowSums(adjacency.class) *
        log(area.class / area.landscape * adjacency.class / rowSums(adjacency.class)),
        na.rm = TRUE) / (2 * log(number.landscape))) * 100",
    label = "contag.landscape", layer = "cover"
  )',
IJI = 'mosaik(...) |>
  msr_adjacency(type = "paired", layer = "cover") |>
  msr_number(scale = "landscape", layer = "cover") |>
  msr(
    equation = "-sum(
        adjacency.class[upper.tri(adjacency.class)] /
          sum(adjacency.class[upper.tri(adjacency.class)]) *
        log(adjacency.class[upper.tri(adjacency.class)] /
          sum(adjacency.class[upper.tri(adjacency.class)])),
        na.rm = TRUE) / log(number.landscape * (number.landscape - 1) / 2) * 100",
    label = "iji.landscape", layer = "cover"
  )',
NP = paste0(HAB, '
  msr_number(scale = "class", layer = "habitat") |>
  msr(
    equation = "number.class",
    label = "np.class", layer = "habitat"
  )'),
PD = paste0(HAB, '
  msr_number(scale = "class", layer = "habitat") |>
  msr_area(scale = "landscape", unit = "map", layer = "habitat") |>
  msr(
    equation = "number.class / area.landscape * 10000 * 100",
    label = "pd.class", layer = "habitat"
  )'),
DIVISION = paste0(HAB, '
  msr_area(scale = "patch", layer = "habitat") |>
  msr_area(scale = "landscape", layer = "habitat") |>
  msr(
    equation = "1 - sum((area.patch / area.landscape)^2, na.rm = TRUE)",
    label = "division.landscape", layer = "habitat"
  )'),
SPLIT = paste0(HAB, '
  msr_area(scale = "patch", layer = "habitat") |>
  msr_area(scale = "landscape", layer = "habitat") |>
  msr(
    equation = "area.landscape^2 / sum(area.patch^2, na.rm = TRUE)",
    label = "split.landscape", layer = "habitat"
  )'),
MESH = paste0(HAB, '
  msr_area(scale = "patch", unit = "map", layer = "habitat") |>
  msr_area(scale = "landscape", unit = "map", layer = "habitat") |>
  msr(
    equation = "sum(area.patch^2, na.rm = TRUE) / area.landscape / 10000",
    label = "mesh.landscape", layer = "habitat"
  )'),
COHESION = paste0(HAB, '
  msr_perimeter(scale = "patch", layer = "habitat") |>
  msr_area(scale = "patch", layer = "habitat") |>
  msr_area(scale = "landscape", layer = "habitat") |>
  msr(
    equation = "(1 - sum(perimeter.patch, na.rm = TRUE) /
      sum(perimeter.patch * sqrt(area.patch), na.rm = TRUE)) /
      (1 - 1 / sqrt(area.landscape)) * 100",
    label = "cohesion.landscape", layer = "habitat"
  )'),
ENN = paste0(HAB, '
  msr_distance(routing = "straight", layer = "habitat") |>
  msr(
    equation = "min(distance.patch)",
    label = "enn.patch", layer = "habitat"
  )'),
PROX = paste0(HAB, '
  msr_area(scale = "patch", unit = "map", layer = "habitat") |>
  msr_distance(routing = "straight", layer = "habitat") |>
  msr(
    # 10: the search radius, in metres
    equation = "sum((area.patch / distance.patch^2)[distance.patch <= 10],
      na.rm = TRUE)",
    label = "prox.patch", layer = "habitat"
  )'),
CONNECT = paste0(HAB, '
  msr_distance(routing = "straight", layer = "habitat") |>
  msr_number(scale = "class", layer = "habitat") |>
  # 10: the distance at which two patches count as joined, in metres
  msr(
    equation = "sum(distance.patch <= 10)",
    label = "joined.patch", layer = "habitat"
  ) |>
  msr(
    equation = "sum(joined.patch, na.rm = TRUE) /
      sum(number.class * (number.class - 1)) * 100",
    label = "connect.landscape", layer = "habitat"
  )'),
PR = 'mosaik(...) |>
  msr_number(scale = "landscape", layer = "cover") |>
  msr(
    equation = "number.landscape",
    label = "pr.landscape", layer = "cover"
  )',
PRD = 'mosaik(...) |>
  msr_number(scale = "landscape", layer = "cover") |>
  msr_area(scale = "landscape", unit = "map", layer = "cover") |>
  msr(
    equation = "number.landscape / area.landscape * 10000 * 100",
    label = "prd.landscape", layer = "cover"
  )',
RPR = 'mosaik(...) |>
  msr_number(scale = "landscape", layer = "cover") |>
  msr(
    # 12: the number of classes the map could have
    equation = "number.landscape / 12 * 100",
    label = "rpr.landscape", layer = "cover"
  )',
SHDI = 'mosaik(...) |>
  msr_area(scale = "class", layer = "cover") |>
  msr_area(scale = "landscape", layer = "cover") |>
  msr(
    equation = "-sum(area.class / area.landscape *
      log(area.class / area.landscape))",
    label = "shdi.landscape", layer = "cover"
  )',
SIDI = 'mosaik(...) |>
  msr_area(scale = "class", layer = "cover") |>
  msr_area(scale = "landscape", layer = "cover") |>
  msr(
    equation = "1 - sum((area.class / area.landscape)^2)",
    label = "sidi.landscape", layer = "cover"
  )',
MSIDI = 'mosaik(...) |>
  msr_area(scale = "class", layer = "cover") |>
  msr_area(scale = "landscape", layer = "cover") |>
  msr(
    equation = "-log(sum((area.class / area.landscape)^2))",
    label = "msidi.landscape", layer = "cover"
  )',
SHEI = 'mosaik(...) |>
  msr_area(scale = "class", layer = "cover") |>
  msr_area(scale = "landscape", layer = "cover") |>
  msr_number(scale = "landscape", layer = "cover") |>
  msr(
    equation = "-sum(area.class / area.landscape *
      log(area.class / area.landscape)) / log(number.landscape)",
    label = "shei.landscape", layer = "cover"
  )',
SIEI = 'mosaik(...) |>
  msr_area(scale = "class", layer = "cover") |>
  msr_area(scale = "landscape", layer = "cover") |>
  msr_number(scale = "landscape", layer = "cover") |>
  msr(
    equation = "(1 - sum((area.class / area.landscape)^2)) /
      (1 - 1 / number.landscape)",
    label = "siei.landscape", layer = "cover"
  )',
MSIEI = 'mosaik(...) |>
  msr_area(scale = "class", layer = "cover") |>
  msr_area(scale = "landscape", layer = "cover") |>
  msr_number(scale = "landscape", layer = "cover") |>
  msr(
    equation = "-log(sum((area.class / area.landscape)^2)) /
      log(number.landscape)",
    label = "msiei.landscape", layer = "cover"
  )'
)

# PROX needs the areas of the other patches, which no variable holds under
# the group rule; to do (see handoff)
codes$PROX <- NA_character_
