# Build the example dataset `landscape`: a 60 x 56 map whose structures are
# placed deliberately, so each recipe in the cookbooks finds something known.
#
# Run from the package root:
#   Rscript --vanilla data-raw/landscape.R          # preview only
#   Rscript --vanilla data-raw/landscape.R save     # also write data/landscape.rda
#
# Coordinates are cell indices, x = 1..60 from the left, y = 1..56 from the
# bottom; only the finished map gets its cell size of 10 m. The class IDs follow the earlier dataset:
#
#   1 river             11 arable land        21 intensive grassland
#  24 extensive grassland 27 fallow, shrub and clear-cuts
#  31 settlement        35 road               41 wetland
#  44 orchard           47 forest and hedgerows
#
# `canopy` is the height of the vegetation in metres. It follows the land cover
# only in part: within the forest it varies with the age of the stand, falls
# towards old edges and drops in small gaps, and outside the forest it shows
# hedgerows and single trees that the land cover does not.

pkgload::load_all(".", quiet = TRUE)
nx <- 60; ny <- 56

cover <- matrix(NA_integer_, nrow = ny, ncol = nx)   # [row, col], row 1 = top
at <- function(x, y) cbind(ny + 1 - y, x)            # (x, y) -> [row, col]
paint <- function(m, x, y, value) {
  cells <- expand.grid(x = x, y = y)
  m[at(cells$x, cells$y)] <- value
  m
}
field <- function(seed, hurst = 0.7) {
  f <- mosaik(extent = c(0, nx, 0, ny), res = 1) |>
    drw_texture(type = "fbm", hurst = hurst, seed = seed) |>
    mdf_scale(layer = "values", range = c(-1, 1))
  matrix(msk_pull(f, "values"), nrow = ny, ncol = nx, byrow = TRUE)
}
wobble <- field(seed = 5)                            # roughens outlines
# an ellipse whose outline is pushed in and out by `wobble`
blob <- function(m, cx, cy, rx, ry, value, rough = 0.25) {
  g <- expand.grid(x = seq_len(nx), y = seq_len(ny))
  r <- sqrt(((g$x - cx) / rx)^2 + ((g$y - cy) / ry)^2)
  hit <- r <= 1 + rough * wobble[at(g$x, g$y)]
  m[at(g$x[hit], g$y[hit])] <- value
  m
}
# a band of `width` cells from (x0, y0) to (x1, y1) that meanders sideways
path <- function(m, x0, y0, x1, y1, value, width = 2, bend = 1.5) {
  n <- max(abs(x1 - x0), abs(y1 - y0)) + 1
  s <- seq(0, 1, length.out = n)
  xs <- x0 + (x1 - x0) * s
  ys <- y0 + (y1 - y0) * s
  off <- bend * sin(s * 2 * pi)
  horizontal <- abs(x1 - x0) >= abs(y1 - y0)
  px <- round(if (horizontal) xs else xs + off)
  py <- round(if (horizontal) ys + off else ys)
  for (i in seq_len(n)) {
    if (horizontal) {
      m <- paint(m, px[i], py[i] + 0:(width - 1), value)
    } else {
      m <- paint(m, px[i] + 0:(width - 1), py[i], value)
    }
    # a step in both directions would leave two cells touching only at a
    # corner; fill the corner so the path stays connected with four neighbours
    if (i > 1 && px[i] != px[i - 1] && py[i] != py[i - 1]) {
      m <- paint(m, px[i], py[i - 1] + if (horizontal) 0:(width - 1) else 0, value)
    }
  }
  m
}

# open land: fields of six classes ---------------------------------------------
base <- mosaik(extent = c(0, nx, 0, ny), res = 1) |>
  drw_tessellation(type = "voronoi", n = 14, seed = 3)
tiles <- matrix(msk_pull(base, "values"), nrow = ny, ncol = nx, byrow = TRUE)
open <- c(11, 21, 24, 27, 31, 44)
cover[] <- open[(tiles + 2) %% length(open) + 1]
# B lies in one large meadow, so it borders few classes, unlike A and C
cover <- blob(cover, 48, 39, 15, 17, 21L, rough = 0.05)

# a wetland (41) with two irregular islands of open land -----------------------
cover <- blob(cover, 50, 11, 7, 6, 41L, rough = 0.2)
cover[at(c(47, 48, 48, 49, 47), c(8, 8, 9, 9, 7))] <- 24L
cover[at(c(52, 53, 53, 54, 53, 52), c(12, 12, 13, 13, 14, 11))] <- 24L

# bocage in the south-west: small irregular fields of arable land and
# grassland, bounded by hedgerows that follow the field boundaries -------------
small <- mosaik(extent = c(0, nx, 0, ny), res = 1) |>
  drw_tessellation(type = "voronoi", n = 90, seed = 11)
fields <- matrix(msk_pull(small, "values"), nrow = ny, ncol = nx, byrow = TRUE)
bocage <- blob(matrix(0L, ny, nx), 8, 14, 8, 13, 1L, rough = 0.2) == 1
cover[bocage] <- c(11L, 21L, 24L)[fields[bocage] %% 3 + 1]
land <- cover                                        # what lies under the forest

# forest (47) -------------------------------------------------------------------
# A: large patch in the north-west, with a fresh clear-cut notch from the
# north, two small indentations in the south, an old irregular glade and an
# irregular clearing made when the road was built
cover <- blob(cover, 13, 40, 9, 10, 47L)
cover <- blob(cover, 12, 29, 6, 4, 47L, rough = 0.15) # south lobe, core of its own
cover <- paint(cover, 8, 25:27, NA)                  # small indentations
cover <- paint(cover, 17:18, 25:26, NA)
clearing <- matrix(FALSE, ny, nx)                    # beside the road, fresh
clearing[at(c(7, 8, 8, 9, 7, 8), c(36, 36, 37, 37, 35, 35))] <- TRUE
glade <- matrix(FALSE, ny, nx)                       # an old windthrow
glade[at(c(17, 18, 18, 19, 19, 20, 18, 17, 19), c(39, 39, 40, 40, 41, 41, 41, 38, 39))] <- TRUE
cover[glade] <- 27L                                  # fallow and shrub
# B: patch in the north-east with an old square glade, a dead-end branch to
# the south and an irregular loop to the east
cover <- blob(cover, 47, 42, 7, 7, 47L)
cover <- paint(cover, 48:49, 43:44, 24L)
cover <- path(cover, 46, 36, 46, 26, 47L, bend = 1)  # branch
cover <- path(cover, 52, 47, 57, 47, 47L, width = 1, bend = 0.7)   # loop
cover <- path(cover, 57, 47, 57, 37, 47L, width = 1, bend = 1)
cover <- path(cover, 57, 37, 52, 37, 47L, width = 1, bend = 0.7)
# C: core in the south-west, joined to A by a corridor
cover <- blob(cover, 20, 14, 6, 5, 47L, rough = 0.15)
cover <- path(cover, 14, 31, 19, 18, 47L)
# hedgerows along the field boundaries of the bocage. The bocage is decaying:
# only some field boundaries still carry a hedge, and those that do have gaps.
# No hedge touches core C, so none merges into it
right <- cbind(fields[, -1], fields[, nx])
up <- rbind(fields[1, ], fields[-ny, ])
other <- ifelse(fields != right, right, up)
edge_cell <- bocage & (fields != right | fields != up)
pair <- paste(pmin(fields, other), pmax(fields, other))
set.seed(4)
pairs <- unique(pair[edge_cell])
standing <- sample(pairs, round(length(pairs) * 0.4))
hedge <- edge_cell & pair %in% standing & runif(length(fields)) > 0.08
near_c <- blob(matrix(0L, ny, nx), 20, 14, 8, 7, 1L, rough = 0.15) == 1
cover[hedge & !near_c & cover != 47] <- 47L
# fresh clear-cuts: the notch in A and the clearing beside the road. They are
# cut out after the forest edge was measured (see canopy), so they have no
# mantle yet
fresh <- clearing
fresh[at(rep(10:15, each = 10), rep(44:53, 6))] <- TRUE
fresh <- fresh & cover == 47
cover[fresh] <- 27L
# corridor from A to B, interrupted where it crosses the river (drawn later)
cover <- path(cover, 21, 42, 41, 42, 47L, bend = 0.6)
# small patches: three islets, a small irregular patch, a thin strip without
# core, two patches touching only at a corner, and one patch at the grid edge
cover[at(c(30, 27, 42), c(51, 25, 30))] <- 47L
cover[at(c(52, 53, 53, 54, 52, 54), c(20, 20, 21, 21, 19, 22))] <- 47L
cover <- path(cover, 16, 5, 28, 5, 47L, bend = 0.8)
cover[at(c(22, 23, 21, 22, 23), c(51, 51, 52, 52, 52))] <- 47L
cover[at(c(24, 25, 24, 25, 26), c(53, 53, 54, 54, 54))] <- 47L
cover <- blob(cover, 55, 1, 5, 3, 47L, rough = 0.2)

# river (1) and a road (35) running beside it ----------------------------------
# the road follows the river from the southern edge for a third of the map,
# first one cell away from it, then two, then swings off in an S towards the
# north-west, crossing forest A on its way to the western edge
rx <- round(36 + 2 * sin(seq_len(ny) / 8))
for (y in seq_len(ny)) cover[at(rx[y], y)] <- 1L
road <- function(m, x, y) { m[at(x, y)] <- 35L; m }
for (y in 1:19) {
  gap <- if (y <= 9) 1 else 2
  # where the river shifts sideways, the road takes both columns, so it stays
  # connected with four neighbours too
  xs <- if (y < 19) range(rx[y], rx[y + 1]) else rx[y]
  cover <- road(cover, seq(xs[1], xs[length(xs)]) - gap - 1, y)
}
# the S: a cubic Bezier curve, sampled densely and joined with four neighbours
p <- rbind(c(rx[19] - 3, 19), c(36, 34), c(12, 26), c(1, 44))
t <- seq(0, 1, length.out = 400)
bz <- (1 - t)^3 %o% p[1, ] + 3 * (1 - t)^2 * t %o% p[2, ] +
  3 * (1 - t) * t^2 %o% p[3, ] + t^3 %o% p[4, ]
bx <- round(bz[, 1]); by <- round(bz[, 2])
keep <- c(TRUE, diff(bx) != 0 | diff(by) != 0)
bx <- bx[keep]; by <- by[keep]
for (i in seq_along(bx)) {
  cover <- road(cover, bx[i], by[i])
  if (i > 1 && bx[i] != bx[i - 1] && by[i] != by[i - 1]) {
    cover <- road(cover, bx[i], by[i - 1])
  }
}

# whatever was cut out of the forest returns to the land beneath it
cover[is.na(cover)] <- land[is.na(cover)]

# canopy height -----------------------------------------------------------------
# each class has a typical height, varying smoothly within it
height <- c("1" = 0, "11" = 0.3, "21" = 0.3, "24" = 0.5, "27" = 1.5,
            "31" = 1, "35" = 0, "41" = 1, "44" = 5, "47" = 25)
vary <- field(seed = 7, hurst = 0.8)
spread <- c("1" = 0, "11" = 0.2, "21" = 0.2, "24" = 0.4, "27" = 1.2,
            "31" = 1, "35" = 0, "41" = 0.8, "44" = 1, "47" = 3)
canopy <- height[as.character(cover)] + spread[as.character(cover)] * vary
# gardens and trees between the houses of the settlement
garden <- cover == 31 & field(seed = 8, hurst = 0.6) > 0.3
canopy[garden] <- 8 + 2 * vary[garden]
# a younger stand: C, the forest in the south-west
young <- blob(matrix(0L, ny, nx), 20, 14, 6, 5, 1L, rough = 0.15) == 1 & cover == 47
canopy[young] <- 13 + 2 * vary[young]
# old edges have grown a mantle, so the trees get lower towards the edge.
# The edge is measured on the forest as it was before the fresh clear-cuts and
# the road were cut through it, so the trees along them still stand at full
# height right up to the new edge: the edge is sharp
# The road cut through a forest that was closed before, so near the road the
# forest is measured as the closing of the forest covers it: A and its south
# lobe stand tall right up to the road, as one forest the road divides
before <- (cover == 47) | fresh | (cover == 35)
square <- msk_struct("square", width = 5)
edge_of <- mosaik(extent = c(0, nx, 0, ny), res = 1,
                  vals = list(f = as.vector(t(before) * 1),
                              road = as.vector(t(cover == 35) * 1))) |>
  mdf_dilate(struct = square, layer = "f", add = "closed") |>
  mdf_erode(struct = square, layer = "closed") |>
  mdf_distance(source = "foreground", layer = "road", add = "to_road") |>
  mdf_filter(f == 1 | (closed == 1 & to_road <= 3), add = "then") |>
  mdf_distance(source = "background", layer = "then", add = "d")
d <- matrix(msk_pull(edge_of, "d"), nrow = ny, ncol = nx, byrow = TRUE)
forest <- cover == 47
drop <- ifelse(forest, 1 - 0.6 * exp(-(d - 1) / 2), 0)
# averaged over the forest cells around each cell: where a strip two cells
# wide steps sideways it is three wide for one cell, and that cell would
# otherwise stand out much higher than the strip around it
pad <- function(m) rbind(0, cbind(0, m, 0), 0)
sum3 <- function(m) {
  p <- pad(m)
  Reduce(`+`, lapply(0:2, function(i) lapply(0:2, function(j)
    p[1:ny + i, 1:nx + j])) |> unlist(recursive = FALSE))
}
drop <- sum3(drop) / pmax(sum3(forest * 1), 1)
# but the outermost row keeps the height of the mantle: averaging would make
# a cell in an inner corner taller and one on an outer corner lower than the
# rest of the edge
rim <- forest & d < 1.5
drop[rim] <- 0.4
# where the corridor from A joins C, C has no edge: its trees run at their
# height into the corridor. For C the edge is measured with the corridor
# widened, so C does not drop in front of it
corridor <- path(matrix(0L, ny, nx), 14, 31, 19, 18, 1L) == 1
joined <- mosaik(extent = c(0, nx, 0, ny), res = 1,
                 vals = list(then = msk_pull(edge_of, "then"),
                             corridor = as.vector(t(corridor) * 1))) |>
  mdf_dilate(struct = square, layer = "corridor", add = "wide") |>
  mdf_filter(then == 1 | wide == 1, add = "joined") |>
  mdf_distance(source = "background", layer = "joined", add = "d")
dj <- matrix(msk_pull(joined, "d"), nrow = ny, ncol = nx, byrow = TRUE)
drop[young] <- pmax(0.4, (1 - 0.6 * exp(-(floor(dj) - 1) / 2)))[young]
canopy[forest] <- canopy[forest] * drop[forest]
# hedgerows
hedgerow <- hedge & !near_c & cover == 47
canopy[hedgerow] <- 7 + vary[hedgerow]
# the old windthrow in A: lying stems, a few survivors and regrowth, uneven
canopy[glade] <- 3.5 + 2.5 * matrix(field(seed = 9, hurst = 0.3), ny, nx)[glade]
# single trees in the grassland and fields, away from the forest
set.seed(6)
open_land <- which(cover %in% c(11, 21, 24) & !bocage &
                     matrix(msk_pull(mosaik(extent = c(0, nx, 0, ny), res = 1,
                                            vals = list(f = as.vector(t(forest) * 1))) |>
                                       mdf_distance(source = "foreground", layer = "f", add = "o"),
                                     "o"), nrow = ny, ncol = nx, byrow = TRUE) > 3)
trees <- sample(open_land, 12)
canopy[trees] <- 12 + 2 * vary[trees]
# just harvested, so bare
cut <- fresh & cover != 35
canopy[cut] <- 0.2
canopy <- pmax(round(canopy, 1), 0)

# the map is built in cells; each cell is 10 m wide, so the map is 600 x 560 m
landscape <- mosaik(extent = c(0, nx * 10, 0, ny * 10), res = 10,
                    vals = list(cover = as.vector(t(cover)),
                                canopy = as.vector(t(canopy))))

# preview and the counts the cookbooks rely on
check <- landscape |>
  mdf_filter(cover == 47, add = "forest") |>
  mdf_fill(layer = "forest", add = "filled")
cols <- c("1" = "#2b6cb0", "11" = "#f3e9c6", "21" = "#e3d39a", "24" = "#cdbb7a",
          "27" = "#c8a27a", "31" = "#9a8f86", "35" = "#3a3a3a", "41" = "#9ac27c",
          "44" = "#c9d9a0", "47" = "#1f5f2e")
png("inst/design/cookbook/landscape_preview.png", width = 1200, height = 560)
msk_vis(check, .layer("cover", colours = cols), .layer("canopy"))
dev.off()
print(table(msk_pull(landscape, "cover")))
fg <- ifelse(msk_pull(check, "forest") == 1, 1, NA)
for (conn in c(4L, 8L)) {
  cc <- msk_pull(mdf_componentise(msk_add(check, fg = fg), connectivity = conn,
                                  layer = "fg", add = "cc"), "cc")
  cat("forest patches,", conn, "neighbours:", length(unique(na.omit(cc))),
      " sizes:", sort(table(cc), decreasing = TRUE), "\n")
}
cat("hole cells in forest:",
    sum(msk_pull(check, "filled") - msk_pull(check, "forest")), "\n")
# edge contrast: how many different classes each forest patch borders
cm <- matrix(msk_pull(landscape, "cover"), nrow = ny, ncol = nx, byrow = TRUE)
pm <- matrix(cc, nrow = ny, ncol = nx, byrow = TRUE)
for (id in sort(unique(na.omit(as.vector(pm))))) {
  cells <- which(pm == id, arr.ind = TRUE)
  nb <- unique(unlist(lapply(seq_len(nrow(cells)), function(i) {
    r <- cells[i, 1] + -1:1; c <- cells[i, 2] + -1:1
    as.vector(cm[r[r >= 1 & r <= ny], c[c >= 1 & c <= nx]])
  })))
  nb <- setdiff(nb, 47)
  cat(sprintf("patch %2d  %3d cells  borders %d classes: %s\n", id, nrow(cells),
              length(nb), paste(sort(nb), collapse = " ")))
}

if ("save" %in% commandArgs(trailingOnly = TRUE)) {
  # class labels, written directly: the constructor has no argument for them
  labels <- c("1" = "river", "11" = "arable land", "21" = "intensive grassland",
              "24" = "extensive grassland", "27" = "fallow, shrub and clear-cuts",
              "31" = "settlement", "35" = "road", "41" = "wetland",
              "44" = "orchard", "47" = "forest and hedgerows")
  gid <- sort(unique(msk_pull(landscape, "cover")))
  # and a colour per class, which msk_vis uses wherever the land cover is drawn
  landscape@categories$cover <- list(gid = gid, val = unname(labels[as.character(gid)]),
                                     colour = unname(cols[as.character(gid)]))
  methods::validObject(landscape)
  save(landscape, file = "data/landscape.rda", compress = "xz")
}
