# Prototype of the groups-and-focus design (inst/design/groups_and_focus.md):
# primitives per class of a layer, variables <name>.<self|others|all>[_layer],
# the label's level sets the focus. Compared with today's code.
suppressMessages(devtools::load_all("D:/rdev/mosaik", quiet = TRUE))
suppressMessages(library(landscapemetrics))
landscape <- mosaik::landscape
source("inst/design/groups_prototype_today.R")   # today's metric code, as in the app

# ---- the layers ----------------------------------------------------------
d <- msk_dims(landscape)
m <- landscape |>
  mdf_filter(cover == 47, add = "habitat") |>
  mdf_componentise(connectivity = 8L, layer = "habitat", add = "patch") |>
  mdf_distance(source = "background", method = "manhattan", layer = "habitat",
               add = "depth") |>
  mdf_filter(depth >= 2, add = "core") |>
  mdf_componentise(connectivity = 8L, layer = "core", add = "coreid") |>
  msk_add(x = rep(seq_len(d[1]), d[2]), y = rep(rev(seq_len(d[2])), each = d[1]))
L <- function(nm) msk_pull(m, nm)
group_layers <- c("patch")          # layers whose classes are groups

# ---- primitives per class of a layer --------------------------------------
cls <- function(v) sort(unique(v[!is.na(v)]))
prim <- function(layer) {
  v <- L(layer)
  if (layer %in% group_layers) v[is.na(v)] <- 0     # outside the groups
  g <- setdiff(cls(v), if (layer %in% group_layers) 0)
  vv <- as.numeric(v)
  cnt <- countCellValuesCpp(vals = vv, nrow = d[2], ncol = d[1])
  ed  <- countCellEdgesCpp(vals = vv, nrow = d[2], ncol = d[1])
  adj <- countCellAdjacenciesCpp(vals = vv, nrow = d[2], ncol = d[1],
                                 doublecount = TRUE, eightconn = FALSE)
  allv <- cls(vv); dimnames(adj) <- list(allv, allv)
  out <- list(gid = g,
              area = cnt$cells[match(g, cnt$value)],
              perimeter = (ed$edgesX + ed$edgesY)[match(g, ed$value)],
              adjacency = adj[as.character(g), as.character(g), drop = FALSE])
  # straight distance between every two classes (nearest cells)
  D <- matrix(Inf, length(g), length(g), dimnames = list(g, g))
  if (layer %in% group_layers) for (i in seq_along(g)) {
    src <- as.numeric(v == g[i])
    dd <- sqrt(distanceCpp(vals = src, nrow = d[2], ncol = d[1], method = "euclidean"))
    for (j in seq_along(g)) if (i != j) D[i, j] <- min(dd[v == g[j]])
  }
  out$distance <- D
  # classes the map border cuts
  grid <- matrix(v, nrow = d[2], ncol = d[1], byrow = TRUE)
  out$border <- g %in% c(grid[1, ], grid[d[2], ], grid[, 1], grid[, d[1]])
  out
}
store <- list(cover = prim("cover"), patch = prim("patch"),
              habitat = prim("habitat"))
cmat <- matrix(0.5, length(store$cover$gid), length(store$cover$gid)); diag(cmat) <- 0
dimnames(cmat) <- list(store$cover$gid, store$cover$gid)
dis <- msr_dissimilarity(landscape, contrast = cmat, scale = "class", layer = "cover")
store$cover$dissimilarity <- msk_categories(dis, "cover")$dissimilarity

# ---- the evaluator ----------------------------------------------------------
proto_msr <- function(equation, label, layer, exclude_cut = TRUE) {
  lab <- regmatches(label, regexec("^([^.]+)\\.(class|landscape)$", label))[[1]]
  e <- parse(text = equation)
  S <- store[[layer]]
  pop <- seq_along(S$gid)
  if (exclude_cut && layer %in% group_layers) pop <- pop[!S$border]
  lv <- L(layer)
  resolve <- function(focus_k) {
    env <- list()
    for (v in all.vars(e)) {
      if (exists(v, envir = baseenv(), inherits = FALSE)) next
      p <- regmatches(v, regexec("^([A-Za-z][A-Za-z0-9]*)\\.(self|others|all)(_(.+))?$", v))[[1]]
      if (!length(p)) stop("bad variable ", v)
      nm <- p[2]; fo <- p[3]; ly <- if (nzchar(p[5])) p[5] else layer
      if (is.null(focus_k) && fo != "all") stop(v, ": only .all at landscape level")
      T <- store[[ly]]
      if (!is.null(T[[nm]])) {                      # a stored value
        val <- T[[nm]]
        p_idx <- if (ly == layer) pop else seq_along(T$gid)
        k <- if (!is.null(focus_k)) match(S$gid[focus_k], T$gid)
        if (is.matrix(val)) {
          env[[v]] <- if (is.null(focus_k)) val[p_idx, p_idx, drop = FALSE] else
            switch(fo, self = val[k, k], others = val[k, setdiff(p_idx, k)],
                   all = val[k, p_idx])
        } else {
          env[[v]] <- if (is.null(focus_k)) val[p_idx] else
            switch(fo, self = val[k], others = val[setdiff(p_idx, k)],
                   all = val[p_idx])
        }
      } else {                                      # a layer: cell values
        cv <- L(nm)
        inpop <- !is.na(lv) & lv %in% S$gid[pop]
        env[[v]] <- if (is.null(focus_k)) cv[inpop] else
          switch(fo, self = cv[!is.na(lv) & lv == S$gid[focus_k]],
                 others = cv[inpop & lv != S$gid[focus_k]], all = cv[inpop])
      }
    }
    env
  }
  if (lab[3] == "landscape") {
    res <- eval(e, resolve(NULL))
  } else {
    res <- rep(NA_real_, length(S$gid))
    for (k in pop) res[k] <- eval(e, resolve(k))
  }
  store[[layer]][[lab[2]]] <<- res
  res
}

# ---- the metrics in the new notation ---------------------------------------
minp <- function(a) "4 * floor(sqrt(A)) + 2 * (A > floor(sqrt(A))^2) + 2 * (A - floor(sqrt(A))^2 > floor(sqrt(A)))"
mp <- function(A) gsub("A", A, minp())
maxadj <- "2 * floor(sqrt(A)) * (floor(sqrt(A)) - 1) + (A > floor(sqrt(A))^2) * (2 * (A - floor(sqrt(A))^2) - 1) - (A - floor(sqrt(A))^2 > floor(sqrt(A)))"
new <- list(
  AREA   = list("area.self / 10000", "areaha.class", "patch"),
  PERIM  = list("perimeter.self", "perim.class", "patch"),
  GYRATE = list("mean(sqrt((x.self - mean(x.self))^2 + (y.self - mean(y.self))^2))", "gyrate.class", "patch"),
  CA     = list("area.self / 10000", "ca.class", "cover"),
  PLAND  = list("area.self / sum(area.all) * 100", "pland.class", "cover"),
  LPI    = list("max(area.all) / sum(area.all_cover) * 100", "lpi.landscape", "patch"),
  TE     = list("perimeter.self", "te.class", "cover"),
  ED     = list("perimeter.self / sum(area.all) * 10000", "ed.class", "cover"),
  PARA   = list("perimeter.self / area.self", "para.class", "patch"),
  SHAPE  = list(paste0("perimeter.self / (", mp("area.self"), ")"), "shape.class", "patch"),
  FRAC   = list("ifelse(area.self == 1, 1, 2 * log(0.25 * perimeter.self) / log(area.self))", "frac.class", "patch"),
  PAFRAC = list("2 / (cov(log(area.all), log(perimeter.all)) / var(log(perimeter.all)))", "pafrac.landscape", "patch"),
  LSI    = list(paste0("perimeter.self / (", mp("area.self"), ")"), "lsi.class", "cover"),
  CORE   = list("sum(core.self) / 10000", "coreha.class", "patch"),
  NCORE  = list("length(unique(na.omit(coreid.self)))", "ncore.class", "patch"),
  CAI    = list("mean(core.self) * 100", "cai.class", "patch"),
  TCA    = list("sum(core.self) / 10000", "tca.class", "habitat"),
  CPLAND = list("sum(core.self) / sum(area.all) * 100", "cpland.class", "habitat"),
  CWED   = list("dissimilarity.self / sum(area.all) * 10000", "cwed.class", "cover"),
  TECI   = list("dissimilarity.self / perimeter.self * 100", "teci.class", "cover"),
  PLADJ  = list("sum(diag(adjacency.all)) / sum(adjacency.all) * 100", "pladj.landscape", "cover"),
  AI     = list(paste0("adjacency.self / 2 / (", gsub("A", "area.self", maxadj), ") * 100"), "ai.class", "cover"),
  CONTAG = list("(1 + sum(area.all / sum(area.all) * adjacency.all / rowSums(adjacency.all) * log(area.all / sum(area.all) * adjacency.all / rowSums(adjacency.all)), na.rm = TRUE) / (2 * log(length(area.all)))) * 100", "contag.landscape", "cover"),
  IJI    = list("-sum(adjacency.all[upper.tri(adjacency.all)] / sum(adjacency.all[upper.tri(adjacency.all)]) * log(adjacency.all[upper.tri(adjacency.all)] / sum(adjacency.all[upper.tri(adjacency.all)])), na.rm = TRUE) / log(length(area.all) * (length(area.all) - 1) / 2) * 100", "iji.landscape", "cover"),
  NP     = list("length(area.all)", "np.landscape", "patch"),
  PD     = list("length(area.all) / sum(area.all_cover) * 10000 * 100", "pd.landscape", "patch"),
  DIVISION = list("1 - sum((area.all / sum(area.all_cover))^2)", "division.landscape", "patch"),
  SPLIT  = list("sum(area.all_cover)^2 / sum(area.all^2)", "split.landscape", "patch"),
  MESH   = list("sum(area.all^2) / sum(area.all_cover) / 10000", "mesh.landscape", "patch"),
  COHESION = list("(1 - sum(perimeter.all) / sum(perimeter.all * sqrt(area.all))) / (1 - 1 / sqrt(sum(area.all_cover))) * 100", "cohesion.landscape", "patch"),
  ENN    = list("min(distance.others)", "enn.class", "patch"),
  PROX   = list("sum(area.others / distance.others^2 * (distance.others <= 10))", "prox.class", "patch"),
  CONNECT = list("sum(distance.others <= 10)", "joined.class", "patch"),
  PR     = list("length(area.all)", "pr.landscape", "cover"),
  PRD    = list("length(area.all) / sum(area.all) * 10000 * 100", "prd.landscape", "cover"),
  RPR    = list("length(area.all) / 12 * 100", "rpr.landscape", "cover"),
  SHDI   = list("-sum(area.all / sum(area.all) * log(area.all / sum(area.all)))", "shdi.landscape", "cover"),
  SIDI   = list("1 - sum((area.all / sum(area.all))^2)", "sidi.landscape", "cover"),
  MSIDI  = list("-log(sum((area.all / sum(area.all))^2))", "msidi.landscape", "cover"),
  SHEI   = list("-sum(area.all / sum(area.all) * log(area.all / sum(area.all))) / log(length(area.all))", "shei.landscape", "cover"),
  SIEI   = list("(1 - sum((area.all / sum(area.all))^2)) / (1 - 1 / length(area.all))", "siei.landscape", "cover"),
  MSIEI  = list("-log(sum((area.all / sum(area.all))^2)) / log(length(area.all))", "msiei.landscape", "cover")
)

# ---- today's values ------------------------------------------------------------
habitat <- 47; my_contrast <- cmat
today <- function(id) {
  if (is.na(codes[[id]])) return(NULL)
  code <- gsub("mosaik(...)", "landscape", codes[[id]], fixed = TRUE)
  mm <- eval(parse(text = code))
  lab <- sub('.*label = "([^".]+)\\.([a-z]+)", layer = "([^"]+)"\\s*\\)\\s*$', "\\1 \\2 \\3", code)
  p <- strsplit(lab, " ")[[1]]
  tbl <- switch(p[2], patch = mm@patches[[p[3]]], class = mm@categories[[p[3]]],
                landscape = mm@global[[p[3]]])
  list(v = tbl[[p[1]]], gid = tbl$gid, patch = tbl$patch, layer = p[3], scale = p[2])
}
old_prox <- c(6.69517, 6.52466, 48.4707, 106.114, 57.8438, 24.5556, 3.10769, 15.552,
              3.62295, NA, NA, 27.1993, 54.2628, 7.50994, 52.4659, 17.7933, NA)

cat(sprintf("%-9s %-8s %s\n", "metric", "result", "new | today"))
for (id in names(new)) {
  q <- new[[id]]
  r <- tryCatch(proto_msr(q[[1]], q[[2]], q[[3]]), error = function(e) e)
  if (inherits(r, "error")) { cat(sprintf("%-9s ERROR    %s\n", id, conditionMessage(r))); next }
  if (id == "CONNECT") {
    r <- proto_msr("sum(joined.all) / (length(joined.all) * (length(joined.all) - 1)) * 100",
                   "connect.landscape", "patch")
  }
  t <- if (id == "PROX") list(v = old_prox) else today(id)
  if (id == "CONNECT") t <- today("CONNECT")
  # pick the comparable values
  if (!is.null(t$scale) && t$scale == "class" && q[[3]] %in% c("cover", "habitat")) {
    tv <- t$v[match(store[[q[[3]]]]$gid, t$gid)]
  } else tv <- t$v
  if (q[[3]] %in% c("cover", "habitat") && grepl("\\.class$", q[[2]]) && id %in% c("TCA", "CPLAND")) {
    r <- r[store$habitat$gid == 1]; tv <- t$v[t$gid == 1]
  }
  ok <- length(r) == length(tv) && isTRUE(all.equal(as.numeric(r), as.numeric(tv), tolerance = 1e-6))
  cat(sprintf("%-9s %-8s %s | %s\n", id, if (ok) "SAME" else "DIFF",
              paste(signif(head(r, 6), 5), collapse = " "),
              paste(signif(head(tv, 6), 5), collapse = " ")))
}

cat("\nENN  new:  ", round(store$patch$enn, 3), "\n")
cat("ENN  today:", round(today("ENN")$v, 3), "\n")
cat("PROX new:  ", round(store$patch$prox, 3), "\n")
cat("PROX old:  ", round(old_prox, 3), "\n")
cat("cut patches:", store$patch$gid[store$patch$border], "\n")
