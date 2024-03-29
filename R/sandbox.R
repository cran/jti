# path <- "/home/mads/Documents/phd/publications/unity_propagation/src/message"
# easypackages::libraries(
#   "tidyverse", # dplyr, ggplot2, readr, stringr, purrr, tidyr
#   "fs",
#   "glue"  
# ) 

# source(glue::glue("{path}/utils.R"))

# RNGkind("L'Ecuyer-CMRG")
# set.seed(300718)
# nc    <- 1L
# kfold <- 10L

# x      <- "chess"
# dat_   <- readRDS(glue::glue("{path}/../data/prc/{x}.rds"))
# dat    <- dat_[[1]]
# cls    <- dat_[[2]]
# g      <- readRDS(glue::glue("{path}/graphs/graph_chess.rds"))

# ne_max <- ncol(dat)-1L
# nobs   <- nrow(dat)
# perm_obs <- sample(1:nobs, nobs)

# cat("\n")

# for (ne in 29:ne_max) {
#   cv <- cv_jt(
#     data       = dat,
#     graph      = g,
#     class_var  = cls,
#     dat_name   = x,
#     perm       = perm_obs,
#     nevidence  = ne,
#     kfold      = kfold,
#     ncores     = nc
#   )
# }
  
# size_mb <- function(x) {
#   format(object.size(x), units = "Mb", standard = "auto", digits = 1L)
# }

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                 MUNIN: 1041-1397-80592
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# l    <- readRDS("../../../../sandbox/r/bns/munin.rds")
# cpts <- jti::bnfit_to_cpts(l)
# cl   <- jti::cpt_list(cpts)
# cp   <- compile(cl)
# .map_lgl(cp$charge$C, function(x) inherits(x, "sparta_unity")) |> sum()

# microbenchmark::microbenchmark(
#   jt(cp, propagate = "collect"),
#   times = 1
# )

# 35s (UP) vs 37s (OLD)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                         LINK: 724-1125-14211
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# l    <- readRDS("../../../../sandbox/r/bns/link.rds")
# cpts <- bnfit_to_cpts(l)
# cl   <- cpt_list(cpts)
# cp   <- compile(cl)

# microbenchmark::microbenchmark(
#   jt(cp, propagate = "collect"),
#   times = 1
# )

# set.seed(70)
# ev <- sapply(dim_names(cl), function(x) x[sample(1:length(x), 1)])[sample(1:700, 50)]



# .map_lgl(cp$charge$C, function(x) inherits(x, "sparta_unity")) |> sum()

# microbenchmark::microbenchmark(
#   jt(cp, propagate = "collect"),
#   jt(cp, propagate = "collect", unity_msg = FALSE),
#   times = 1
# )

# 42s (UP) vs 70s (OLD)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                        DIABETES: 413-602-429409
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# l    <- readRDS("../../../../sandbox/r/bns/diabetes.rds")
# cpts <- bnfit_to_cpts(l)
# cl   <- cpt_list(cpts)

# microbenchmark::microbenchmark(
#   triangulate(cl, mpd_based = TRUE),
#   triangulate(cl, mpd_based = FALSE),
#   times = 2
# )

# tt <- triangulate(cl, mpd_based = TRUE)
# tf <- triangulate(cl, mpd_based = FALSE)

# jt_nbinary_ops(tt) / jt_nbinary_ops(tf)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                      ANDES: 223-338-1157
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# l    <- readRDS("../../../../sandbox/r/bns/andes.rds")
# cpts <- bnfit_to_cpts(l)
# cl   <- cpt_list(cpts)
# cp   <- compile(cl)
# .map_lgl(cp$charge$C, function(x) inherits(x, "sparta_unity")) |> sum()

# microbenchmark::microbenchmark(
#   jt(cp, propagate = "collect"),
#   times = 10
# )

# 0.82s (UP) vs 0.94s (OLD)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                      MILDEW: 35-46-540150
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# l    <- readRDS("../../../../sandbox/r/bns/link.rds")
# cpts <- bnfit_to_cpts(l)
# cl   <- cpt_list(cpts)
# tt   <- triangulate(cl)

# p   <- "/home/mads/Documents/phd/publications/evidence_triangulation/src/tri_fail_ev_lai_1.rds"
# tri <- readRDS(p)
# ct  <- tri$junction_tree_collect + t(tri$junction_tree_collect)
# g   <- igraph::graph_from_adjacency_matrix(ct, "undirected")
# plot(g, vertex.size = 1)
# jt_nbinary_ops(tri)


# ev <- replicate(500, sample(names(cl), sample(1:200)))

# jt_nbinary_ops(tt, ev, nc = 3) |> sum()

# microbenchmark::microbenchmark(
#   jt_nbinary_ops(tt, ev, nc = 3),
#   jt_nbinary_ops(tt, ev, nc = 1),
#   times = 2
# )

# .map_lgl(cp$charge$C, function(x) inherits(x, "sparta_unity")) |> sum()

# microbenchmark::microbenchmark(
#   jt(cp, propagate = "collect"),
#   times = 15
# )

# 0.280 (UP) vs 0.432 (OLD)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                      HAILFINDER: 56-66-2656
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# l    <- readRDS("../../../../sandbox/r/bns/hailfinder.rds")
# cpts <- bnfit_to_cpts(l)
# cl   <- cpt_list(cpts)
# cp   <- jti::compile(cl)

# microbenchmark::microbenchmark(
#   jt(cp, propagate = "collect"),
#   times = 15
# )

# 0.027 (UP) vs 0.017 (OLD)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                      BARLEY: 48-84-114005
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# l <- readRDS(url("https://www.bnlearn.com/bnrepository/diabetes/diabetes.rds"))
# cpts <- bnfit_to_cpts(l)
# cl   <- cpt_list(cpts)

# mp  <- mpd(cl)
# fps <- unique(unlist(mp$primes_int[which(mp$flawed)]))
# fps <- colnames(mp$graph)[fps]
# nfps <- length(fps)

# set.seed(9)

# pe <- structure(runif(nfps), names = names(fps)) # rep(.1, nfps)

# evidence <- unique(lapply(1:1000, function(x) {
#   len <- sample(3:100, 1)
#   sample(fps, len, prob = pe)
# }))



# tris <- list(
#   triangulate(cl, tri = "min_fill", mpd_based = FALSE),
#   triangulate(cl, tri = "min_sp", mpd_based = FALSE),
#   triangulate(cl, tri = "min_ssp", mpd_based = FALSE),
#   triangulate(cl, tri = "min_lsp", mpd_based = FALSE),
#   triangulate(cl, tri = "min_lssp", mpd_based = FALSE),
#   triangulate(cl, tri = "min_elsp", pmf_evidence = pe, mpd_based = FALSE),
#   triangulate(cl, tri = "min_elssp", pmf_evidence = pe, mpd_based = FALSE)
# )

# (s <- structure(sapply(tris, function(x) jt_nbinary_ops(x, evidence) |> sum()), names = 1:7))
# sort(s)
# which.min(s)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                      CHILD: 20-25-230
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# l    <- readRDS("../../../../sandbox/r/bns/child.rds")
# cpts <- bnfit_to_cpts(l)
# cl   <- cpt_list(cpts)
# cp   <- jti::compile(cl)
# triangulate(cl)

# .map_lgl(cp$charge$C, function(x) inherits(x, "sparta_unity")) |> sum()

# none unities


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                      HEPAR2: 70-123-1453
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# l    <- readRDS("../../../../sandbox/r/bns/hepar2.rds")
# cpts <- bnfit_to_cpts(l)
# cl   <- cpt_list(cpts)
# cp   <- jti::compile(cl)

# .map_lgl(cp$charge$C, function(x) inherits(x, "sparta_unity")) |> sum()

# microbenchmark::microbenchmark(
#   jt(cp, propagate = "collect"),
#   times = 15
# )

# 30ms (UP) vs 23.5ms (OLD)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                      INSURANCE: 27-52-984
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# l    <- readRDS("../../../../sandbox/r/bns/insurance.rds")
# cpts <- bnfit_to_cpts(l)
# cl   <- cpt_list(cpts)
# triangulate(cl)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                      CANCER: 5-4-10
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## l    <- readRDS("../../../../sandbox/r/bns/cancer.rds")
## cpts <- bnfit_to_cpts(l)
## cl   <- cpt_list(cpts)
## cp   <- compile(cl, save_graph = TRUE)
## j    <- jt(cp)
## .map_dbl(j$charge$C, function(x) sum(x))


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                      ASIA
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# l <- readRDS(url("https://www.bnlearn.com/bnrepository/asia/asia.rds"))
# cpts <- bnfit_to_cpts(l)
# cl   <- cpt_list(cpts)

# pe <- c(smoke = .7, either = .4)

# tris <- list(
#   triangulate(cl, tri = "min_fill", mpd_based = FALSE),
#   triangulate(cl, tri = "min_sp", mpd_based = FALSE),
#   triangulate(cl, tri = "min_ssp", mpd_based = FALSE),
#   triangulate(cl, tri = "min_lsp", mpd_based = FALSE),
#   triangulate(cl, tri = "min_lssp", mpd_based = FALSE),
#   triangulate(cl, tri = "min_elsp", pmf_evidence = pe, mpd_based = FALSE),
#   triangulate(cl, tri = "min_elssp", pmf_evidence = pe, mpd_based = FALSE)
# )

# sapply(tris, function(x) jt_nbinary_ops(x, list(names(pe)), root = 1))

# plot(get_graph(cl))


# # min_fill based:
#  "min_fill"
#  "min_rfill"

# # min_sp based:
#  "min_sp"    # clique statespace   
#  "min_ssp"   # separator statespace   
#  "min_lsp"   # log clique statespace   
#  "min_lssp"  # log separator statespace 
#  "min_elsp"  # expected log clique statespace
#  "min_elssp" # expected log separtor statespace

# # usage:
# l <- readRDS(url("https://www.bnlearn.com/bnrepository/asia/asia.rds"))
# cpts <- bnfit_to_cpts(l)
# cl   <- cpt_list(cpts)
# pe <- c(smoke = .7, either = .4, dysp = 0.2)
# N <- length(names(cl))
# t <- triangulate(cl, tri = "min_elssp", pmf_evidence = pe, mpd_based = TRUE, perm = TRUE)
# jt_nbinary_ops(t, list(names(pe)))


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                      EARTHQUAKE: 5-4-10
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## l    <- readRDS("../../../../sandbox/r/bns/earthquake.rds")
## cpts <- bnfit_to_cpts(l)
## cl   <- cpt_list(cpts)
## cp   <- compile(cl, save_graph = TRUE)
## j    <- jt(cp)
## .map_dbl(j$charge$C, function(x) sum(x))



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                      SACHS: 11-17-178
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## l <- readRDS("../../../../sandbox/r/bns/sachs.rds")
## cpts <- bnfit_to_cpts(l)
## cl   <- cpt_list(cpts)
## cp   <- compile(cl, save_graph = TRUE)
## j    <- jt(cp)
## .map_dbl(j$charge$C, function(x) sum(x))


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                      ALARM: 37-46-509
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## l    <- readRDS("../../../../sandbox/r/bns/alarm.rds")
## cpts <- bnfit_to_cpts(l)
## cl   <- cpt_list(cpts)
## cp   <- jti::compile(cl)
## j    <- jt(cp)
## .map_dbl(j$charge$C, function(x) sum(x))


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                      WIN95PTS: 76-112-574
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## l    <- readRDS("../../../../sandbox/r/bns/win95pts.rds")
## cpts <- bnfit_to_cpts(l)
## cl   <- cpt_list(cpts)
## cp   <- jti::compile(cl)
## j    <- jt(cp)
## .map_dbl(j$charge$C, function(x) sum(x))
## size_mb(j)


# library(Rcpp)
# Rcpp::sourceCpp("../src/rooted_junction_tree.cpp")
# A <- matrix(0, nrow = 3, ncol = 3)
# A[1, 2] <- 1L
# A[2, 1] <- 1L
# A[2, 3] <- 1L
# A[3, 2] <- 1L

# B <- root_clique_tree(A, 1)
# B

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#        
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# library(ess)
# g <- ess::fit_graph(derma, q = 0, sparse_qic = TRUE)
# cl1 <- cpt_list(derma, ess::as_igraph(g))

# plot(g)



# cp1 <- compile(cl1)
# j1  <- jt(cp1)
# query_belief(j1, "ES")



# pl  <- pot_list(derma, ess::as_igraph(g))
# cpl <- compile(pl)
# jl  <- jt(cpl) 
# query_belief(j, "ES")


# cl <- cpt_list(derma, ess::as_igraph(g))
# cp <- compile(cl)
# j  <- jt(cp)
# query_belief(j, "ES")


# microbenchmark::microbenchmark(
#   pot_list(derma, ess::as_igraph(g)),
#   cpt_list(derma, ess::as_igraph(g)),
#   times = 5,
#   unit = "ms"
# )



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#        
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# library(ess)
# g <- ess::fit_graph(derma, sparse_qic = TRUE)
# cpt_list(derma, as_igraph(g))



