# size_mb <- function(x) {
#   format(object.size(x), units = "Mb", standard = "auto", digits = 1L)
# }

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                 MUNIN: 1041-1397-80592
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# l    <- readRDS("../../../../sandbox/r/bns/munin.rds")
# cpts <- jti::bnfit_to_cpts(l)
# cl   <- jti::cpt_list(cpts)

# # mean: 1.5 sec
# microbenchmark::microbenchmark(
#   xx <- Reduce(sparta::mult, cl[1:20]),
#   times = 2
# )


# k <- 5
# x <- letters[1:k]
# X <- array(
#   1L,
#   rep(k,k),
#   structure(replicate(k, x, FALSE), names = LETTERS[1:k])
# )

# SX <- sparta::as_sparta(X)

# # CRAN: 4.23 ms
# # DEVTOOLS: 8.08 ms
# # SOURCE BUILD: 4.29 ms
# microbenchmark::microbenchmark(
#   for (j in 1:k) sparta::marg(SX, names(SX)[1:j]),
#   times = 1000
# )

## size_mb(cp_munin) # 5.8 Mb

# tictoc::tic()
# j <- jt(cp_munin)
# tictoc::toc() # 130 -> 58 -> 24

## size_mb(j) # 151 MB
## .map_dbl(j$charge$C, sum)

# library(stringr)
# library(igraph)
# library(graph)

# g <- get_graph(cl)
# h <- igraph::igraph.to.graphNEL(g)
# plot(h)


# parents  <- get.edgelist(g)[, 1]
# children <- get.edgelist(g)[, 2]
# leaves   <- unique(children[which(!(children %in% parents))])
# leaves   <- leaves[!stringr::str_detect(leaves, "DUMMY")]

# # grouping leave nodes:
# grps_chr <- .map_chr(str_split(leaves, "_"), function(s) {
#   paste0(s[1], "_", s[2])
# })
# grps <- split(leaves, grps_chr)

# set.seed(3)
# evars <- unname(.map_chr(grps, function(x) x[sample(1:length(x), 1L)]))

# e <- .map_chr(evars, function(x) {
#   dd <- dim_names(cl)[[x]]
#   sample(dd, 1L)
# })

# v <- names(e)

# # TODO: Deduce evidence_nodes from evidence parameter?
# cp1 <- compile(cl, evidence = e, tri = "evidence", evidence_nodes = v)
# cp2 <- compile(cl, evidence = e, tri = "min_fill")

# microbenchmark::microbenchmark(
#   j1 <- jt(cp1), # 55
#   j2 <- jt(cp2), # 55
#   times = 1
# )

# sum(.map_lgl(get_cliques(j1), function(x) "sort" %in% x))
# sum(.map_lgl(get_cliques(j2), function(x) "sort" %in% x))


# tictoc::tic()
# cp_munin <- jti::compile(cl, r, tri = "minimal")
# tictoc::toc() # 12 sec

# j <- jt(cp_munin)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                         LINK: 724-1125-14211
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# library(jti)

# l    <- readRDS("../../../../sandbox/r/bns/link.rds")
# cpts <- bnfit_to_cpts(l)
# cl   <- cpt_list(cpts)

# cp   <- compile(cl, tri = "min_fill")

# tictoc::tic()
# j  <- jt(cp)
# tictoc::toc()


# 113 : jti 0.7.0, sparta 0.7.2


# Try with this jti and CRAN sparta



# # 324 -> 260 -> 92
# size_mb(j) # 1,6Gb (4Gb max when running)

# tictoc::tic()
# j  <- jt(cp, e)
# tictoc::toc()

# size_mb(j) # 181.1Gb (4Gb max when running)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                        DIABETES: 413-602-429409
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# l    <- readRDS("../../../../sandbox/r/bns/diabetes.rds")
# cpts <- bnfit_to_cpts(l)
# cl   <- cpt_list(cpts)

# dim_names(cl)

# e <- c(
#   dm_1   = "0_00_kg_m2",
#   dm_3   = "0_00_kg_m2",
#   dm_4   = "0_00_kg_m2",
#   foto_1 = "0_00_kg_m2",
#   straaling_4 = "300___349_MJ_m2"
# )


# set.seed(65)
# e <- unlist(lapply(dim_names(cl)[sample(1:413, 1)], function(x) x[1]))
# v <- names(e)

# cp1 <- compile(cl, evidence = e, tri = "evidence", evidence_nodes = v)
# cp2 <- compile(cl, evidence = e, tri = "min_fill")

# j <- jt(cp2)

# idx <- which(.map_lgl(j$cliques, function(x) v %in% x))

# j$cliques[[idx[1]]]
# names(j$charge$C[[idx[1]]])

# microbenchmark::microbenchmark(
#   jt(cp1),
#   jt(cp2),
#   times = 1
# )

# tt <- triangulate(cl)

# tictoc::tic()
# cp   <- jti::compile(cl)
# tictoc::toc() # 33 s -> 2s

# size_mb(cp) # 5 mb

# e    <- .map_chr(attr(cl, "dim_names"), `[[`, 1L)[sample(1:400, 8)]
# tictoc::tic()
# j    <- jt(cp, e)
# tictoc::toc() # 23.5 s -> 21.3

## size_mb(j) # 126 mb

## .map_dbl(j$charge$C, function(x) sum(x))


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                      MILDEW: 35-46-540150
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# l    <- readRDS("../../../../sandbox/r/bns/mildew.rds")
# cpts <- bnfit_to_cpts(l)
# cl   <- cpt_list(cpts)

# ####### TODO: SHOW THIS EXAMPLE!!!!
# prod(.map_int(sparta::dim_names(cp$charge$C$C1), length))
# ncol(cp$charge$C$C1) * 4
# ####### TODO: SHOW THIS EXAMPLE!!!!

# prod(.map_int(sparta::dim_names(cp$charge$C$C2), length))
# ncol(cp$charge$C$C2) * 4

# # TODO:
# # 1) Somehow highlight this benchmark
# # 2) Show the memory saving of the CR
# # 3) Remember Sorens idea

# cp <- compile(cl, tri = "min_fill")
# gr <- gRain::grain(gg)

# microbenchmark::microbenchmark(
#   j <- jt(cp, propagate = "full"),
#   g <- gRbase::compile(gr, propagate = TRUE),
#   times = 10
# )

# # TODO: ALSO SEE THIS...
# prod(.map_int(sparta::dim_names(j$charge$C$C2), length))
# ncol(j$charge$C$C2) * 4

# MORALE::: FÅ VARIABLE, HVER MED MEGET STORT STATESPACE!!!


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                      HAILFINDER: 56-66-2656
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# l    <- readRDS("../../../../sandbox/r/bns/hailfinder.rds")
# cpts <- bnfit_to_cpts(l)
# cl   <- cpt_list(cpts)


## tictoc::tic()
## cp   <- jti::compile(cl, opt = "min_sp")
## tictoc::toc()

## size_mb(cp)

## j <- jt(cp, propagate = "full")
## .map_dbl(j$charge$C, function(x) sum(x))


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                      BARLEY: 48-84-114005
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# library(jti)

# l    <- readRDS("../../../../sandbox/r/bns/barley.rds")
# cpts <- bnfit_to_cpts(l)
# cl   <- cpt_list(cpts)
# cp   <- compile(cl, tri = "min_fill")

# tictoc::tic()
# jt(cp)
# tictoc::toc()


# current jti, current sparta; 105.301 sec elapsed

# CRAN jti and CRAN sparta; 45.067 sec elapsed

# plot(get_graph(cl))

# e <- c(
#   s2528  = "x0_5",
#   sort   = "Baracuda",
#   ngodnt = "x_45"
# )

# v <- names(e)

# # TODO: Deduce evidence_nodes from evidence parameter?
# cp1 <- compile(cl, evidence = e, tri = "evidence", evidence_nodes = v)
# cp2 <- compile(cl, evidence = e, tri = "min_fill")

# # NOTE: Barley is hard for jti/sparta since it is "completely dense",
# # with tables having >7e+10 elements.

# microbenchmark::microbenchmark(
#   j1 <- jt(cp1),
#   j2 <- jt(cp2),
#   times = 1
# )

# sum(.map_lgl(get_cliques(j1), function(x) "sort" %in% x))
# sum(.map_lgl(get_cliques(j2), function(x) "sort" %in% x))

# jt(cp2)

## tictoc::tic()
## cp   <- jti::compile(cl)
## tictoc::toc()

## tictoc::tic()
## j    <- jt(cp)
## tictoc::toc()  # 118 -> 98

## .map_dbl(j$charge$C, function(x) sum(x))


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                      CHILD: 20-25-230
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# l    <- readRDS("../../../../sandbox/r/bns/child.rds")
# cpts <- bnfit_to_cpts(l)
# cl   <- cpt_list(cpts)
# cp   <- jti::compile(cl, save_graph = TRUE)
# j    <- jt(cp)

# m    <- mpd(cl)
# g    <- m$graph
# h    <- g[m$primes_int[[5]], m$primes_int[[5]]]
# library(dplyr)
# igraph::graph_from_adjacency_matrix(h, "undirected") %>%
#   plot(vertex.size = 10, vertex.label = NA)

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                      HEPAR2: 70-123-1453
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# l    <- readRDS("../../../../sandbox/r/bns/hepar2.rds")
# cpts <- bnfit_to_cpts(l)
# cl   <- cpt_list(cpts)

# plot(get_graph(cl), vertex.size = 1)


# dim_names(cl)

# e <- c(
#   fatigue    = "present",
#   gallstones = "present",
#   PBC        = "absent",
#   bleeding   = "present",
#   Steatosis  = "absent",
#   ChHepatitis = "persistent"
#   # Cirrhosis  = "absent"
# )
# v <- names(e)

# cp1 <- compile(cl, evidence = e, tri = "evidence", evidence_nodes = v)
# cp2 <- compile(cl, evidence = e, tri = "min_fill")

# microbenchmark::microbenchmark(
#   j1 <- jt(cp1),
#   j2 <- jt(cp2),
#   times = 25
# )

# sum(.map_lgl(get_cliques(j1), function(x) "fatigue" %in% x))
# sum(.map_lgl(get_cliques(j2), function(x) "fatigue" %in% x))

# sum(.map_lgl(get_cliques(j1), function(x) "gallstones" %in% x))
# sum(.map_lgl(get_cliques(j2), function(x) "gallstones" %in% x))

# sum(.map_lgl(get_cliques(j1), function(x) "bleeding" %in% x))
# sum(.map_lgl(get_cliques(j2), function(x) "bleeding" %in% x))

# sum(.map_lgl(get_cliques(j1), function(x) "Steatosis" %in% x))
# sum(.map_lgl(get_cliques(j2), function(x) "Steatosis" %in% x))

# sum(.map_lgl(get_cliques(j1), function(x) "ChHepatitis" %in% x))
# sum(.map_lgl(get_cliques(j2), function(x) "ChHepatitis" %in% x))

# jt(cp2)
 

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                      INSURANCE: 27-52-984
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# l    <- readRDS("../../../../sandbox/r/bns/insurance.rds")
# cpts <- bnfit_to_cpts(l)
# cl   <- cpt_list(cpts)

# alpha <- sample(1:length(attr(cl, "dim_names")), length(attr(cl, "dim_names")))
# cp   <- jti::compile(cl, tri = "alpha", alpha = alpha)

# cp   <- jti::compile(cl, tri = "min_fill")
# e    <- .map_chr(attr(cl, "dim_names"), `[[`, 1L)

# j    <- jt(cp, e)
# .map_dbl(j$charge$C, function(x) sum(x))


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
# l    <- readRDS("../../../../sandbox/r/bns/asia.rds")
# cpts <- bnfit_to_cpts(l)
# cl   <- cpt_list(cpts)
# triangulate(cl, tri = "min_rfill")$alpha
# triangulate(cl)

# cp    <- compile(cl, tri = "alpha", alpha = sample(names(cl), length(names(cl))))
# jt   <- jt(cp, propagate = "no")
# tt   <- triangulate(cl)
# jt_nbinary_ops(tt)
# jt_nbinary_ops(jt)

# microbenchmark::microbenchmark(
#   compile(cl, tri = "evidence", evidence_nodes = character(0), .mpd = TRUE),
#   compile(cl, tri = "evidence", evidence_nodes = character(0)),
#   times = 1
# )

# plot(get_graph(cl))

# # Evidence triangulation
# obj <-  new_min_fill_triang(g)
# elim_game(obj)

# cp1   <- compile(cl, c(lung = "yes"), tri = "evidence", evidence_nodes = "lung")
# cp1$cliques

# cp2   <- compile(cl, tri = "evidence", evidence_nodes = "smoke")
# cp2$cliques

# cp3   <- compile(cl, tri = "evidence", evidence_nodes = "bronc")
# cp3$cliques


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
#                      ANDES: 223-338-1157
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# l    <- readRDS("../../../../sandbox/r/bns/andes.rds")
# cpts <- bnfit_to_cpts(l)
# cl   <- cpt_list(cpts)

# m    <- mpd(cl)
# g    <- m$graph
# h    <- g[m$primes_int[[100]], m$primes_int[[100]]]
# library(dplyr)
# igraph::graph_from_adjacency_matrix(h, "undirected") %>%
#   plot(vertex.size = 1, vertex.label = NA)

## cp   <- jti::compile(cl)
## j    <- jt(cp)
## size_mb(cp)

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


# make_all_evidence_data_and_pmf2 <- function(dat, mrf, nsims, mar = TRUE) {
#   # mrf = named vector
#   nr       <- nrow(dat)
#   dat_miss <- dat
#   for (en in names(mrf)) {
#     col <- dat_miss[, en]
#     na_idx <- sample(c(FALSE, TRUE), nr, TRUE, c(1-mrf[en], mrf[en]))
#     col[na_idx] <- NA
#     dat_miss[, en] <- col
#   }
#   dat_miss
# }

        
# # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# easypackages::libraries(
#   "tidyverse", # dplyr, ggplot2, readr, stringr, purrr, tidyr
#   "fs",
#   "glue",
#   "igraph"
# ) 

# nsims   <- 1500

# net <- c(
#   1,10,
#   2,10,
#   10,3,
#   10,9,
#   4,3,
#   4,5,
#   5,6,
#   9,8,
#   6,7,
#   8,7
# ) %>%
#   as.character() %>%
#   graph(directed = TRUE)

# # BN
# plot(net)

# moral_graph <- jti:::moralize_igraph(net, jti:::parents_igraph(net))
# plot(moral_graph)

# nodes_net <- igraph::V(net)$name
# lvls_net  <- structure(sample(3:9, length(nodes_net), TRUE), names = nodes_net)
# dat       <- jti::sim_data_from_bn(net, lvls_net, nsims)

# # EVIDENCE PMF:
# evars <- c("8", "5")# sample(nodes_net, 4)
# mrf_ <- structure(runif(length(evars), 0.2, .6), names = evars)

# nodes_mrf   <- names(mrf_)
# dat_miss    <- make_all_evidence_data_and_pmf2(dat, mrf = mrf_, nsims)
# cl          <- cpt_list(dat_miss, net)
# t_evidence  <- triangulate(cl, tri = "evidence2", pmf_evidence = mrf_)
