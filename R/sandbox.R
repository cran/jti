# size_mb <- function(x) {
#   format(object.size(x), units = "Mb", standard = "auto", digits = 1L)
# }

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                 MUNIN: 1041-1397-80592
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# l    <- readRDS("../../../../sandbox/r/bns/munin.rds")
# cpts <- jti::bnfit_to_cpts(l)

# tictoc::tic()
# cl <- jti::cpt_list(cpts)
# cp_munin <- jti::compile(cl)
# j  <- jt(cp_munin)
# tictoc::toc() # 130 -> 58 -> 24

# 75
# 68

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                         LINK: 724-1125-14211
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# l    <- readRDS("../../../../sandbox/r/bns/link.rds")
# cpts <- bnfit_to_cpts(l)
# cl   <- cpt_list(cpts)
# cp   <- compile(cl, tri = "min_fill")

# .map_lgl(cp$charge$C, function(x) inherits(x, "sparta_unity")) |> sum()
# .map_dbl(cp$charge$C, sparta::table_size) |> max()

# tictoc::tic()
# j   <- jt(cp, propagate = "full") # 92 sec.
# tictoc::toc()

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

# library(dplyr)
# l    <- readRDS(url("https://www.bnlearn.com/bnrepository/link/link.rds"))
# cl   <- l %>%
#   # lapply(unclass) %>%
#   bnfit_to_cpts() %>%
#   # lapply(function(x) x + .1) %>%
#   cpt_list()

# m <- mpd(cl)
# i <- which.max(sapply(m$primes_int, length))
# es <- colnames(m$graph)[m$primes_int[[i]]]


# set.seed(1)


# pmf <- structure(runif(5, .6, 1), names = sample(es, 5))
# pmf <- structure(runif(1, .6, 1), names = sample(es, 1))

# pmf <- structure(runif(13, 1, 1), names = sample(es, 13))

# t0 <- triangulate(cl, tri = "min_sfill")
# t2 <- triangulate(cl, tri = "min_fill")

# .map_int(t0$cliques, length) |> max()
# .map_int(t2$cliques, length) |> max()


# res <- unlist(parallel::mclapply(mc.cores = 3, X = 1:500, FUN = function(a) {
#   print(a)
  
#   ev <- lapply(1:500, function(k) {
#     sample(es, 12)
#   })

#   nemf <- sapply(ev, function(x) {
#     jt_nbinary_ops(t0, x)
#   }) |> sum()

#   nmf <- sapply(ev, function(x) {
#     jt_nbinary_ops(t2, x)
#   }) |> sum()

#   nemf / nmf
  
# }))


# TODO: Are we doing the right cals in the benchmark?
# make new_triang constructers that accepts a permutation
# of the moral graph!


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
# l <- readRDS(url("https://www.bnlearn.com/bnrepository/asia/asia.rds"))
# cpts <- bnfit_to_cpts(l)
# cl   <- cpt_list(cpts)
# triangulate(cl, tri = "min_fill", perm = sample(1:8, 8))

# plot(get_graph(cl))

# e <- c(
#   either = "no",
#   lung   = "no",
#   tub    = "yes",
#   xray   = "no",
#   dysp   = "no",
#   bronc  = "yes",
#   smoke  = "yes",
#   asia   = "yes"
# )

# cp <- compile(cl, e)
# j  <- jt(cp, propagate = "no")

# # j <- send_messages(j)
# # parents(j)
# # leaves(j)

# query_belief(j, c("bronc"))

# mpe(j)
# j$charge$C

# h  <- jt(cp, flow = "sum", propagate = "no")
# plot(h)

# sparta::marg(j$charge$C$C3, c("either", "lung"))
# sparta::marg(j$charge$C$C4, c("smoke", "lung"))

# query_belief(j, c("bronc"), "joint")

# has_inconsistencies(cp)
# query_evidence(j)

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
#            TEST PERMUTATION AND RANDOM TIES
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# library(igraph)
# library(dplyr)

# # make net
# net <- c(
#   "a", "b",
#   "a", "k",
#   "a", "j",
#   "b", "c",
#   "c", "d",
#   "d", "e",
#   "e", "f",
#   "g", "f",
#   "h", "g",
#   "i", "h",
#   "i", "d",
#   "i", "g",
#   "j", "i",
#   "k", "j"
# ) %>%
#   as.character() %>%
#   graph(directed = TRUE)

# m <- moralize_igraph(net, parents_igraph(net))

# # # visualize
# par(mfrow = c(1, 1))
# plot(net); plot(m)

# # # make dummy data
# d <- replicate(
#   igraph::vcount(net),
#   as.character(sample(1:2, 1000, replace = TRUE))
# ) |> as.data.frame()
# colnames(d) <- V(net)$name
# head(d)


# # # inspection of elimination game
# cl   <- cpt_list(d, net)
# perm <- sample(1:ncol(d), ncol(d))
# perm

# triangulate(cl, tri = "min_fill", perm = perm)
