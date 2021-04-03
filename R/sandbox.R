# size_mb <- function(x) {
#   format(object.size(x), units = "Mb", standard = "auto", digits = 1L)
# }

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                 MUNIN: 1041-1397-80592
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# l    <- readRDS("../../../../sandbox/r/bns/munin.rds")
# cpts <- jti::bnfit_to_cpts(l)
# cl   <- cpt_list(cpts)

## lapply(cp_munin$charge$C, names)

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
# grps_chr <- .map_chr(str_split(ll, "_"), function(s) {
#   paste0(s[1], "_", s[2])
# })
# grps <- split(leaves, grps_chr)
# set.seed(3007)
# evars <- unname(.map_chr(grps, function(x) x[sample(1:length(x), 1L)]))

# e <- .map_chr(evars, function(x) {
#   dd <- dim_names(cl)[[x]]
#   sample(dd, 1L)
# })


# tictoc::tic()
# cp_munin <- jti::compile(cl, r, tri = "minimal")
# tictoc::toc() # 12 sec

# j <- jt(cp_munin)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                         LINK: 724-1125-14211
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# l    <- readRDS("../../../../sandbox/r/bns/link.rds")
# cpts <- bnfit_to_cpts(l)
# cl   <- cpt_list(cpts)

# tictoc::tic()
# j  <- jt(cp)
# tictoc::toc()

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
## l    <- readRDS("../../../../sandbox/r/bns/mildew.rds")
## cpts <- bnfit_to_cpts(l)
## cl   <- cpt_list(cpts)
## cp   <- compile(cl, save_graph = TRUE)

## tictoc::tic()
## j    <- jt(cp, propagate = "full")
## tictoc::toc() # 2.5 -> 2

## .map_dbl(j$charge$C, function(x) sum(x))
## size_mb(j)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                      HAILFINDER: 56-66-2656
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## l    <- readRDS("../../../../sandbox/r/bns/hailfinder.rds")
## cpts <- bnfit_to_cpts(l)
## cl   <- cpt_list(cpts)

## tictoc::tic()
## cp   <- jti::compile(cl, opt = "min_sp")
## tictoc::toc()

## size_mb(cp)

## j <- jt(cp, propagate = "full")
## .map_dbl(j$charge$C, function(x) sum(x))


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                      BARLEY: 48-84-114005
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## l    <- readRDS("../../../../sandbox/r/bns/barley.rds")
## cpts <- bnfit_to_cpts(l)
## cl   <- cpt_list(cpts)

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


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                      HEPAR2: 70-123-1453
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## l    <- readRDS("../../../../sandbox/r/bns/hepar2.rds")
## cpts <- bnfit_to_cpts(l)
## cl   <- cpt_list(cpts)
## cp   <- jti::compile(cl) 
## j    <- jt(cp)
## size_mb(j) 
 

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
# triangulate(cl)
# cp   <- compile(cl, save_graph = TRUE)
# j    <- jt(cp)
# .map_dbl(j$charge$C, function(x) sum(x))


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
## l    <- readRDS("../../../../sandbox/r/bns/andes.rds")
## cpts <- bnfit_to_cpts(l)
## cl   <- cpt_list(cpts)
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
