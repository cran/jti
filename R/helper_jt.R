leaves_jt <- function(x) {
  # x: rooted tree structure of a junctions tree (jt$schedule$collect$tree)
  which(colSums(x) == 0L)
}

parents_jt <- function(x, lvs) {
  # x:   rooted tree structure of a junctions tree (jt$schedule$collect$tree)
  # lvs: leaves of the junction tree
  par <- vector("list", length = length(lvs))
  for (i in seq_along(lvs)) {
    pari <- which(x[lvs[i], ] == 1L)
    par[[i]] <- pari
  }
  return(par)
}

valid_evidence <- function(dim_names, e) {
  lookup <- mapply(match, e, dim_names[names(e)])
  if (anyNA(lookup)) {
    return(FALSE)
  } else {
    return(TRUE)
  }
}

has_root_node <- function(x) UseMethod("has_root_node")

has_root_node.jt <- function(x) attr(x, "root_node") != ""

## new_schedule <- function(cliques) {
##   nc <- length(cliques)
##   clique_graph <- matrix(0L, nc, nc)
##   coll_tree   <- clique_graph
##   dist_tree   <- clique_graph

##   # Alg. 4.8 - Probabilistic Expert Systems (p.55)
##   for (i in seq_along(cliques)[-1L]) {
##     for (j in 1:(i-1L)) {
##       Ci   <- cliques[[i]]
##       Cj   <- cliques[[j]]
##       Hi_1 <- unlist(cliques[1:(i-1L)])
##       Si   <- intersect(Ci, Hi_1)
##       if (all(Si %in% Cj)) {
##         clique_graph[i, j] <- 1L
##         clique_graph[j, i] <- 1L
##         is_new_directed_edge <- !neq_empt_int(which(coll_tree[i, ] == 1L))
##         if (is_new_directed_edge) {
##           coll_tree[i, j] <- 1L
##           dist_tree[j, i] <- 1L
##         }
##       }
##     }
##   }

##   coll_lvs <- leaves_jt(coll_tree)
##   dist_lvs <- leaves_jt(dist_tree)

##   attr(coll_tree, "leaves")  <- coll_lvs
##   attr(dist_tree, "leaves")  <- dist_lvs

##   attr(coll_tree, "parents") <- parents_jt(coll_tree, coll_lvs)
##   attr(dist_tree, "parents") <- parents_jt(dist_tree, dist_lvs)

##   collect    <- list(cliques = cliques, tree = coll_tree)
##   distribute <- list(cliques = cliques, tree = dist_tree)
  
##   return(
##     list(
##       collect = collect ,
##       distribute = distribute,
##       clique_graph = clique_graph,
##       clique_root = "C1"
##     )
##   )
  
## }

new_schedule2 <- function(cliques_chr, cliques_int, root_node) {
  
  # mcs promise that the root_node lives in clique one
  jrn <- if (root_node != "") 1L else 0L
  
  rjt <- rooted_junction_tree(cliques_int, jrn)

  coll_tree    <- rjt$collect
  dist_tree    <- rjt$distribute
  clique_root  <- rjt$clique_root
  clique_graph <- coll_tree + dist_tree


  coll_lvs <- leaves_jt(coll_tree)
  dist_lvs <- leaves_jt(dist_tree)

  attr(coll_tree, "leaves")  <- coll_lvs
  attr(dist_tree, "leaves")  <- dist_lvs

  attr(coll_tree, "parents") <- parents_jt(coll_tree, coll_lvs)
  attr(dist_tree, "parents") <- parents_jt(dist_tree, dist_lvs)

  collect    <- list(cliques = cliques_chr, tree = coll_tree)
  distribute <- list(cliques = cliques_chr, tree = dist_tree)
  
  return(
    list(
      collect      = collect ,
      distribute   = distribute,
      clique_graph = clique_graph,
      clique_root  = paste("C", clique_root, sep = "") # TODO: Just return the index
    )
  )
}


prune_jt <- function(jt) {

  direction <- attr(jt, "direction")
  x <- if (direction == "collect") jt$schedule$collect else jt$schedule$distribute

  if (identical(x, "full")) {
    stop("The junction tree has already been propagated in this direction!")  
  }

  leaves <- attr(x$tree, "leaves")
  pars   <- attr(x$tree, "parents")

  
  if (length(leaves) == ncol(x$tree)) { # If all nodes left are singletons in distribute
    x$cliques <- NULL
  } else {
    x$cliques <- x$cliques[-leaves]
    x$tree    <- x$tree[-leaves, -leaves]   
  }
  
  has_arrived_at_root <- length(x$cliques) < 2L
  if (has_arrived_at_root) {
    if (direction == "collect") {
      jt$schedule$collect    <- "full"
      attr(jt, "direction")  <- "distribute"

      # Normalize clique_root
      cr <- attr(jt, "clique_root")
      probability_of_evidence <- sum(sparta::vals(jt$charge$C[[cr]]))
      attr(jt, "probability_of_evidence") <- probability_of_evidence
      jt$charge$C[[cr]] <- sparta::normalize(jt$charge$C[[cr]])

    } else {
      jt$schedule$distribute <- "full"
      attr(jt, "direction")  <- "full"
    }
    return(jt)
  }
  
  attr(x$tree, "leaves")  <- leaves_jt(x$tree)
  attr(x$tree, "parents") <- parents_jt(x$tree, attr(x$tree, "leaves"))

  if (direction == "collect") {
    jt$schedule$collect <- list(cliques = x$cliques, tree = x$tree)
  } else {
    jt$schedule$distribute <- list(cliques = x$cliques, tree = x$tree)
  }
  
  return(jt)
}

set_evidence_jt <- function(charge, cliques, evidence) {
  for (k in seq_along(charge$C)) {
    Ck <- names(charge$C[[k]])
    for (i in seq_along(evidence)) {
      e     <- evidence[i]
      e_var <- names(e)
      e_val <- unname(e)
      if (e_var %in% Ck) {
        m <- try(sparta::slice(charge$C[[k]], e), silent = TRUE)
        if (inherits(m, "try-error")) {
          stop(
            "The evidence leads to a degenerate distribution ",
            "since the evidence was never observed in one or ",
            "more of the clique potentials.",
            call. = FALSE
          )
        }
        charge$C[[k]] <- m
      }
    }
  }
  return(charge)
}


## new_schedule_grain <- function(grain_obj) {
##   # grain_obj = jt_gr$rip

##   cliques <- grain_obj$cliques
  
##   nc <- length(cliques)
##   clique_graph <- matrix(0L, nc, nc)
##   coll_tree   <- clique_graph
##   dist_tree   <- clique_graph

##   for (k in seq_along(grain_obj$childList)) {
##     # parent = k
##     for (i in grain_obj$childList[[k]]) {
##       # i = child
##       coll_tree[i, k] <- 1L
##       dist_tree[k, i] <- 1L
##     }
##   }

##   coll_lvs <- leaves_jt(coll_tree)
##   dist_lvs <- leaves_jt(dist_tree)

##   attr(coll_tree, "leaves")  <- coll_lvs
##   attr(dist_tree, "leaves")  <- dist_lvs

##   attr(coll_tree, "parents") <- parents_jt(coll_tree, coll_lvs)
##   attr(dist_tree, "parents") <- parents_jt(dist_tree, dist_lvs)

##   collect    <- list(cliques = cliques, tree = coll_tree)
##   distribute <- list(cliques = cliques, tree = dist_tree)
  
##   return(
##     list(
##       collect = collect ,
##       distribute = distribute,
##       clique_graph = clique_graph,
##       clique_root = "C1"
##     )
##   )
## }


new_jt <- function(x, evidence = NULL, flow = "sum") {
  # x: a charge object returned from compile
  #  - a list with the charge and the cliques

  charge  <- x$charge
  cliques <- x$cliques

  if (!is.null(evidence)) charge <- set_evidence_jt(charge, cliques, evidence)

  # schedule  <- new_schedule_grain(grain_obj)
  # schedule  <- new_schedule(cliques)
  schedule  <- new_schedule2(cliques, attr(x, "cliques_int"), attr(x, "root_node"))
  attr(x, "cliques_int") <- NULL

  jt <- list(
    schedule = schedule[1:2], # collect and distribute
    charge   = charge,
    cliques  = cliques,
    clique_graph = schedule$clique_graph
  )
  
  class(jt)               <- c("jt", class(jt))
  attr(jt, "direction")   <- "collect" # collect, distribute or full
  attr(jt, "lookup")      <- attr(x, "lookup")
  attr(jt, "flow")        <- flow
  attr(jt, "root_node")   <- attr(x, "root_node")
  attr(jt, "clique_root") <- schedule$clique_root
  
  if (flow == "max") {
    # most probable explanation
    all_vars <- names(attr(x, "dim_names"))
    attr(jt, "mpe") <- structure(
      vector("character", length = length(all_vars)),
      names = all_vars
    )
  }
  return(jt)
}

send_messages <- function(jt, flow = "sum") {

  direction <- attr(jt, "direction")
  if (direction == "full") {
    message("The junction tree is already propagated fully. jt is returned")
    return(jt)
  }

  x   <- if (direction == "collect") jt$schedule$collect else jt$schedule$distribute
  lvs <- attr(x$tree, "leaves")
  par <- attr(x$tree, "parents")
  
  for (k in seq_along(lvs)) {

    lvs_k <- lvs[k]
    par_k <- par[[k]]
    
    for (pk in par_k) {
      
      C_lvs_k <- x$cliques[[lvs_k]]
      C_par_k <- x$cliques[[pk]]
      Sk      <- intersect(C_lvs_k, C_par_k)

      # TODO: Whats the point in naming them Ck etc.? Just use the index?
      C_lvs_k_name <- names(x$cliques)[lvs_k]
      C_par_k_name <- names(x$cliques)[pk]
      
      message_k_names <- setdiff(C_lvs_k, Sk)
      
      if (direction == "collect") {
        message_k <- sparta::marg(jt$charge$C[[C_lvs_k_name]], message_k_names, attr(jt, "flow"))
        jt$charge$C[[C_par_k_name]] <- sparta::mult(jt$charge$C[[C_par_k_name]], message_k)
        jt$charge$C[[C_lvs_k_name]] <- sparta::div(jt$charge$C[[C_lvs_k_name]], message_k)  
      }

      if (direction == "distribute") {
        if (attr(jt, "flow") == "max") {

          ## Find the max cell and change the potential before sending the information:
          max_idx  <- sparta::which_max_idx(jt$charge$C[[C_lvs_k_name]])
          max_cell <- sparta::which_max_cell(jt$charge$C[[C_lvs_k_name]])
          max_mat  <- jt$charge$C[[C_lvs_k_name]][, max_idx, drop = FALSE]
          max_val  <- attr(jt$charge$C[[C_lvs_k_name]], "vals")[max_idx]
          max_dn   <- attr(jt$charge$C[[C_lvs_k_name]], "dim_names")

          jt$charge$C[[C_lvs_k_name]] <- sparta::sparta_struct(
            max_mat,
            max_val,
            max_dn
          )
          attr(jt, "mpe")[names(max_cell)] <- max_cell
        }

        # Send the message
        message_k <- sparta::marg(jt$charge$C[[C_lvs_k_name]], message_k_names, attr(jt, "flow"))
        jt$charge$C[[C_par_k_name]] <- sparta::mult(jt$charge$C[[C_par_k_name]], message_k)
        jt$charge$S[[paste("S", pk, sep = "")]] <- message_k
        
        if (attr(jt, "flow") == "max") {
          ## Record the max cell for the parent potential
          max_cell <- sparta::which_max_cell(jt$charge$C[[C_par_k_name]])
          attr(jt, "mpe")[names(max_cell)] <- max_cell
        }
      }
    }
  }
  prune_jt(jt)
}