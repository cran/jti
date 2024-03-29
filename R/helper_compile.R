as_undirected_igraph <- function(g) igraph::as.undirected(g)

parents_igraph <- function(g) {
  # g: dag
  # out: named list where a name is the child and elements are parents
  A  <- igraph::as_adjacency_matrix(g, sparse = FALSE)
  cn <- colnames(A)
  if (is.null(cn)) {
    stop("The vertices in the igraph object must have names")
  }

  par_lst <- structure(vector("list", ncol(A)), names = cn)
  for (name in cn) {
    z <- which(A[, name] == 1L)
    par_lst[[name]] <- names(z)
  }

  par_lst
}

parents_cpt_list <- function(x) {
  parents <- structure(lapply(seq_along(x), function(i) {
    child   <- names(x)[i]
    setdiff(names(x[[i]]), child)
  }), names  = names(x))
}


graph_from_cpt_list <- function(x) {
  g <- igraph::make_empty_graph(n = length(x))
  g <- igraph::set_vertex_attr(g, "label", value = names(x))
  g <- igraph::set_vertex_attr(g, "name", value = names(x))
  edges <- lapply(seq_along(x), function(i) {
    child   <- names(x)[i]
    parents <- setdiff(names(x[[i]]), child)
    as.matrix(expand.grid(parents, child, stringsAsFactors = FALSE))
  })
  edges <- do.call(rbind, edges)
  for (k in 1:nrow(edges)) {
    g <- igraph::add_edges(g, unname(edges[k, ]))
  }
  g
}

graph_from_pot_list <- function(x, cliques, nodes) {
  g <- igraph::make_empty_graph(n = length(nodes), directed = FALSE)
  g <- igraph::set_vertex_attr(g, "label", value = nodes)
  g <- igraph::set_vertex_attr(g, "name", value = nodes)
  edges <- lapply(cliques, function(clique) {
    if (length(clique) == 1L) return(NULL)
    t(utils::combn(clique, 2, simplify = TRUE))
  })
  edges <- Filter(neq_null, edges)
  edges <- unique(do.call(rbind, edges))
  for (k in 1:nrow(edges)) {
    g <- igraph::add_edges(g, unname(edges[k, ]))
  }
  g
}


moralize_igraph <- function(g, parents) {
  g <- igraph::as.undirected(g)
  for (p in parents) {
    if (length(p) > 1) {
      pairs <- utils::combn(p, 2, simplify = FALSE)
      for (ps in pairs) {
        if (!igraph::are_adjacent(g, ps[1], ps[2])) {
          g <- g + igraph::edge(ps[1], ps[2]) 
        }
      }
    }
  }
  g
}

add_joint_vars_igraph <- function(g, joint_vars) {
  if (length(joint_vars) < 1) return(g)
  pairs <- utils::combn(joint_vars, 2,  simplify = FALSE)
  for (ps in pairs) {
    if (!igraph::are_adjacent(g, ps[1], ps[2])) {
      g <- g + igraph::edge(ps[1], ps[2]) 
    }
  }
  g
}

triangulate_igraph <- function(g) {
  igraph::is.chordal(g, fillin = FALSE, newgraph = TRUE)$newgraph
}

construct_cliques <- function(adj) {
  rip_ <- rip(adj, check = FALSE)
  structure(rip_$C, names = paste("C", 1:length(rip_$C), sep = ""))
}

construct_cliques_int <- function(adj_mat) {
  # cliques_int is needed to construct the junction tree in new_jt -> new_schedule
  dimnames(adj_mat) <- lapply(dimnames(adj_mat), function(x) 1:nrow(adj_mat))
  adj_lst_int       <- as_adj_lst(adj_mat)
  cliques_int       <- rip(adj_lst_int)$C
  lapply(cliques_int, as.integer)
}
