% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helper_binary_ops.R
\name{jt_nbinary_ops}
\alias{jt_nbinary_ops}
\alias{jt_nbinary_ops.triangulation}
\title{Number of Binary Operations}
\usage{
jt_nbinary_ops(x, evidence = list(), root = NULL, nc = 1)

\method{jt_nbinary_ops}{triangulation}(x, evidence = list(), root = NULL, nc = 1)
}
\arguments{
\item{x}{A junction tree object or an object returned from
the triangulation function}

\item{evidence}{List of character vectors with evidence nodes}

\item{root}{Integer specifying the root node in the junction tree}

\item{nc}{Integer. The number of cores to be used in parallel}
}
\description{
Number of binary operations needed to propagate in a junction tree
given evidence, using the Lauritzen-Spiegelhalter scheme
}
