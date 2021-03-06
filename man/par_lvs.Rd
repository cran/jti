% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/api_jt.R
\name{leaves}
\alias{leaves}
\alias{leaves.jt}
\alias{parents}
\alias{parents.jt}
\title{Query Parents or Leaves in a Junction Tree}
\usage{
leaves(jt)

\method{leaves}{jt}(jt)

parents(jt)

\method{parents}{jt}(jt)
}
\arguments{
\item{jt}{A junction tree object, \code{jt}.}
}
\description{
Return the clique indices of current parents or leaves
in a junction tree
}
\examples{
# See example 6 in the help page for the jt function
}
\seealso{
\code{\link{jt}}, \code{\link{get_cliques}}
}
