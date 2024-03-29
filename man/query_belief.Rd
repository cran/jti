% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/api_query.R
\name{query_belief}
\alias{query_belief}
\alias{query_belief.jt}
\title{Query probabilities}
\usage{
query_belief(x, nodes, type = "marginal")

\method{query_belief}{jt}(x, nodes, type = "marginal")
}
\arguments{
\item{x}{A junction tree object, \code{jt}.}

\item{nodes}{The nodes for which the probability is desired}

\item{type}{Either 'marginal' or 'joint'}
}
\description{
Get probabilities from a junction tree object
}
\examples{
# See the 'jt' function
}
\seealso{
\code{\link{jt}}, \code{\link{mpe}}
}
