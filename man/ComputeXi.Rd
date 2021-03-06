% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/computexi.R
\name{ComputeXi}
\alias{ComputeXi}
\title{Auxiliary function. Computes combinatorial matrix.}
\usage{
ComputeXi(adj, directed, selfloops, regular = FALSE)
}
\arguments{
\item{adj}{adjacency matrix}

\item{directed}{boolean, whether the model is for a directed network}

\item{selfloops}{boolean, whether the model contains selfloops}

\item{regular}{boolean. Is the combinatorial matrix computed for configuration model or for regular gnp model? default FALSE.}
}
\value{
combinatorial matrix
}
\description{
Combinatorial matrix computed according to soft
configuration model or 'regular' gnp model.
}
\examples{
data('adj_karate')
xi = ComputeXi(adj_karate, directed = FALSE, selfloops = FALSE)

}
