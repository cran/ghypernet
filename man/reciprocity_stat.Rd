% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/endogenousstats.R
\name{reciprocity_stat}
\alias{reciprocity_stat}
\title{Calculate weighted reciprocity change statistics for multi-edge graphs.}
\usage{
reciprocity_stat(graph, nodes = NULL, zero_values = NULL)
}
\arguments{
\item{graph}{A graph adjacency matrix or an edgelist.  The edgelist needs to
have 3 columns: a sender vector, a target vector and an edgecount vector.}

\item{nodes}{optional character/factor vector. If an edgelist is provided,
you have to provide a list of unique identifiers of your nodes in the graph.
This is because in the edgelist, isolates are usually not recorded.
If you do not specify isolates in your nodes object, they are excluded
from the analysis (falsifies data).}

\item{zero_values}{optional numeric value. Use this to substitute zero-values
in your reciprocity change statistic matrix. Zero values in the predictors
are recognized in the gHypEG regression as structural zeros. To ensure this
does not happen, please recode your zero-values in all your predictors,
ideally using a dummy variable fitting an optimal value for the zeroes.}
}
\value{
Reciprocity change statistic matrix.
}
\description{
The function takes either an edgelist or an adjacency matrix and returns an
adjacency matrix with the reciprocity change statistic. This reciprocity
matrix can then be used as a predictor in the gHypEG regression.
}
\examples{
recip_stat <- reciprocity_stat(adj_karate)
recip_stat_dummy <- get_zero_dummy(recip_stat, name = 'reciprocity')
nrm(w=recip_stat_dummy, adj_karate, directed = FALSE, selfloops = FALSE)

}
\seealso{
\code{\link{sharedPartner_stat}} or \code{\link{homophily_stat}}
}
\author{
LB, GC
}
