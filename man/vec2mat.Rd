% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/matrixoperations.R
\name{vec2mat}
\alias{vec2mat}
\title{Auxiliary function, produces matrix from vector}
\usage{
vec2mat(vec, directed, selfloops, n)
}
\arguments{
\item{vec}{vector to be put in matrix form}

\item{directed}{a boolean argument specifying whether object is directed or not.}

\item{selfloops}{a boolean argument specifying whether the model should
incorporate selfloops.}

\item{n}{vector. if length(n)==1, n is the number of vertices. If length(n)==3
first element is number of vertices, second and third elements are number of
vertices for row and column of bipartite matrix.}
}
\value{
matrix nxn generated from vector.
}
\description{
The number of elements of vec are the number of non-zero elements in the
adjacency matrix.
It performs the opposite operation of `mat2vec.ix`.
}
\examples{
data('adj_karate')
ix <- mat2vec.ix(adj_karate, FALSE, FALSE)
vec <- adj_karate[ix]
vec2mat(vec, FALSE, FALSE, nrow(adj_karate))

}
