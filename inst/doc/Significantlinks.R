## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## -----------------------------------------------------------------------------
library(ghypernet)
library(igraph)
data("adj_karate")
data("vertexlabels")
karate <- graph_from_adjacency_matrix(adjmatrix = adj_karate, mode='undirected', weighted=TRUE)
V(karate)$color <- vertexlabels

## -----------------------------------------------------------------------------
mod <- scm(adj_karate,directed = F, selfloops = F)

blockModel <- bccm(adj = adj_karate, labels = vertexlabels, directed = F, selfloops = F, homophily = F)

## -----------------------------------------------------------------------------
# obtain significance matrix
signmat <- link_significance(adj_karate, mod, under = FALSE)

# filter adjacency matrix
adjfiltered <- adj_karate
adjfiltered[signmat>(1/mod$m)] <- 0
adjfiltered[signmat<(1/mod$m) & adj_karate==0] <- 1
diag(adjfiltered) <- 0

## -----------------------------------------------------------------------------
adjcolor <- adj_karate
adjcolor[adj_karate>0] <- 2
adjcolor[signmat<(1/mod$m)] <- 1
diag(adjcolor) <- 0
gfiltered <- graph_from_adjacency_matrix(adjfiltered, mode = 'upper')
g <-  graph_from_adjacency_matrix(adjcolor, mode = 'upper', weighted = 'color')
E(g)$color[E(g)$color==1] <- "red"
E(g)$color[E(g)$color==2] <- "black"

V(gfiltered)$color <- V(g)$color <- vertexlabels
plot(karate)
plot(g)
plot(gfiltered)

## -----------------------------------------------------------------------------
signmat <- link_significance(adj_karate, blockModel, under=FALSE)

# filter adjacency matrix
adjfiltered <- adj_karate
adjfiltered[signmat>(1/mod$m)] <- 0
adjfiltered[signmat<(1/mod$m) & adj_karate==0] <- 1
diag(adjfiltered) <- 0

## -----------------------------------------------------------------------------
adjcolor <- adj_karate
adjcolor[adj_karate>0] <- 2
adjcolor[signmat<(1/mod$m)] <- 1
diag(adjcolor) <- 0
gfiltered <- graph_from_adjacency_matrix(adjfiltered, mode = 'upper')
g <-  graph_from_adjacency_matrix(adjcolor, mode = 'upper', weighted = 'color')
E(g)$color[E(g)$color==1] <- "red"
E(g)$color[E(g)$color==2] <- "black"

V(gfiltered)$color <- V(g)$color <- V(karate)$color
plot(karate)
plot(g)
plot(gfiltered)

