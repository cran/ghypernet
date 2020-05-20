## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
library(ghypernet)

data("adj_karate")
directed <- FALSE
selfloops <- FALSE
print(adj_karate[1:10,1:10])

## -----------------------------------------------------------------------------
(regularmodel <- regularm(graph = adj_karate, directed = directed, selfloops = selfloops))

## -----------------------------------------------------------------------------
(confmodel <- scm(graph = adj_karate, directed = directed, selfloops = selfloops))

## -----------------------------------------------------------------------------
random_graphs_rm <- rghype(nsamples = 10, model = regularmodel, seed = 123)
random_graphs_scm <- rghype(nsamples = 10, model = confmodel, seed = 123)

## -----------------------------------------------------------------------------
AIC(regularmodel)
AIC(confmodel)

# difference in AICs, high value means more complex model is better
AIC(regularmodel) - AIC(confmodel)

## -----------------------------------------------------------------------------
# Generate regular models and configuration models for random graphs
regularmodels <- lapply(X = random_graphs_rm, FUN = regularm, directed=directed, selfloops=selfloops)
confmodels <- lapply(X = random_graphs_rm, FUN = scm, directed=directed, selfloops=selfloops)
# Compute AICs
AIC_regularmodels <- sapply(X = regularmodels,FUN = AIC)
AIC_confmodels <- sapply(X = confmodels,FUN = AIC)
# differences in AIC, high value means more complex model is better
AIC_regularmodels - AIC_confmodels

## -----------------------------------------------------------------------------
conf.test(graph = adj_karate, directed = directed, selfloops = selfloops, seed=123)

## -----------------------------------------------------------------------------
tests <- lapply(X = random_graphs_rm, FUN = conf.test, directed = directed, selfloops = selfloops, seed=123)
sapply(X = tests, FUN = function(x) x$p.value)

## -----------------------------------------------------------------------------
data("vertexlabels")
(blockmodel <- bccm(adj = adj_karate, labels = vertexlabels, directed = directed, selfloops = selfloops))
print(blockmodel$blockOmega)

## -----------------------------------------------------------------------------
(blockmodel_2 <- bccm(adj = adj_karate, labels = vertexlabels, directed = directed, selfloops = selfloops, homophily = TRUE))

## -----------------------------------------------------------------------------
lr.test(nullmodel = confmodel, altmodel = blockmodel, seed = 123)
lr.test(nullmodel = confmodel, altmodel = blockmodel_2, seed = 123)

## -----------------------------------------------------------------------------
AIC(confmodel)
AIC(blockmodel_2)
AIC(blockmodel)

## -----------------------------------------------------------------------------
lr.test(nullmodel = blockmodel_2, altmodel = blockmodel, seed=123)

## -----------------------------------------------------------------------------
# First generate random sample from blockmodel
random_graphs_bccm2 <- rghype(nsamples = 100, model = blockmodel_2, seed = 123)
# Generate the two models for random graphs
blockmodels <- lapply(X = random_graphs_bccm2, FUN = bccm, labels = vertexlabels, directed=directed, selfloops=selfloops)
blockmodel_2s <- lapply(X = random_graphs_bccm2, FUN = bccm, labels = vertexlabels, directed=directed, selfloops=selfloops, homophily = TRUE)
# Compute AICs
AIC_blockmodels <- sapply(X = blockmodels, FUN = AIC)
AIC_blockmodel_2s <- sapply(X = blockmodel_2s, FUN = AIC)
# mean difference in AIC, high value means more complex model is better
summary(AIC_blockmodel_2s - AIC_blockmodels)


## -----------------------------------------------------------------------------
fullmodel <- ghype(graph = adj_karate, directed = directed, selfloops = selfloops, unbiased = FALSE)

lr.test(nullmodel = blockmodel_2, altmodel = fullmodel, seed = 123)
gof.test(model = blockmodel_2, seed = 123)

