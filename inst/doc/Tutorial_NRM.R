## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
library(ghypernet)
library(texreg, quietly = TRUE) # for regression tables
library(ggplot2) # for plotting
library(ggraph) #for network plots using ggplot2

## ---- eval=FALSE, echo=FALSE--------------------------------------------------
#  data("swissParliament_network", package = "ghypernet")

## -----------------------------------------------------------------------------
el <- adj2el(cospons_mat, directed = TRUE)

## -----------------------------------------------------------------------------
summary(el$edgecount)

## -----------------------------------------------------------------------------
nodes <- colnames(cospons_mat)
adj_mat <- el2adj(el, nodes = nodes)

## -----------------------------------------------------------------------------
identical(cospons_mat, adj_mat)

## -----------------------------------------------------------------------------
identical(rownames(cospons_mat), dt$idMP)

## -----------------------------------------------------------------------------
dt_unsorted <- dt[order(dt$firstName),]
identical(rownames(cospons_mat), dt_unsorted$idMP)

## -----------------------------------------------------------------------------
dtsorted <- data.frame(idMP = rownames(cospons_mat)) 
dtsorted <- dplyr::left_join(dtsorted, dt_unsorted, by = "idMP") 
identical(dt$idMP, dtsorted$idMP)

## -----------------------------------------------------------------------------
dim(cospons_mat) == dim(onlinesim_mat)

## -----------------------------------------------------------------------------
table(rownames(cospons_mat) == rownames(onlinesim_mat))

## -----------------------------------------------------------------------------
recip_cospons <- reciprocity_stat(cospons_mat)
recip_cospons[1:5, 1:3]

## ---- out.width='100%', fig.align='center', fig.cap='Figure 1: Triadic closure: (a) undirected triangle, (b) transitive triplet, (c) edge-wise shared partners. Source: Brandenberger et al. [-@brandenberger2019quantifying].', echo=FALSE----
knitr::include_graphics('images/tikz_nweffects.pdf')

## -----------------------------------------------------------------------------
shp_cospons_unweighted <- sharedPartner_stat(cospons_mat, directed = TRUE, weighted = FALSE)
shp_cospons_unweighted[1:5, 1:3]
shp_cospons_weighted <- sharedPartner_stat(cospons_mat, directed = TRUE)
shp_cospons_weighted[1:5, 1:3]

## -----------------------------------------------------------------------------
shp_cospons_incoming <- sharedPartner_stat(cospons_mat, directed = TRUE,
                                          triad.type = 'directed.incoming')
shp_cospons_incoming[1:5, 1:3]
shp_cospons_outgoing <- sharedPartner_stat(cospons_mat, directed = TRUE,
                                          triad.type = 'directed.outgoing')
shp_cospons_outgoing[1:5, 1:3]

## -----------------------------------------------------------------------------
canton_homophilymat <- homophily_stat(dt$canton, type = 'categorical', 
                                      nodes = dt$idMP)
canton_homophilymat[1:5, 1:3]

## -----------------------------------------------------------------------------
canton_BE_homophilymat <- homophily_stat(dt$canton, type = 'categorical', 
                                      nodes = dt$idMP, these.categories.only = 'Bern')

## -----------------------------------------------------------------------------
canton_BEZH_homophilymat <- homophily_stat(dt$canton, type = 'categorical', 
                                      nodes = dt$idMP, 
                                      these.categories.only = c('Bern', 'Zuerich'))

## -----------------------------------------------------------------------------
party_homophilymat <- homophily_stat(dt$party, type = 'categorical', nodes = dt$idMP)
parlgroup_homophilymat <- homophily_stat(dt$parlGroup, type = 'categorical', nodes = dt$idMP)
gender_homophilymat <- homophily_stat(dt$gender, type = 'categorical', nodes = dt$idMP)

## -----------------------------------------------------------------------------
dt$age <- 2019 - as.numeric(format(as.Date(dt$birthdate, format = '%d.%m.%Y'), "%Y"))
age_absdiffmat <- homophily_stat(dt$age,  type = 'absdiff', nodes = dt$idMP)
age_absdiffmat[1:5, 1:3]

## -----------------------------------------------------------------------------
head(dtcommittee)

## -----------------------------------------------------------------------------
## This is just one potential way of accomplishing this! 
identical(as.character(dtcommittee$idMP), rownames(cospons_mat))
shared_committee <- matrix(0, nrow = nrow(cospons_mat), ncol = ncol(cospons_mat))
rownames(shared_committee) <- rownames(cospons_mat)
colnames(shared_committee) <- colnames(cospons_mat)
for(i in 1:nrow(shared_committee)){
  for(j in 1:ncol(shared_committee)){
    committees_i <- unlist(strsplit(as.character(dtcommittee$committee_names[i]), ";"))
    committees_j <- unlist(strsplit(as.character(dtcommittee$committee_names[j]), ";"))
    shared_committee[i, j] <- length(intersect(committees_i, committees_j))
  }
}
shared_committee[1:5, 1:3]

## -----------------------------------------------------------------------------
dt$degree <- rowSums(cospons_mat) + colSums(cospons_mat)
degreemat <- cospons_mat
for(i in 1:nrow(cospons_mat)){
  for(j in 1:ncol(cospons_mat)){
    degreemat[i, j] <- sum(dt$degree[i], dt$degree[j])
  }
}

## -----------------------------------------------------------------------------
age_activity_mat <- matrix(rep(dt$age, ncol(cospons_mat)),
                                nrow = nrow(cospons_mat), byrow = FALSE)
svp_activity_mat <- matrix(rep(dt$party, ncol(cospons_mat)),
                           nrow = nrow(cospons_mat), byrow = FALSE)
svp_activity_mat <- ifelse(svp_activity_mat == 'SVP', exp(1), 1)

## -----------------------------------------------------------------------------
age_popularity_mat <- matrix(rep(dt$age, ncol(cospons_mat)),
                                nrow = nrow(cospons_mat), byrow = TRUE)
svp_popularity_mat <- matrix(rep(dt$party, ncol(cospons_mat)),
                           nrow = nrow(cospons_mat), byrow = TRUE)
svp_popularity_mat <- ifelse(svp_popularity_mat == 'SVP', exp(1), 1)

## -----------------------------------------------------------------------------
recip_cospons <- get_zero_dummy(recip_cospons, name = 'reciprocity')
age_absdiffmat <- get_zero_dummy(age_absdiffmat, name = 'age')
shared_committee <- get_zero_dummy(shared_committee, name = 'committee')

## ---- eval=FALSE,echo=TRUE----------------------------------------------------
#  fit <- nrm(adj = cospons_mat, w = recip_cospons,
#             directed = TRUE, selfloops = FALSE, regular = FALSE)

## -----------------------------------------------------------------------------
nfit1 <- nrm(adj = cospons_mat, 
                      w = list(same_canton = canton_homophilymat), 
                      directed = TRUE)
summary(nfit1)

## -----------------------------------------------------------------------------
nfit1 <- nrm(adj = cospons_mat, 
                      w = list(same_canton = canton_homophilymat), 
                      directed = TRUE,
                      init = c(0.208))
summary(nfit1)

## -----------------------------------------------------------------------------
texreg::screenreg(nfit1)

## -----------------------------------------------------------------------------
nfit2 <- nrm(adj = cospons_mat, 
             w = c(
                   recip_cospons,
                   list(party = party_homophilymat,
                        canton = canton_homophilymat,
                        gender = gender_homophilymat),
                   age_absdiffmat,
                   shared_committee,
                   list(online_similarity = onlinesim_mat)
             ), 
             directed = TRUE,
             init = c(.1,-.9, 1.2, .2, .2, 0, 0,0, -.2,-.1))

## -----------------------------------------------------------------------------
screenreg(nfit2, 
          groups = list('Endogenous' = 1:2, 
                     'Homophily' = c(3:7), 
                     'Exogenous' = c(8:10)))

## -----------------------------------------------------------------------------
nfit3 <- nrm(adj = cospons_mat, 
              w = c(
                get_zero_dummy(degreemat, name = 'degree'),
                     recip_cospons,
                list(party = party_homophilymat,
                     svp_in = svp_popularity_mat, 
                     svp_out = svp_activity_mat,
                     canton = canton_homophilymat, 
                     gender = gender_homophilymat),
                     age_absdiffmat,
                list(agein = age_popularity_mat,
                     ageout = age_activity_mat),
                     shared_committee,
                list(online_similarity = onlinesim_mat)
                ), 
              directed = TRUE, regular = TRUE,
              init = c(1,0,0,0, 0.1, 0.5, 0, 0, .1, 0,0, 0,0, .1, .01))
summary(nfit3)

## -----------------------------------------------------------------------------
screenreg(list(nfit2, nfit3), 
          custom.model.names = c('with degree correction', 'without deg. cor.'))

## -----------------------------------------------------------------------------
fullfit <- ghype(graph = cospons_mat, directed = TRUE, selfloops = FALSE)

## -----------------------------------------------------------------------------
nfit2omega <- data.frame(omega = as.vector(nfit2$omega), 
                         cosponsfull = as.vector(cospons_mat), 
                         age_absdiff = as.vector(age_absdiffmat$age), 
                         sameparty = as.vector(party_homophilymat))
nfit2omega[nfit2omega == 0] <- NA
nfit2omega <- na.omit(nfit2omega)

## ---- fig.height=4, fig.width=7-----------------------------------------------
ggplot(nfit2omega, aes(x = age_absdiff, y = omega, color = factor(sameparty)))+
  geom_point(alpha = .1) +
  geom_smooth() + theme(legend.position = 'bottom') + 
  scale_color_manual("", values = c('#E41A1C', '#377EB8'), labels = c('Between parties', 'Within party'))+
  xlab("Age difference") + ylab("Tie propensities")+
  ggtitle('Model (2): Marginal effects of age difference')

## -----------------------------------------------------------------------------
simnw <- rghype(nsamples = 1, model = nfit2, seed = 1253)

## ---- fig.height=5, fig.width=5-----------------------------------------------
ggraph(graph = simnw, layout = 'stress') +
  geom_edge_link(aes(filter = weight>5, alpha=weight)) +
  geom_node_point(aes(colour = dt$parlGroup), size=10*apply(simnw,1,sum)/max(apply(simnw,1,sum))) +
  scale_colour_manual("", values = c('orange', 'yellow', 'blue', 'green', 'grey',
                                    'darkblue', 'red', 'darkgreen', 'purple')) +
  theme(legend.position = 'bottom') + coord_fixed() + theme_graph()

