## ---- include = FALSE----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup---------------------------------------------------------------
library(clusterProb)

## ----data, echo=TRUE-----------------------------------------------------
dat <- clusterProb::unc_breast_data
gene_dat <- as.matrix(dat[7:1671 , 5:341])
gene_dat <- apply(gene_dat, 2, as.numeric)
rownames(gene_dat) <- as.character(dat[7:1671, 3])

gene_dat = scale(gene_dat)
t_dat = t(gene_dat)

## ----clustering, echo=TRUE-----------------------------------------------
hdat = hclust(dist(t_dat), method = "ward.D2")

## ----cluster_prob, echo=TRUE---------------------------------------------
prob_result <- cluster_prob(t_dat, hdat, 5, 200, 50)
heatmap_prob(prob_result, -2, 2)

