## ---- include = FALSE----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup---------------------------------------------------------------
library(clusterProb)

## ----data, echo=TRUE-----------------------------------------------------
# Read data
dat  <- clusterProb::pan_cancer_data
num_dat <- dat[ , 7:187]
rownames(num_dat) <- dat$TCGA_ID

## ----clustering, echo=TRUE-----------------------------------------------
# Clustering
cor_dat <- cor(t(num_dat), method = "spearman")
cor_dat = 1 - cor_dat
cor_dat = as.dist(cor_dat)
hdat <- hclust(cor_dat, method = "ward.D2")

## ----clusterProb, echo=TRUE----------------------------------------------
prob_result <- cluster_prob(num_dat, hdat, 8, 500, 50)
heatmap_prob(prob_result, -2.5, 2.5)

