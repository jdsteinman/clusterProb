#' Confidence for Clustering Results
#'
#' Calculates probabilities for cluster membership using random forest
#'
#' @param data A matrix or matrix-like object. The `cluster_prob` function treats the
#' rows as the observations and the columns as the features.
#' @param hdat An object of class 'hclust'
#' @param nclust number of clusters
#' @param nfake number of fake features
#' @param cutoff maximum number of important features
#'
#' @details This function uses a random forest model to produce a confidence measure for clustering results.
#' trains a random forest model on the data to predict the clustering labels.
#' The model uses the features (columns) from the input data to predict the cluster membership of each observation.
#'
#' @return a list containing the following information:
#'     \item{select_dat}{A matrix containing the data with only important features}
#'     \item{probability}{A matrix containing the cluster probabilities}
#'     \item{global_imp}{global features importances}
#'     \item{class_imp}{class-wise feature importance}
#'     \item{hdat}{clustering data from `hclust`}
#'     \item{nclust}{number of clusters}
#'
#' @examples
#' ### Example 1: Simulated Data
#'
#' ## Create Data
#' dat <- matrix(data = rnorm(20*20), nrow = 20)
#' dat <- dat / max(dat)   # max value = 1
#' rownames(dat) <- letters[1:20]
#' colnames(dat) <- LETTERS[1:20]
#' dat[1:5, ] = dat[1:5, ] + 1    # first cluster
#' dat[5:10, ] = dat[5:10, ] - 1  # second cluster
#'
#' ## Clustering
#' hdat <- hclust(dist(dat), method = "ward.D2")
#'
#' ## cluster_prob
#' prob_result <- cluster_prob(dat, hdat, nclust = 3, nfake = 10, cutoff = 21)
#'
#' ## heatmap_prob
#' ht <- heatmap_prob(prob_result, -2, 2)
#'
#' ### Example 2: Car Models
#' car_dat <- mtcars
#' car_dat <- scale(car_dat)
#' hdat <- hclust(dist(car_dat))
#' prob_result <- cluster_prob(car_dat, hdat, 3, 10, 10)
#'
#' @import ComplexHeatmap RColorBrewer circlize
#'
#' @import caret randomForest
#' @importFrom stats cutree rnorm
#'
#' @export

cluster_prob = function(data, hdat, nclust, nfake = NULL, cutoff = NULL){
  # Format hdat ===================================================================
  cluster_ID <- as.factor(cutree(hdat, nclust))

  # Generate Fake Data ============================================================
  if(is.null(nfake) == TRUE){
    nfake = round(ncol(data)*0.25)
  }

  fake_dat <- matrix(data = rnorm(nfake*nrow(data)), nrow = nrow(data), ncol = nfake)
  prefix   <- "fake"
  suffix   <- 1:nfake
  fake_names <- paste(prefix, suffix, sep = "_")
  rownames(fake_dat) <- rownames(data)
  colnames(fake_dat) <- fake_names
  new_dat  <- cbind(data, fake_dat) # attach fake data


  # Set model control =============================================================
  my_control = caret::trainControl(
    method = "cv", ## cross validation
    number = 5,   ## 5-fold
    summaryFunction = multiClassSummary,
    verboseIter = FALSE)

  # Train Model ===================================================================
  model = caret::train(x = new_dat, y = cluster_ID,
                       method = "rf",
                       importance = TRUE,
                       tuneLength = 15,
                       ntree = 500,
                       trControl = my_control)
  probability = t(model$finalModel$votes)   # stores votes as confidence measure

  # Importance ===================================================================
  global_imp   = as.matrix(model$finalModel$importance[, "MeanDecreaseAccuracy"]) # matrix to save names
  global_order = order(global_imp, decreasing = TRUE)
  global_imp   = as.matrix(global_imp[global_order, ])
  fake_index   = vector(length = nfake)

  for(i in 1:nfake){
    fake_index[i] = which(rownames(global_imp) == paste(prefix, i, sep = "_"))
  }
  fake_index = sort(fake_index, decreasing = FALSE)

  # Feature Selection ============================================================
  if(is.null(cutoff == TRUE)){
    if(ncol(data) <= 50){
      cutoff = ncol(data)
    } else{
      cutoff = 50
    }
  }

  if(min(fake_index) <= cutoff){
    cut = min(fake_index) - 1
    select_dat = new_dat[ , global_order[1:cut]]
    global_imp = global_imp[global_order[1:cut]]
    class_imp  = as.matrix(model$finalModel$importance[global_order[1:cut], 1:nclust])
  } else {
    cut = cutoff
    select_dat = new_dat[ , global_order[1:cut]]
    global_imp = global_imp[global_order[1:cut]]
    class_imp  = as.matrix(model$finalModel$importance[global_order[1:cut], 1:nclust])
  }
  cluster_labels = paste("cluster", 1:nclust, sep = "_")
  colnames(class_imp) = cluster_labels

  # Output =======================================================================
  output = list(select_dat  = select_dat,
                probability = probability,
                global_imp  = global_imp,
                class_imp   = class_imp,
                hdat        = hdat,
                nclust      = nclust)
  return(output)
}
