#' Draws heat map with cluster probabilities
#'
#' @param result an object produced by 'cluster_prob'
#'
#' @param low_col lower expression bound
#' @param high_col upper expression bound
#'
#' @return an object of class `HeatmapList` (see \link[pck=ComplexHeatmap]{HeatmapList-class})
#'
#' @seealso
#' \link[pck=ComplexHeatmap]{HeatmapList-class}
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
#' ## cluster_prob
#' prob_result <- cluster_prob(dat, hdat, nclust = 3, nfake = 10, cutoff = 21)
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
#' @export
heatmap_prob = function(result, low_col = -1, high_col = 1){

  # colors ====================================================================
  ht_color = circlize::colorRamp2(c(low_col, 0, high_col), c("blue", "white", "red")) #colors for gene expression
  prob_colors = RColorBrewer::brewer.pal(result$nclust, 'Set3') # potential problem: max 12 disting colors
  prob_col_list = vector("list", result$nclust)
  for(i in 1:result$nclust){
    prob_col_list[[i]] = circlize::colorRamp2(c(0, 1), c("white", prob_colors[i]))
  }

  # new class imp ========================================================
  shifted_imp  = as.matrix(result$class_imp)
  if(min(shifted_imp) < 0){
    shifted_imp = shifted_imp + abs(min(shifted_imp))
    print("shift = TRUE")
  }

  w = 30  #heatmap width

  # Draw heat maps
  ha = ComplexHeatmap::rowAnnotation(#Text = anno_text(rownames(shifted_imp)),
    Importance = ComplexHeatmap::anno_barplot(shifted_imp,
                                              axis_param = list(direction = "reverse"),
                                              gp = grid::gpar(col = prob_colors, fill = prob_colors),
                                              height = grid::unit(20, "cm"), width = grid::unit(4, "cm")
    )
  )

  ht = ComplexHeatmap::Heatmap(t(result$select_dat), name = "Gene Expression", column_title = "Cluster_prob Heat Map" ,
                               col = ht_color, left_annotation = ha,
                               cluster_columns = FALSE, column_order = result$hdat$order,
                               clustering_distance_rows = "euclidean", clustering_method_rows = "ward.D2",
                               show_row_dend = FALSE, show_row_names = TRUE, row_names_side = "right",
                               show_column_names = FALSE,
                               width = grid::unit(w, "cm"), height = grid::unit(20, "cm"))
  ht_list = ht

  for (i in 1:result$nclust) {
    temp_ht = ComplexHeatmap::Heatmap(t(result$probability[i,]), name = paste("Cluster ", i),
                                      col = prob_col_list[[i]],
                                      cluster_columns = FALSE, show_column_names = FALSE,
                                      width = grid::unit(w, "cm"), height = grid::unit(0.5, "cm"))
    ht_list = ht_list %v% temp_ht
  }

  ComplexHeatmap::draw(ht_list)
}
