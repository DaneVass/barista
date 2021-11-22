#' Sample to sample distances
#'
#' Plot sample distances between barcode sets / samples
#'
#'
#' @param counts.obj matrix/dataframe containing raw or normalised counts
#'
#'
#' @return Returns a heatmap of sample distances using desired clustering
#' @export
#' @examples
#' test.mat <- matrix(rnorm(10*10,mean=0,sd=2), 10, 10)
#' plotDistances(test.mat)
#'
plotDistances <- function(counts.obj,
                          method = "euclidean",
                          upper = T,
                          name = "Sample Distances"){

  # generate correlation matrix
  sample.distances <- stats::dist(base::t(counts), method = method)
  sample.matrix <- base::as.matrix(sample.distances)

  # take only upper triangle if specified
  if(base::isTRUE(upper)){
    sample.matrix <- barcodeR::get_upper_tri(sample.matrix)
  }

  # Melt the correlation matrix
  melted_dist <- reshape2::melt(sample.matrix, na.rm = TRUE)

  # Create a ggheatmap
  name <- paste(name, "-", method)
  ggheatmap <- ggplot2::ggplot(melted_dist, aes(Var2, Var1, fill = value)) +
    ggplot2::geom_tile(color = "white") +
    viridis::scale_fill_viridis(option = "inferno") +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 1, size = 6, hjust = 1),
                   axis.text.y = ggplot2::element_text(size = 6)) +
    ggplot2::coord_fixed() +
    ggplot2::labs(x = "", y = "", title = name)

  # Print the heatmap
  print(ggheatmap)
}

#' Sample to sample correlation
#'
#' Plot sample correlation between barcode sets / samples
#'
#'
#' @param counts.obj matrix/dataframe containing raw or normalised counts
#'
#'
#' @return Returns a heatmap of sample distances using desired clustering
#' @export
#' @examples
#' test.mat <- matrix(rnorm(10*10,mean=0,sd=2), 10, 10)
#' plotCorrelation(test.mat, method = "pearson")
#'
plotCorrelation <- function(counts,
                            method = "pearson",
                            upper = T,
                            clustered = T,
                            name = "Sample Correlation"){

  # generate correlation matrix
  cormat <- base::round(stats::cor(counts, method = method),2)
  melted_cormat <- reshape2::melt(cormat)

  # cluster matrix if specified
  if(base::isTRUE(clustered)){
    cormat <- barcodeR::cluster_cormat(cormat)
  }

  # take only upper triangle if specified
  if(base::isTRUE(upper)){
    cormat <- barcodeR::get_upper_tri(cormat)
  }

  # Melt the correlation matrix
  melted_cormat <- reshape2::melt(cormat, na.rm = TRUE)

  # Create a ggheatmap
  name <- paste(name, "-", method)
  ggheatmap <- ggplot2::ggplot(melted_cormat, aes(Var2, Var1, fill = value)) +
    ggplot2::geom_tile(color = "white") +
    viridis::scale_fill_viridis() +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 1, size = 6, hjust = 1),
                   axis.text.y = ggplot2::element_text(size = 6)) +
    ggplot2::coord_fixed() +
    ggplot2::labs(x = "", y = "", title = name)

  # Print the heatmap
  print(ggheatmap)
}


#' Cluster correlation matrix
#'
#' cluster a correlation matrix using hierarchical clustering
#'
#'
#' @param cormat matrix of correlation values
#'
#'
#' @return Returns a matrix of correlation values with columns and rows hierarchically clustered
#' @export
#' @examples
#' test.mat <- matrix(rnorm(10*10,mean=0,sd=2), 10, 10)
#' cluster_cormat(test.mat)
#'
cluster_cormat <- function(cormat){
  dd <- stats::as.dist((1-cormat)/2)
  hc <- stats::hclust(dd)
  cormat <- cormat[hc$order, hc$order]
  return(cormat)
}

#' Get upper triangle
#'
#' change lower triangle values in a mirrored matrix to NA
#'
#'
#' @param cormat matrix of correlation values
#'
#'
#' @return Returns a matrix of correlation values with lower triangle values changed to NA
#' @export
#' @examples
#' test.mat <- matrix(rnorm(10*10,mean=0,sd=2), 10, 10)
#' get_upper_tri(test.mat)
#'
get_upper_tri <- function(cormat){
  cormat[base::lower.tri(cormat)] <- NA
  return(cormat)
}
