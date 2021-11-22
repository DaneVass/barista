#' plotScatter
#'
#' Generate a simple scatterplot between two or more sets of sample counts.
#'
#' @param obj DGEList object containing samples to be compared
#' @param samps character vector of samplenames of length 2 or more. more then 2 samples will create a pairwise scatter plot.
#' @param title desired name of output plot
#' @param rug boolean. include geom_rug density information on the axes?
#' @param trendline boolean. include linear trendline using stat_smooth()?
#' @param stats boolean. include model stats in the plot subtitle?
#' @param trans From ggplot2. For continuous scales, the name of a transformation object or the object itself. Built-in transformations include "asn", "atanh", "boxcox", "date", "exp", "hms", "identity", "log", "log10", "log1p", "log2", "logit", "modulus", "probability", "probit", "pseudo_log", "reciprocal", "reverse", "sqrt" and "time".
#'
#' @return Returns a scatterplot between two or more sets of sample counts.
#' @export
#' @examples
#' plotScatter(obj, samps = c("Samp_1", "Samp_2"), title = "Scatter plot", rug = T, trendline = T, stats = T, trans = NULL)

plotRegression <- function(obj = NULL, samps = NULL, title = "Scatter plot", rug = T, trendline = T, stats = T, trans = NULL) {
  if(base::is.null(obj) | class(obj) != "DGEList" | class(obj) != "data.frame" | class(obj) != "matrix"){
    stop("please supply a valid DGEList object, data frame or matrix of counts.")
  }

  if(base::is.null(samps) | length(samps) < 2){
    stop("please supply two or more samplenames from DGEList object")
  }

  if(length(samps) == 2){
    p <- ggplot2::ggplot(counts(), ggplot2::aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) +
      ggplot2::geom_point() +
      ggplot2::theme_bw() +
      ggplot2::ggtitle(title)
  } else {
    p <- ggplot2::ggplot(fit$model, ggplot2::aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) +
      ggplot2::geom_point() +
      ggplot2::theme_bw() +
      ggplot2::ggtitle(title)
  }

  if(!base::is.null(trans)){
    p <- p + ggplot2::scale_x_continuous(trans = trans)
    p <- p + ggplot2::scale_y_continuous(trans = trans)
  }

  if(base::isTRUE(rug)){
    p <- p + ggplot2::geom_rug(col=rgb(.5,0,0,alpha=.2))
  }

  if(base::isTRUE(trendline)){
    p <- p + ggplot2::stat_smooth(method = "lm", col = "blue")
  }

  if(base::isTRUE(stats)){
    p <- p + ggplot2::ggtitle(title, subtitle = base::paste("Adj R2 = ", base::signif(base::summary(fit)$adj.r.squared, 5), " ",
                                                            " Intercept = ", base::signif(fit$coef[[1]], 5), " ",
                                                            " Slope = ", base::signif(fit$coef[[2]], 5), " ",
                                                            " P = ", base::summary(fit)$coef[2,4]))
  }
  return(p)
}
