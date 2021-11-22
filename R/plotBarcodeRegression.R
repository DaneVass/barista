#' plotRegression
#'
#' Generate a linear regression scatterplot from a model fit using lm() between two sets of sample counts.
#'
#' @param dge DGEList object containing grouping variable to fit linear model
#' @param group sample group in dge$samples$group to plot correlation for. must be paired
#' @param fit fit of linear model between two samples generated using lm
#' @param title desired name of output plot
#' @param rug boolean. include geom_rug density information on the axes?
#' @param trendline boolean. include linear trendline using stat_smooth()?
#' @param stats boolean. include model stats in the plot subtitle?
#' @param trans From ggplot2. For continuous scales, the name of a transformation object or the object itself. Built-in transformations include "asn", "atanh", "boxcox", "date", "exp", "hms", "identity", "log", "log10", "log1p", "log2", "logit", "modulus", "probability", "probit", "pseudo_log", "reciprocal", "reverse", "sqrt" and "time".
#'
#' @return Returns a scatterplot represented by proportion of total pool
#' @export
#' @examples
#' plotRegression(fit, title = "Proportional Bubble Plot", rug = T)

plotRegression <- function(dge, group, title = paste("Regression plot:", group), rug = T, trendline = T, stats = T, trans = NULL) {
  dat <- as.data.frame(dge$counts[,dge$samples$group == as.character(group)])
  if( dim(dat)[2] != 2){
    print("must be two replicates for group.")
    print(paste("only found", dim(dat)[2]))
    stop()
  }

  # generate lm fit
  fit <- lm(dat[,1] ~ dat[,2])


  # plot fit data
  p <- ggplot2::ggplot(fit$model, ggplot2::aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) +
    ggplot2::geom_point() +
    ggplot2::theme_bw() +
    ggplot2::ggtitle(title) +
    labs(y = colnames(dat)[1], x = colnames(dat)[2])

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
