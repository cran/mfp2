#' Plot response functions from a fitted `mfp2` object
#'
#' Plots the partial linear predictors with confidence limits 
#' against the selected covariate(s) of interest. 
#'
#' @param model fitted `mfp2` model.
#' @param terms character vector with variable names to be plotted.
#' @param partial_only a logical value indicating whether only the partial 
#' predictor (component) is drawn (`TRUE`), or also component-plus-residual 
#' (`FALSE`, the default). Only used if `type = "terms"`. See below for details. 
#' @param type,ref,terms_seq arguments of [predict.mfp2()]. Only 
#' `type = "terms"` and `type = "contrasts"` are supported by this function. 
#' @param alpha `alpha` argument of [predict.mfp2()].
#' @param shape,size_points,color_points `ggplot2` properties of drawn 
#' data points.
#' @param color_line,linetype,linewidth `ggplot2` properties of line for 
#' partial predictor.
#' @param color_fill,alpha_fill `ggplot2` properties of ribbon for confidence
#' interval.
#' @param ... used in alias `plot_mfp` to pass arguments. 
#' 
#' @details 
#' The confidence limits of the partial linear predictors or contrasts are
#' obtained from the variance–covariance matrix of the final fitted model, 
#' which takes into account the uncertainty in estimating the model parameters 
#' but not the FP powers. This can lead to narrow confidence intervals. A simple
#' way to obtain more realistic confidence intervals within the FP is by using
#' bootstrap, which is not currently implemented. See Royston and Sauerbrei (2008)
#' chapter 4.9.2 for guidance on conducting bootstrapping within 
#' the FP class.
#' 
#' The component-plus-residual, is the partial linear predictor plus residuals, 
#' where deviance residuals are used in generalized linear regression models, 
#' while martingale residuals are used in Cox models, as done in Stata mfp program.
#' This kind of plot is only available if `type = "terms"`.
#' @examples
#'
#' # Gaussian response
#' data("prostate")
#' x = as.matrix(prostate[,2:8])
#' y = as.numeric(prostate$lpsa)
#' # default interface
#' fit = mfp2(x, y, verbose = FALSE)
#' fracplot(fit) # generate plots
#'
#' @return 
#' A list of `ggplot2` plot objects, one for each term requested. Can be 
#' drawn as individual plots or facetted / combined easily using e.g. 
#' `patchwork::wrap_plots` and further customized. 
#' 
#' @seealso 
#' [predict.mfp2()]
#' 
#' @import ggplot2
#' @importFrom ggplot2 .data
#' @export
fracplot <- function(model, 
                     terms = NULL, 
                     partial_only = FALSE, 
                     type = c("terms","contrasts"),
                     ref = NULL,
                     terms_seq = c("data", "equidistant"),
                     alpha = 0.05,
                     color_points = "#AAAAAA",
                     color_line = "#000000", 
                     color_fill = "#000000",
                     shape = 1,
                     size_points = 1,
                     linetype = "solid", 
                     linewidth = 1,
                     alpha_fill = 0.1) {
  
  # Does a check if ggplot2 is available
  # It should be as it is in the imports section but in CRAN checks some systems don't have it!
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package \"ggplot2\" needed for this function to work. Please install it.",
         call. = FALSE
    )
  }
  # assert that the object must be mfp2
  if (!inherits(model, "mfp2")) 
    stop("The model entered is not an mfp2 object.", call. = FALSE)
  
 # check if ref is a list (only if it is not NULL)
  if (!is.null(ref) && !is.list(ref))
    stop("ref must be a list", call. = FALSE)
  
 # set defaults depending on type
  type <- match.arg(type)
  if (type == "contrasts") {
    partial_only <- TRUE
    if (missing(terms_seq))
      terms_seq <- "equidistant"
  }
  
  terms_seq <- match.arg(terms_seq)
  
  pred <- predict(model, 
                  type = type, 
                  terms = terms,
                  ref = ref,
                  terms_seq = terms_seq,
                  alpha = alpha)
  # for points also need the data point predictions
  pred_data <- pred
  if (!partial_only && terms_seq != "data") {
    pred_data <- predict(model,
                         type = "terms", 
                         terms = terms, 
                         terms_seq = "data")
    }
  
  if (length(pred) == 0) {
    warning("i Variables specified in terms not used in final model.")
    return()
  }
  
  ylab <- "Partial Predictor"
  
  if (!partial_only) {
    # compute residuals to plot points
    # for glm, deviance residuals are required
    # while for cox martingale residuals
    resid <- if (model$family_string == "cox") {
      model$residuals
    } else {
      residuals.glm(model, type = "deviance")
    }
    # add residuals to the data
    #pred_data <- lapply(pred_data, function(v) transform(v, resid = resid))
    pred_data <- lapply(pred_data, function(v) {
      v$resid <- resid
    v})
    
    # y label for the plot
    ylab <- "Partial Predictor + residuals"
  }
  
  # Preallocate the list for plots
  plots <- setNames(vector("list", length(names(pred))), names(pred))  
  
  # in the calls to ggplot2::aes use the .data pronoun to avoid notes 
  # generated by R CMD CHECK about missing bindings of global variables
  for (v in names(pred)) {
    plots[[v]] <- ggplot2::ggplot(data = pred[[v]], 
                         ggplot2::aes(x = .data$variable, 
                                      y = .data$value)) + 
      ggplot2::geom_line(linewidth = linewidth,
                         linetype = linetype,
                         color = color_line) + 
      ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$lower, 
                                        ymax = .data$upper), 
                           alpha = alpha_fill) + 
      ggplot2::ggtitle(sprintf("FP%s(%s)%s", 
                               ifelse(model$fp_terms[v, "acd"], "1", ""),
                               paste0(model$fp_powers[[v]], collapse = ", "),
                               ifelse(model$fp_terms[v, "acd"], ":ACD", "")
      )) +
      ggplot2::xlab(v) + ggplot2::ylab(ylab) +
      ggplot2::theme_bw()
    
    if (!partial_only) 
      plots[[v]] <- plots[[v]] + 
        ggplot2::geom_point(data = pred_data[[v]], 
                            ggplot2::aes(y = .data$value + .data$resid),
                            color = color_points,
                            size = size_points, 
                            shape = shape)
  }
  
  plots
}

#' @describeIn fracplot Alias for fracplot.
plot_mfp <- function(...) {
    fracplot(...)
}
