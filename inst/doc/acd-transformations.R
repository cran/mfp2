## ----setup, include=FALSE-------------------------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  message = FALSE,
  warning = FALSE,
  fig.width = 7,
  fig.height = 5,
  out.width = "100%",
  fig.align = "center",
  dev = "png",
  dpi = 96
)

oldwidth <- getOption("width")
options(width = 100)

library(mfp2)
has_ggplot2 <- requireNamespace("ggplot2", quietly = TRUE)
has_geepack <- requireNamespace("geepack", quietly = TRUE)
has_patchwork <- requireNamespace("patchwork", quietly = TRUE)
library(survival)

## ----acd-model-table, echo=FALSE, results='asis'--------------------------------------------------
acd_models <- data.frame(
  Model = c("M1", "M2", "M3", "M4", "M5", "M6"),
  Notation = c(
    "$FP1(p_1, p_2)$",
    "$FP1(p_1, \\cdot)$",
    "$FP1(\\cdot, p_2)$",
    "$FP1(1, \\cdot)$",
    "$FP1(\\cdot, 1)$",
    "$FP1(\\cdot, \\cdot)$"
  ),
  Function = c(
    "$\\beta_1 h_{p_1}(x+s) + \\beta_2 h_{p_2}(A_x)$",
    "$\\beta_1 h_{p_1}(x+s)$",
    "$\\beta_2 h_{p_2}(A_x)$",
    "$\\beta_1(x+s)$",
    "$\\beta_2 A_x$",
    "Null model"
  ),
  Interpretation = c(
    "Most complex allowed ACD-extended function",
    "Standard FP1 function of x",
    "FP1 function of ACD(x)",
    "Linear function of x",
    "Linear function of ACD(x)",
    "x omitted"
  ),
  check.names = FALSE
)

knitr::kable(
  acd_models,
  align = c("l", "l", "l", "l"),
  escape = FALSE,
  caption = "Six candidate model families considered by the FSPA."
)

## ----simulate-acd---------------------------------------------------------------------------------
set.seed(54)

n <- 500
x <- matrix(rnorm(n), ncol = 1, dimnames = list(NULL, "x"))

sigmoid <- function(x) {
  1 / (1 + exp(-1.7 * x))
}

y <- as.numeric(20 * sigmoid(x) + rnorm(n, mean = 0, sd = 0.5))

## ----fit-acd-matrix-------------------------------------------------------------------------------
fit_acd_matrix <- mfp2(
  x, y,
  acd_vars = "x",
  verbose = FALSE
)

fit_acd_matrix$fp_terms

## ----fit-acd-formula------------------------------------------------------------------------------
dat <- data.frame(y = y, x = as.numeric(x))

fit_acd_formula <- mfp2(
  y ~ fp(x, acd = TRUE),
  data = dat,
  verbose = FALSE
)

fit_acd_formula$fp_terms

## ----fit-standard-mfp-----------------------------------------------------------------------------
fit_standard <- mfp2(
  y ~ fp(x, acd = FALSE),
  data = dat,
  verbose = FALSE
)

fit_standard$fp_terms

## ----plot-acd, eval=has_patchwork && has_ggplot2, fig.cap="ACD model compared with a standard MFP model.", fig.width=8, fig.height=4----
p1 <- plot(fit_acd_matrix, terms = "x")[[1]] +
  ggplot2::ylab("Partial predictor + residuals") +
  ggplot2::ggtitle("ACD model")

p2 <- plot(fit_standard, terms = "x")[[1]] +
  ggplot2::ylab("Partial predictor + residuals") +
  ggplot2::ggtitle("Standard MFP model")

patchwork::wrap_plots(p1, p2, ncol = 2)

## ----cleanup, include=FALSE---------------------------------------------------
options(width = oldwidth)

