## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 6,
  fig.height = 4.5,
  fig.align = "center",
  out.width = "100%"
)
options(width = 80)
has_ggplot2 <- requireNamespace("ggplot2", quietly = TRUE)
has_geepack <- requireNamespace("geepack", quietly = TRUE)

## ----packages-and-data--------------------------------------------------------
library(mfp2)

data("prostate", package = "mfp2")

## ----prostate-zero-summary----------------------------------------------------
nrow(prostate)
table(pgg45_zero = prostate$pgg45 == 0)
round(mean(prostate$pgg45 == 0), 2)
length(unique(prostate$pgg45[prostate$pgg45 > 0]))  # distinct positive values

## ----univariable-saz-fit------------------------------------------------------
fit_pgg45 <- mfp2(
  lpsa ~ fp(pgg45, df = 4, select = 0.05, alpha = 0.05, spike = TRUE),
  data = prostate,
  family = "gaussian",
  criterion = "pvalue",
  verbose = FALSE
)

print(fit_pgg45, detailed_settings = FALSE, notes = FALSE)

## ----univariable-saz-predictions----------------------------------------------
newdat <- data.frame(pgg45 = c(0, 4, 20, 45, 60, 100))
predict(fit_pgg45, newdata = newdat, type = "link")

## ----univariable-saz-plot, eval=has_ggplot2, fig.cap="Binary-only SAZ model for pgg45 (paper, Figure 1)."----
p1 <- plot(fit_pgg45, terms = "pgg45", partial_only = TRUE)
p1[["pgg45"]]

## ----select-one---------------------------------------------------------------
fit_forced <- mfp2(
  lpsa ~ fp(pgg45, select = 1, alpha = 0.05, spike = TRUE),
  data = prostate,
  family = gaussian(),
  criterion = "pvalue",
  verbose = FALSE
)
print(fit_forced, detailed_settings = FALSE, notes = FALSE)

## ----keep---------------------------------------------------------------------
fit_keep <- mfp2(
  lpsa ~ fp(pgg45, spike = TRUE) + fp(age) +
    fp(cavol) + fp(weight) + fp(bph) + fp(cp) + svi,
  data = prostate,
  family = gaussian(),
  criterion = "pvalue",
  select = 0.05,
  alpha = 0.05,
  keep = "pgg45",
  verbose = FALSE
)
print(fit_keep, detailed_settings = FALSE, notes = FALSE)

## ----force-max-fp-------------------------------------------------------------
fit_maxfp <- mfp2(
  lpsa ~ fp(pgg45, spike = TRUE, force_max_fp = TRUE),
  data = prostate,
  family = gaussian(),
  criterion = "pvalue",
  verbose = FALSE
)
print(fit_maxfp, detailed_settings = FALSE, notes = FALSE)

## ----multivariable-saz-fit----------------------------------------------------
fit_multi <- mfp2(
  lpsa ~ fp(pgg45, spike = TRUE) + fp(age) +
    fp(cavol) + fp(weight) + fp(bph) + fp(cp) + svi,
  data = prostate,
  family = gaussian(),
  criterion = "pvalue",
  select = 0.05,
  alpha = 0.05,
  verbose = FALSE
)
print(fit_multi, detailed_settings = FALSE, notes = FALSE)

## ----aic-saz-fit--------------------------------------------------------------
fit_pgg45_aic <- mfp2(
  lpsa ~ fp(pgg45, spike = TRUE),
  data = prostate,
  family = "gaussian",
  criterion = "aic",
  verbose = FALSE
)
print(fit_pgg45_aic, detailed_settings = FALSE, notes = FALSE)

## ----aic-saz-predictions------------------------------------------------------
newdat <- data.frame(pgg45 = c(0, 4, 25, 100))
predict(fit_pgg45_aic, newdata = newdat, type = "link")

## ----aic-saz-plot, eval=has_ggplot2, fig.height=5, fig.cap="Continuous-only SAZ model for pgg45 selected by AIC (paper, Figure 2)."----
p2 <- plot(fit_pgg45_aic, terms = "pgg45", partial_only = TRUE,
           show_grid = TRUE)
p2[["pgg45"]]

## ----bic-saz-fit--------------------------------------------------------------
fit_pgg45_bic <- mfp2(
  lpsa ~ fp(pgg45, spike = TRUE),
  data = prostate,
  family = "gaussian",
  criterion = "bic",
  verbose = FALSE
)
print(fit_pgg45_bic, detailed_settings = FALSE, notes = FALSE)

## ----simulate-both-components, eval=FALSE-------------------------------------
# set.seed(2107)
# n <- 5000
# 
# is_zero <- rbinom(n, size = 1, prob = 0.30) == 1
# x_saz <- ifelse(is_zero, 0, runif(n, min = 1, max = 12))
# 
# # True FP2 function for positive values; 0 at x = 0
# positive <- x_saz > 0
# f_saz <- numeric(n)
# f_saz[positive] <- 4.0 * sqrt(x_saz[positive]) - 0.9 * x_saz[positive]
# 
# alpha <- 1.5  # intercept
# gamma <- 2.0  # binary-indicator effect
# z_saz <- as.numeric(x_saz == 0)
# y <- alpha + gamma * z_saz + f_saz + rnorm(n, sd = 1)
# 
# dat_bincont <- data.frame(y = y, x_saz = x_saz)

## ----simulated-saz-fit, eval=FALSE--------------------------------------------
# fit_bincont <- mfp2(
#   y ~ fp(x_saz, center = FALSE, spike = TRUE),
#   data = dat_bincont,
#   family = gaussian(),
#   criterion = "pvalue",
#   min_saz_prop = 0.10,
#   verbose = FALSE
# )
# print(fit_bincont, detailed_settings = FALSE, notes = FALSE)
# 
# predict(fit_bincont, newdata = data.frame(x_saz = c(0, 1, 5, 7, 12)))
# 
# p3 <- plot(fit_bincont, terms = "x_saz", partial_only = TRUE,
#            show_grid = TRUE)
# p3[["x_saz"]]

## ----simulate-multiple-saz, eval=FALSE----------------------------------------
# set.seed(2107)
# n <- 10000
# p <- 6
# rho <- 0.20
# mu <- rep(3, p)
# sigma <- 0.60
# Sigma <- matrix(rho * sigma^2, nrow = p, ncol = p)
# diag(Sigma) <- sigma^2
# 
# repeat {
#   X <- MASS::mvrnorm(n = n, mu = mu, Sigma = Sigma)
#   if (all(X > 0)) break
# }
# colnames(X) <- paste0("x", 1:6)
# x1 <- X[, 1]; x2 <- X[, 2]; x3 <- X[, 3]
# x4 <- X[, 4]; x5 <- X[, 5]; x6 <- X[, 6]
# 
# # Introduce zeros in the three SAZ covariates
# x1[rbinom(n, 1, 0.30) == 1] <- 0
# x2[rbinom(n, 1, 0.25) == 1] <- 0
# x3[rbinom(n, 1, 0.35) == 1] <- 0
# 
# # FP transformation of positive values; 0 at x = 0
# positive_fp <- function(x, power) {
#   out <- numeric(length(x))
#   pos <- x > 0
#   out[pos] <- if (power == 0) log(x[pos]) else x[pos]^power
#   out
# }
# 
# # Linear predictor (paper, Equation 15)
# eta <- 0.5 +
#   20.00 * positive_fp(x1, -2) - 1.40 * (x1 == 0) +   # FP1 + binary
#   15.00 * positive_fp(x2, -2) +
#    0.25 * positive_fp(x2, 2) +                         # FP2 only
#    1.40 * (x3 == 0) +                                  # binary only
#    4.00 * log(x4)                                      # ordinary FP1
# 
# y <- eta + rnorm(n, mean = 0, sd = 0.80)
# simdat <- data.frame(y, x1, x2, x3, x4, x5, x6)
# 
# colMeans(simdat[c("x1", "x2", "x3")] == 0)  # zero proportions

## ----multiple-saz-fit, eval=FALSE---------------------------------------------
# fit_sim <- mfp2(
#   y ~ fp(x1, spike = TRUE, center = FALSE) +
#     fp(x2, spike = TRUE, center = FALSE) +
#     fp(x3, spike = TRUE, center = FALSE) +
#     fp(x4, center = FALSE) +
#     fp(x5, center = FALSE) +
#     fp(x6, center = FALSE),
#   data = simdat,
#   family = "gaussian",
#   criterion = "bic",
#   verbose = FALSE
# )
# print(fit_sim, detailed_settings = FALSE, notes = FALSE)

## ----matrix-fits--------------------------------------------------------------
y <- prostate$lpsa
x <- as.matrix(prostate[, "pgg45", drop = FALSE])  # drop = FALSE keeps a matrix
x_multi <- as.matrix(
  prostate[, c("pgg45", "age", "cavol", "weight", "bph", "cp", "svi")]
)

fit_pgg45_mat <- mfp2(x, y, family = "gaussian", criterion = "pvalue",
                      select = 0.05, alpha = 0.05, spike_vars = "pgg45",
                      verbose = FALSE)

fit_multi_mat <- mfp2(x_multi, y, family = "gaussian", criterion = "pvalue",
                      select = 0.05, alpha = 0.05, spike_vars = "pgg45",
                      verbose = FALSE)

fit_aic_mat <- mfp2(x, y, family = "gaussian", criterion = "aic",
                    spike_vars = "pgg45", verbose = FALSE)

print(fit_pgg45_mat, detailed_settings = FALSE, notes = FALSE)

