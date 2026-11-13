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
library(ggplot2)
has_patchwork <- requireNamespace("patchwork", quietly = TRUE)
library(survival)

has_rms <- requireNamespace("rms", quietly = TRUE) &&
  utils::packageVersion("rms") >= "7.0-0"

## ----eval=FALSE-----------------------------------------------------------------------------------
# # install the package
# install.packages("mfp2")
# 
# # load the package
# library(mfp2)

## ----prostate-data, out.lines=10------------------------------------------------------------------
data("prostate")
head(prostate, 6)

## ----prostate-formula-fit-------------------------------------------------------------------------
fit_formula <- mfp2(
  lpsa ~ fp(age) + svi + fp(pgg45) + fp(cavol) +
    fp(weight) + fp(bph) + fp(cp),
  data = prostate,
  family = gaussian(),
  verbose = FALSE
)
summary(fit_formula)

## ----prostate-matrix-fit--------------------------------------------------------------------------
x <- as.matrix(prostate[
  , c("age", "svi", "pgg45", "cavol", "weight", "bph", "cp")
])
y <- prostate$lpsa

fit_matrix <- mfp2(x, y, family = gaussian(), verbose = FALSE)
summary(fit_matrix)

## ----inspect-shift-scale--------------------------------------------------------------------------
fit_formula$transformations[, c("shift", "scale"), drop = FALSE]

## ----manual-shift-scale-formula, eval=FALSE-------------------------------------------------------
# fit_manual <- mfp2(
#   lpsa ~
#     fp(age, shift = 0, scale = 10) +
#     fp(pgg45, shift = 1, scale = 100) +
#     fp(cavol, shift = 0, scale = 10) +
#     fp(weight, shift = 0, scale = 100) +
#     fp(bph, shift = 1, scale = 1) +
#     fp(cp, shift = 1, scale = 10) +
#     svi,
#   data = prostate,
#   family = gaussian()
# )

## ----subset-formula-examples, eval=FALSE----------------------------------------------------------
# keep_rows <- seq_len(nrow(prostate)) %% 2L == 0L
# 
# mfp2(lpsa ~ fp(age) + fp(cavol), data = prostate, subset = age >= 60)
# mfp2(lpsa ~ fp(age) + fp(cavol), data = prostate, subset = keep_rows)
# mfp2(lpsa ~ fp(age) + fp(cavol), data = prostate, subset = c(10, 3, 8, 15, 20))

## ----predictor-df-settings, eval=FALSE------------------------------------------------------------
# mfp2(
#   lpsa ~ fp(age, df = 2) + svi + fp(pgg45) + fp(cavol, df = 4) +
#     fp(weight) + fp(bph) + fp(cp, df = 1),
#   data = prostate,
#   family = gaussian()
# )

## ----default-predictor-df, eval=FALSE-------------------------------------------------------------
# mfp2(
#   lpsa ~ age + svi + pgg45 + cavol + weight + bph + cp,
#   data = prostate,
#   family = gaussian()
# )

## ----retained-predictors, eval=FALSE--------------------------------------------------------------
# mfp2(
#   lpsa ~ fp(age, select = 1) + fp(svi, df = 1, select = 1) +
#     fp(pgg45) + fp(cavol) + fp(weight) + fp(bph) + fp(cp),
#   data = prostate,
#   family = gaussian(),
#   criterion = "pvalue"
# )

## ----aic-selection, eval=FALSE--------------------------------------------------------------------
# mfp2(
#   lpsa ~ fp(age) + svi + fp(pgg45) + fp(cavol) +
#     fp(weight) + fp(bph) + fp(cp),
#   data = prostate,
#   family = gaussian(),
#   criterion = "aic",
#   keep = c("age", "svi")
# )

## ----gaussian-ftest, eval=FALSE-------------------------------------------------------------------
# mfp2(
#   lpsa ~ fp(age) + svi + fp(pgg45) + fp(cavol) +
#     fp(weight) + fp(bph) + fp(cp),
#   data = prostate,
#   family = gaussian(),
#   criterion = "pvalue",
#   ftest = TRUE
# )

## ----candidate-powers, eval=FALSE-----------------------------------------------------------------
# custom_powers <- list(age = c(0.5, 1), cavol = c(0, 1))
# 
# mfp2(
#   lpsa ~ fp(age) + svi + fp(pgg45) + fp(cavol) +
#     fp(weight) + fp(bph) + fp(cp),
#   data = prostate,
#   family = gaussian(),
#   powers = custom_powers
# )

## ----candidate-powers-in-formula, eval=FALSE------------------------------------------------------
# mfp2(
#   lpsa ~ fp(age, powers = c(0.5, 1)) + svi + fp(pgg45) +
#     fp(cavol, powers = c(0, 1)) + fp(weight) + fp(bph) + fp(cp),
#   data = prostate,
#   family = gaussian()
# )

## ----algorithm-output, out.lines=103--------------------------------------------------------------
fit <- mfp2(
  lpsa ~ fp(age) + svi + fp(pgg45) + fp(cavol) +
    fp(weight) + fp(bph) + fp(cp),
  data = prostate,
  family = gaussian(),
  criterion = "pvalue",
  select = 0.05,
  alpha = 0.05,
  ftest = TRUE
)

## ----gaussian-methods, eval=TRUE------------------------------------------------------------------
fit <- mfp2(
  lpsa ~ fp(age) + svi + fp(pgg45) + fp(cavol) +
    fp(weight) + fp(bph) + fp(cp),
  data = prostate,
  family = gaussian(),
  verbose = FALSE
)

class(fit)

## ----print-gaussian-model, eval=TRUE--------------------------------------------------------------
print(fit)

## ----compact-model-report, eval=TRUE--------------------------------------------------------------
print(fit, detailed_settings = FALSE)

## ----extract-selection-details, eval=TRUE---------------------------------------------------------
get_selected_variable_names(fit)
fit$fp_terms
fit$transformations

## ----summarise-gaussian-model, eval=TRUE----------------------------------------------------------
summary(fit)

## ----detailed-model-summary, eval=TRUE------------------------------------------------------------
summary(fit, formulas = TRUE, basis = TRUE)
summary(fit, raw = TRUE)

## ----extract-model-coefficients, eval=TRUE--------------------------------------------------------
coef(fit)
vcov(fit)
confint(fit)

## ----gaussian-prediction, eval=TRUE---------------------------------------------------------------
new_observations <- prostate[1:5, , drop = FALSE]

pred <- predict(
  fit,
  newdata = new_observations,
  type = "response",
  se.fit = TRUE
)

pred$fit
pred$se.fit

## ----predict-selected-term, eval=TRUE-------------------------------------------------------------
predict(
  fit,
  newdata = new_observations,
  type = "terms",
  terms = "cavol",
  terms_seq = "data",
  add_intercept = FALSE
)

## ----plot-cavol, eval=TRUE, fig.cap="Prostate data. Fitted function for cavol with pointwise 95% confidence intervals and component-plus-residual points.", fig.width=7, fig.height=5, out.width="100%", fig.align="center"----
plots <- plot(fit, terms = "cavol")
plots[["cavol"]]

## ----plot-cavol-function-only, eval=TRUE----------------------------------------------------------
curves <- plot(fit, terms = "cavol", partial_only = TRUE)
curves[["cavol"]]

## ----plot-cavol-contrasts, eval=TRUE, fig.cap="Prostate data. Differences in fitted mean lpsa relative to cavol = 30, with pointwise 95% confidence intervals.", fig.width=7, fig.height=5, out.width="100%", fig.align="center"----
contrasts <- plot(
  fit,
  terms = "cavol",
  type = "contrasts",
  ref = list(cavol = 30)
)

contrasts[["cavol"]]

## ----plot-evaluation-values, eval=has_patchwork, fig.cap="Art data. The same fitted function evaluated at equally spaced values (left) and observed values (right). Residual points remain at the observed values in both panels.", fig.width=8, fig.height=5, out.width="100%", fig.align="center"----
data("art")

fit2 <- mfp2(
  y ~ fp(x5),
  data = art,
  family = gaussian(),
  verbose = FALSE
)

observed <- plot(fit2, terms = "x5", terms_seq = "data")
equidistant <- plot(fit2, terms = "x5", terms_seq = "equidistant")

patchwork::wrap_plots(
  equidistant[["x5"]] + ggplot2::ggtitle("Equally spaced values"),
  observed[["x5"]] + ggplot2::ggtitle("Observed values"),
  ncol = 2
)

## ----ordinal-score-example, eval=FALSE------------------------------------------------------------
# dat$severity_score <- NA_real_
# dat$severity_score[dat$severity == "mild"] <- 0
# dat$severity_score[dat$severity == "moderate"] <- 1
# dat$severity_score[dat$severity == "severe"] <- 2
# 
# fit <- mfp2(
#   y ~ fp(severity_score, df = 1) + fp(age),
#   data = dat
# )

## ----categorical-formula-example, eval=TRUE-------------------------------------------------------
data("art")
dat <- art

dat$x4 <- ordered(dat$x4, levels = c(1, 2, 3))
dat$x9 <- factor(dat$x9, levels = c(1, 2, 3))

fit_formula <- mfp2(
  y ~ fp(x1) + x2 + fp(x3) + x4 + fp(x5) + fp(x6) +
    fp(x7) + x8 + x9 + fp(x10),
  data = dat,
  family = gaussian(),
  keep = c("x4", "x9"),
  verbose = FALSE
)

fit_formula

## ----categorical-matrix-example, eval=FALSE-------------------------------------------------------
# fit_matrix <- mfp2(
#   x,
#   y,
#   family = gaussian(),
#   term_groups = list(stage = c("stageII", "stageIII")),
#   keep = "stage"
# )

## ----logistic-fit, eval=TRUE----------------------------------------------------------------------
data("pima")
head(pima)

fit <- mfp2(
  y ~ fp(pregnant) + fp(glucose) + fp(diastolic) + fp(triceps) +
    fp(insulin) + fp(bmi) + fp(diabetes) + fp(age),
  data = pima,
  family = binomial(),
  verbose = FALSE
)

fit

## ----logistic-functions, eval=has_patchwork, fig.cap="Pima data. Selected predictor contributions on the log-odds scale, with pointwise 95% confidence intervals and component-plus-residual points.", fig.width=7, fig.height=6, out.width="100%", fig.align="center"----
plots <- plot(fit)
patchwork::wrap_plots(plots, ncol = 2)

## ----logistic-probabilities, eval=TRUE------------------------------------------------------------
new_observations <- pima[1:5, , drop = FALSE]

predict(fit, newdata = new_observations, type = "response")

## ----cox-fit, eval=TRUE, warning=FALSE------------------------------------------------------------
library(survival)
data("gbsg")

# Specify the grade order and cumulative contrasts.
grade_levels <- sort(unique(as.character(gbsg$grade)))
gbsg$grade_ord <- ordered(gbsg$grade, levels = grade_levels)

contrasts(gbsg$grade_ord) <- cbind(
  grade_1 = c(0, 1, 1),
  grade_2 = c(0, 0, 1)
)

fit1 <- mfp2(
  Surv(rectime, censrec) ~
    fp(age) + meno + fp(size) + fp(nodes) + fp(pgr) + fp(er) +
    grade_ord + hormon,
  data = gbsg,
  family = "cox",
  ties = "breslow",
  keep = "hormon",
  control = coxph.control(iter.max = 50),
  verbose = FALSE
)

fit1

## ----cox-monotonic-fit, eval=TRUE-----------------------------------------------------------------
gbsg$enodes <- exp(-0.12 * gbsg$nodes)

fit2 <- mfp2(
  Surv(rectime, censrec) ~
    fp(age) + meno + fp(size) +
    fp(enodes, df = 2, powers = c(0.5, 1, 2, 3)) +
    fp(pgr) + fp(er) + grade_ord + hormon,
  data = gbsg,
  family = "cox",
  ties = "breslow",
  keep = "hormon",
  control = coxph.control(iter.max = 50),
  verbose = FALSE
)

fit2

## ----stratified-cox-fit, eval=TRUE----------------------------------------------------------------
fit3 <- mfp2(
  Surv(rectime, censrec) ~
    fp(age) + meno + fp(size) +
    fp(enodes, df = 2, powers = c(0.5, 1, 2, 3)) +
    fp(pgr) + fp(er) + grade_ord + strata(hormon),
  data = gbsg,
  family = "cox",
  ties = "breslow",
  control = coxph.control(iter.max = 50),
  verbose = FALSE
)

fit3

## ----multinomial-data-----------------------------------------------------------------------------
set.seed(1101)
n_multinom <- 500

dat_multinom <- data.frame(
  x = runif(n_multinom, 0.5, 4),
  group = factor(sample(c("A", "B"), n_multinom, replace = TRUE))
)

eta2 <- 0.2 + 5 * log(dat_multinom$x) +
  0.4 * (dat_multinom$group == "B")
eta3 <- -0.3 + 0.4 * dat_multinom$x -
  0.4 * (dat_multinom$group == "B")
probabilities <- cbind(class1 = 1, class2 = exp(eta2), class3 = exp(eta3))
probabilities <- probabilities / rowSums(probabilities)

dat_multinom$outcome <- factor(
  vapply(seq_len(n_multinom), function(i) {
    sample(colnames(probabilities), size = 1, prob = probabilities[i, ])
  }, character(1)),
  levels = colnames(probabilities)
)

head(dat_multinom)

## ----multinomial-fit------------------------------------------------------------------------------
fit_multinom <- mfp2(
  outcome ~ fp(x, df = 2) + group,
  data = dat_multinom,
  family = multinomial_family(reference = "class1"),
  keep = c("x", "group"),
  verbose = FALSE
)

## ----multinomial-prediction-----------------------------------------------------------------------
predict(
  fit_multinom,
  newdata = dat_multinom[1:6, c("x", "group")],
  type = "response"
)

## ----multinomial-plot-----------------------------------------------------------------------------
multinom_plots <- plot(fit_multinom, terms = "x", partial_only = TRUE)
print(multinom_plots[[1]])

## ----ordinal-data, eval=has_rms-------------------------------------------------------------------
set.seed(1102)
n_ordinal <- 500

dat_ordinal <- data.frame(
  x = runif(n_ordinal, 0.5, 4),
  group = factor(sample(c("A", "B"), n_ordinal, replace = TRUE))
)
latent <- 2 * (dat_ordinal$x - 2)^2 - 2 +
  0.5 * (dat_ordinal$group == "B") + rlogis(n_ordinal)

dat_ordinal$outcome <- ordered(
  cut(
    latent,
    breaks = c(-Inf, 0.2, 1.3, Inf),
    labels = c("low", "medium", "high")
  ),
  levels = c("low", "medium", "high")
)

head(dat_ordinal)

## ----ordinal-fit, eval=has_rms--------------------------------------------------------------------
fit_ordinal <- mfp2(
  outcome ~ fp(x, df = 4) + group,
  data = dat_ordinal,
  family = ordinal_family(),
  keep = c("x", "group"),
  verbose = FALSE
)

## ----ordinal-prediction, eval=has_rms-------------------------------------------------------------
predict(
  fit_ordinal,
  newdata = dat_ordinal[1:5, c("x", "group")],
  type = "response",
  se.fit = FALSE
)

## ----ordinal-plot, eval=has_rms-------------------------------------------------------------------
ordinal_plots <- plot(fit_ordinal, terms = "x", partial_only = TRUE)
print(ordinal_plots[[1]])

## ----gee-data-------------------------------------------------------------------------------------
set.seed(1103)
n_individuals <- 50
n_visits <- 4

dat_gee <- data.frame(
  id = rep(seq_len(n_individuals), each = n_visits),
  visit = rep(seq_len(n_visits), times = n_individuals),
  x = runif(n_individuals * n_visits, 0.5, 4),
  group = factor(rep(
    sample(c("A", "B"), n_individuals, replace = TRUE),
    each = n_visits
  ))
)

individual_effect <- rep(rnorm(n_individuals, sd = 0.6), each = n_visits)
dat_gee$outcome <- 0.8 * log(dat_gee$x) +
  0.4 * (dat_gee$group == "B") + individual_effect +
  rnorm(nrow(dat_gee), sd = 0.5)

## ----gee-fit--------------------------------------------------------------------------------------
fit_gee <- mfp2(
  outcome ~ fp(x, df = 2) + group,
  data = dat_gee,
  family = gee_family(family = gaussian(), corstr = "exchangeable"),
  id = id,
  waves = visit,
  criterion = "aic",
  keep = c("x", "group"),
  verbose = FALSE
)

## ----gee-prediction-------------------------------------------------------------------------------
predict(
  fit_gee,
  newdata = dat_gee[1:5, c("x", "group")],
  type = "response",
  se.fit = TRUE
)

## ----gee-plot-------------------------------------------------------------------------------------
gee_plots <- plot(fit_gee, terms = "x", partial_only = TRUE)
print(gee_plots[[1]])

## ----cleanup, include=FALSE---------------------------------------------------
options(width = oldwidth)

