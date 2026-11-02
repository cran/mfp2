## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 9,
  fig.height = 4.5,
  fig.align = "center",
  out.width = "100%"
)
options(width = 80)
has_ggplot2 <- requireNamespace("ggplot2", quietly = TRUE)
has_geepack <- requireNamespace("geepack", quietly = TRUE)
has_patchwork <- requireNamespace("patchwork", quietly = TRUE)

## ----packages-----------------------------------------------------------------
library(mfp2)
library(survival)

## ----data---------------------------------------------------------------------
data("advanced_prostate_cancer", package = "mfp2")
advanced_prostate_cancer$rx <- factor(
  advanced_prostate_cancer$rx,
  levels = c(0, 1), labels = c("Placebo", "Treatment")
)

table(advanced_prostate_cancer$rx)
tapply(advanced_prostate_cancer$age, advanced_prostate_cancer$rx, range)

## ----overall-treatment-effect-------------------------------------------------
survdiff(Surv(survtime, cens) ~ rx, data = advanced_prostate_cancer)

cox <- coxph(Surv(survtime, cens) ~ rx, data = advanced_prostate_cancer)
summary(cox)$conf.int

## ----unadjusted-fit-----------------------------------------------------------
fit_unadjusted <- mfpi(
  Surv(survtime, cens) ~ rx + fp(age),
  data = advanced_prostate_cancer,
  family = "cox",
  group_var = "rx",
  interaction_vars = "age",
  interaction_forms = c(age = "fp2"),
  flex = "flex1",
  ties = "breslow"
)
print(fit_unadjusted)

## ----unadjusted-plot, eval=has_patchwork && has_ggplot2, fig.cap="Unadjusted analysis: fitted age functions by group (left) and the Treatment-minus-Placebo log hazard ratio with 95% CI (right); paper, Figure 2."----
p_unadj <- plot(fit_unadjusted, terms = "age", plot_type = "both",
                show_title = FALSE)

## ----unadjusted-predict-------------------------------------------------------
age_pred <- predict(
  fit_unadjusted,
  newdata = data.frame(age = c(50, 60, 70)),
  terms = "age",
  level = 0.95,
  type = "both"
)
age_pred

## ----adjusted-fit-------------------------------------------------------------
fit_adjusted <- mfpi(
  Surv(survtime, cens) ~ rx + fp(age) + fp(wt) + fp(sbp) +
    fp(dbp) + fp(sz) + fp(ap) + fp(hg) + fp(sg) +
    pf + hx + bm + stage + ekg,
  data = advanced_prostate_cancer,
  family = "cox",
  group_var = "rx",
  interaction_vars = "age",
  interaction_forms = c(age = "fp2"),
  flex = "flex1",
  ties = "breslow"
)
print(fit_adjusted)

## ----adjusted-plot, eval=has_patchwork && has_ggplot2, fig.cap="Adjusted analysis: fitted age functions by group (left) and the Treatment-minus-Placebo log hazard ratio with 95% CI (right); paper, Figure 3."----
p_adj <- plot(fit_adjusted, terms = "age", plot_type = "both",
              show_title = FALSE)

## ----aic-fit------------------------------------------------------------------
fit_aic <- mfpi(
  Surv(survtime, cens) ~
    rx + fp(age) + fp(wt) + fp(sbp) + fp(dbp) +
    fp(sz) + fp(ap) + fp(hg) + fp(sg) +
    pf + hx + bm + stage + ekg,
  data = advanced_prostate_cancer,
  family = "cox",
  group_var = "rx",
  interaction_vars = "age",
  interaction_forms = c(age = "fp2"),
  flex = "flex1",
  ties = "breslow",
  criterion = "aic",
  min_improvement = 2
)
print(fit_aic)

## ----multiple-fit-------------------------------------------------------------
fit_multiple_interactions <- mfpi(
  Surv(survtime, cens) ~
    rx + fp(age) + fp(wt) + fp(sbp) + fp(dbp) +
    fp(sz) + fp(ap) + fp(hg) + fp(sg) +
    pf + hx + bm + stage + ekg,
  data = advanced_prostate_cancer,
  family = "cox",
  group_var = "rx",
  interaction_vars = c("age", "wt", "sbp", "dbp", "sz", "ap", "hg"),
  interaction_forms = c(age = "fp2"),
  flex = "flex1",
  ties = "breslow",
  p_adjust_method = "holm",
  verbose = FALSE
)
print(fit_multiple_interactions)

## ----fixed-linear-fit---------------------------------------------------------
fit_fixed_linear <- mfpi(
  Surv(survtime, cens) ~
    rx +
    fp(age, df = 1) + fp(wt, df = 1) + fp(sbp, df = 1) +
    fp(dbp, df = 1) + fp(sz, df = 1) + fp(ap, df = 1) +
    fp(hg, df = 1) + fp(sg, df = 1) +
    pf + hx + bm + stage + ekg,
  data = advanced_prostate_cancer,
  family = "cox",
  group_var = "rx",
  interaction_vars = "age",
  interaction_forms = c(age = "fp2"),
  flex = "flex1",
  ties = "breslow",
  criterion = "pvalue",
  keep = c("age", "wt", "sbp", "dbp", "sz", "ap", "hg", "sg",
           "pf", "hx", "bm", "stage", "ekg")
)
print(fit_fixed_linear)

## ----xy-data------------------------------------------------------------------
y <- with(advanced_prostate_cancer, Surv(survtime, cens))
x <- advanced_prostate_cancer[, c(
  "rx", "age", "wt", "sbp", "dbp", "sz", "ap", "hg", "sg",
  "pf", "hx", "bm", "stage", "ekg"
)]

## ----xy-unadjusted------------------------------------------------------------
fit_unadjusted_xy <- mfpi(
  x = x[, c("rx", "age")],
  y = y,
  family = "cox",
  group_var = "rx",
  interaction_vars = "age",
  interaction_forms = c(age = "fp2"),
  flex = "flex1",
  ties = "breslow"
)
print(fit_unadjusted_xy)

## ----xy-other-analyses, eval=FALSE--------------------------------------------
# fit_adjusted_xy <- mfpi(
#   x = x, y = y, family = "cox",
#   group_var = "rx", interaction_vars = "age",
#   interaction_forms = c(age = "fp2"), flex = "flex1", ties = "breslow"
# )
# 
# fit_aic_xy <- mfpi(
#   x = x, y = y, family = "cox",
#   group_var = "rx", interaction_vars = "age",
#   interaction_forms = c(age = "fp2"), flex = "flex1", ties = "breslow",
#   criterion = "aic", min_improvement = 2
# )
# 
# fit_multiple_xy <- mfpi(
#   x = x, y = y, family = "cox",
#   group_var = "rx",
#   interaction_vars = c("age", "wt", "sbp", "dbp", "sz", "ap", "hg"),
#   interaction_forms = c(age = "fp2"), flex = "flex1", ties = "breslow",
#   p_adjust_method = "holm", verbose = FALSE
# )

## ----xy-fixed-linear, eval=FALSE----------------------------------------------
# fit_fixed_linear_xy <- mfpi(
#   x = x, y = y, family = "cox",
#   group_var = "rx", interaction_vars = "age",
#   interaction_forms = c(age = "fp2"), flex = "flex1", ties = "breslow",
#   criterion = "pvalue",
#   df = c(age = 1, wt = 1, sbp = 1, dbp = 1,
#          sz = 1, ap = 1, hg = 1, sg = 1),
#   keep = c("age", "wt", "sbp", "dbp", "sz", "ap", "hg", "sg",
#            "pf", "hx", "bm", "stage", "ekg")
# )

