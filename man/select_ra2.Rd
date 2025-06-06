% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mfp_step.R
\name{select_ra2}
\alias{select_ra2}
\title{Function selection procedure based on closed testing procedure}
\usage{
select_ra2(
  x,
  xi,
  keep,
  degree,
  acdx,
  y,
  powers_current,
  powers,
  criterion,
  ftest,
  select,
  alpha,
  family,
  ...
)
}
\arguments{
\item{x}{an input matrix of dimensions nobs x nvars. Does not contain
intercept, but columns are already expanded into dummy variables as
necessary. Data are assumed to be shifted and scaled.}

\item{xi}{a character string indicating the name of the current variable
of interest, for which the best fractional polynomial transformation is
to be estimated in the current step.}

\item{keep}{a character vector with names of variables to be kept
in the model.}

\item{degree}{integer > 0 giving the degree for the FP transformation.}

\item{acdx}{a logical vector of length nvars indicating continuous variables
to undergo the approximate cumulative distribution (ACD) transformation.}

\item{y}{a vector for the response variable or a \code{Surv} object.}

\item{powers_current}{a list of length equal to the number of variables,
indicating the fp powers to be used in the current step for all variables
(except \code{xi}).}

\item{powers}{a named list of numeric values that sets the permitted FP
powers for each covariate.}

\item{criterion}{a character string defining the criterion used to select
variables and FP models of different degrees.}

\item{ftest}{a logical indicating the use of the F-test for Gaussian models.}

\item{select}{a numeric value indicating the significance level
for backward elimination of \code{xi}.}

\item{alpha}{a numeric value indicating the significance level
for tests between FP models of different degrees for \code{xi}.}

\item{family}{a character string representing a family object.}

\item{...}{passed to fitting functions.}
}
\value{
A list with several components:
\itemize{
\item \code{keep}: logical indicating if \code{xi} is forced into model.
\item \code{acd}: logical indicating if an ACD transformation was applied for \code{xi},
i.e. \code{FALSE} in this case.
\item \code{powers}: (best) fp powers investigated in step, indexing \code{metrics}.
Always starts with highest power, then null, then linear, then FP in
increasing degree (e.g. FP2, null, linear, FP1).
\item \code{power_best}: a numeric vector with the best power found. The returned
best power may be \code{NA}, indicating the variable has been removed from the
model.
\item \code{metrics}: a matrix with performance indices for all models investigated.
Same number of rows as, and indexed by, \code{powers}.
\item \code{model_best}: row index of best model in \code{metrics}.
\item \code{pvalue}: p-value for comparison of linear and null model.
\item \code{statistic}: test statistic used, depends on \code{ftest}.
}
}
\description{
Used in \code{\link[=find_best_fp_step]{find_best_fp_step()}} when \code{criterion = "pvalue"}.
For parameter explanations, see \code{\link[=find_best_fp_step]{find_best_fp_step()}}. All parameters
captured by \code{...} are passed on to \code{\link[=fit_model]{fit_model()}}.
}
\details{
In case \code{criterion = "pvalue"} the function selection procedure as outlined
in Chapters 4 and 6 of Royston and Sauerbrei (2008) is used.
\itemize{
\item \emph{Step 1}: test the best FPm function against a null model at level
\code{select} with 2m df. If not significant, the variable is excluded.
Otherwise continue with step 2.
\item \emph{Step 2}: test the best FPm versus a linear model at level \code{alpha}
with 2m - 1 df. If not significant, use a linear model.
Otherwise continue with step 3.
\item \emph{Step 3}: test the best FPm versus the best FP1 at
level \code{alpha} with 2m - 2 df. If not significant, use the best FP1 model.
Otherwise, repeat this step for all remaining higher order FPs until
FPm-1, which is tested at level \code{alpha} with 2 df against FPm.
If the final test is not significant, use a FPm-1 model, otherwise use FPm.
}

Note that the "best" FPx model used in each step is given by the model using
a FPx transformation for the variable of interest and having the highest
likelihood of all such models given the current powers for all other
variables, as outlined in Section 4.8 of Royston and Sauerbrei (2008).
These best FPx models are computed in \code{\link[=find_best_fpm_step]{find_best_fpm_step()}}.

When a variable is forced into the model by including it in \code{keep}, then
this function will not exclude it from the model (by setting its power to
\code{NA}), but will only choose its functional form.
}
\references{
Royston, P. and Sauerbrei, W., 2008. \emph{Multivariable Model - Building:
A Pragmatic Approach to Regression Anaylsis based on Fractional Polynomials
for Modelling Continuous Variables. John Wiley & Sons.}
}
\seealso{
\code{\link[=select_ra2_acd]{select_ra2_acd()}}
}
