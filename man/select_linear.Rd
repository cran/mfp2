% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mfp_step.R
\name{select_linear}
\alias{select_linear}
\title{Helper to select between null and linear term for a single variable}
\usage{
select_linear(
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

\item{degree}{not used.}

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
\item \code{acd}: logical indicating if an ACD transformation was applied for \code{xi}.
\item \code{powers}: fp powers investigated in step, indexing \code{metrics}.
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
To be used in \code{\link[=find_best_fp_step]{find_best_fp_step()}}. Only used if \code{df = 1} for a variable.
Handles all criteria for selection.
For parameter explanations, see \code{\link[=find_best_fp_step]{find_best_fp_step()}}. All parameters
captured by \code{...} are passed on to \code{\link[=fit_model]{fit_model()}}.
}
\details{
This function assesses a single variable of interest \code{xi} regarding its
functional form in the current working model as indicated by
\code{powers_current}, with the choice between a excluding \code{xi} ("null model") and
including a linear term ("linear fp") for \code{xi}.

Note that this function handles an ACD transformation for \code{xi} as well.

When a variable is forced into the model by including it in \code{keep}, then
this function will not exclude it from the model (by setting its power to
\code{NA}), but will only choose the linear model.
}
