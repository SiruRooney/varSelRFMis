% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/varselrfmisMIfun4varsel.R
\name{varSelRF_mis}
\alias{varSelRF_mis}
\title{varSelRF_mis}
\usage{
varSelRF_mis(
  xdata,
  Class,
  c.sd = 1,
  mtryFactor = 1,
  ntree = 5000,
  ntreeIterat = 2000,
  vars.drop.num = NULL,
  vars.drop.frac = 0.2,
  whole.range = TRUE,
  recompute.var.imp = FALSE,
  verbose = FALSE,
  returnFirstForest = TRUE,
  fitted.rf = NULL,
  keep.forest = FALSE,
  caseweights
)
}
\arguments{
\item{xdata:}{predictor matrix.}

\item{class:}{the response. Must be a factor}

\item{c.sd:}{The factor that multiplies the sd. to decide on stopping the iterations or choosing the final solution.}

\item{mtryFactor:}{The multiplication factor}

\item{ntree:}{The number of trees to use for the first forest; same as ntree for randomForest}

\item{ntreeIterat:}{The number of trees to use (ntree of randomForest) for all additional forests}

\item{vars.drop.num:}{The number of variables to exclude at each iteration.}

\item{vars.drop.frac:}{The fraction of variables, from those in the previous forest, to exclude at each iteration.}

\item{whole.range:}{If TRUE continue dropping variables until a forest with only two variables is built, and choose the best model from the complete series of models. If FALSE, stop the iterations if the current OOB error becomes larger than the initial OOB error}

\item{recompute.var.imp:}{If TRUE recompute variable importances at each new iteration.}

\item{verbose:}{Give more information about what is being done}

\item{returnFirstForest:}{If TRUE the random forest from the complete set of variables is returned.}

\item{fitted.rf:}{An (optional) object of class randomForest previously fitted.}

\item{keep.forest:}{Same argument as in randomForest function.}

\item{caseweights:}{each weight corresponds to the individual.}
}
\value{
a list names "varSelRF"
}
\description{
FUNCTION: Using random forests to select important variables. Mainly written by Diaz
}
