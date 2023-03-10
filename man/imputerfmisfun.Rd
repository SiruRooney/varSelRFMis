% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/varselrfmisMIfun.R
\name{imputerfmisfun}
\alias{imputerfmisfun}
\title{imputerfmisfun}
\usage{
imputerfmisfun(
  predimp.w,
  predy.w,
  misindcom.tw,
  dectnalist,
  weight,
  itvs,
  optmethod,
  dummytran,
  bicglmnet,
  bdiag_m,
  indmisfun,
  ridgselectedvars
)
}
\arguments{
\item{predimp.w:}{A dataframe of SNPs after multiple imputation at last iteration.}

\item{predy.w:}{A vector of phenotype after multiple imputation at last iteration.}

\item{misindcom.tw:}{A missing indicator matrix for a dataframe of SNPs and phenotype after multiple imputation.}

\item{dectnalist:}{A vector which indicates whether SNPs and phenotype are missing.}

\item{weight:}{A vector indicating each weight in every individual.}

\item{itvs:}{Count for the number of the loop.}

\item{dummytran:}{A function which converts each SNP vector into the dummy variable}

\item{bicglmnet:}{A function which chooses an optimal tunning parameter using the model criterion BIC.}

\item{bdiag_m:}{A function which generates a block diagonal matrix.}

\item{indmisfun:}{A function which indicates the value of a SNP ("0", "1" or "2").}

\item{ridgselectedvars:}{A list which contains selected important variables in a series of models for missing SNPs, corresponding missing indicators and phenotype.}
}
\value{
a list named "ridgimp" containing the results from fitting the ridge regression models
}
\description{
FUNCTION: To apply glmnet function to choose the optimal tunning parameter in the constructed models.
}
