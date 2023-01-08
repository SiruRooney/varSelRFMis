
#'gensetf
#'@description Function name: genset.f. It is used to generate all possible conbinations of factor levels in an experiment.Written by Guoqi Qian, 17/11/2002.
#'@param gvec: gvec: vector consisting of numbers of levels of all factors.
#'
#'@return a matrix containing all possible combination
#'@export

gensetf<-function(gvec){

  nvar <- length(gvec)
  genset <- matrix(1, nvar, prod(gvec))
  genset[nvar, 1:gvec[nvar]] <- 1:gvec[nvar]
  for(i in 2:nvar) {
    temp <- prod(gvec[(nvar - i + 2):nvar])
    for(j in 2:gvec[nvar - i + 1]) {
      genset[nvar - i + 1, ((j - 1) * temp + 1):(j * temp)] <- j
      genset[(nvar - i + 2):nvar, ((j - 1) * temp + 1):(j * temp)] <- genset[(nvar - i + 2):nvar,
                                                                             1:temp]
    }
  }
  return(genset)
}
