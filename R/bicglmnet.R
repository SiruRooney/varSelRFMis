#' bdiag_m
#' @description FUNCTION:Fast version of block diag. Copyright (C) 2016 Martin Maechler, ETH Zurich
#' @param lamt: A list of block matrices.
#' 
#' @import Matrix
#' @return a matrix
#' @export
bdiag_m <- function(lmat) {
#FUNCTION:Fast version of block diag.
  ## Copyright (C) 2016 Martin Maechler, ETH Zurich
#Argments:
#lamt: A list of block matrices. 
  
  if(!length(lmat)) return(new("dgCMatrix"))
  stopifnot(is.list(lmat), is.matrix(lmat[[1]]),
            (k <- (d <- dim(lmat[[1]]))[1]) == d[2], # k x k
            all(vapply(lmat, dim, integer(2)) == k)) # all of them
  N <- length(lmat)
  if(N * k > .Machine$integer.max)
    stop("resulting matrix too large; would be  M x M, with M=", N*k)
  M <- as.integer(N * k)
  ## result: an   M x M  matrix
  new("dgCMatrix", Dim = c(M,M),
      ## 'i :' maybe there's a faster way (w/o matrix indexing), but elegant?
      i = as.vector(matrix(0L:(M-1L), nrow=k)[, rep(seq_len(N), each=k)]),
      p = k * 0L:M,
      x = as.double(unlist(lmat, recursive=FALSE, use.names=FALSE)))
}

#' bicglmnet
#' @description FUNCTION: Defined BIC for binomial distribution and trinomial distribution.
#' @param x: Predictor matrix.
#' @param y: Reponse.
#' @param indmisfunction: To indicate which value of one SNP belongs to.
#' @param bdiag_M: Fast version of block diag.
#' @import Matrix
#' @import glmnet
#' @import psych
#'
#' @return a list contained the optimal ridge regression model
#' @export
bicglmnet<-function(x,y,weight,indmisfun,bdiag_m){
#FUNCTION: Defined BIC for binomial distribution and trinomial distribution.
#Argments:
#x: Predictor matrix.
#y: Reponse.
#indmisfunction: To indicate which value of one SNP belongs to.
#bdiag_M: Fast version of block diag.
  
  
  tre=glmnet(x=x,y=y,weights=weight,alpha=0,family="multinomial")
  lam=tre$lambda
  ppro=predict(object=tre,newx=x,type="response") #dim(ppro)=nrow*(2/3)*length(lambda)
  
  if(length(levels(y))==2){
    #binomial distribution
    wcl=weight*ppro[,-1,]*(1-ppro[,-1,])
    

    trhat=foreach(ilam=1:length(lam),.combine='c')%dopar%{
      library(Matrix)#bdiag_m
      library(psych)
     xwx=t(x)%*%diag(wcl[,ilam])%*%x
     lami=diag(lam[ilam],NCOL(x))
     trhat=tr(solve(xwx+NROW(x)*lami)%*%xwx)
     return(trhat)
    }

    #trhat is a vector of length(lambda)

  }else{
    z=matrix(0,2*NROW(x),2*NCOL(x))
    z[seq(1,2*NROW(x),by=2),1:NCOL(x)]=x
    z[seq(2,2*NROW(x),by=2),-(1:NCOL(x))]=x
    
    #multinomial distribution
    wcl=array(0,dim=c(NROW(x),4,length(lam)))
    wcl[,1,]=weight*ppro[,2,]*(1-ppro[,2,])
    wcl[,2,]=-weight*ppro[,2,]*ppro[,3,];wcl[,3,]=-weight*ppro[,2,]*ppro[,3,];
    wcl[,4,]=weight*ppro[,3,]*(1-ppro[,3,])
    

    trhat=foreach(ilam=1:length(lam),.combine='c')%dopar%{
      library(Matrix)#bdiag_m
      library(psych)
      cat("ilam=",ilam,"\n")
      wclist=lapply(as.list(1:dim(wcl)[1]),function(nr,wcc){
        
        return(matrix(wcc[nr,],nrow=2,ncol=2))
      },wcl[,,ilam])
      #wcliat is a list of nrow(x), each is a matrix of 2 times 2
      W=as.matrix(bdiag_m(wclist))
      
      lami=diag(lam[ilam],2*NCOL(x))
      
      zwz=t(z)%*%W%*%z
      
      trhat=tr(solve(zwz+NROW(x)*lami)%*%zwz)
      return(trhat)
    }
  
    #W is a list of length lambda
    #trhat is a vector of length(lambda)
  }  
  logL=tre$nulldev*(1-tre$dev.ratio)
  
  absbeta=lapply(tre$beta,function(tb){
    bet=apply(tb,2,function(tbb){return(tbb%*%tbb)})
    return(bet)
  })
  #absbeta is a list of length levels, each is a vector of length(lambda)
  absbeta=matrix(unlist(absbeta),nrow=length(lam))
  absbeta=rowSums(absbeta)
 
  bic=logL+lam*absbeta+log(NROW(x))*trhat/NROW(x)
  return(list(bic=bic,lam.min=tre$lambda[which.min(bic)],fit=tre))
}


