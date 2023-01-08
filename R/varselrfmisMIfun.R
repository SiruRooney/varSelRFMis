#' detctnafun
#' @description FUNCTION: To generate a list of which includes a vector for indicating whether SNPs and phenotype are missing and a matrix for indicating what components in genomic matrix and phenotype are missing.
#' @param predt: dataframe of SNPs
#' @param resps: phenotype factor
#' @return a list which contains "detna" indicating which individual have missing values and "misind" indicating which entry is missing in the dataset
#' @export
detctnafun<-function(predt,reps){

  detna=list()
  detna$zsnp=(!complete.cases(t(predt)))
  detna$y=any(is.na(reps))
  
  misind=list()
  misind$zsnp=matrix(0,nrow=NROW(predt),ncol=NCOL(predt))
  misind$zsnp[,(detna$zsnp)][is.na(predt[,(detna$zsnp)])]=1
  misind$y=rep(0,length(reps))
  if(detna$y) misind$y[is.na(reps)]=1
  
  misind$zsnp=lapply(as.list(1:dim(misind$zsnp)[2]),function(x,y){
    return(as.factor(y[,x]))
  },misind$zsnp)
  
  names(misind$zsnp)=paste(rep("ind",length(detna$zsnp)),colnames(predt),sep="")
 
  misind$zsnp=as.data.frame(misind$zsnp)
  
  misind$y=as.factor(misind$y)
  names(misind)[2]="indy"
  return(list(detna=detna,misind=misind))
  
}

#' dummytran 
#' @description FUNCTION: To convert each SNP vector into the dummy variable
#' @param cl: a vector of SNP
#' @return a matrix which contains the transformed dummy data
#' @export

dummytran <- function(cl){
  cl=as.factor(cl)
  n <- length(cl)
  x <- matrix(0, n, length(levels(cl)))
  x[(1L:n) + n * (as.integer(cl) - 1L)] <- 1
  dimnames(x) <- list(names(cl), levels(cl))
  
  return(x[,-1])
  
}

#'rfmisinitfun
#'@description FUNCTION: To initially impute missing values which are present in the genomic data and phenotype.
#'@param predt: dataframe of SNPs.
#'@param resps: phenotype factor.
#'@param misind: a list of length 2. The first part is a matrix which is composed of "1" and "0" indicating whether each component is missing. The second part is a vector of indicating which component in the phenotype is missing.
#'@param detna: a list of length 2. The first vector indicates which SNP contain missing values, and the second vector indicates whether phenotype is missing.
#'@import randomForest
#'@return a list which contains initialized "predimp" and "predy"
#'@export

rfmisinitfun<-function(predt,resps,misind,detna){

  impc=lapply(predt[,detna$zsnp],function(co){
    wy=(is.na(co))
    yry=co[!is.na(co)]
    co[is.na(co)]=sample(yry,size=sum(wy),replace=TRUE)
    return(co)
  })
  
  impy=sample(levels(resps),size=sum(is.na(resps)),replace=TRUE)
  
  predt[,detna$zsnp]=impc
  resps[misind$indy=="1"]=impy
  #require(randomForest)
  for (i.initm in which(detna$zsnp==TRUE)){
    cat("i.initm",i.initm,"\n")
    
    ridsnp=randomForest(predt[(misind$zsnp[,i.initm]!="1"),-i.initm],as.factor(predt[,i.initm][misind$zsnp[,i.initm]!="1"]))
    
    cc=predict(ridsnp,predt[(misind$zsnp[,i.initm]=="1"),-i.initm])
    predt[misind$zsnp[,i.initm]=="1",i.initm]=as.character(as.numeric(cc)-1)
    
  }
  
  if(detna$y==TRUE){
    ridsnp=randomForest(predt[(misind$indy!="1"),],resps[misind$indy!="1"])
    cc=predict(ridsnp,predt[(misind$indy=="1"),])
    resps[misind$indy=="1"]=cc}
 
  return(list(predimp=predt,predy=resps))    
}


#'priowinitfun
#'@description FUNCTION: To initially add some weights to individuals which contain few missing values.  
#'@param predt: dataframe of SNPs.
#'@param resps: phenotype factor.
#'@param misind: a list of length 2. The first part is a matrix which is composed of "1" and "0" indicating whether each component is missing. The second part is a vector of indicating which component in the phenotype is missing.
#'@param complet.ind: a vector which indicates whether each individual contain missing values.
#' @import doParallel
#' @import parallel
#' @import iterators
#'
#'@return a list containing "predimp", "predy", "misidimp" which would be fitted by the models for missingness mechanism and "weight".
#'@export

priowinitfun<-function(predt,resps,misind,complet.ind){

  submisind=cbind(misind$zsnp[!complet.ind,],indy=misind$indy[!complet.ind])
  subpredy=cbind(predt[!complet.ind,],predy=resps[!complet.ind])
 
  
 # library(doParallel)
  cl<-makeCluster(30)
  registerDoParallel(cl)
  mulpredt=foreach(irow=1:NROW(submisind),.combine=rbind)%dopar%{
    cat("irow=",irow,"\n")
    temprety=data.frame(rep(subpredy[irow,][[1]],12))
    for(ip in 2:length(subpredy[irow,])){
      #cat("ip=",ip,"\n")
      temprety=cbind(temprety,rep(subpredy[irow,][[ip]],12))
    }
    
    names(temprety)=names(subpredy)

    tempmisid=data.frame(rep(submisind[irow,][[1]],12))
    for(ip in 2:length(submisind[irow,])){
      tempmisid=cbind(tempmisid,rep(submisind[irow,][[ip]],12))
    }
    
    names(tempmisid)=names(submisind)

    if(NCOL(cbind(predt,resps)[,submisind[irow,]=="1"])!=1){
      ttemp=lapply(cbind(predt,resps)[,submisind[irow,]=="1"],function(ip){
        sample(ip[!is.na(ip)],size=12,replace=TRUE)})
    }else{
      ttemp=sample(cbind(predt,resps)[,submisind[irow,]=="1"][!is.na(cbind(predt,resps)[,submisind[irow,]=="1"])],
                   size=12,replace=TRUE)
    }
    

    temprety[,submisind[irow,]=="1"]=data.frame(ttemp)
    
    return(data.frame(temprety,tempmisid,weight=rep(1/12,12)))
  }
  
  stopCluster(cl)
  predtimp=rbind(mulpredt[,1:dim(predt)[2]],cbind(predt[complet.ind,]))
  weight=append(mulpredt$weight,rep(1,dim(predt[complet.ind,])[1]))
  repimp=unlist(list(mulpredt$predy,resps[complet.ind]))
  misidimp=rbind(mulpredt[(dim(predt)[2]+2):(2*dim(predt)[2]+2)],cbind(misind$zsnp[complet.ind,],indy=misind$indy[complet.ind]))
  return(list(predimp=predtimp,predy=repimp,misidimp=misidimp,weight=weight))  
}

#'vsrfmisfun
#'@description FUNCTION: To select important SNPs using the random forests approach introduced by Diaz
#'@param predyimp: The dataframe contains SNP and phenotype data after multiple imputation.
#'@param predymisind: The dataframe contains SNP data, phenotype and missing indicator matrix which only contains SNPs and phenotype with missing values.
#'@param dectnalist: A vector which indicates whether SNPs and phenotype are missing.
#'@param na.phe: A missing indicator of phenotype.
#'@param weight: A vector indicating each weight in every individual.
#'@param varSelRF_mis: A function for selecting important variables.
#' @import doParallel
#' @import parallel
#' @import iterators
#' @import glmnet
#' @import ranger
#' @return a list containing the results of selected variables using random forest
#' @export
vsrfmisfun<-function(predyimp,predymisind,dectnalist,na.phe,weight,varSelRF_mis){

  #require(doParallel)
  cl<-makeCluster(30)
  registerDoParallel(cl)
  if(na.phe==T){
    cat("na.phe",na.phe,"\n")
    rfvsparall=foreach(imis=iter(1:length(dectnalist),checkFunc=function(ii){dectnalist[ii]==T}),.combine=cbind,.packages="ranger")%dopar%{
      cat("imis",imis,"\n")
      dectnatemp=dectnalist
      if(imis!=1){
        dectnatemp[1:(imis-1)]=FALSE 
      }
      cat("dim(dectnatemp)",dim(predyimp[,!dectnatemp]),"\n");
      rfvs=varSelRF_mis(predyimp[,!dectnatemp],predyimp[,imis],ntree=100,c.sd=1,whole.range=FALSE,caseweights=weight)#why here I choose ntree=10, instead of 100. reference EM algorithm 
      return(rfvs)
    } 
  }else{
    cat("na.phe",na.phe,"\n")
    dectnalist2=identity(dectnalist)
    dectnalist2[length(dectnalist)]=TRUE
    rfvsparall=foreach(imis=iter(1:length(dectnalist2),checkFunc=function(ii){dectnalist2[ii]==T}),.combine=cbind,.packages="ranger")%dopar%{
      cat("imis",imis,"\n")
      dectnatemp=dectnalist2
      if(imis!=1){
        dectnatemp[1:(imis-1)]=FALSE 
      }
      cat("dim(dectnatemp)",dim(predyimp[,!dectnatemp]),"\n");
      rfvs=varSelRF_mis(predyimp[,!dectnatemp],predyimp[,imis],ntree=100,c.sd=1,whole.range=FALSE,caseweights=weight)#why here I choose ntree=10, instead of 100. reference EM algorithm 
      
      return(rfvs)
    } 
    
  }
  
  rfvsparallmisind=foreach(imis=1:length(dectnalist[dectnalist==T]),.combine=cbind,.packages="ranger")%dopar%{
    cat("imis",imis,"\n") 
    cat("dim(dectnatemp)",dim(predymisind[,1:(imis+length(dectnalist)-1)]),"\n");
    rfvs=varSelRF_mis(predymisind[,1:(imis+length(dectnalist)-1)],predymisind[,imis+length(dectnalist)],ntree=100,c.sd=1,whole.range=FALSE,caseweights=weight)#why here I choose ntree=10, instead of 100. reference EM algorithm
    
    return(rfvs)
  }
  stopCluster(cl)
  
  #rfvsparall is a matrix of size 11*length(miscovariates+misphenotype);
  #rfvsparallmisind is a matrix of size 11*length(miscovariates+misphenotype)
  return(list(rfvsparall=rfvsparall,rfvsparallmisind=rfvsparallmisind))
}
 
#'imputerfmisfun
#'@description FUNCTION: To apply glmnet function to choose the optimal tunning parameter in the constructed models.
#'@param predimp.w: A dataframe of SNPs after multiple imputation at last iteration.
#'@param predy.w: A vector of phenotype after multiple imputation at last iteration.
#'@param misindcom.tw: A missing indicator matrix for a dataframe of SNPs and phenotype after multiple imputation.
#'@param dectnalist: A vector which indicates whether SNPs and phenotype are missing.
#'@param weight: A vector indicating each weight in every individual.
#'@param itvs: Count for the number of the loop.
#'@param dummytran: A function which converts each SNP vector into the dummy variable
#'@param bicglmnet: A function which chooses an optimal tunning parameter using the model criterion BIC.
#'@param bdiag_m: A function which generates a block diagonal matrix.
#'@param indmisfun: A function which indicates the value of a SNP ("0", "1" or "2").
#'@param ridgselectedvars: A list which contains selected important variables in a series of models for missing SNPs, corresponding missing indicators and phenotype.
#' @import doParallel
#' @import parallel
#' @import iterators
#' @import glmnet
#' @return a list named "ridgimp" containing the results from fitting the ridge regression models
#' @export

imputerfmisfun<-function(predimp.w,predy.w,misindcom.tw,dectnalist,weight,itvs,optmethod,dummytran,bicglmnet,bdiag_m,indmisfun,ridgselectedvars){

  predyimp.w=data.frame(predimp.w,predy=predy.w)
  
  predymisind.w=data.frame(predimp.w,predy=predy.w,misindcom.tw)
 
  initimputecom.w=list(predyimp.w=predyimp.w,predymisind.w=predymisind.w)
  
  
  #require(doParallel)
  cl<-makeCluster(2)
  registerDoParallel(cl)
  ridgimp=foreach(imis2=1:length(ridgselectedvars),.packages="doParallel")%dopar%{
    
    ridgimpute.part=foreach(imis=1:(ifelse((dectnalist[length(dectnalist)]==FALSE)&(imis2==1),length(dectnalist[dectnalist])+1,length(dectnalist[dectnalist]))),.packages="glmnet")%do%{
      cat("imis",imis)
      cat("imis2=",imis2,"\n")
       
      selectedvars.dum=apply(initimputecom.w[[imis2]][,ridgselectedvars[[imis2]][,imis]$selected.vars],2,dummytran)
      
      if(is.list(selectedvars.dum)==T){
        selectedvars.dum=matrix(unlist(selectedvars.dum),nrow=NROW(predimp.w)) 
      }else{
        selectedvars.dum=matrix(selectedvars.dum,nrow=NROW(predimp.w))
      }
      
  
      if(imis2==1){
        cat("when imis2 is",imis2,"\n")
  
        if(imis==length(dectnalist[dectnalist])+1){
          selectedresp=initimputecom.w$predyimp.w[,c(dectnalist[-length(dectnalist)],TRUE)][,imis]
        }else{
          selectedresp=initimputecom.w$predyimp.w[,dectnalist][,imis]
        }
        
      }else{
        cat("when imis2 is",imis2,"\n")

        selectedresp=initimputecom.w$predymisind.w[,imis+NCOL(initimputecom.w$predyimp.w)]
      }
      
      if(optmethod=="bic"){
        if(itvs==1){
          cat("itvs",itvs,"\n")
          ridgmis=bicglmnet(x=selectedvars.dum,y=selectedresp,weight,indmisfun,bdiag_m)
          #ridgmis=cv.glmnet(x=selectedvars.dum,y=selectedresp,weights=weight,alpha=0,nfolds=5,family="multinomial",type.measure="class",parallel=TRUE)
        }else{
          cat("itvs",itvs,"\n")
          ridgmis=bicglmnet(x=selectedvars.dum,y=selectedresp,weight,indmisfun,bdiag_m)
          #ridgmis=cv.glmnet(x=selectedvars.dum,y=selectedresp,weights=weight,alpha=0,family="multinomial",type.measure="class",parallel=TRUE)
        }
      }else{
        if(itvs==1){
          cat("itvs",itvs,"\n")
          ridgmis=cv.glmnet(x=selectedvars.dum,y=selectedresp,weights=weight,alpha=0,nfolds=5,family="multinomial",type.measure="class",parallel=TRUE)
        }else{
          cat("itvs",itvs,"\n")
          ridgmis=cv.glmnet(x=selectedvars.dum,y=selectedresp,weights=weight,alpha=0,family="multinomial",type.measure="class",parallel=TRUE)
        }
      }


      #ridgmis is a list of size 2, each maintains 5 classes named "cv.glmnet"
      return(ridgmis)
    }
    
    return(ridgimpute.part)  
  }
  stopCluster(cl)
  
  return(ridgimp)
}

#' jimpmisfun
#'@description FUNCTION: To multiple impute missing values in SNPs and phenotype. 
#'@param misindcom: A missing indicator matrix for a dataframe of SNPs and phenotype after multiple imputation.
#'@param xycat: a vector indicating the number of levels in each SNPs and phenotype.
#'@param dectnalist: A vector which indicates whether SNPs and phenotype are missing. 
#'@param ncol.predyimp: The number of SNPs and phenotype.
#'@param predymisind: The dataframe contains SNP data, phenotype and missing indicator matrix which only contains SNPs and phenotype with missing values.
#'@param complet.ind: A vector which indicates whether each individual contain missing values.
#'@param ridgimp: A list of ridge models with the optimal tunning parameter.
#'@param ridgselectedvars: A list which contains selected important variables in a series of models for missing SNPs, corresponding missing indicators and phenotype.
#'@param gensetf: A function for generating all possible conbinations of factor levels in an experiment 
#'@param dummytran: A function which converts each SNP vector into the dummy variable.
#'@param indmisfun: A function which indicates the value of a SNP ("0", "1" or "2").
#'@param mulimpmisfun: A function which generates multiple imputation dataset in each individual with missing values.
#' @import doParallel
#' @import parallel
#' @import iterators
#'
#' @return a dataframe named "mdatamis" comprised of multiple imputed data and the observed data with the corresponding weights
#' @export


jimpmisfun<-function(misindcom,xycat,dectnalist,ncol.predyimp,predymisind,complet.ind,ridgimp,ridgselectedvars,optmethod,gensetf,dummytran,indmisfun,mulimpmisfun){

  submisind=misindcom[(complet.ind==F),]
  subpredymis=predymisind[(complet.ind==F),]
  
  cl<-makeCluster(35)
  registerDoParallel(cl)
  mdatamis=foreach(irow=1:NROW(submisind),.combine=rbind,.packages="doParallel")%dopar%{
    cat("irow",irow,"\n")

    multidat=mulimpmisfun(submisind[irow,],xycat,dectnalist,ncol.predyimp,subpredymis[irow,],ridgimp,ridgselectedvars,optmethod,gensetf,dummytran,indmisfun) 
    
    return(multidat)
  }
  stopCluster(cl)
  
  mdatamis=mdatamis[mdatamis$weight!=0,]
  mdatamis=rbind(mdatamis,data.frame(predymisind[(complet.ind==T),],weight=rep(1,dim(predymisind[(complet.ind==T),])[1])))
  return(mdatamis)
}

#' mulimpmisfun
#' @description FUNCTION: To generate multiple imputation dataset in each individual with missing values.
#' @param submisrow: One row of "misindcom" matrix.
#' @param xycat: a vector indicating the number of levels in each SNPs and phenotype.
#' @param dectnalist: A vector which indicates whether SNPs and phenotype are missing. 
#' @param ncol.predyimp: The number of SNPs and phenotype.
#' @param xyindrow: One row of multiple imputed dataframe of SNPs and phenotype.
#' @param ridgimp: A list of ridge models with the optimal tunning parameter.
#' @param ridgselectedvars:  A list which contains selected important variables in a series of models for missing SNPs, corresponding missing indicators and phenotype.
#' @param gensetf: A function for generating all possible conbinations of factor levels in an experiment.
#' @param dummytran: A function which converts each SNP vector into the dummy variable.
#' @param indmisfun: A function which indicates the value of a SNP ("0", "1" or "2").
#' @import doParallel
#' @import parallel
#' @import iterators
#' @import glmnet
#' @return a dataframe named "mulimp" which contains the multiple imputed data with the corresponding weights
#' @export

mulimpmisfun<-function(submisrow,xycat,dectnalist,ncol.predyimp,xyindrow,ridgimp,ridgselectedvars,optmethod,gensetf,dummytran,indmisfun){

  nsim<-prod(xycat[submisrow==1])
  nmis=length(submisrow[submisrow==1])
  
  if(nmis==1){
    temp<-1:xycat[submisrow==1]
  }else{
    temp<-t(gensetf(xycat[submisrow==1]))
  }
  #xyindsim<-matrix(xyindrow,ncol=length(xyindrow),nrow=nsim,dimnames=list(NULL,names(xyindrow)),byrow=T)
  xyindsim=data.frame(rep(xyindrow[[1]],nsim))
  for(ixy in 2:length(xyindrow)){xyindsim=cbind(xyindsim,rep(xyindrow[[ixy]],nsim))}
  names(xyindsim)=names(xyindrow)
  
  temp2=data.frame(lapply(data.frame(temp-1),as.factor))
  xyindsim[,1:ncol.predyimp][,submisrow==1]=temp2
  
  snporind=foreach(indlist=1:length(ridgimp),.packages="doParallel")%dopar%{
    cat("indlist",indlist,"\n")
    cond.part=foreach(iconmis3=1:length(ridgimp[[indlist]]),.packages="glmnet")%dopar%{
      cat("iconmis3",iconmis3,"\n")
      
      newxx=xyindsim[,ridgselectedvars[[indlist]][,iconmis3]$selected.vars]
      if(NCOL(newxx)>1){
        newxxdum=lapply(newxx,dummytran)
      }else{
        newxxdum=dummytran(newxx)
      }
      
      if(is.list(newxxdum)){
        newxxdum=matrix(unlist(newxxdum),nrow=dim(xyindsim)[1])  
      }else{
        newxxdum=matrix(newxxdum,nrow=dim(xyindsim)[1])
      }
      
      if(optmethod=="bic"){
        if(NCOL(newxxdum)>1){
          cond.part=predict(object=ridgimp[[indlist]][[iconmis3]]$fit,newx=newxxdum,s=ridgimp[[indlist]][[iconmis3]]$lam.min,type="response")
        }else{
          cond.part=predict(object=ridgimp[[indlist]][[iconmis3]]$fit,newx=newxxdum,s=ridgimp[[indlist]][[iconmis3]]$lam.min,type="response")
        }
      }else{
        if(NCOL(newxxdum)>1){
          cond.part=predict(object=ridgimp[[indlist]][[iconmis3]],newx=newxxdum,s="lambda.min",type="response")
        }else{
          cond.part=predict(object=ridgimp[[indlist]][[iconmis3]],newx=newxxdum,s="lambda.min",type="response") 
        }
      }

      
      
      return(cond.part)
    }
    
    return(cond.part)
  }
  
  if(dectnalist[length(dectnalist)]==T){
    dectnalist2=dectnalist
  }else{
    dectnalist2=dectnalist
    dectnalist2[length(dectnalist)]=TRUE
  }
  
  poit=list(poitxy=xyindsim[,1:ncol.predyimp][,dectnalist2],poitmis=xyindsim[,-(1:ncol.predyimp)])
  a=lapply(poit$poitxy,indmisfun)  
  b=lapply(poit$poitmis,indmisfun)
  
  
  ja=mapply(function(as,ss){print(dim(ss));print(dim(as));return(rowSums(ss[,,1]*as))},a,snporind[[1]])
  jb=mapply(function(bs,ss){print(dim(ss));print(dim(bs));return(rowSums(ss[,,1]*bs))},b,snporind[[2]])
  jab=apply(cbind(ja,jb),1,prod)
  jall=jab/sum(jab)
  
  
  if(length(jall)>=10){
    sort10=sort(jall,decreasing = TRUE)[1:10] 
    jall10=jall[which(jall%in%sort10)]
    jall10=jall10/sum(jall10)
    xyindsim10=xyindsim[which(jall%in%sort10),]
  }else{
    jall10=jall
    xyindsim10=xyindsim
  }
  mulimp=data.frame(xyindsim10,weight=jall10)
  
  return(mulimp)
}

#' varselrfmisfun
#' @description FUNCTION: To select important SNPs with non-ignorable missing values using random forests.
#' @param predt: dataframe of SNPs.
#' @param resps: phenotype factor.
#' @param misind: A list of length 2. The first part is a matrix which is composed of "1" and "0" indicating whether each component is missing. The second part is a vector of indicating which component in the phenotype is missing.
#' @param detna: A list of length 2. The first vector indicates which SNP contain missing values, and the second vector indicates whether phenotype is missing.
#' @param iter.varsel: The number of iterations.
#' @param initmethod: "priorweight" or "init". "priorweight" is used for initializing SNPs with few missing values. If there are no SNPs with few missing values, choose "init".
#' @param optmethod: "bic" or "cv". It is used for choosing the optimal tunning parameter.
#' @param dummytran: A function which converts each SNP vector into the dummy variable.
#' @param vsrfmisfun: A function which selects important SNPs using the random forests approach introduced by Diaz
#' @param jimpmisfun: A function to multiple impute missing values in SNPs and phenotype. 
#' @param mulimpmisfun: A function which generates multiple imputation dataset in each individual with missing values.
#' @param imputerfmisfun: A function to choose the optimal tunning parameter.
#' @param varSelRF_mis: A function for selecting important variables.
#' @param gensetf: A function for generating all possible conbinations of factor levels in an experiment.
#' @param bicglmnet: A function which chooses an optimal tunning parameter using the model criterion BIC.
#' @param bdiag_m: A function which generates a block diagonal matrix.
#'@import doParallel
#' @return a list containing the final selected result
#' @export

varselrfmisfun<-function(predt,resps,misind,detna,iter.varsel=10,initmethod="priorweight",optmethod="bic",dummytran,vsrfmisfun,jimpmisfun,mulimpmisfun,
                         imputerfmisfun,varSelRF_mis,gensetf,bicglmnet,bdiag_m){

  complet.ind=complete.cases(cbind(predt,resps))
  if(!is.data.frame(class(predt))){predt=as.data.frame(predt)}
  
  if(!is.factor(class(resps))){resps=as.factor(resps)}
  
  if(any(!is.factor(apply(predt,2,class)))) {
    predt=lapply(predt,as.factor);
    predt=as.data.frame(predt)}
  
  
   
  misindcom=data.frame(misind$zsnp,indy=misind$indy)

  ncol.predyimp<-dim(predt)[2]+1
  dectnalist=unlist(detna)
  xycat=c(unlist(lapply(predt,function(cc){return(length(levels(cc)))})),phe=length(levels(resps)))
  predymisind=data.frame(predt,predy=resps,misindcom[,dectnalist==T])
  
  indmisfun<-function(ct){
    
    ct=as.factor(ct)
    indn<-length(ct)
    indx<-matrix(0,indn,length(levels(ct)))
    indx[(1L:indn) +indn*(as.integer(ct) - 1L)] <- 1
    dimnames(indx) <- list(names(ct), levels(ct))
    
    return(indx)
  }
  
  
  
  record.selvars=matrix(NA,nrow=ncol.predyimp,ncol=iter.varsel)
  ridgselectedv.list=NULL;
  
  if(initmethod=="priorweight"){
    rfmisinit<-priowinitfun(predt,resps,misind,complet.ind) 
  }else{
    rfmisinit<-rfmisinitfun(predt,resps,misind,detna);
  }

  for (itvs in 1:iter.varsel){
    print("################################################################")
    
    cat("iter.varsel",itvs,"\n")
    print("################################################################")
    
    if(itvs==1){
      predimp.w=rfmisinit$predimp;
      predy.w=rfmisinit$predy;
      if(initmethod=="priorweight"){
        misindcom.tw=rfmisinit$misidimp[,dectnalist==T];
        weight=rfmisinit$weight  
      }else{
        misindcom.tw=misindcom[,dectnalist==T];
        weight=rep(1,NROW(predimp.w))
      }

      predyimp.w=data.frame(predimp.w,predy=predy.w)
      predymisind.w=data.frame(predimp.w,predy=predy.w,misindcom.tw)
    }else{
      
      predimp.w=jimpmis[,1:(ncol.predyimp-1)]
      predy.w=jimpmis$predy
      misindcom.tw=jimpmis[,-c(1:ncol.predyimp,dim(jimpmis)[2])]
      weight=jimpmis$weight
      predyimp.w=jimpmis[,1:ncol.predyimp]
      predymisind.w=jimpmis[,-dim(jimpmis)[2]]
      
    }
    
    ridgselectedvars<-vsrfmisfun(predyimp.w,predymisind.w,dectnalist,detna$y,weight,varSelRF_mis)
    ridgselectedv.list=cbind(ridgselectedv.list,ridgselectedvars)
    if(detna$y){
      print(ridgselectedvars$rfvsparall[,length(dectnalist[dectnalist==T])]$selected.vars)
      record.selvars[,itvs][1:length(ridgselectedvars$rfvsparall[,length(dectnalist[dectnalist==T])]$selected.vars)]=ridgselectedvars$rfvsparall[,length(dectnalist[dectnalist==TRUE])]$selected.vars
    }else{
      print(ridgselectedvars$rfvsparall[,length(dectnalist[dectnalist])+1]$selected.vars)
      record.selvars[,itvs][1:length(ridgselectedvars$rfvsparall[,length(dectnalist[dectnalist])+1]$selected.vars)]=ridgselectedvars$rfvsparall[,length(dectnalist[dectnalist])+1]$selected.vars
    }
    
    
    ridgimp<-imputerfmisfun(predimp.w,predy.w,misindcom.tw,dectnalist,weight,itvs,optmethod,dummytran,bicglmnet,bdiag_m,indmisfun,ridgselectedvars);
    
    #gibsiterf<-gibiterfmisfun(misindcom,dectnalist,ncol.predyimp,predymisind,complet.ind,ridgimp,ridgselectedvars,
    #                          dummytran,indmisfun,conditprobgibbs,conditproblayer,gibsamprfmisfun,updatedatset);
    
    jimpmis<-jimpmisfun(misindcom,xycat,dectnalist,ncol.predyimp,predymisind,complet.ind,ridgimp,ridgselectedvars,optmethod,gensetf,dummytran,indmisfun,mulimpmisfun)
    
  }
  
  
  return(list(record.selvars=record.selvars,
              ridgselectedv.list=ridgselectedv.list,
              ridgimp=ridgimp,
              multidat=jimpmis,
              weight=jimpmis$weight)) 
}
