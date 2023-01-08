#'varSelRF_mis
#'@description FUNCTION: Using random forests to select important variables. Mainly written by Diaz
#'@param xdata: predictor matrix.
#'@param class: the response. Must be a factor
#'@param c.sd: The factor that multiplies the sd. to decide on stopping the iterations or choosing the final solution.
#'@param mtryFactor: The multiplication factor
#'@param ntree: The number of trees to use for the first forest; same as ntree for randomForest
#'@param  ntreeIterat: The number of trees to use (ntree of randomForest) for all additional forests
#'@param vars.drop.num: The number of variables to exclude at each iteration.
#'@param vars.drop.frac: The fraction of variables, from those in the previous forest, to exclude at each iteration.
#'@param whole.range: If TRUE continue dropping variables until a forest with only two variables is built, and choose the best model from the complete series of models. If FALSE, stop the iterations if the current OOB error becomes larger than the initial OOB error
#'@param recompute.var.imp: If TRUE recompute variable importances at each new iteration.
#'@param verbose: Give more information about what is being done
#'@param returnFirstForest: If TRUE the random forest from the complete set of variables is returned.
#'@param fitted.rf: An (optional) object of class randomForest previously fitted.
#'@param keep.forest: Same argument as in randomForest function.
#'@param caseweights: each weight corresponds to the individual.
#'@import ranger
#'@import randomForest
#'@import varSelRF
#' @return a list names "varSelRF"
#' @export


varSelRF_mis <- function(xdata, Class,
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
                     caseweights) {
#FUNCTION: Using random forests to select important variables.
#Argments:
#xdata:predictor matrix.
#class: the response.
#Other variables: Same as varSelRF() function and ranger() function in R packages (VarSelRF and ranger).
  
  if(!is.factor(Class))
    stop("Class should be a factor")
  if( (is.null(vars.drop.num) & is.null(vars.drop.frac)) |
      (!is.null(vars.drop.num) & !is.null(vars.drop.frac)))
    stop("One (and only one) of vars.drop.frac and vars.drop.num must be NULL and the other set")
  
  max.num.steps <- dim(xdata)[2]
  num.subjects <- dim(xdata)[1]
  
  if(is.null(colnames(xdata)))
    colnames(xdata) <- paste("v", 1:dim(xdata)[2], sep ="")
  
  ##oversize the vectors; will prune later.
  n.vars <- vars <- OOB.rf <- OOB.sd <- rep(NA, max.num.steps)
  
  oobError <- function(rf) { ## should not be exported in the namespace.
    ## The out of bag error
    
    ooo <- rf$confusion.matrix
    s.ooo <- sum(ooo)
    diag(ooo) <- 0
    sum(ooo)/s.ooo
  }
  
  
  if(!is.null(fitted.rf)) {
   # if(ncol(fitted.rf$importance) < 2)
   #    stop("The fitted rf was not fitted with importance = TRUE")
    n.ntree <- fitted.rf$num.tree
    mtry <- fitted.rf$mtry
    n.mtryFactor <- mtry/sqrt(ncol(xdata))
    if((n.ntree != ntree) | (n.mtryFactor != mtryFactor))
      warning("Using as ntree and mtry the parameters obtained from fitted.rf",
              immediate.= TRUE)
    ntree <- n.ntree
    mtryFactor <- n.mtryFactor
    rm(n.ntree, n.mtryFactor)
    rf <- fitted.rf
  } else {
    mtry <- floor(sqrt(ncol(xdata)) * mtryFactor)
    data2=data.frame(phe=Class,xdata)
    rf<-ranger(phe~.,data=data2,importance="impurity",case.weights=caseweights,
               num.trees=ntree,mtry=mtry,write.forest=keep.forest)
    #rf <- randomForest(x = xdata, y = Class,
    #                   ntree = ntree, mtry = mtry,
    #                   importance = TRUE,
    #                   keep.forest = keep.forest)
  }
  
  if(returnFirstForest)
    FirstForest <- rf
  else
    FirstForest <- NULL
  m.iterated.ob.error <- m.initial.ob.error <- oobError(rf)
  sd.iterated.ob.error <- sd.initial.ob.error <-
    sqrt(m.iterated.ob.error * (1 - m.iterated.ob.error) *
           (1/num.subjects))
  
  if(verbose) {
    print(paste("Initial OOB error: mean = ",
                round(m.initial.ob.error, 4),
                "; sd = ", round(sd.initial.ob.error, 4), sep = ""))
  }
  
  importances <- importance(rf)#importance() in ranger
  selected.vars <- order(importances, decreasing = TRUE)
  ordered.importances <- importances[selected.vars]
  
  initialImportances <- importances
  initialOrderedImportances <- ordered.importances
  
  j <- 1
  n.vars[j] <- dim(xdata)[2] 
  vars[j] <- paste(colnames(xdata), collapse = " + ")
  OOB.rf[j] <- m.iterated.ob.error
  OOB.sd[j] <- sd.iterated.ob.error
  
  var.simplify <- TRUE
  
  while(var.simplify) {
    if (verbose){
      print("gc inside loop of varSelRF")
      print(gc())
    } else {
      gc()
    }
    
    last.rf <- rf 
    last.vars <- selected.vars
    previous.m.error <- m.iterated.ob.error
    previous.sd.error <- sd.iterated.ob.error
    
    ## If this is left as only
    ## "if(length(selected.vars) <= 2) var.simplify <- FALSE"
    ## under the if((length(selected.vars) < 2) | (any(selected.vars < 1))),
    ## as it used to be, then we fit a 2 model var, which might be
    ## better or as good as others, but as we never re-enter,
    ## we cannot return it, even if we see it in the history.
    
    ## Alternatively, we cannot just convert
    ## "if((length(selected.vars) < 2) | (any(selected.vars < 1))"
    ## to the <= 2, as we then would re-enter many times because
    ## of the way selected.vars <- selected.vars[1:2] when
    ## num.vars < (vars.drop + 2)
    
    ## This way, we enter just to set last.rf, last.vars and
    ## we bail out
    
    if(length(selected.vars) <= 2) {
      var.simplify <- FALSE
      break
    }
    
    
    if(recompute.var.imp & (j > 1)) {
      importances <- importance(rf)
      tmp.order <- order(importances, decreasing = TRUE)
      selected.vars <- selected.vars[tmp.order]
      ordered.importances <- importances[tmp.order]
    }
    
    num.vars <- length(selected.vars)
    
    if(is.null(vars.drop.num))
      vars.drop <- round(num.vars * vars.drop.frac)
    else vars.drop <- vars.drop.num
    
    if(num.vars >= (vars.drop + 2)) {
      ## prevent infinite looping when num.vars is, say, 3 and
      ## vars.drop.frac < 0.17
      ## We must drop a variable for sure, since > 2 and we
      ## are still simplifying
      ## Alternatively, use ceiling instead of round, above.
      if(vars.drop == 0) {
        vars.drop <- 1
        if( (num.vars - vars.drop) < 1)
          stop("vars.drop = 0 and num.vars -vars.drop < 1!")
      }
      selected.vars <- selected.vars[1: (num.vars - vars.drop)]
      ordered.importances <- ordered.importances[1: (num.vars - vars.drop)]
    }
    else {
      selected.vars <- selected.vars[1:2]
      ordered.importances <- ordered.importances[1:2]
    }
    
    ## couldn't we eliminate the following?
    if((length(selected.vars) < 2) | (any(selected.vars < 1))) {
      var.simplify <- FALSE
      break
    }
    
    
    
    mtry <- floor(sqrt(length(selected.vars)) * mtryFactor)
    if(mtry > length(selected.vars)) mtry <- length(selected.vars)
    
    if(recompute.var.imp){
      data2=data.frame(phe=Class,xdata[,selected.vars])
      rf<-ranger(phe~.,data=data2,importance="impurity",case.weights=caseweights,
                 num.trees=ntree,mtry=mtry,write.forest=keep.forest)
      # rf <- randomForest(x = xdata[, selected.vars],
      #                     y = Class, importance= TRUE,
      #                    ntree = ntree, mtry = mtry,
      #                     keep.forest = keep.forest) 
    }else{
      data2=data.frame(phe=Class,xdata[,selected.vars])
      rf<-ranger(phe~.,data=data2,importance="impurity",case.weights=caseweights,
                 num.trees=ntreeIterat,mtry=mtry,write.forest=keep.forest)
     # rf <- randomForest(x = xdata[, selected.vars],
     #                     y = Class, importance= FALSE,
      #                   ntree = ntreeIterat, mtry = mtry,
      #                   keep.forest = keep.forest)
    }
    
    m.iterated.ob.error <- oobError(rf)
    sd.iterated.ob.error <-
      sqrt(m.iterated.ob.error * (1 - m.iterated.ob.error) *
             (1/num.subjects))
    
    if(verbose) {
      print(paste("..... iteration ", j, "; OOB error: mean = ",
                  round(m.iterated.ob.error, 4),
                  "; sd = ", round(sd.iterated.ob.error, 4),
                  "; num. vars = ", length(selected.vars), 
                  sep = ""))
    }
    j <- j + 1
    
    
    n.vars[j] <- length(selected.vars)
    vars[j] <- paste(colnames(xdata)[selected.vars],
                     collapse = " + ")
    OOB.rf[j] <- m.iterated.ob.error
    OOB.sd[j] <- sd.iterated.ob.error
    
    
    if(!whole.range &
       (
         (m.iterated.ob.error >
          (m.initial.ob.error + c.sd*sd.initial.ob.error))
         |
         (m.iterated.ob.error >
          (previous.m.error + c.sd*previous.sd.error)))
    )
      var.simplify <- FALSE
  }
  
  if (!whole.range) {
    if(!is.null(colnames(xdata)))
      selected.vars <- sort(colnames(xdata)[last.vars])
    else
      selected.vars <- last.vars
    
    out <- list(selec.history = data.frame(
      Number.Variables = n.vars,
      Vars.in.Forest = vars,
      OOB = OOB.rf,
      sd.OOB = OOB.sd)[1:j,],
      rf.model = last.rf,
      selected.vars = selected.vars,
      selected.model =  paste(selected.vars, collapse = " + "),
      best.model.nvars = length(selected.vars),
      initialImportances = initialImportances,
      initialOrderedImportances = initialOrderedImportances,
      ntree = ntree,
      ntreeIterat = ntreeIterat,
      mtryFactor = mtryFactor,
      #                    mtry = mtry,
      firstForest = FirstForest)
    class(out) <- "varSelRF"
    return(out)
  }
  else {
    ## Prune the too long vectors created at begin.
    ## not needed above, because we select the 1:j rows
    ## of the return matrix selec.history.
    n.vars <- n.vars[1:j]
    vars <- vars[1:j]
    OOB.rf<- OOB.rf[1:j]
    OOB.sd <- OOB.sd[1:j]
    min.oob.ci <- min(OOB.rf) + c.sd * OOB.sd[which.min(OOB.rf)]
    best.pos <-
      which(OOB.rf <= min.oob.ci)[which.min(n.vars[which(OOB.rf <= min.oob.ci)])]
    
    selected.vars <- sort(unlist(strsplit(vars[best.pos],
                                          " + ", fixed = TRUE)))
    out <- list(selec.history = data.frame(
      Number.Variables = n.vars,
      Vars.in.Forest = vars,
      OOB = OOB.rf,
      sd.OOB = OOB.sd),
      rf.model = NA,
      selected.vars = selected.vars,
      selected.model = paste(selected.vars, collapse = " + "),
      best.model.nvars = n.vars[best.pos],
      initialImportances = initialImportances,
      initialOrderedImportances = initialOrderedImportances,
      ntree = ntree,
      ntreeIterat = ntreeIterat,
      mtryFactor = mtryFactor,
      ##mtry = mtry,
      firstForest = FirstForest)
    class(out) <- "varSelRF"
    return(out)
  }
}
