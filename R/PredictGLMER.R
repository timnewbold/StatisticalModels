PredictGLMER <- function(model,data,se.fit=FALSE,seMultiplier = 1.96,randEffs = FALSE){
  
  # stopifnot((class(model)[1] == "lmerMod") | (class(model)[1] == "glmerMod"))
  
  # Get names of model factors from prediction data frame
  model.factors <- names(which(sapply(data,is.factor)))
  
  # For each factor, check that levels in prediction data frame match levels
  # in model data frame
  invisible(sapply(X = model.factors,FUN = function(fac){
    stopifnot(all(levels(data[[fac]]) == levels(model[[fac]])))
  }))
  
  mm<-model.matrix(terms(model),data)
  
  if(ncol(mm)>length(lme4::fixef(model))){
    mm <- mm[,-which(!(names(mm[1,]) %in% names(
      lme4::fixef(model))))]
  }
  
  y <- mm %*% lme4::fixef(model)
  
  if (se.fit){
    pvar1 <- diag(mm %*% base::tcrossprod(as.matrix(vcov(model)),mm))
    
    if(randEffs){
      cat('WARNING: May not work for models with random slopes\n')
      pvar1 <- pvar1 + sum(unlist(lapply(X = VarCorr(sr1$model),
                                         FUN = function(x) return(x[1]))))
    }
    
    yplus <- y + seMultiplier * sqrt(pvar1)
    yminus <- y - seMultiplier * sqrt(pvar1)
    
    
    
    
    return(data.frame(y=y,yplus=yplus,yminus=yminus))
    
  } else {
    return(data.frame(y=y))
  }
  
}