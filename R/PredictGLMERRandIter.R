PredictGLMERRandIter <- function(model,data,nIters=1000){
  
  # stopifnot((class(model)[1] == "lmerMod") | (class(model)[1] == "glmerMod"))
  
  # Get names of model factors from prediction data frame
  model.factors <- names(which(sapply(data,is.factor)))
  
  # For each factor, check that levels in prediction data frame match levels
  # in model data frame
  invisible(sapply(X = model.factors,FUN = function(fac){
    stopifnot(all(levels(data[[fac]]) == levels(model@frame[[fac]])))
  }))
  
  preds <- sapply(X = 1:nIters,FUN = function(i){
    
    mm<-model.matrix(terms(model),data)
    
    if(ncol(mm)>length(lme4::fixef(model))){
      mm <- mm[,-which(!(names(mm[1,]) %in% names(
        lme4::fixef(model))))]
    }
    
    fe.draw <- mvrnorm(n = 1,mu = fixef(model),Sigma = vcov(model))
    
    y <- mm %*% fe.draw
    
  })
  
  return(preds)
}