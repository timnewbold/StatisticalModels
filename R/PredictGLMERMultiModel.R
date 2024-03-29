PredictGLMERMultiModel <- function(models,data,nIters=1000){
  
  # Get names of model factors from prediction data frame
  model.factors <- names(which(sapply(data,is.factor)))
  
  # For each candidate model and for each factor, 
  # check that levels in prediction data frame match levels
  # in model data frame
  invisible(lapply(X = models,FUN = function(mod){
    invisible(sapply(X = model.factors,FUN = function(fac){
      stopifnot(all(levels(data[[fac]]) == levels(mod@frame[[fac]])))
    }))
  }))
  
  # Get AIC values of all models
  aics <- unlist(lapply(models,AIC))
  
  # Get differences in model AICs compared to best-fitting model
  aic.diffs <- aics - min(aics)
  
  # Calculate AIC weights for each model
  aic.weights <- exp(-0.5 * aic.diffs)/sum(exp(-0.5 * aic.diffs))
  
  # Iterate to specified number of draws
  preds <- sapply(X = 1:nIters,FUN = function(i){
    
    # Select model at random, weighted by AIC weight
    model <- models[[sample(
      x = 1:length(models),size = 1,
      replace = FALSE,prob = aic.weights)]]
    
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