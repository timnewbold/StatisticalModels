PredictGLMER <- function(model,data,se.fit=FALSE,seMultiplier = 1.96){
  
  stopifnot((class(model)=="lmerMod") | (class(model)=="glmerMod"))
  
  mm<-model.matrix(terms(model),data)
  y <- mm %*% fixef(model)
  
  if (se.fit){
    pvar1 <- diag(mm %*% base::tcrossprod(as.matrix(vcov(model)),mm))
    yplus <- y + seMultiplier * sqrt(pvar1)
    yminus <- y - seMultiplier * sqrt(pvar1)
    
    
    return(data.frame(y=y,yplus=yplus,yminus=yminus))
    
  } else {
    return(data.frame(y=y))
  }
  
}