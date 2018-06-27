PredictGLMER <- function(model,data,se.fit=FALSE,seMultiplier = 1.96){
  
  # stopifnot((class(model)[1] == "lmerMod") | (class(model)[1] == "glmerMod"))
  
  mm<-model.matrix(terms(model),data)
  
  if(ncol(mm)>length(lme4::fixef(model))){
    mm <- mm[,-which(!(names(mm[1,]) %in% names(lme4::fixef(model$model))))]
  }
  
  y <- mm %*% lme4::fixef(model)
  
  if (se.fit){
    pvar1 <- diag(mm %*% base::tcrossprod(as.matrix(vcov(model)),mm))
    yplus <- y + seMultiplier * sqrt(pvar1)
    yminus <- y - seMultiplier * sqrt(pvar1)
    
    
    return(data.frame(y=y,yplus=yplus,yminus=yminus))
    
  } else {
    return(data.frame(y=y))
  }
  
}