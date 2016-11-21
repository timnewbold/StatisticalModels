GLMEROverdispersion<-function(model){
  stopifnot((class(model)=="lmerMod") | (class(model)=="glmerMod"))
  
  vpars <- function(m) {
    nrow(m)*(nrow(m)+1)/2
  }
  model.df <- sum(sapply(VarCorr(model),vpars))+length(fixef(model))
  resid.df <- nrow(model.frame(model))-model.df
  
  resids<-residuals(model,type="pearson")
  resid.dev<-sum(residuals(model,type="pearson")^2)
  
  ratio<-resid.dev/resid.df
  
  p<-pchisq(resid.dev,df=resid.df,lower.tail=FALSE)
  
  return(list(residDev=resid.dev,residDF=resid.df,ratio=ratio,P.ChiSq=p))
}