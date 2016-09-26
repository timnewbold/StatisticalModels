
R2GLMER<-function(model){
  
  stopifnot(class(model) == "lmerMod")
  
  # Calculated following Nakagawa & Schielzeth (2013) MEE
  VarF <- var(as.vector(fixef(model) %*% t(getME(model,"X"))))
  # Conditional R2 for the full model
  cond.r2<-(VarF + sum(unlist(lapply(VarCorr(model),function(x) return(x[1])))))/
    (VarF + sum(unlist(lapply(VarCorr(model),function(x) return(x[1])))) + 
       (attr(VarCorr(model), "sc")^2))
  
  # Marginal R2 for the fixed effects
  marg.r2<-VarF/(VarF + 
                   sum(unlist(lapply(VarCorr(model),function(x) return(x[1])))) +
                   (attr(VarCorr(model), "sc")^2))
  
  return(list(conditional=cond.r2,marginal=marg.r2))
  
}



