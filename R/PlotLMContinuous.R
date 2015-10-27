PlotLMContinuous <- function(model,terms){
  
  stopifnot(is.LM(model))
  
  allTerms<-names(attr(lossModel$model$terms,"dataClasses"))[-1]
  
  newdat<-data.frame(id = 1:100)
  
  for (t in allTerms){
    if (t %in% terms){
      newdat[,t] <- seq(from=min(model$data[,t]),to=max(model$data[,t]),length.out = 100)
    } else {
      newdat[,t] <- 
    }
  }
  
  
  # newdat[,strsplit(harvestModel$final.call,'~')[[1]][1]]<-0
  
  
  if (length(terms)==1){
    
    
    
  } else if (length(terms)==2){
    stop("Error: plotting of two-way interactions between continuous effects not yet supported")
  } else {
    stop("Error: this function can only plot single effects or two-way interactions")
  }
  
  return(newdat)
  
}