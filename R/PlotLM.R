PlotLMContinuous <- function(model,terms){
  
  if (length(terms)==1){
    
    newdat<-data.frame(id = 1:100)
    newdat[,terms[1]] <- seq(from=min(model$data[,terms[1]]),to=max(model$data[,terms[1]]),length.out = 100)
    # newdat[,strsplit(harvestModel$final.call,'~')[[1]][1]]<-0
    
    
  } else if (length(terms)==2){
    stop("Error: plotting of two-way interactions between continuous effects not yet supported")
  } else {
    stop("Error: this function can only plot single effects or two-way interactions")
  }
  
  return(newdat)
  
}