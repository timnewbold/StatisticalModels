.Log <- function(...) {
  if(getOption('SMVerbose', TRUE)) cat(...)
}

.SetPar<-function(parList){
  
  stopifnot(is.list(parList))
  
  for (p in names(parList)){
    eval(parse(text=gsub("x",p,"par(x=parList$x)")))
  }
  
}

.ConstructCall<-function(responseVar,fixedStruct,randomStruct){
  return(paste(responseVar,"~",fixedStruct,"+",randomStruct,sep=""))
}