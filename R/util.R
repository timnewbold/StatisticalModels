.Log <- function(...) {
  if(getOption('SMverbose', TRUE)) cat(...)
}

.SetPar<-function(parList){
  
  stopifnot(is.list(parList))
  
  for (p in names(parList)){
    eval(parse(text=gsub("x",p,"par(x=parList$x)")))
  }
  
}