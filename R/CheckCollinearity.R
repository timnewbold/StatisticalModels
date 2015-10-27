CheckCollinearity<-function(data,allTerms){
  
  results<-data.frame(matrix(nrow=length(allTerms),ncol=length(allTerms)),row.names=allTerms)
  names(results)<-allTerms
  
  for (t1 in allTerms){
    for (t2 in allTerms){
      eval(substitute(c1<-class(data$x),list(x=t1)))
      eval(substitute(c2<-class(data$x),list(x=t2)))
      
      if (t1 != t2){
        if ((c1=="factor") & (c2=="factor")){
          
        } else if (c1=="factor") {
          eval(substitute(r.squared<-summary(lm(data$y~data$x))$r.squared,list(x=t1,y=t2)))
          eval(substitute(results[x,y]<-r.squared,list(x=t1,y=t2)))
        } else {
          eval(substitute(r.squared<-summary(lm(data$y~data$x))$r.squared,list(x=t2,y=t1)))
          eval(substitute(results[x,y]<-r.squared,list(x=t2,y=t1)))
        }
        
      }
      
      
    }
  }
  
  return(results)
  
}


