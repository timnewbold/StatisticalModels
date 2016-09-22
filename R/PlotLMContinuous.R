PlotLMContinuous <- function(model,terms,se.fit=FALSE,seMultiplier=1.96,
                             xlim=NULL,ylim=NULL,xlab=NULL,ylab=NULL,
                             params=list(),
                             byFactor=NULL,
                             yTransform=function(x){return(x)},
                             xTransform=function(x){return(x)},
                             log=""){
  
  stopifnot(is.LM(model))
  
  allTerms<-names(attr(model$model$terms,"dataClasses"))[-1]
  allTerms<-gsub("poly[(]","",gsub(", [0-9]{1,2}[)]","",allTerms))
  
  nd<-data.frame(id = 1:100)
  
  for (t in allTerms){
    if (t %in% terms){
      nd[,t] <- seq(from=min(model$data[,t]),to=max(model$data[,t]),length.out = 100)
    } else {
      if(is.numeric(model$data[,t])){
        nd[,t] <- rep(median(model$data[,t]),100)
      } else if (is.factor(model$data[,t])){
        nd[,t] <- factor(rep(levels(model$data[,t])[1],100),
                             levels=levels(model$data[,t]))
      }
      
    }
  }
  
  newdat<-list()
  if (is.null(byFactor)){
    newdat[[1]]<-nd
  } else {
    i<-1
    for (l in levels(model$data[,byFactor])){
      newdat[[i]]<-nd
      nd[,byFactor]<-l
      i <- i+1
    }
  }
  
  if (length(terms)==1){
    if (is.null(xlim)){
      xlim <- xTransform(range(newdat[[1]][,terms[1]]))
    }
    if (is.null(ylim)){
      preds <- predict(model$model,se.fit=se.fit)
      if (se.fit){
        y.minus <- preds$fit - seMultiplier * preds$se.fit
        y.plus <- preds$fit + seMultiplier * preds$se.fit
      }
      
      if (se.fit){
        ylim <- yTransform(c(min(y.minus),
                  max(y.plus)))
      } else {
        ylim <- yTransform(range(preds))
      }
    }
    
    if (is.null(xlab)){
      xlab <- terms[1]
    }
    if (is.null(ylab)){
      ylab <- strsplit(model$final.call,'~')[[1]][1]
    }
    
    par(mar=c(3,3,0.5,0.5),
        mgp=c(1.6,0.4,0),
        tck=-0.01,
        las=1)
    
    .SetPar(params)
    
    if (se.fit){
      plot(9e99,9e99,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,log=log)
    } else {
      plot(9e99,9e99,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,log=log)
    }
    
    for (nd in newdat){
      preds <- predict(model$model,se.fit=se.fit,newdata=nd)
      if (se.fit){
        y <- preds$fit
        y.minus <- preds$fit - seMultiplier * preds$se.fit
        y.plus <- preds$fit + seMultiplier * preds$se.fit
      } else {
        y <- preds
      }
      
      y <- yTransform(y)
      if (se.fit) {
        y.minus <- yTransform(y.minus)
        y.plus <- yTransform(y.plus)
      }
      
      
      
      if(se.fit){
        
        X.Vec <- xTransform(c(nd[,terms[1]], max(nd[,terms[1]]), 
                   rev(nd[,terms[1]]), min(nd[,terms[1]])))
        Y.Vec <- c(y.minus, tail(y.plus, 1), rev(y.plus), (y.minus)[1])
        polygon(X.Vec,Y.Vec,col="#00000033",border = NA)
        points(xTransform(nd[,terms[1]]),y,type="l")
        
        
      } else {
        
        
        points(xTransform(nd[,terms[1]]),y,type="l")
      }
      
    }
    
    
  } else if (length(terms)==2){
    stop("Error: plotting of two-way interactions between continuous effects not yet supported")
  } else {
    stop("Error: this function can only plot single effects or two-way interactions")
  }
  
}