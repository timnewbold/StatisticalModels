PlotLMContinuous <- function(model,terms,se.fit=FALSE,seMultiplier=1.96,
                             line.col="#000000",
                             xlim=NULL,ylim=NULL,zlim=NULL,
                             xlab=NULL,ylab=NULL,zlab=NULL,
                             params=list(),
                             byFactor=NULL,
                             yTransform=function(x){return(x)},
                             xTransform=function(x){return(x)},
                             log="",exclude.too.far=0.2){
  
  stopifnot(is.LM(model))
  
  allTerms<-names(attr(model$model$terms,"dataClasses"))[-1]
  allTerms<-gsub("poly[(]","",gsub(", [0-9]{1,2}[)]","",allTerms))
  
  if (length(terms)==1){
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
        polygon(X.Vec,Y.Vec,col=paste0(line.col,"33"),border = NA)
        points(xTransform(nd[,terms[1]]),y,type="l",col=line.col)
        
        
      } else {
        
        
        points(xTransform(nd[,terms[1]]),y,type="l")
      }
      
    }
    
    
  } else if (length(terms)==2){
    if (is.null(xlab)){
      xlab <- terms[1]
    }
    if (is.null(ylab)){
      ylab <- terms[2]
    }
    if (is.null(zlab)){
      zlab <- strsplit(model$final.call,'~')[[1]][1]
    }
    
    e1<-seq(from=min(model$data[,terms[1]],na.rm=TRUE),to=max(model$data[,terms[1]],na.rm=TRUE),length=30)
    e2<-seq(from=min(model$data[,terms[2]],na.rm=TRUE),to=max(model$data[,terms[2]],na.rm=TRUE),length=30)
    
    z<-array(dim=c(30,30))
    zu<-array(dim=c(30,30))
    zl<-array(dim=c(30,30))
    
    gx<-rep(e1,30)
    gy<-rep(e2,rep(30,30))
    
    tf<-matrix(exclude.too.far(gx,gy,model$data[,terms[1]],model$data[,terms[2]],exclude.too.far),
               byrow=F,ncol=30)
    
    for (x in 1:length(e1)){
      for (y in 1:length(e2)){
        
        v1<-e1[x]
        v2<-e2[y]
        nd<-data.frame(v1,v2)
        names(nd)<-terms
        
        for (t in allTerms){
          if (!(t %in% terms)){
            if(is.numeric(model$data[,t])){
              nd[,t] <- median(model$data[,t])
            } else if (is.factor(model$data[,t])){
              nd[,t] <- factor(levels(model$data[,t])[1],
                               levels=levels(model$data[,t]))
            }
            
          }
        }
        
        
        if(tf[x,y]==FALSE){
          preds <- predict(model$model,se.fit=se.fit,newdata=nd)
          if (se.fit){
            y.mean <- preds$fit
            y.minus <- preds$fit - seMultiplier * preds$se.fit
            y.plus <- preds$fit + seMultiplier * preds$se.fit
          } else {
            y.mean <- preds
          }
          
          z[x,y] <- yTransform(y.mean)
          if (se.fit) {
            zl[x,y] <- yTransform(y.minus)
            zu[x,y] <- yTransform(y.plus)
          }
        }
      }
    }
    
    z[(zu-zl)>(quantile(z,0.975,na.rm=TRUE)-min(quantile(z,0.025,na.rm=TRUE)))]<-NA
    z[z > max(predict(model$model))]<-NA
    z[z < min(predict(model$model))]<-NA
    
    #     ry.cols<-colorRampPalette(c("red","orange","yellow","green"))
    #     colour<-ry.cols(64)
    colour<-brewer.pal(n=11,name="RdYlBu")
    nrz <- nrow(z)
    ncz <- ncol(z)
    zfacet <- z[-1, -1] + z[-1, -ncz] + z[-nrz, -1] + z[-nrz, -ncz]
    facetcol <- cut(zfacet, 11)
    if (se.fit){
      plims=c(min(zl,na.rm=TRUE),max(zu,na.rm=TRUE))
    } else {
      plims=c(min(z,na.rm=TRUE),max(z,na.rm=TRUE))
    }
    
    if(!is.null(zlim)){
      plims<-zlim
    }
    
    if (se.fit){
      persp(e1,e2,zl,theta=45,phi=35,xlab=NA,ylab=NA,zlab=NA,zlim=plims,
            col=NA,border="#CCCCCC44",axes=FALSE)
      par(new=T)
    }
    persp(e1,e2,z,theta=45,phi=35,xlab=paste("\n",xlab,sep=""),
          ylab=paste("\n",ylab,sep=""),zlab=paste("\n\n",zlab,sep=""),
          zlim=plims,ticktype="detailed",nticks=6,col=colour[facetcol])
    if(se.fit){
      par(new=T)
      persp(e1,e2,zu,theta=45,phi=35,xlab=NA,ylab=NA,zlab=NA,zlim=plims,
            col=NA,border="#CCCCCC44",axes=FALSE)
      
    }
    
    
  } else {
    stop("Error: this function can only plot single effects or two-way interactions")
  }
  
}