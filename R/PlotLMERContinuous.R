
PlotLMERContinuous<-function(model,data,effects,otherContEffects=character(0),
                          otherFactors=character(0),xlab,ylab,
                          byFactor=NULL,byContEffect=NULL,
                          zlab=NULL,outDir=NULL,logLink="n",plotRug=FALSE,
                          seMultiplier=1.96,
                          params=list(),add=FALSE,xlim=NULL,ylim=NULL,zlim=NULL,
                          line.cols=NULL,line.types=NULL,plotUncertainty=TRUE,
                          nPanels=1,main=NULL,yDiv=1,transformX=FALSE){  
  
  if ((length(effects)>1) && (!is.null(byFactor))){
    stop("If plotting by a factor is specified, only one effect can be plotted")
  }
  
  if ((!is.null(byFactor)) && (!is.null(byContEffect))){
    stop("You cannot specify both a factor and a continuous effect to divide the results")
  }
  
  if ((length(effects)==1) && (is.null(byFactor)) && (is.null(byContEffect))){
    eval(substitute(newdat<-data.frame(seq(from=min(data$x,na.rm=TRUE),to=max(data$x,na.rm=TRUE),length.out=100)),list(x=effects)))
    names(newdat)<-effects
    for (o in otherContEffects){
      newdat[,o] <- median(data[,o])
    }
    
    # Add a column to newdat for each of otherFactors
    for(i in 1:length(otherFactors)) {
      col <- names(otherFactors)[i]
      newdat[,col] <- factor(otherFactors[i], levels=levels(data[,col]))
    }
    
    
    newdat[,names(model@frame)[1]] <- 0
    
    mm<-model.matrix(terms(model),newdat)
    pvar1 <- diag(mm %*% base::tcrossprod(as.matrix(vcov(model)),mm))
    newdat$y<-mm %*% fixef(model)
    
    if (logLink=="e"){
      newdat$yplus<-exp(newdat$y+seMultiplier*sqrt(pvar1))
      newdat$yminus<-exp(newdat$y-seMultiplier*sqrt(pvar1))
      newdat$y<-exp(newdat$y)
    } else if (logLink=="10"){
      newdat$yplus<-10^(newdat$y+seMultiplier*sqrt(pvar1))
      newdat$yminus<-10^(newdat$y-seMultiplier*sqrt(pvar1))
      newdat$y<-10^(newdat$y)
    } else if (logLink=="n") {
      newdat$yplus<-(newdat$y+seMultiplier*sqrt(pvar1))
      newdat$yminus<-(newdat$y-seMultiplier*sqrt(pvar1))
      newdat$y<-(newdat$y)
    } else if (logLink=="b"){
      newdat$yplus<-1/(1+exp(-(newdat$y+seMultiplier*sqrt(pvar1))))
      newdat$yminus<-1/(1+exp(-(newdat$y-seMultiplier*sqrt(pvar1))))
      newdat$y<-1/(1+exp(-(newdat$y)))
    } else {
      stop("Error: the specified log link is not supported")
    }
    
    if (!is.null(outDir)){
      png(paste(outDir,"/",effects," effect_",names(model@frame)[1],".png",
                sep=""),width=22.86/2.54,height=12.57/2.54,units="in",res=1200)
    }
    
    par(mgp=c(2.5,1,0))
    par(mar=c(4,4,1,1))
    par(las=1)
    par(cex.lab=1.5)
    
    for (p in names(params)){
      eval(parse(text=gsub("x",p,"par(x=params$x)")))
    }
    
    if (!transformX){
      xlims <- c(min(newdat[,effects]),max(newdat[,effects]))
    } else {
      xlims <- c(exp(min(newdat[,effects])),exp(max(newdat[,effects])))
    }
    if (!is.null(xlim)){
      xlims <- xlim
    }
    
    
    if (!add){
      if (is.null(ylim)){
        if(!transformX){
          plot(-1e5,-1e5,xlim=xlims,
               ylim=c(min(newdat$yminus/yDiv),max(newdat$yplus/yDiv)),
               xlab=xlab,ylab=ylab,main=main)
        } else {
          plot(9e99,9e99,xlim=xlims,
               ylim=c(min(newdat$yminus/yDiv),max(newdat$yplus/yDiv)),
               xlab=xlab,ylab=ylab,
               log="x",main=main)
        }
        
      } else {
        if (!transformX){
          plot(-1e5,-1e5,xlim=xlims,
               ylim=ylim,
               xlab=xlab,ylab=ylab,main=main)
        } else {
          plot(9e99,9e99,xlim=xlims,
               ylim=ylim,
               xlab=xlab,ylab=ylab,main=main,
               log="x")
        }
        
      }
      
    }
    
    line.col<-ifelse(is.null(line.cols),"#000000",line.cols[1])
    
    if(plotUncertainty){
      if (transformX){
        X.Vec <- exp(c(newdat[,effects], max(newdat[,effects]), 
                       rev(newdat[,effects]), min(newdat[,effects])))
      } else {
        X.Vec <- c(newdat[,effects], max(newdat[,effects]), 
                   rev(newdat[,effects]), min(newdat[,effects]))
      }
      
      Y.Vec <- c(newdat$yminus/yDiv, tail(newdat$yplus/yDiv, 1), rev(newdat$yplus/yDiv), (newdat$yminus/yDiv)[1])
      polygon(X.Vec, Y.Vec, col = paste(line.col,"33",sep=""), border = NA)
      
    }
    
    if(transformX){
      points(exp(newdat[,effects]),
             newdat$y/yDiv,type="l",col=line.col)
    } else {
      points(newdat[,effects],
             newdat$y/yDiv,type="l",col=line.col)
    }
    
    if(plotRug){
      if(transformX){
        rug(eval(substitute(exp(data$x),list(x=effects[1]))))
      } else {
        rug(eval(substitute(data$x,list(x=effects[1]))))
      }
      
    }
    
    if (!is.null(outDir)){
      dev.off()
    }
  } else if ((length(effects)==1) && (!is.null(byFactor))){
    
    if((!is.null(line.cols)) & (length(line.cols)!=length(levels(model@frame[,byFactor])))){
      stop("Specified colours must be of the same number as the number of factor levels")
    }
    if((!is.null(line.types)) & (length(line.types)!=length(levels(model@frame[,byFactor])))){
      stop("Specified line types must be of the same number as the number of factor levels")
    }
    
    
    eval(substitute(newdat<-data.frame(seq(from=min(data$x,na.rm=TRUE),to=max(data$x,na.rm=TRUE),length.out=100)),list(x=effects)))
    names(newdat)<-effects
    
    # Add a column to newdat for each of otherFactors
    for(i in 1:length(otherFactors)) {
      col <- names(otherFactors)[i]
      newdat[,col] <- factor(otherFactors[i], levels=levels(data[,col]))
    }
    
    
    newdat[,names(model@frame)[1]] <- 0
    
    y<-data.frame(x=newdat[,1])
    yplus<-data.frame(x=newdat[,1])
    yminus<-data.frame(x=newdat[,1])
    
    for (l in levels(model@frame[,byFactor])){
      newdat[,byFactor]<-factor(l,levels=levels(model@frame[,byFactor]))
      
      for (o in otherContEffects){
        newdat[,o] <- median(data[,o][(data[,byFactor]==l)])
      }
      
      mm<-model.matrix(terms(model),newdat)
      pvar1 <- diag(mm %*% base::tcrossprod(as.matrix(vcov(model)),mm))
      y[,l]<-mm %*% fixef(model)
      
      if (logLink=="e"){
        yplus[,l]<-exp(y[,l]+seMultiplier*sqrt(pvar1))
        yminus[,l]<-exp(y[,l]-seMultiplier*sqrt(pvar1))
        y[,l]<-exp(y[,l])
      } else if (logLink=="10"){
        yplus[,l]<-10^(y[,l]+seMultiplier*sqrt(pvar1))
        yminus[,l]<-10^(y[,l]-seMultiplier*sqrt(pvar1))
        y[,l]<-10^(y[,l])
      } else if (logLink=="b"){
        yplus[,l]<-1/(1+exp(-(y[,l]+seMultiplier*sqrt(pvar1))))
        yminus[,l]<-1/(1+exp(-(y[,l]-seMultiplier*sqrt(pvar1))))
        y[,l]<-1/(1+exp(-(y[,l])))
      } else if (logLink=="n") {
        yplus[,l]<-(y[,l]+seMultiplier*sqrt(pvar1))
        yminus[,l]<-(y[,l]-seMultiplier*sqrt(pvar1))
      } else {
        stop("Error: the specified log link is not supported")
      }
      
      try(y[y$x>quantile(data[data[,byFactor]==l,][,effects[1]],0.975),][,l]<-NA,silent=TRUE)
      try(yplus[yplus$x>quantile(data[data[,byFactor]==l,][,effects[1]],0.975),][,l]<-NA,silent=TRUE)
      try(yminus[yminus$x>quantile(data[data[,byFactor]==l,][,effects[1]],0.975),][,l]<-NA,silent=TRUE)
      try(y[y$x<quantile(data[data[,byFactor]==l,][,effects[1]],0.025),][,l]<-NA,silent=TRUE)
      try(yplus[yplus$x<quantile(data[data[,byFactor]==l,][,effects[1]],0.025),][,l]<-NA,silent=TRUE)
      try(yminus[yminus$x<quantile(data[data[,byFactor]==l,][,effects[1]],0.025),][,l]<-NA,silent=TRUE)
      
    }
    
    if (!is.null(outDir)){
      png(paste(outDir,"/",effects," effect_",names(model@frame)[1],".png",
                sep=""),width=22.86/2.54,height=12.57/2.54,units="in",res=1200)
    }
    
    par(mgp=c(2.5,1,0))
    par(mar=c(4,4,1,1))
    par(las=1)
    par(cex.lab=1.5)
    if (nPanels>1) par(mfrow=c(1,nPanels))
    
    for (p in names(params)){
      eval(parse(text=gsub("x",p,"par(x=params$x)")))
    }
    
    allLevels<-1:length(levels(model@frame[,byFactor]))
    
    nEffectsPerPlot<-ceiling(length(allLevels)/nPanels)
    effectsPerPlot<-split(allLevels,ceiling(allLevels/nEffectsPerPlot))
    
    if(!is.null(line.cols)){
      cols<-split(line.cols,ceiling(allLevels/nEffectsPerPlot))
    }
    if(!is.null(line.types)){
      line.types<-split(line.types,ceiling(allLevels/nEffectsPerPlot))
    }
    
    for (p in 1:nPanels){
      
      if (!add){
        options(scipen=0)
        
        if(!transformX){
          xlims<-c(min(y$x),max(y$x))
        } else {
          xlims<-c(exp(min(y$x)),exp(max(y$x)))
        }
        if(!is.null(xlim)){
          xlims<-xlim
        }
        
        if (is.null(ylim)){
          if(!transformX){
            plot(-99999,-99999,xlim=xlims,
                 ylim=c(min((yminus/yDiv)[,-1],na.rm=TRUE),
                        max((yplus/yDiv)[,-1],na.rm=TRUE)),
                 xlab=xlab,ylab=ylab,main=main)
          } else {
            plot(9e99,9e99,xlim=xlims,
                 ylim=c(min((yminus/yDiv)[,-1],na.rm=TRUE),
                        max((yplus/yDiv)[,-1],na.rm=TRUE)),
                 xlab=xlab,ylab=ylab,log="x",main=main)
          }
          
        } else {
          if(!transformX){
            plot(-99999,-99999,xlim=xlims,
                 ylim=ylim,
                 xlab=xlab,ylab=ylab,main=main)
          } else {
            plot(9e99,9e99,xlim=xlims,
                 ylim=ylim,
                 xlab=xlab,ylab=ylab,log="x",main=main)
          }
          
        }
        
      }
      
      if (is.null(line.cols)){
        cols<-list()
        cols[1][[1]]<-rep("#000000",length(levels(model@frame[,byFactor])))
      } else {
      }
      
      
      
      o<-c(o,setdiff(1:length(levels(model@frame[,byFactor])),o))
      
      if (plotUncertainty){
        i<-1
        for (l in levels(model@frame[,byFactor])[effectsPerPlot[p][[1]]]){
          y.temp<-y[,c('x',l)]
          y.temp<-na.omit(y.temp)
          yplus.temp<-yplus[,c('x',l)]
          yplus.temp<-na.omit(yplus.temp)
          yminus.temp<-yminus[,c('x',l)]
          yminus.temp<-na.omit(yminus.temp)
          if(transformX){
            X.Vec <- exp(c(y.temp$x, max(y.temp$x), rev(y.temp$x), min(y.temp$x)))
          } else {
            X.Vec <- c(y.temp$x, max(y.temp$x), rev(y.temp$x), min(y.temp$x))
          }
          
          Y.Vec <- c((yminus.temp/yDiv)[,l], tail((yplus.temp/yDiv)[,l], 1), rev((yplus.temp/yDiv)[,l]), (yminus.temp/yDiv)[,l][1])
          polygon(X.Vec, Y.Vec, col = paste(cols[p][[1]][i],"33",sep=""), border = NA)
          i<-i+1
        }
        
      }
      
      i<-1
      for (l in levels(model@frame[,byFactor])[effectsPerPlot[p][[1]]]){
        if(transformX){
          points(exp(y$x),(y/yDiv)[,l],type="l",col=cols[p][[1]][i],lty=ifelse(is.null(line.types),1,line.types[p][[1]][i]))
          if (plotRug){
            rug(exp(data[data[,byFactor]==l,][,effects[1]]),col=cols[p][[1]][i])
          }
        } else {
          points(y$x,(y/yDiv)[,l],type="l",col=cols[p][[1]][i],lty=ifelse(is.null(line.types),1,line.types[p][[1]][i]))
          if (plotRug){
            rug(data[data[,byFactor]==l,][,effects[1]],col=cols[p][[1]][i])
          }
        }
        
        i<-i+1
      }
      
    }
    
  } else if ((length(effects)==1) && (!is.null(byContEffect))){
    
    if((!is.null(line.cols)) & (length(line.cols)!=3)){
      stop("Specified colours must be of the same number as the number of factor levels")
    }
    if((!is.null(line.types)) & (length(line.types)!=3)){
      stop("Specified line types must be of the same number as the number of factor levels")
    }
    
    
    eval(substitute(newdat<-data.frame(seq(from=min(data$x,na.rm=TRUE),to=max(data$x,na.rm=TRUE),length.out=100)),
                    list(x=effects)))
    names(newdat)<-effects
    
    # Add a column to newdat for each of otherFactors
    for(i in 1:length(otherFactors)) {
      col <- names(otherFactors)[i]
      newdat[,col] <- factor(otherFactors[i], levels=levels(data[,col]))
    }
    
    newdat[,names(model@frame)[1]] <- 0
    
    y<-data.frame(x=newdat[,1])
    yplus<-data.frame(x=newdat[,1])
    yminus<-data.frame(x=newdat[,1])
    
    byContEffectVals <- c(quantile(data[,byContEffect],probs = 0.025,na.rm=TRUE),
                          median(data[,byContEffect],na.rm=TRUE),
                          quantile(data[,byContEffect],probs = 0.975,na.rm=TRUE))
    
    binWidth <- 0.1
    prbs <- c(binWidth,0.5-(0.5*binWidth),0.5+(0.5*binWidth),1-binWidth)
    brks <- c(-Inf,quantile(data[,byContEffect],
                     probs = prbs,na.rm = TRUE),Inf)
    brks <- brks + seq_along(brks) * 10 * .Machine$double.eps
    data$ByContEffectTerciles <- cut(data[,byContEffect],
                                            breaks = brks,
                                            labels = c(1,"NA1",2,"NA2",3))
    emptyBin <- any(unlist(lapply(as.list(c(1,2,3)),function(x){
      return(length(which(data$ByContEffectTerciles==x))==0)
      })))
    while(emptyBin){
      binWidth <- binWidth + 0.01
      prbs <- c(binWidth,0.5-(0.5*binWidth),0.5+(0.5*binWidth),1-binWidth)
      brks <- c(-Inf,quantile(data[,byContEffect],
                              probs = prbs,na.rm = TRUE),Inf) + 
        seq_along(brks) * 10 * .Machine$double.eps
      data$ByContEffectTerciles <- cut(data[,byContEffect],
                                       breaks = brks,
                                       labels = c(1,"NA1",2,"NA2",3))
      emptyBin <- any(unlist(lapply(as.list(c(1,2,3)),function(x){
        return(length(which(data$ByContEffectTerciles==x))==0)
      })))
    }
    
    for (i in 1:length(byContEffectVals)){
      newdat[,byContEffect] <- byContEffectVals[i]
      
      for (o in otherContEffects){
        newdat[,o] <- median(data[,o][(data$ByContEffectTerciles==i)],na.rm=TRUE)
      }
      
      mm<-model.matrix(terms(model),newdat)
      pvar1 <- diag(mm %*% base::tcrossprod(as.matrix(vcov(model)),mm))
      y[,paste(i)]<-mm %*% fixef(model)
      
      if (logLink=="e"){
        yplus[,paste(i)]<-exp(y[,paste(i)]+seMultiplier*sqrt(pvar1))
        yminus[,paste(i)]<-exp(y[,paste(i)]-seMultiplier*sqrt(pvar1))
        y[,paste(i)]<-exp(y[,paste(i)])
      } else if (logLink=="10"){
        yplus[,paste(i)]<-10^(y[,paste(i)]+seMultiplier*sqrt(pvar1))
        yminus[,paste(i)]<-10^(y[,paste(i)]-seMultiplier*sqrt(pvar1))
        y[,paste(i)]<-10^(y[,paste(i)])
      } else if (logLink=="b"){
        yplus[,paste(i)]<-1/(1+exp(-(y[,paste(i)]+seMultiplier*sqrt(pvar1))))
        yminus[,paste(i)]<-1/(1+exp(-(y[,paste(i)]-seMultiplier*sqrt(pvar1))))
        y[,paste(i)]<-1/(1+exp(-(y[,paste(i)])))
      } else if (logLink=="n") {
        yplus[,paste(i)]<-(y[,paste(i)]+seMultiplier*sqrt(pvar1))
        yminus[,paste(i)]<-(y[,paste(i)]-seMultiplier*sqrt(pvar1))
      } else {
        stop("Error: the specified log link is not supported")
      }
      
      try(y[y$x>quantile(data[data$ByContEffectTerciles==i,][,effects[1]],0.975),][,paste(i)]<-NA,silent=TRUE)
      try(yplus[yplus$x>quantile(data[data$ByContEffectTerciles==i,][,effects[1]],0.975),][,paste(i)]<-NA,silent=TRUE)
      try(yminus[yminus$x>quantile(data[data$ByContEffectTerciles==i,][,effects[1]],0.975),][,paste(i)]<-NA,silent=TRUE)
      try(y[y$x<quantile(data[data$ByContEffectTerciles==i,][,effects[1]],0.025),][,paste(i)]<-NA,silent=TRUE)
      try(yplus[yplus$x<quantile(data[data$ByContEffectTerciles==i,][,effects[1]],0.025),][,paste(i)]<-NA,silent=TRUE)
      try(yminus[yminus$x<quantile(data[data$ByContEffectTerciles==i,][,effects[1]],0.025),][,paste(i)]<-NA,silent=TRUE)
      
    }
    
    if (!is.null(outDir)){
      png(paste(outDir,"/",effects," effect_",names(model@frame)[1],".png",
                sep=""),width=22.86/2.54,height=12.57/2.54,units="in",res=1200)
    }
    
    par(mgp=c(2.5,1,0))
    par(mar=c(4,4,1,1))
    par(las=1)
    par(cex.lab=1.5)
    
    for (p in names(params)){
      eval(parse(text=gsub("x",p,"par(x=params$x)")))
    }
     
    if (!add){
      options(scipen=0)
      
      if(!transformX){
        xlims<-c(min(y$x),max(y$x))
      } else {
        xlims<-c(exp(min(y$x)),exp(max(y$x)))
      }
      if(!is.null(xlim)){
        xlims<-xlim
      }
      
      if (is.null(ylim)){
        if(!transformX){
          plot(-99999,-99999,xlim=xlims,
               ylim=c(min((yminus/yDiv)[,-1],na.rm=TRUE),
                      max((yplus/yDiv)[,-1],na.rm=TRUE)),
               xlab=xlab,ylab=ylab,main=main)
        } else {
          plot(9e99,9e99,xlim=xlims,
               ylim=c(min((yminus/yDiv)[,-1],na.rm=TRUE),
                      max((yplus/yDiv)[,-1],na.rm=TRUE)),
               xlab=xlab,ylab=ylab,log="x",main=main)
        }
        
      } else {
        if(!transformX){
          plot(-99999,-99999,xlim=xlims,
               ylim=ylim,
               xlab=xlab,ylab=ylab,main=main)
        } else {
          plot(9e99,9e99,xlim=xlims,
               ylim=ylim,
               xlab=xlab,ylab=ylab,log="x",main=main)
        }
        
      }
      
    }
    
    if (is.null(line.cols)){
      cols<-rep("#000000",length(byContEffectVals))
    } else {
      cols <- line.cols
    }
    
    if (plotUncertainty){
      for (i in 1:length(byContEffectVals)){
        y.temp<-y[,c('x',paste(i))]
        y.temp<-na.omit(y.temp)
        yplus.temp<-yplus[,c('x',paste(i))]
        yplus.temp<-na.omit(yplus.temp)
        yminus.temp<-yminus[,c('x',paste(i))]
        yminus.temp<-na.omit(yminus.temp)
        if(transformX){
          X.Vec <- exp(c(y.temp$x, max(y.temp$x), rev(y.temp$x), min(y.temp$x)))
        } else {
          X.Vec <- c(y.temp$x, max(y.temp$x), rev(y.temp$x), min(y.temp$x))
        }
        
        Y.Vec <- c((yminus.temp/yDiv)[,paste(i)], tail((yplus.temp/yDiv)[,paste(i)], 1), rev((yplus.temp/yDiv)[,paste(i)]), (yminus.temp/yDiv)[,paste(i)][1])
        polygon(X.Vec, Y.Vec, col = paste(cols[i],"33",sep=""), border = NA)
        i<-i+1
      }
      
    }
    
    i<-1
    for (i in 1:length(byContEffectVals)){
      if(transformX){
        points(exp(y$x),(y/yDiv)[,paste(i)],type="l",col=cols[i],lty=ifelse(is.null(line.types),1,line.types[i]))
#         if (plotRug){
#           rug(exp(data[data[,byFactor]==l,][,effects[1]]),col=cols[p][[1]][i])
#         }
      } else {
        points(y$x,(y/yDiv)[,paste(i)],type="l",col=cols[i],lty=ifelse(is.null(line.types),1,line.types[i]))
#         if (plotRug){
#           rug(data[data[,byFactor]==l,][,effects[1]],col=cols[p][[1]][i])
#         }
      }
      
    }
      
    
  } else if (length(effects==2)){
    
    if (!is.null(outDir)){
      png(paste(outDir,"/",effects[1]," and ",effects[2]," effect_",names(model@frame)[1],".png",
                sep=""),width=22.86/2.54,height=12.57/2.54,units="in",res=1200)
    }
    
    par(mar=c(2,2,2,2))
    par(mgp=c(2.5,1,0))
    par(las=0)
    par(cex.lab=1.5)
    
    for (p in names(params)){
      eval(parse(text=gsub("x",p,"par(x=params$x)")))
    }
    
    
    eval(substitute(e1<-seq(from=min(data$x,na.rm=TRUE),to=max(data$x,na.rm=TRUE),length=30),list(x=effects[1])))
    eval(substitute(e2<-seq(from=min(data$x,na.rm=TRUE),to=max(data$x,na.rm=TRUE),length=30),list(x=effects[2])))
    z<-array(dim=c(30,30))
    zu<-array(dim=c(30,30))
    zl<-array(dim=c(30,30))
    
    gx<-rep(e1,30)
    gy<-rep(e2,rep(30,30))
    eval(substitute(tf<-matrix(exclude.too.far(gx,gy,data$x,data$y,0.2),
                               byrow=F,ncol=30),list(x=effects[1],y=effects[2])))
    
    for (x in 1:length(e1)){
      for (y in 1:length(e2)){
        v1<-e1[x]
        v2<-e2[y]
        newdat<-data.frame(v1,v2)
        names(newdat)<-effects
        # Add a column to newdat for each other continuous effect
        for (o in otherContEffects){
          newdat[,o] <- median(data[,o])
        }
        # Add a column to newdat for each other factor
        for(i in 1:length(otherFactors)) {
          col <- names(otherFactors)[i]
          newdat[,col] <- factor(otherFactors[i], levels=levels(data[,col]))
        }
        newdat[,names(model@frame)[1]] <- 0
        
        mm<-model.matrix(terms(model),newdat)
        pvar1 <- diag(mm %*% base::tcrossprod(as.matrix(vcov(model)),mm))
        if(tf[x,y]==FALSE){
          pred<-mm %*% fixef(model)
          if (logLink=="e"){
            z[x,y]<-exp(pred)
            zu[x,y]<-exp(pred+seMultiplier*sqrt(pvar1))
            zl[x,y]<-exp(pred-seMultiplier*sqrt(pvar1))
          } else if (logLink=="10"){
            z[x,y]<-10^(pred)
            zu[x,y]<-10^(pred+seMultiplier*sqrt(pvar1))
            zl[x,y]<-10^(pred-seMultiplier*sqrt(pvar1))
          } else if (logLink=="b"){
            z[x,y]<-1/(1+exp(-(pred)))
            zu[x,y]<-1/(1+exp(-(pred+seMultiplier*sqrt(pvar1))))
            zl[x,y]<-1/(1+exp(-(pred-seMultiplier*sqrt(pvar1))))
          } else if (logLink=="n") {
            z[x,y]<-mm %*% fixef(model)
            zu[x,y]<-z[x,y]+seMultiplier*sqrt(pvar1)
            zl[x,y]<-z[x,y]-seMultiplier*sqrt(pvar1)
          } else {
            stop("Error: the specified log link is not supported")
          }
        } else {
          z[x,y]<-NA
          zu[x,y]<-NA
          zl[x,y]<-NA
        }
      }
    }
    
    z[(zu-zl)>(quantile(z,0.975,na.rm=TRUE)-min(quantile(z,0.025,na.rm=TRUE)))]<-NA
    
    #     ry.cols<-colorRampPalette(c("red","orange","yellow","green"))
    #     colour<-ry.cols(64)
    colour<-brewer.pal(n=11,name="RdYlBu")
    nrz <- nrow(z)
    ncz <- ncol(z)
    zfacet <- z[-1, -1] + z[-1, -ncz] + z[-nrz, -1] + z[-nrz, -ncz]
    facetcol <- cut(zfacet, 11)
    if (plotUncertainty){
      plims=c(min(zl,na.rm=TRUE),max(zu,na.rm=TRUE))
    } else {
      plims=c(min(z,na.rm=TRUE),max(z,na.rm=TRUE))
    }
    
    if(!is.null(zlim)){
      plims<-zlim
    }
    
    if (plotUncertainty){
      persp(e1,e2,zl,theta=45,phi=35,xlab=NA,ylab=NA,zlab=NA,zlim=plims,
            col=NA,border="#CCCCCC44",axes=FALSE)
      par(new=T)
    }
    persp(e1,e2,z,theta=45,phi=35,xlab=paste("\n",xlab,sep=""),
          ylab=paste("\n",ylab,sep=""),zlab=paste("\n\n",zlab,sep=""),
          zlim=plims,ticktype="detailed",nticks=6,col=colour[facetcol])
    if(plotUncertainty){
      par(new=T)
      persp(e1,e2,zu,theta=45,phi=35,xlab=NA,ylab=NA,zlab=NA,zlim=plims,
            col=NA,border="#CCCCCC44",axes=FALSE)
      
    }
    
    
    
    if (!is.null(outDir)){
      dev.off()
    }
    
  } else {
    stop("Error: a maximum of 2 effects can be specified")
  }
  
  
}