
PlotGLMERFactor<-function(model,data,responseVar,seMultiplier=1.96,
                     logLink="n",catEffects=NULL,
                     contEffects=list(),contEffectsLabels=NULL,
                     otherCatEffects=list(),
                     otherContEffects=character(0),
                     xtext.srt=0,ylim=NA,yaxp=NULL,order=NULL,rescale=NULL,
                     errbar.cols=NULL,pt.pch=NULL,
                     errbar.lty=1,
                     params=list(),add=FALSE,offset=0,
                     plotLabels=TRUE,cex.txt=NULL,pt.cex=1,
                     pt.bg="white",main=NULL,type="percent"){
  
  if (!is.null(contEffectsLabels)){
    if (length(contEffects)!=length(contEffectsLabels)){
      stop("Labels for continuous effects must be of the same number as the effects")
    }
  }
  
  labels<-character(0)
  coef.labels<-character(0)
  # Get the names of the labels and coefficient names for each factor
  # If align is specified, then use the full set of labels defined above
  for (e in catEffects){
    names<-levels(model@frame[,e])
    labels<-c(labels,names)
    coef.labels<-c(coef.labels,paste(e,names,sep=""))
  }
  
  # Get coefficient and standard error estimates from the model
  o<-match(tolower(coef.labels),tolower(names(fixef(model))))
  y<-fixef(model)[o]
  yplus<-y+se.fixef(model)[o]*seMultiplier
  yminus<-y-se.fixef(model)[o]*seMultiplier
  
  # For each categorical effect, get the reference level
  for (e in catEffects){
    ref.name<-paste(e,levels(model@frame[,e])[1],sep="")
    o<-match(tolower(ref.name),tolower(coef.labels))
    y[o]<-0
    yplus[o]<-0
    yminus[o]<-0
    
  }
  
  if (!is.null(order)){
    labels<-labels[order]
    coef.labels<-coef.labels[order]
    y<-y[order]
    yplus<-yplus[order]
    yminus<-yminus[order]
  }
  
  if(!is.null(rescale)){
    y<-y+rescale
    yplus<-yplus+rescale
    yminus<-yminus+rescale
  }
  
  intercept<-fixef(model)['(Intercept)']
  if (logLink=="e"){
    y<-(exp(y+ifelse(type=="percent",0,intercept)))
    yplus<-(exp(yplus+ifelse(type=="percent",0,intercept)))
    yminus<-(exp(yminus+ifelse(type=="percent",0,intercept)))
  } else if (logLink=="10") {
    y<-(10^(y+ifelse(type=="percent",0,intercept)))
    yplus<-(10^(yplus+ifelse(type=="percent",0,intercept)))
    yminus<-(10^(yminus+ifelse(type=="percent",0,intercept)))
  } else if (logLink=="inv10") {
    y<-(1/(10^(y)))
    yplus<-(1/(10^(yplus+ifelse(type=="percent",0,intercept))))
    yminus<-(1/(10^(yminus+ifelse(type=="percent",0,intercept))))
  } else if (logLink=="b"){
    y<-(1/(1+exp(-(intercept+y))))
    yplus<-(1/(1+exp(-(intercept+yplus))))
    yminus<-(1/(1+exp(-(intercept+yminus))))
  } else if (logLink=="n"){
    
  } else {
    stop("Error: the specified log link is not supported")
  }
  
  if((type=="percent") & (logLink!="n")){
    if(logLink=="b"){
      y<-((y/(1/(1+exp(-(intercept)))))*100)-100
      yplus<-((yplus/(1/(1+exp(-(intercept)))))*100)-100
      yminus<-((yminus/(1/(1+exp(-(intercept)))))*100)-100
    } else {
      y<-y*100-100
      yplus<-yplus*100-100
      yminus<-yminus*100-100
    }
    
    
  }
  
  predRange<-max(yplus,na.rm=T)-min(yminus,na.rm=T)
  if (all(is.na(ylim))){
    plotLims<-c(min(yminus,na.rm=T)-0.35*predRange,max(yplus,na.rm=T)+0.05*predRange)
  } else {
    plotLims<-ylim
    predRange<-ylim[2]-ylim[1]
  }
  
  par(mar=c(0.2,3.5,0.2,0.2))
  par(cex.lab=1)
  par(cex.axis=0.7)
  cex.pt<-0.5
  par(mgp=c(2,1,0))
  cex.txt<-ifelse(is.null(cex.txt),0.75,cex.txt)
  par(las=1)
  
  for (p in names(params)){
    eval(parse(text=gsub("x",p,"par(x=params$x)")))
  }
  
  if (responseVar != ""){
    if (logLink=="n"){
      ylabel=paste("Relative",responseVar)
    } else {
      if (type!="response"){
        ylabel=paste(responseVar,"difference (%)")
      } else {
        ylabel=paste(responseVar)
      }
      
    }
    
  } else {
    ylabel<-""
  }
  
  plot.cols<-"#000000"
  
  if (!is.null(errbar.cols)){
    plot.cols<-errbar.cols
  }
  
  if(!add){
    plot(-99999,-99999,xlim=c(0.5,length(labels)+0.5),ylim=plotLims,xlab=NA,ylab=
           ylabel,xaxt="n",bty="n",main=main,yaxp=yaxp)
  }
  
  for (i in 1:length(catEffects)){
    if ((i!=length(catEffects)) | (length(contEffects)>0)){
      abline(v=max(grep(catEffects[i],coef.labels))+0.5,lwd=1,
             col="dark grey")
    }
    
  }
  if (any(contEffects>0)){
    if (length(contEffects)>0){
      for (i in 1:length(contEffects)){
        if (i != length(contEffects)){
          if (is.null(contEffectsLabels)){
            abline(v=which(labels==names(contEffects)[i])+1.5,lwd=1,
                   col="dark grey")
          } else {
            abline(v=which(labels==contEffectsLabels[i])+1.5,lwd=1,
                   col="dark grey")
          }
          
        }
      }
      
    }
    
  }
  
  if (!add){
    if ((plotLabels) & (length(which(plot.cols=="black"))>0)){
      text(which(plot.cols=="black"),plotLims[1]+0.05*predRange,
           labels[which(plot.cols=="black")],cex=cex.txt,srt=xtext.srt,xpd=TRUE)
      
    }
    
    if(logLink=="n"){
      abline(h=0,col="dark grey")
    } else {
      abline(h=0,col="dark grey")
    }
    
  }
  
  if (length(contEffects>0)){
    if ((logLink=="10") | (logLink=="e")){
      text((which(labels=="")[1]:(which(labels=="")[1]+
                                    length(contEffects)*3))+offset,
           yminus[which(labels=="")[1]:(which(labels=="")[1]+
                                          length(contEffects)*3)]-(predRange/13),
           rep(c("L","M","H"),length(contEffects)),
           col="black",cex=cex.txt,srt=45)
    } else if (logLink=="inv10") {
      text((which(labels=="")[1]:(which(labels=="")[1]+
                                    length(contEffects)*3))+offset,
           yplus[which(labels=="")[1]:(which(labels=="")[1]+
                                         length(contEffects)*3)]-(predRange/13),
           rep(c("L","M","H"),length(contEffects)),
           col="black",cex=cex.txt,srt=45)
    } else {
      text((which(labels=="")[1]:(which(labels=="")[1]+
                                    length(contEffects)*3))+offset,
           yminus[which(labels=="")[1]:(which(labels=="")[1]+
                                          length(contEffects)*3)]-(predRange/13),
           rep(c("L","M","H"),length(contEffects)),
           col="black",cex=cex.txt,srt=45)
    }
    
  }
  
  if(is.null(pt.pch)){
    pt.pch<-16
  }
  
  if(is.null(yaxp)){
    errbar((1:length(labels))+offset,y,yplus,yminus,col=plot.cols,errbar.col=plot.cols,
           add=T,pch=pt.pch,cex=pt.cex,lty=errbar.lty)
  } else {
    errbar((1:length(labels))+offset,y,yplus,yminus,col=plot.cols,errbar.col=plot.cols,
           add=T,pch=pt.pch,cex=pt.cex,lty=errbar.lty)
  }
  
  if (pt.bg == "match"){
    points((1:length(labels))+offset,y,col=plot.cols,bg=plot.cols,pch=pt.pch,cex=pt.cex)
  } else {
    points((1:length(labels))+offset,y,col=plot.cols,bg=pt.bg,pch=pt.pch,cex=pt.cex)
  }
  
}