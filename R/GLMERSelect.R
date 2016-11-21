GLMERSelect<-function(modelData,responseVar,fitFamily,fixedFactors=
                         character(0),fixedTerms=list(),randomStruct,
                       fixedInteractions=character(0),fitInteractions=FALSE,
                      verbose=FALSE,saveVars=character(0),
                      optimizer="bobyqa",maxIters=10000){
  
  contEffectNames<-names(fixedTerms)
  
  randomTerms <- strsplit(randomStruct,'[+]')[[1]]
  randomTerms <- unlist(sapply(randomTerms,function(x) strsplit(x,'[|]')))
  randomTerms <- randomTerms[(randomTerms!="(1")]
  randomTerms <- gsub("[)]","",randomTerms)
  
  allTerms <- unique(c(randomTerms,fixedFactors,names(fixedTerms),responseVar,saveVars))
  
  modelData<-subset(modelData,select=c(allTerms))
  modelData<-na.omit(modelData)
  
  cat<-sapply(modelData,is.factor)
  modelData[cat]<-lapply(modelData[cat],factor)
  
  for (fe in fixedFactors){
    eval(substitute(modelData$x<-factor(modelData$x),list(x=fe)))
  }
  
  results<-list(fixef=character(),AIC=numeric())
  
  allTerms<-character(0)
  fixedStruct<-""
  for (i in 1:length(fixedFactors)){
    fixedStruct<-paste(fixedStruct,fixedFactors[i],sep="")
    allTerms<-c(allTerms,fixedFactors[i])
    if ((i != length(fixedFactors)) | (length(fixedTerms)>0) | 
          ((length(fixedTerms)==0) & (
            length(fixedInteractions)>0))){
      fixedStruct<-paste(fixedStruct,"+",sep="")
    }
  }
  if (length(fixedTerms)>0){
    for (i in 1:length(fixedTerms)){
      term<-paste("poly(",names(fixedTerms)[i],
                  ",",fixedTerms[i],")",sep="")
      fixedStruct<-paste(fixedStruct,term,sep="")
      allTerms<-c(allTerms,term)
      if ((i != length(fixedTerms)) | (length(fixedInteractions)>0)){
        fixedStruct<-paste(fixedStruct,"+",sep="")
      }
    }
  }
  
  if (fitInteractions){
    fixedStruct<-paste(fixedStruct,"+")
    
    mainTerms<-allTerms
    
    for (i in 1:(length(mainTerms)-1)){
      for (j in (i+1):length(mainTerms)){
        term<-paste(mainTerms[i],mainTerms[j],sep=":")
        fixedStruct<-paste(fixedStruct,term)
        allTerms<-c(allTerms,term)
      }
    }
    
  }
  
  if (length(fixedInteractions)>0){
    for (i in 1:length(fixedInteractions)){
      fixedStruct<-paste(fixedStruct,fixedInteractions[i],sep="")
      allTerms<-c(allTerms,fixedInteractions[i])
      if (i != length(fixedInteractions)){
        fixedStruct<-paste(fixedStruct,"+",sep="")
      }
    }
  }
  
  randomStruct<-gsub(" ","",randomStruct)
  
  call.old<-paste(responseVar,"~",allTerms[1],sep="")
  if (length(allTerms)>1){
    for (t in 2:length(allTerms)){
      call.old<-paste(call.old,"+",allTerms[t],sep="")
    }
  }
  call.old<-paste(call.old,"+",randomStruct,sep="")
  
  stats<-data.frame(allTerms)
  names(stats)<-"terms"
  polyTerms2 <- paste(stats$terms[grep("poly[(].{1,},2[)]", 
                                       stats$terms)])
  polyInters <- which(grepl(":", polyTerms2))
  if (length(polyInters) > 0) 
    polyTerms2 <- polyTerms2[-polyInters]
  polyTerms2 <- gsub(",2", ",1", polyTerms2)
  stats <- rbind(stats, data.frame(terms = polyTerms2))
  polyTerms3 <- paste(stats$terms[grep("poly[(].{1,},3[)]", 
                                       stats$terms)])
  polyInters <- which(grepl(":", polyTerms3))
  if (length(polyInters) > 0) 
    polyTerms3 <- polyTerms3[-polyInters]
  polyTerms3 <- gsub(",3", ",2", polyTerms3)
  stats <- rbind(stats, data.frame(terms = polyTerms3))
  polyTerms3 <- gsub(",2", ",1", polyTerms3)
  stats <- rbind(stats, data.frame(terms = polyTerms3))
  stats$terms<-paste(stats$terms)
  stats$ChiSq<-NA
  stats$Df<-NA
  stats$P<-NA
  stats$dAIC<-NA
  
  iter<-1
  
  
  repeat {
    
    print(paste("Performing round ",iter," of interaction-term removal",sep=""))
    
    if (verbose) print(call.old)
     
    if (fitFamily=="gaussian"){
      mOld<-lmer(call.old,data=modelData,REML=FALSE,
                 lmerControl(optimizer = optimizer,optCtrl=list(maxfun=maxIters)))
    } else {
      mOld<-glmer(call.old,family=fitFamily,data=modelData,
                  control=glmerControl(optimizer = optimizer,optCtrl=list(maxfun=maxIters)))
    }
    
    if (iter==1) iTerms<-allTerms[grep(":",allTerms)]
    
    iTerms<-gsub(" ","",iTerms)
    
    pVals<-numeric()
    Chis<-numeric()
    Dfs<-character()
    dAICs<-numeric()
    
    for (t in iTerms){
      
      t1<-gsub("[(]","[(]",t)
      t2<-gsub("[)]","[)]",t1)
      t3<-paste(t2,"[+]",sep="")
      
      call.new<-gsub(t3,"",call.old)
      
      if (fitFamily=="gaussian"){
        mNew<-lmer(call.new,data=modelData,REML=FALSE,
                   lmerControl(optimizer = optimizer,optCtrl=list(maxfun=maxIters)))
      } else {
        mNew<-glmer(call.new,family=fitFamily,data=modelData,
                    control=glmerControl(optimizer = optimizer,optCtrl=list(maxfun=maxIters)))
      }
      
      
      pVals<-c(pVals,anova(mOld,mNew)$Pr[2])
      Chis<-c(Chis,anova(mOld,mNew)$Chisq[2])
      Dfs<-c(Dfs,paste(anova(mOld,mNew)$'Chi Df'[2],",",anova(mOld,mNew)$Df[2]))
      dAICs<-c(dAICs,AIC(mNew)-AIC(mOld))
    }
    
    if (verbose){
      print(iTerms)
      print(pVals)
    }
    
    print(paste(length(which(pVals>0.05))," interaction terms have P-values >0.05",sep=""))
    
    if (length(which(pVals>0.05))==0) break
    
    dropI<-iTerms[order(pVals)[length(order(pVals))]]
    
    stats$ChiSq[which(stats$terms==dropI)]<-Chis[order(pVals)[length(order(pVals))]]
    stats$Df[which(stats$terms==dropI)]<-Dfs[order(pVals)[length(order(pVals))]]
    stats$P[which(stats$terms==dropI)]<-pVals[order(pVals)[length(order(pVals))]]
    stats$dAIC[which(stats$terms==dropI)]<-dAICs[order(pVals)[length(order(pVals))]]
    
    print(paste("Dropping ",dropI,sep=""))
    
    t1<-gsub("[(]","[(]",dropI)
    t2<-gsub("[)]","[)]",t1)
    t3<-paste(t2,"[+]",sep="")
    
    call.old<-gsub(t3,"",call.old)
    
    iTerms<-iTerms[-order(pVals)[length(order(pVals))]]
    allTerms<-allTerms[-which(allTerms==dropI)]
    
    iter<-iter+1 
    
  }
  
  stats$ChiSq[na.omit(match(iTerms,stats$terms))]<-Chis
  stats$Df[na.omit(match(iTerms,stats$terms))]<-Dfs
  stats$P[na.omit(match(iTerms,stats$terms))]<-pVals
  stats$dAIC[na.omit(match(iTerms,stats$terms))]<-dAICs
  
  
  for (t in iTerms){
    t1<-gsub("[(]","[(]",t)
    t2<-gsub("[)]","[)]",t1)
    t3<-paste(t2,"[+]",sep="")
    
    call.old<-gsub(t3,"",call.old)
  }
  
  itersRemaining<-which(allTerms %in% iTerms)
  if (length(itersRemaining)>0) allTerms<-allTerms[-itersRemaining]
  
  iter<-1
  
  repeat {
    
    print(paste("Performing round ",iter," of main-effect removal",sep=""))
    
    if (verbose) print(call.old)
    
    if (fitFamily=="gaussian"){
      mOld<-lmer(call.old,data=modelData,REML=FALSE,
                 lmerControl(optimizer = optimizer,optCtrl=list(maxfun=maxIters)))
    } else {
      mOld<-glmer(call.old,family=fitFamily,data=modelData,
                  control=glmerControl(optimizer = optimizer,optCtrl=list(maxfun=maxIters)))
    }
    
    
    
    mTerms<-allTerms
    
    mTerms<-gsub(" ","",mTerms)
    
    pVals<-numeric()
    Chis<-numeric()
    Dfs<-character()
    dAICs<-numeric()
    
    for (t in mTerms){
      if ((grepl("poly",t)) & grepl(",2",t)){
        t1<-gsub("[(]","[(]",t)
        t2<-gsub("[)]","[)]",t1)
        t3<-gsub(",2",",1",t)
        
        call.new<-gsub(t2,t3,call.old)
        
      } else {
        t1<-gsub("[(]","[(]",t)
        t2<-gsub("[)]","[)]",t1)
        t3<-paste(t2,"[+]",sep="")
        
        call.new<-sub(t3,"",call.old)
      }
      
      if (fitFamily=="gaussian"){
        mNew<-lmer(call.new,data=modelData,REML=FALSE,
                   lmerControl(optimizer = optimizer,optCtrl=list(maxfun=maxIters)))
      } else {
        mNew<-glmer(call.new,family=fitFamily,data=modelData,
                    control=glmerControl(optimizer = optimizer,optCtrl=list(maxfun=maxIters)))
      }
      
      
      pVals<-c(pVals,anova(mOld,mNew)$Pr[2])
      Chis<-c(Chis,anova(mOld,mNew)$Chisq[2])
      Dfs<-c(Dfs,paste(anova(mOld,mNew)$'Chi Df'[2],",",anova(mOld,mNew)$Df[2]))
      dAICs<-c(dAICs,AIC(mNew)-AIC(mOld))
    }
    
    print(paste(length(which(pVals>0.05))," candidate main effects have P-values >0.05",sep=""))
    
    if (verbose){
      print(mTerms)
      print(pVals)
    }
    
    if (length(which(pVals>0.05))==0) break
    
    dropM<-mTerms[order(pVals)[length(order(pVals))]]
    
    stats$ChiSq[which(stats$terms==dropM)]<-Chis[order(pVals)[length(order(pVals))]]
    stats$Df[which(stats$terms==dropM)]<-Dfs[order(pVals)[length(order(pVals))]]
    stats$P[which(stats$terms==dropM)]<-pVals[order(pVals)[length(order(pVals))]]
    stats$dAIC[which(stats$terms==dropM)]<-dAICs[order(pVals)[length(order(pVals))]]
    
    if ((grepl("poly", dropM)) & grepl(",3", dropM)) {
      print(paste("Simplifying ", dropM, sep = ""))
      d1 <- gsub("[(]", "[(]", dropM)
      d2 <- gsub("[)]", "[)]", d1)
      d3 <- gsub(",3", ",2", dropM)
      call.old <- gsub(d2, d3, call.old)
      mTerms <- gsub(d2, d3, mTerms)
      allTerms <- gsub(d2, d3, allTerms)
    } else if ((grepl("poly",dropM)) & grepl(",2",dropM)){
      print(paste("Simplifying ",dropM,sep=""))
      
      d1<-gsub("[(]","[(]",dropM)
      d2<-gsub("[)]","[)]",d1)
      d3<-gsub(",2",",1",dropM)
      
      call.old<-gsub(d2,d3,call.old)
      
      mTerms<-gsub(d2,d3,mTerms)
      allTerms<-gsub(d2,d3,allTerms)
      
    } else {
      print(paste("Dropping ",dropM,sep=""))
      t1<-gsub("[(]","[(]",dropM)
      t2<-gsub("[)]","[)]",t1)
      t3<-paste(t2,"[+]",sep="")
      
      
      call.old<-gsub(t3,"",call.old)
      
      mTerms<-mTerms[-order(pVals)[length(order(pVals))]]
      allTerms<-allTerms[-which(allTerms==dropM)]
      
    }
    
    
    
    iter<-iter+1 
    
    
  }
  
  stats$ChiSq[na.omit(match(mTerms,stats$terms))]<-Chis
  stats$Df[na.omit(match(mTerms,stats$terms))]<-Dfs
  stats$P[na.omit(match(mTerms,stats$terms))]<-pVals
  stats$dAIC[na.omit(match(mTerms,stats$terms))]<-dAICs
  
  fixedStruct<-""
  sig.terms<-stats[stats$P<0.05,]
  sig.terms<-na.omit(sig.terms)
  if (dim(sig.terms)[1]>0){
    sig.terms<-paste(sig.terms$terms)
    sig.inter<-sig.terms[grepl(":",sig.terms)]
    inter.mains<-unique(unlist(strsplit(sig.inter,":")))
    sig.terms<-c(sig.terms,inter.mains[!(inter.mains %in% sig.terms)])
    for (t in names(fixedTerms)){
      mainMatches<-sig.terms[which((grepl(paste("poly[(]",t,",[0-9]{1}[)]",sep=""),sig.terms)) & 
                                     !(grepl(":",sig.terms)))]
      mainMatchPosits<-which((grepl(paste("poly[(]",t,",[0-9]{1}[)]",sep=""),sig.terms)) & 
                               !(grepl(":",sig.terms)))
      if (length(mainMatches)>1){
        sig.terms<-sig.terms[-mainMatchPosits[order(mainMatches,decreasing=TRUE)][-1]]
      }
    }
    for (i in 1:length(sig.terms)){
      fixedStruct<-paste(fixedStruct,sig.terms[i],sep="")
      if (i != length(sig.terms)) fixedStruct<-paste(fixedStruct,"+",sep="")
    }
    call.best<-.ConstructCall(responseVar,fixedStruct,randomStruct)
    if (verbose) print(call.best)
    if (fitFamily=="gaussian"){
      mBest<-lmer(call.best,data=modelData,REML=TRUE,
                  lmerControl(optimizer = optimizer,optCtrl=list(maxfun=maxIters)))
    } else {
      mBest<-glmer(call.best,family=fitFamily,data=modelData,
                   control=glmerControl(optimizer = optimizer,optCtrl=list(maxfun=maxIters)))
    }
    return(list(model=mBest,data=modelData,stats=stats,final.call=call.best))
  } else {
    print("Warning: all terms were dropped from the model")
    return(list(model=NULL,data=modelData,stats=stats,final.call=NULL))
  }
  
}