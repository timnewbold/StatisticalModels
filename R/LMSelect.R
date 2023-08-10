
LMSelect <- function(modelData,responseVar,fitFamily,factors=
                       character(0),contEffects=list(),
                     interactions=character(0),
                     allInteractions=FALSE,
                     alpha=0.05,
                     saveVars=character(0)){
  
  contEffectNames<-names(contEffects)
  
  if ((length(interactions)>0) & (allInteractions)){
    stop("Error: specifying particular interactions and all two-way interactions will not work!")
  }
  
  model.data<-subset(modelData,select=c(factors,names(contEffects),
                                        responseVar,saveVars))
  
  
  model.data<-na.omit(model.data)
  cat<-sapply(model.data,is.factor)
  model.data[cat]<-lapply(model.data[cat],factor)
  
  
  for (fe in factors){
    eval(substitute(model.data$x<-factor(model.data$x),list(x=fe)))
  }
  
  allTerms<-character(0)
  allTerms<-c(allTerms,factors)
  allTerms<-c(allTerms,paste("poly(",names(contEffects),",",contEffects,")",sep=""))
  allTerms<-c(allTerms,interactions)
  if (allInteractions){
    mainTerms<-allTerms
    
    for (i in 1:(length(mainTerms)-1)){
      for (j in (i+1):length(mainTerms)){
        term<-paste(mainTerms[i],mainTerms[j],sep=":")
        allTerms<-c(allTerms,term)
      }
    }
  }
  
  call.old<-paste(responseVar,"~",paste(allTerms,collapse="+"),sep="")
  
  stats<-data.frame(terms=allTerms)
  
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
  if (fitFamily=="gaussian"){
    stats$F<-NA
  } else {
    stats$Deviance<-NA
  }
  stats$Df<-NA
  stats$P<-NA
  if (fitFamily != "quasipoisson"){
    stats$dAIC<-NA
  }
  
  iter<-1
  
  repeat {
    
    .Log(paste("Performing round ",iter," of interaction-term removal\n",sep=""))
    
    .Log(paste(call.old,"\n",sep=""))
    
    if (fitFamily=="gaussian"){
      mOld<-lm(call.old,data=model.data)
    } else {
      mOld<-glm(call.old,family=fitFamily,data=model.data)
    }
    
    if (iter==1) iTerms<-allTerms[grep(":",allTerms)]
    
    iTerms<-gsub(" ","",iTerms)
    
    pVals<-numeric()
    if (fitFamily=="gaussian"){
      Fs<-numeric()
    } else {
      Devs <- numeric()
    }
    Dfs<-character()
    if (fitFamily != "quasipoisson"){
      dAICs<-numeric()
    }
    
    for (t in iTerms){
      t1<-gsub("[(]","[(]",t)
      t2<-gsub("[)]","[)]",t1)
      if (t != tail(allTerms,1)){
        t3<-paste(t2,"[+]",sep="")
      } else {
        t3<-paste("[+]",t2,sep="")
      }
      
      call.new<-gsub(t3,"",call.old)
      
      if (fitFamily=="gaussian"){
        mNew<-lm(call.new,data=model.data)
      } else {
        mNew<-glm(call.new,family=fitFamily,data=model.data)
      }
      
      if (fitFamily=="gaussian"){
        pVals<-c(pVals,anova(mOld,mNew)$Pr[2])
        Fs<-c(Fs,anova(mOld,mNew)$F[2])
        Dfs<-c(Dfs,paste(
          anova(mOld,mNew)$Df[2],
          -anova(mOld,mNew)$Df[2],sep=","))
        dAICs<-c(dAICs,AIC(mNew)-AIC(mOld))
      } else {
        pVals<-c(pVals,anova(mOld,mNew,test = "Chi")$Pr[2])
        Devs<-c(Devs,-anova(mOld,mNew,test = "Chi")$Deviance[2])
        Dfs<-c(Dfs,paste(
          -anova(mOld,mNew,test = "Chi")$Df[2],
          anova(mOld,mNew,test = "Chi")$Res.Df[2],
          sep=","))
        if (fitFamily != "quasipoisson"){
          dAICs<-c(dAICs,AIC(mNew)-AIC(mOld))
        }
      }
      
      
    }
    
    .Log(paste(length(which(pVals>alpha))," interaction terms have P-values >alpha\n",sep=""))
    
    if (length(which(pVals>alpha))==0) break
    
    dropI<-iTerms[order(pVals)[length(order(pVals))]]
    if (fitFamily=="gaussian"){
      stats$F[which(stats$terms==dropI)]<-Fs[order(pVals)[length(order(pVals))]]
    } else {
      stats$Deviance[which(stats$terms==dropI)]<-Devs[order(pVals)[length(order(pVals))]]
    }
    stats$Df[which(stats$terms==dropI)]<- Dfs[order(pVals)[length(order(pVals))]]
    stats$P[which(stats$terms==dropI)]<-pVals[order(pVals)[length(order(pVals))]]
    if (fitFamily != "quasipoisson"){
      stats$dAIC[which(stats$terms==dropI)]<-dAICs[order(pVals)[length(order(pVals))]]
    }
    
    .Log(paste("Dropping ",dropI,"\n",sep=""))
    
    t1<-gsub("[(]","[(]",dropI)
    t2<-gsub("[)]","[)]",t1)
    if (t != tail(allTerms,1)){
      t3<-paste(t2,"[+]",sep="")
    } else {
      t3<-paste("[+]",t2,sep="")
    }
    
    call.old<-gsub(t3,"",call.old)
    
    iTerms<-iTerms[-order(pVals)[length(order(pVals))]]
    allTerms<-allTerms[-which(allTerms==dropI)]
    
    iter<-iter+1 
    
  }
  
  if (fitFamily=="gaussian"){
    stats$F[na.omit(match(iTerms,stats$terms))]<-Fs
  } else {
    stats$Deviance[na.omit(match(iTerms,stats$terms))]<-Devs
  }
  stats$Df[na.omit(match(iTerms,stats$terms))]<- Dfs
  stats$P[na.omit(match(iTerms,stats$terms))]<-pVals
  if (fitFamily != "quasipoisson"){
    stats$dAIC[na.omit(match(iTerms,stats$terms))]<-dAICs
  }
  
  for (t in iTerms){
    t1<-gsub("[(]","[(]",t)
    t2<-gsub("[)]","[)]",t1)
    if (t != tail(allTerms,1)){
      t3<-paste(t2,"[+]",sep="")
    } else {
      t3<-paste("[+]",t2,sep="")
    }
    
    call.old<-gsub(t3,"",call.old)
  }
  
  itersRemaining<-which(allTerms %in% iTerms)
  if (length(itersRemaining)>0) allTerms<-allTerms[-itersRemaining]
  
  iter<-1
  
  repeat {
    
    .Log(paste("Performing round ",iter," of main-effect removal\n",sep=""))
    
    .Log(paste(call.old,"\n",sep=""))
    
    if (fitFamily=="gaussian"){
      mOld<-lm(call.old,data=model.data)
    } else {
      mOld<-glm(call.old,family=fitFamily,data=model.data)
    }
    
    mTerms<-allTerms
    
    mTerms<-gsub(" ","",mTerms)
    
    pVals<-numeric()
    if(fitFamily=="gaussian"){
      Fs<-numeric()
    } else {
      Devs <- numeric()
    }
    Dfs<-character()
    if (fitFamily != "quasipoisson"){
      dAICs<-numeric()
    }
    
    for (t in mTerms){
      if ((grepl("poly",t)) & grepl(",3",t)){
        t1<-gsub("[(]","[(]",t)
        t2<-gsub("[)]","[)]",t1)
        t3<-gsub(",3",",2",t)
        
        call.new<-gsub(t2,t3,call.old)
        
      }else if ((grepl("poly",t)) & grepl(",2",t)){
        t1<-gsub("[(]","[(]",t)
        t2<-gsub("[)]","[)]",t1)
        t3<-gsub(",2",",1",t)
        
        call.new<-gsub(t2,t3,call.old)
        
      } else {
        t1<-gsub("[(]","[(]",t)
        t2<-gsub("[)]","[)]",t1)
        if (t != tail(allTerms,1)){
          t3<-paste(t2,"[+]",sep="")
        } else {
          t3<-paste("[+]",t2,sep="")
        }
        call.new<-sub(t3,"",call.old)
      }
      
      if (fitFamily=="gaussian"){
        mNew<-lm(call.new,data=model.data)
      } else {
        mNew<-glm(call.new,family=fitFamily,data=model.data)
      }
      
      if (fitFamily=="gaussian"){
        pVals<-c(pVals,anova(mOld,mNew)$Pr[2])
        Fs<-c(Fs,anova(mOld,mNew)$F[2])
        Dfs<-c(Dfs,paste(
          -anova(mOld,mNew)$Df[2],
          anova(mOld,mNew)$Res.Df[2],sep=","))
        dAICs<-c(dAICs,AIC(mNew)-AIC(mOld))
      } else {
        pVals<-c(pVals,anova(mOld,mNew,test = "Chi")$Pr[2])
        Devs<-c(Devs,-anova(mOld,mNew,test = "Chi")$Deviance[2])
        Dfs<-c(Dfs,paste(
          -anova(mOld,mNew,test = "Chi")$Df[2],
          anova(mOld,mNew,test = "Chi")$Res.Df[2],
          sep=","))
        if (fitFamily != "quasipoisson"){
          dAICs<-c(dAICs,AIC(mNew)-AIC(mOld))
        }
      }
    }
    
    .Log(paste(length(which(pVals>alpha))," candidate main effects have P-values > alpha\n",sep=""))
    
    if (length(which(pVals>alpha))==0) break
    
    dropM<-mTerms[order(pVals)[length(order(pVals))]]
    
    if (fitFamily=="gaussian"){
      stats$F[which(stats$terms==dropM)]<-Fs[order(pVals)[length(order(pVals))]]
    } else {
      stats$Deviance[which(stats$terms==dropM)]<-Devs[order(pVals)[length(order(pVals))]]
    }
    stats$Df[which(stats$terms==dropM)]<- Dfs[order(pVals)[length(order(pVals))]]
    stats$P[which(stats$terms==dropM)]<-pVals[order(pVals)[length(order(pVals))]]
    if (fitFamily != "quasipoisson"){
      stats$dAIC[which(stats$terms==dropM)]<-dAICs[order(pVals)[length(order(pVals))]]
    }
    
    if ((grepl("poly", dropM)) & grepl(",3", dropM)) {
      .Log(paste("Simplifying ", dropM, "\n", sep = ""))
      d1 <- gsub("[(]", "[(]", dropM)
      d2 <- gsub("[)]", "[)]", d1)
      d3 <- gsub(",3", ",2", dropM)
      call.old <- gsub(d2, d3, call.old)
      mTerms <- gsub(d2, d3, mTerms)
      allTerms <- gsub(d2, d3, allTerms)
    } else if ((grepl("poly",dropM)) & grepl(",2",dropM)){
      .Log(paste("Simplifying ",dropM,"\n",sep=""))
      
      d1<-gsub("[(]","[(]",dropM)
      d2<-gsub("[)]","[)]",d1)
      d3<-gsub(",2",",1",dropM)
      
      call.old<-gsub(d2,d3,call.old)
      
      mTerms<-gsub(d2,d3,mTerms)
      allTerms<-gsub(d2,d3,allTerms)
      
    } else {
      .Log(paste("Dropping ",dropM,"\n",sep=""))
      t1<-gsub("[(]","[(]",dropM)
      t2<-gsub("[)]","[)]",t1)
      if (t != tail(allTerms,1)){
        t3<-paste(t2,"[+]",sep="")
      } else {
        t3<-paste("[+]",t2,sep="")
      }
      
      call.old<-gsub(t3,"",call.old)
      
      mTerms<-mTerms[-order(pVals)[length(order(pVals))]]
      allTerms<-allTerms[-which(allTerms==dropM)]
      
    }
    
    
    
    iter<-iter+1 
    
    
  }
  
  if (fitFamily=="gaussian"){
    stats$F[na.omit(match(mTerms,stats$terms))]<-Fs
  } else {
    stats$Deviance[na.omit(match(mTerms,stats$terms))]<-Devs
  }
  stats$Df[na.omit(match(mTerms,stats$terms))]<- Dfs
  stats$P[na.omit(match(mTerms,stats$terms))]<-pVals
  if (fitFamily != "quasipoisson"){
    stats$dAIC[na.omit(match(mTerms,stats$terms))]<-dAICs
  }
  
  sig.terms<-stats[stats$P<alpha,]
  sig.terms<-na.omit(sig.terms)
  if (dim(sig.terms)[1]>0){
    sig.terms<-paste(sig.terms$terms)
    sig.inter<-sig.terms[grepl(":",sig.terms)]
    inter.mains<-unique(unlist(strsplit(sig.inter,":")))
    sig.terms<-c(sig.terms,inter.mains[!(inter.mains %in% sig.terms)])
    for (t in names(contEffects)){
      mainMatches<-sig.terms[which((grepl(paste("poly[(]",t,",[0-9]{1}[)]",sep=""),sig.terms)) & 
                                     !(grepl(":",sig.terms)))]
      mainMatchPosits<-which((grepl(paste("poly[(]",t,",[0-9]{1}[)]",sep=""),sig.terms)) & 
                               !(grepl(":",sig.terms)))
      if (length(mainMatches)>1){
        sig.terms<-sig.terms[-mainMatchPosits[order(mainMatches,decreasing=TRUE)][-1]]
      }
    }
    
    call.best<-paste(responseVar,"~",paste(sig.terms,collapse="+"),sep="")
    
    if (fitFamily=="gaussian"){
      mBest<-lm(call.best,data=model.data)
    } else {
      mBest<-glm(call.best,family=fitFamily,data=model.data)
    }
    return(LM(model=mBest,data=model.data,stats=stats,final.call=call.best,
              family=fitFamily))
  } else {
    .Log("Warning: all terms dropped from the model")
    return(LM(model=lm(0~0),data=model.data,stats=stats,final.call="",
              family=fitFamily))
  }
  
  
  
}
