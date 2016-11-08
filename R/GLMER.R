GLMER<-function(modelData,responseVar,fitFamily,fixedStruct,
                      randomStruct,saveVars=character(0),REML=TRUE,
                     optimizer="bobyqa",maxIters=10000){
  
  call.best<-.ConstructCall(responseVar,fixedStruct,randomStruct)
  
  allTerms<-unlist(strsplit(call.best,"[+]"))
  allTerms<-unlist(strsplit(allTerms,"[-]"))
  allTerms<-unlist(strsplit(allTerms,"[~]"))
  allTerms<-unlist(strsplit(allTerms,"[(]"))
  allTerms<-unlist(strsplit(allTerms,"[)]"))
  allTerms<-unlist(strsplit(allTerms,"[|]"))
  allTerms<-unlist(strsplit(allTerms,"[:]"))
  allTerms<-unlist(strsplit(allTerms,"[*]"))
  allTerms<-unique(allTerms)
  allTerms<-gsub(" ","",allTerms)
  allTerms<-allTerms[allTerms!=""]
  allTerms<-allTerms[allTerms!="1"]
  allTerms<-allTerms[allTerms!="0"]
  allTerms<-allTerms[!grepl("poly",allTerms)]
  
  allTerms<-gsub("[,][0-9]","",allTerms)
  
  allTerms<-c(allTerms,saveVars)
  
  allTerms<-unique(allTerms)
  
  modelData<-subset(modelData,select=c(allTerms))
  modelData<-na.omit(modelData)
  
  if (fitFamily=="gaussian"){
    eval(substitute(m<-lmer(cb,data=modelData,lmerControl(optimizer = optimizer,
                                                           optCtrl = list(maxfun=maxIters)),
                            REML=REML),
                    list(cb=call.best)))
  } else {
    eval(substitute(m<-glmer(cb,family=fitFamily,data=modelData,
                             control=glmerControl(optimizer = optimizer,
                                                  optCtrl = list(maxfun=maxIters))),
                    list(cb=call.best)))
  }
  
  
  
  return(list(model=m,data=modelData))
  
  
}