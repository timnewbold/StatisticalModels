GLMERCandidates <- function(modelData,responseVar,fitFamily,candidateFixedStructs,
                            randomStruct,saveVars=character(0),REML=TRUE,
                            optimizer="bobyqa",maxIters=10000,ADMB=FALSE){
  
  # Ensure that models will be fit to same dataset, by removing NAs
  modelData <- na.omit(modelData)
  
  # Find the number of cores on this computer
  cores <- parallel::detectCores()
  
  # Make a cluster with one fewer cores than the computer possesses
  cl <- snow::makeCluster(cores-1)
  
  # Export necessary data to the cluster
  snow::clusterExport(cl = cl,list = c("modelData","responseVar","fitFamily",
                                       "randomStruct","saveVars","REML",
                                       "optimizer","maxIters","ADMB"))
  
  # Load the lme4 package on each node in the cluster
  invisible(snow::clusterCall(cl = cl,fun = function() library(lme4)))
  
  # Iterate over candidate fixed effects and build a model
  models <- snow::parLapply(cl = cl,x = candidateFixedStructs,fun = function(cand){
    
    mod <- StatisticalModels::GLMER(modelData = modelData,
                                    responseVar = responseVar,
                                    fitFamily = fitFamily,
                                    fixedStruct = cand,
                                    randomStruct = randomStruct,
                                    saveVars = saveVars,
                                    REML = REML,
                                    optimizer = optimizer,
                                    maxIters = maxIters,
                                    ADMB = ADMB)
    
    return(mod$model)
    
  })
  
  # Stop the cluster
  snow::stopCluster(cl)
  
  return(models)
  
}