
SpatialAutocorrelationTest <- function(model,ranefGrouping=NULL){
  
  if(is.null(ranefGrouping)){
    stopifnot(is.LM(model))
    
    stopifnot("Latitude" %in% names(model$data))
    stopifnot("Longitude" %in% names(model$data))
    
    ds.nb<-try(dnearneigh(cbind(model$data$Longitude,model$data$Latitude),
                          d1=0.00000001,d2=10),silent=TRUE)
    ds.listw<-try(nb2listw(ds.nb,style = "W",zero.policy = TRUE),silent=TRUE)
    mt<-tryCatch(moran.test(residuals(model$model),ds.listw),silent=TRUE,error=function(e) e, 
                 warning=function(w) w)
    
    if(class(mt)[1]=="htest"){
      
      return(list(i=mt$statistic,p=mt$p.value))
      
    } else {
      cat("Error: spatial autocorrelation test failed\n")
      return(NULL)
    }
  } else {
    cat('CAUTION: function will only behave sensibly if data points each represent a unique spatial location\n')
    
    groups<-character()
    failed<-character()
    moran.i<-numeric()
    moran.p<-numeric()
    
    res <- residuals(model$model)
    
    i=1
    for (grp in unique(model$data[,ranefGrouping])){
      cat(paste("\rProcessing group ",i," of ",
                length(unique(model$data[,ranefGrouping])),sep=""))
      data.sub<-droplevels(model$data[model$data[,ranefGrouping]==grp,])
      
      ds.nb<-try(dnearneigh(cbind(data.sub$Longitude,data.sub$Latitude),
                            d1=0.00000001,d2=10),silent=TRUE)
      ds.listw<-try(nb2listw(ds.nb),silent=TRUE)
      mt<-tryCatch(moran.test(res,ds.listw),silent=TRUE,error=function(e) e, 
                   warning=function(w) w)
      
      if(class(mt)[1]=="htest"){
        if ((!is.na(mt$statistic))){
          groups<-c(groups,grp)
          moran.i<-c(moran.i,mt$statistic)
          moran.p<-c(moran.p,mt$p.value)
        } else {
          failed<-c(failed,ss)
        }
        
      } else {
        failed<-c(failed,ss)
      }
      
      i<-i+1
    }
    
    return(list(studies=studies,I=moran.i,P=moran.p,failed=failed))
    
  }
  
  
  
}