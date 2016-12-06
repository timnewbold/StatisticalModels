
SpatialAutocorrelationTest <- function(model){
  
  stopifnot(is.LM(model))
  
  stopifnot("Latitude" %in% names(model$data))
  stopifnot("Longitude" %in% names(model$data))
  
  ds.nb<-try(dnearneigh(cbind(model$data$Longitude,model$data$Latitude),
                        d1=0.00000001,d2=10),silent=TRUE)
  ds.listw<-try(nb2listw(ds.nb),silent=TRUE)
  mt<-tryCatch(moran.test(residuals(model$model),ds.listw),silent=TRUE,error=function(e) e, 
               warning=function(w) w)
  
  if(class(mt)[1]=="htest"){
    
    return(list(i=mt$statistic,p=mt$p.value))
    
  } else {
    cat("Error: spatial autocorrelation test failed\n")
    return(NULL)
  }
  
}