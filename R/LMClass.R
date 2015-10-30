# LM. An S3 class that represents a linear or generalized linear model

is.LM <- function(x) {
  return (inherits(x, 'LM'))
}

"[<-.LM" <- function(x, i, value) {
  if(!is.LM(x)) stop('Not an LM')
  stop("Can't assign to an LM")
}

"[[<-.LM" <- function(x, i, value) {
  if(!is.LM(x)) stop('Not an LM')
  stop("Can't assign to an LM")
}

'$<-.LM' <- function(object, x, value) {
  if(!is.LM(object)) stop('Not an LM')
  stop("Can't assign to an LM")
}

"names<-.LM" <- function(x, value) {
  if(!is.LM(x)) stop('Not an LM')
  stop("Can't assign to an LM")
}

"length<-.LM" <- function(x, value) {
  if(!is.LM(x)) stop('Not an LM')
  stop("Can't change length of an LM")
}

"levels<-.LM" <- function(x, value) {
  if(!is.LM(x)) stop('Not an LM')
  stop("Can't change levels of an LM")
}

"dim<-.LM" <- function(x, value) {
  if(!is.LM(x)) stop('Not an LM')
  stop("Can't change dim of an LM")
}

print.LM <- function(x, ...) {
  if(!is.LM(x)) stop('Not an LM')
  if(x$family=="gaussian"){
    cat('A linear model\n')
  } else {
    cat('A generalized linear model of family ',x$family,'\n',sep='')
  }
  cat('Final call: ',x$final.call,'\n',sep='')
  invisible(x)
}

summary.LM <- function(object, ...) {
  if(!is.LM(object)) stop('Not an LM')
  
  return(list(R.squared=summary(object$model)$r.squared,
              family=object$family,
              call=object$final.call,
              stats=object$stats))
  
}

plot.LM <- function(x, ...) {
  if(!is.LM(x)) stop('Not an LM')
  plot(x$model)
}

LM<-function(model,data,stats,final.call,family){
  stopifnot("lm"==class(model))
  stopifnot(is.data.frame(data))
  stopifnot(is.data.frame(stats))
  stopifnot(is.character(final.call))
  stopifnot(is.character(family))
  
  self<-list(model=model,data=data,stats=stats,final.call=final.call,family=family)
  
  class(self)<-'LM'
  
  return(self)
}
  
  


