\name{PlotLMERContinuous}
\alias{PlotLMERContinuous}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plot modelled effects of continuous variables with uncertainty}
\description{
Plots modelled effects on a specified response variable of one or two specified continuous variables. Where two continuous effects are specified, a perspective plot is used.
}
\usage{
PlotLMERContinuous(model,data,effects,otherContEffects=character(0),
                          otherFactors=character(0),xlab,ylab,byFactor=NULL,
                          zlab=NULL,outDir=NULL,logLink="n",plotRug=FALSE,
                          seMultiplier=1.96,params=list(),
                          add=FALSE,ylim=NULL,zlim=NULL,line.cols=NULL,
                          line.types=NULL,plotUncertainty=TRUE,nPanels=1,
                          main=NULL,yDiv=1,transformX=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
The model to derive parameter estimates from
}
  \item{data}{
A data frame containing all variables in the model
}
  \item{effects}{
The continuous variables (up to a maximum of 2) to plot
}
  \item{otherContEffects}{
Any continuous variables other than those selected for display that were fitted in the model. Median values of these variables will be used as the reference level.
}
  \item{otherFactors}{
Any factors that were fitted in the model, given as a list where the names are the names of the factors and the values are the corresponding factor levels; shouldn't include any factors specified in byFactor
}
  \item{xlab}{
The x-axis label
}
  \item{ylab}{
The y-axis label
}
  \item{byFactor}{
Plot separate lines for each level in this factor
}
  \item{zlab}{
The z-axis label, for perspective plots showing two continuous effects
}
  \item{outDir}{
Optional: directory to write the graph to. If not specified, then the graph will be written to the native graphics device
}
  \item{logLink}{
The link function used by the model: one of "n" (identity), "e" (log, base e), "10" (log, base 10) or "b" (logit - used in binomial models)
}
  \item{plotRug}{
  Whether to add a "rug" showing the values of the continuous variable represented in the data
}
  \item{seMultiplier}{
The number of standard errors to use for the error bars. Default is 1.96, which shows 95 percent confidence intervals
}
  \item{params}{
Any R graphical parameters, which will overwrite the defaults
}
  \item{add}{
Whether to add the plot to an existing plot; default is FALSE
}
  \item{ylim}{
Fixed y-axis limits. If not specified, then resolved automatically
}
  \item{zlim}{
For plots with two continuous effects (perspective plots), optional limits for the z-axis  
}
  \item{line.types}{
The line types to use (solid, dashed, dotted etc.). Uses the same codes as R's graphical parameter 'lty'  
}
  \item{line.cols}{
The colours to use for the lines
}
  \item{plotUncertainty}{
Whether to include uncertainty estimates in the plot  
}
  \item{nPanels}{
  When byFactor is specified, the number of panels to divide the factor levels across
}
  \item{main}{
A title for the plot  
}
  \item{yDiv}{
Factor to divide y-axis values by (for example to remove scientific notation)  
}
  \item{transformX}{
Whether to back-transform log-transformed x-axis values to identity space. Currently only works with base-e 
}
}
\details{
Code for calculating predicted values and confidence intervals was taken from the GLMM wiki (see references).
}
\references{
http://glmm.wikidot.com/faq
}
\author{
Tim Newbold <Tim.Newbold@unep-wcmc.org>
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
# Load example PREDICTS site-level data for the UK
data(sites)
sites$UseIntensity<-sites$Use_intensity
sites$UI<-factor(paste(sites$LandUse,sites$Use_intensity))
sites$UI<-relevel(sites$UI,ref="Natural Minimal use")

# Load human population density map-header and data for the UK
data(hpd)

# Convert to a SpatialGridDataFrame object
hpd<-SpatialGridDataFrame(grid=GridTopology(
  cellcentre.offset=c(hpdHeader[3,2]+(hpdHeader[5,2]/2),
                      hpdHeader[4,2]++(hpdHeader[5,2]/2)),
  cellsize=c(hpdHeader[5,2],hpdHeader[5,2]),
  cells.dim=c(hpdHeader[1,2],hpdHeader[2,2])),
  data=data.frame(band1=c(t(hpdData))))

# Add human population density estimates to the sites data frame
sites<-AddGridData(gridData=hpd,dataFrame=sites,
                   columnName="HPD")
    
# Log estimates of human population density
sites$LogHPD<-log(sites$HPD+1)

# Example plot of a model of species richness against human population density
m1<-FixedModel(all.data=sites,responseVar="Species_richness",
fitFamily="poisson",fixedStruct="poly(LogHPD,2)",
randomStruct="(1|SS)+(1|SSB)")

PlotLMERContinuous(model=m1$model,data=m1$data,effects="LogHPD",
xlab="Human Population Density",ylab="Species richness",logLink="e",
plotRug=TRUE,line.cols="#ff0000")

# Example plot of a model with both human population density and land use
# Note the requirement here to set the otherFactors parameter
m2<-FixedModel(all.data=sites,responseVar="Species_richness",
fitFamily="poisson",fixedStruct="LandUse+poly(LogHPD,2)",
randomStruct="(1|SS)+(1|SSB)")

PlotLMERContinuous(model=m2$model,data=m2$data,effects="LogHPD",
otherFactors=list(LandUse="Natural"),xlab="Human Population Density",
ylab="Species richness",logLink="e",plotRug=TRUE)

# Example plot of a model with an interaction between human population density
# and land use
m3<-FixedModel(all.data=sites,responseVar="Species_richness",
fitFamily="poisson",
fixedStruct="LandUse+poly(LogHPD,2)+LandUse:poly(LogHPD,2)",
randomStruct="(1|SS)+(1|SSB)")

PlotLMERContinuous(model=m3$model,data=m3$data,effects="LogHPD",
xlab="Human Population Density",ylab="Species richness",logLink="e",
plotRug=TRUE,byFactor="LandUse",seMultiplier=0.5)


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
