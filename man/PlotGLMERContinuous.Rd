\name{PlotGLMERContinuous}
\alias{PlotGLMERContinuous}
\title{Plot modelled effects of continuous variables with uncertainty}
\description{
Plots estimated effects on a specified response variable of one or two specified continuous variables, from a generalized linear mixed-effects model. Where two continuous effects are specified, a perspective plot is used. Alternatively if byContEffect is used to specify the second effect, the first effect is shown at 3 levels of the second effect (2.5th percentile, median and 97.5th percentile).
}
\usage{
PlotGLMERContinuous(model,data,effects,byContEffect=NULL,
                        otherContEffects=character(0),
                        otherFactors=character(0),xlab,ylab,
                        byFactor=NULL,zlab=NULL,outDir=NULL,
                        logLink="n",plotRug=FALSE,seMultiplier=1.96,
                        params=list(),add=FALSE,xlim=NULL,ylim=NULL,
                        zlim=NULL,line.cols=NULL,line.types=NULL,
                        plotUncertainty=TRUE,nPanels=1,main=NULL,
                        yDiv=1,transformX=FALSE)
}
\arguments{
  \item{model}{
The model to derive parameter estimates from
}
  \item{data}{
A data frame containing all variables in the model
}
  \item{effects}{
The continuous variables (up to a maximum of 2) to plot an a continuous axis
}
  \item{byContEffect}{
An  optional interacting continuous variable to plot at three levels (2.5th percentile, 
median and 97.5th percentile)  
}
  \item{otherContEffects}{
Any continuous variables other than those selected for display that were fitted in the model. 
Median values of these variables will be used as the reference level.
}
  \item{otherFactors}{
Any factors that were fitted in the model, given as a list where the names are the names of the factors and the values are the factor levels for which to plot; shouldn't include any factors specified in byFactor
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
A list of R graphical parameters, which will overwrite the defaults
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
  \item{line.cols}{
The colours to use for the lines
}
  \item{line.types}{
The line types to use (solid, dashed, dotted etc.). Uses the same codes as R's graphical parameter 'lty'  
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
Factor to divide y-axis values by (useful for very large numbers, 
which would otherwise be printed in standard form)  
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
Tim Newbold <t.newbold@ucl.ac.uk>
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
# Load example data (site-level effects of land use on biodiversity from the PREDICTS database)
data(PREDICTSSiteData)

# Run a model of species richness as a function of land use, human population density
# and distance to nearest road (with an interaction between human population density
# and road distance)
m1 <- GLMER(modelData = PREDICTSSites,responseVar = "Species_richness",fitFamily = "poisson",fixedStruct = "LandUse+poly(logHPD.rs,2)+poly(logDistRd.rs,2)+poly(logHPD.rs,2):poly(logDistRd.rs,2)",randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)",REML = TRUE)

# Plot the effect of human population density
PlotGLMERContinuous(model = m1$model,data = m1$data,effects = "logHPD.rs",otherContEffects = "logDistRd.rs",otherFactors = list(LandUse="Primary Vegetation"),xlab = "Human Population Density",ylab = "Species richness",logLink = "e",plotRug = TRUE)

# Plot a perspective plot of the effects on species richness of human population density
# and distance to nearest road
PlotGLMERContinuous(model = m1$model,data = m1$data,effects = c("logHPD.rs","logDistRd.rs"),otherFactors = list(LandUse="Primary Vegetation"),xlab = "Human population density",ylab = "Distance to road",zlab = "Species richness",logLink = "e",plotUncertainty = FALSE)

# Run a model of species as a function of land use, human population density, and their interaction
m1 <- GLMER(modelData = PREDICTSSites,responseVar = "Species_richness",fitFamily = "poisson",fixedStruct = "LandUse+poly(logHPD.rs,2)+LandUse:poly(logHPD.rs,2):",randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)",REML = TRUE)

# Plot the effect on species richness of the interaction between human population density and land use
PlotGLMERContinuous(model = m1$model,data = m1$data,effects = "logHPD.rs",otherFactors = list(LandUse="Primary Vegetation"),xlab = "Human population density",ylab = "Species richness",byFactor = "LandUse",logLink = "e",line.cols=c("#E6AB02","#D95F02","#7570B3","#66A61E","#1B9E77","#E7298A"))
legend(0.75,15,c("Primary","Secondary","Plantation","Cropland","Pasture","Urban"),
col=c("#66A61E","#1B9E77","#7570B3","#E6AB02","#D95F02","#E7298A"),bty="n",lty=1)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
