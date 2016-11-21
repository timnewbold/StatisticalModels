\name{PlotGLMERFactor}
\alias{PlotGLMERFactor}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Make an error-bar plot showing model coefficients}
\description{
Plots modelled effects of specified factors on a specified response variable.
}
\usage{
PlotGLMERFactor(model,data,responseVar,seMultiplier=1.96,
                     logLink="n",catEffects=NULL,
                     xtext.srt=0,ylim=NA,yaxp=NULL,order=NULL,rescale=NULL,
                     errbar.cols=NULL,pt.pch=NULL,
                     errbar.lty=1,
                     params=list(),add=FALSE,offset=0,
                     plotLabels=TRUE,cex.txt=NULL,pt.cex=1,
                     pt.bg="white",main=NULL,type="percent")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
The model to derive parameter estimates from
}
  \item{data}{
A data frame containing all variables in the model
}
  \item{responseVar}{
The name of the response variable fitted, as text, to display on the y-axis of the graph
}
  \item{seMultiplier}{
The number of standard errors to use for the error bars. Default is 1.96, which shows 95 percent confidence intervals
}
  \item{logLink}{
The transformation applied to the response variable before or during modelling: one of "n" (identity), "e" (log, base e), "10" (log, base 10) or "b" (logit - used in binomial models)
}
  \item{catEffects}{
The categorical effects (factors) to plot, as a vector  
}
  \item{xtext.srt}{
The number of degrees by which to rotate the x-axis labels
}
  \item{ylim}{
Fixed y-axis limits. If not specified, then resolved automatically
}
  \item{yaxp}{
Specified values to include on the y axis, as used in general R graphical parameters  
}
  \item{order}{
The order to plot the categorical factor levels in. If not specified, then as in the original data frame
}
  \item{rescale}{
The amount by which to rescale y-axis values for the categorical effects by. If not specified, then no rescaling is performed
}
  \item{errbar.cols}{
Colours to use for the error bars
}
  \item{pt.pch}{
The point types (using R's 'pch' values), with which to overwrite the default
}
  \item{errbar.lty}{
The line types to use for the error bars; values as in the general R graphical parameter 'lty'  
}
  \item{params}{
Any R graphical parameters, which will overwrite the defaults  
}
  \item{add}{
Whether to add the plot to an existing plot; default is FALSE
}
  \item{offset}{
The amount by which to offset points on the x axis; default is zero
}
  \item{plotLabels}{
Whether to plot individual labels for each factor level  
}
  \item{cex.txt}{
Text size for additional labels (i.e. text not on the axes or in the axis labels, which should be specified using the params argument)  
}
  \item{pt.cex}{
Point size to use. Default is 1.  
}
  \item{pt.bg}{
The background colour to use for the points (as long as this is compatible with the point type used)  
}
  \item{main}{
The title for the plot  
}
  \item{type}{
Whether to plot as a percentage change ('percent' - the default) or directly as the response variable ('response')  
}
}
\details{
Code for calculating predicted values and confidence intervals for continuous effects was taken from the GLMM wiki (see references).
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

PREDICTSSites$LandUse <- factor(PREDICTSSites$LandUse,levels=c("Primary Vegetation","Secondary Vegetation",
"Plantation forest","Cropland","Pasture","Urban"))

# Run a model of log-transformed total abundance as a function of land use, human population density
# and distance to nearest road (with an interaction between human population density
# and road distance)
m1 <- GLMER(modelData = PREDICTSSites,responseVar = "LogAbund",fitFamily = "gaussian",fixedStruct = "LandUse+poly(logHPD.rs,2)+poly(logDistRd.rs,2)+poly(logHPD.rs,2):poly(logDistRd.rs,2)",randomStruct = "(1|SS)+(1|SSB)",REML = TRUE)

# Plot the effect of land use as an error bar, showing +/- 1 standard error around the mean coefficient estimates
PlotGLMERFactor(model = m1$model,data = m1$data,responseVar = "Total abundance",seMultiplier = 1,
logLink = "e",catEffects = "LandUse")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
