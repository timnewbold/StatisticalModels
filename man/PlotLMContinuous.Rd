\name{PlotLMContinuous}
\alias{PlotLMContinuous}
\title{Plot modelled effects of continuous variables with uncertainty}
\description{
Plots estimated effects on a specified response variable of one or two specified continuous variables, from a linear or generalized linear model. Where two continuous effects are specified, a perspective plot is used. If a factor is specified as 'byFactor', modelled effects of the continuous variable at each factor level are shown.
}
\usage{
PlotLMContinuous(model,terms,se.fit=FALSE,seMultiplier=1.96,
                             xlim=NULL,ylim=NULL,zlim=NULL,
                             xlab=NULL,ylab=NULL,zlab=NULL,
                             params=list(),
                             byFactor=NULL,
                             yTransform=function(x){return(x)},
                             xTransform=function(x){return(x)},
                             log="",exclude.too.far=0.2)
}
\arguments{
  \item{model}{
    The model to derive parameter estimates from
  }
  \item{terms}{
    The continuous variables (up to a maximum of 2) to plot an a continuous axis
  }
  \item{se.fit}{
    Whether to show uncertainty
  }
  \item{xlim}{
  Specified x-axis limits
  }
  \item{ylim}{
  Specified y-axis limits
  }
  \item{zlim}{
  Specified z-axis limits (for perspective plots of 2 continuous effects)
  }
  \item{xlab}{
  The x-axis label
  }
  \item{ylab}{
  The y-axis label
  }
  \item{zlab}{
  The z-axis label (for perspective plots of 2 continuous effects)
  }
  \item{params}{
  A list of R graphical parameters, which will overwrite the defaults
  }
  \item{byFactor}{
  Plot separate lines for each level in this factor
  }
  \item{yTransform}{
  A function to use to back transform y-axis values
  }
  \item{xTransform}{
  A function to use to back transform x-axis values
  }
  \item{log}{
  Whether to log-transform both the x- and y-axis values "xy", one set of values "x" or "y", or neither "" (the default) 
  }
  \item{exclude.too.far}{
  When plotting persepective plots based on 2 continuous variables, takes the 2-dimensional surface described by the two variables and excludes values that are too far from the data used to fit the model. See exclude.too.far function in the mgcv package
  }
  }
\details{
}
\references{
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
