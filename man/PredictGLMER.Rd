\name{PredictGLMER}
\alias{PredictGLMER}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Predicted values for GLMERs}
\description{
Generates predicted values from a generalized linear mixed-effects model and a data frame with values of the explanatory variables
}
\usage{
PredictGLMER(model,data,se.fit=FALSE,seMultiplier = 1.96)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
  A mixed-effects model, of class 'lmerMod' or 'glmerMod'
  }
  \item{data}{
  A data frame containing values of the explanatory variables for which to make predictions
  }
  \item{se.fit}{
  Whether to estimate uncertainty around the predictions (default is False)
  }
  \item{seMultiplier}{
  The multiplier to apply to the uncertainty estimates (default is 1.96, which generates 95% confidence limits)
  }
}
\details{
Code for calculating predicted values and confidence intervals was taken from the GLMM wiki (see references). 
}
\value{

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
m1 <- GLMER(modelData = PREDICTSSites,responseVar = "Species_richness",fitFamily = "poisson",
fixedStruct = "LandUse",randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)",REML = TRUE)

# Generate predicted values for each land use
preds <- PredictGLMER()
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
    