\name{GLMEROverdispersion}
\alias{GLMEROverdispersion}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Test for over-dispersion}
\description{
Tests for over-dispersion in the residuals of a mixed-effects model
}
\usage{
GLMEROverdispersion(model)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
  A mixed-effects model, of class 'lmerMod'
  }
}
\details{
Compares residual deviance and residual degrees of freedom of a mixed-effects model to a chi-sqaure distribution. Adapted from code taken from the GLMM Wiki (see references) 
}
\value{
Reports the residual deviance, residual degrees of freedom and the ratio of these, and also a the P-value of a chi-square test comparing the residual deviance and degrees of freedom to a chi-square distribution
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
fixedStruct = "LandUse",randomStruct = "(1|SS)+(1|SSB)",REML = TRUE)

# Test for overdispersion
GLMEROverdispersion(m1$model)

# Try adding an observation-level random intercept (which can control for overdispersion)
m2 <- GLMER(modelData = PREDICTSSites,responseVar = "Species_richness",fitFamily = "poisson",
fixedStruct = "LandUse",randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)",REML = TRUE)

# Test for overdispersion again
GLMEROverdispersion(m2$model)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
    