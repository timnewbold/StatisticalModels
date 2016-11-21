\name{R2GLMER}
\alias{R2GLMER}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{GLMER Pseudo-R2}
\description{
Calculates psuedo-R2 values for GLMER and LMER models
}
\usage{
R2GLMER(model)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
  A mixed-effects model, of class 'lmerMod' or 'glmerMod'
  }
}
\details{
Calculated conditional (i.e. variance explained by fixed and random effects) and marginal (fixed effects only) R2 values, using the method proposed by Nakagawa & Schielzeth (2013)
}
\value{
conditional: the conditional R2 value, i.e. the variance explained by fixed and random effects
marginal: the marginal R2 value, i.e. the variance explained by the fixed effects
}
\references{
Nakagawa, S. & Schielzeth, H. (2013). A general and simple method for obtaining R2 from generalized linear mixed-effects models. <i> Methods in Ecology & Evolution</i> <b>4</b>: 133-142.
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
# Load example PREDICTS site-level data for the UK
data(sites)

# Build a model
m1<-FixedModel(all.data=sites,responseVar="Species_richness",
fitFamily="poisson",
fixedStruct="LandUse",
randomStruct="(1|SS)+(1|SSB)+(1|SSBS)")

# Perform test for over-dispersion
R2GLMER(m1$model)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
    