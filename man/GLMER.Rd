\name{GLMER}
\alias{GLMER}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Fit a fixed-structure generalized linear mixed-effects model}
\description{
Fit a generalized linear mixed-effects model with a fixed model structure
}
\usage{
GLMER(modelData,responseVar,fitFamily,fixedStruct,
                      randomStruct,saveVars=character(0),REML=TRUE,
                     optimizer="bobyqa",maxIters=10000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{modelData}{
A data frame containing the response variable, and all fixed and random effects in the specified model structure
}
  \item{responseVar}{
The response variable to fit in the model
}
  \item{fitFamily}{
The family to use for the generalized linear mixed effects model
}
  \item{fixedStruct}{
The fixed effects to include, in the format of a glmer model-call   
}
  \item{randomStruct}{
The random effects to include, in the format of a glmer call  
}
  \item{saveVars}{
Any variables in the original data frame to retain in the model data frame for later analysis  
}
  \item{REML}{
Whether to use Restricted Maximum Likelihood for fitting the model. Alternative is simple maximum likelihood. Default is to use REML  
}
  \item{optimizer}{
The GLMER optimizer to use. Options are 'bobyqa' (the default) and 'Nelder_Mead'
}
  \item{maxIters}{
The maximum number of iterations to allow by the optimizer (default is 10,000) 
  }
}
\details{
  }
\value{
model: the model

data: the dataset used in to fit the model, i.e. a subset of the original data frame, containing only the variables fit in the model, variables specified to be saved, and with any rows containing NA values removed
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
}