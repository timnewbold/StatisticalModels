\name{GLMERSelect}
\alias{GLMERSelect}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Backward stepwise selection of GLMER fixed effects}
\description{
Performs backward stepwise selection of fixed effects in a generalized linear mixed-effects model. Tests interaction terms first, and then drops them to test main effects. Main effects that are part of interaction terms will be retained, regardless of their significance as main effects
}
\usage{
GLMERSelect(modelData,responseVar,fitFamily,fixedFactors=
                         character(0),fixedTerms=list(),randomStruct,
                       fixedInteractions=character(0),fitInteractions=FALSE,
                      alpha=0.05,verbose=FALSE,saveVars=character(0),
                      optimizer="bobyqa",maxIters=10000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{modelData}{
A data frame containing the response variable, all fixed effects to be considered, and all terms in the specified random-effects structure
}
  \item{responseVar}{
The response variable to fit in the model
}
  \item{fitFamily}{
The family to use for the generalized linear mixed effects model
}
  \item{fixedFactors}{
The fixed-effect factors to consider in the model, specified as a vector of strings that correspond to the column names in modelData
}
  \item{fixedTerms}{
The fixed-effect continuous variables to consider in the model, specified as a list where the item names correspond to the column names in modelData and the values are integers specifying the maximum complexity of the polynomial term to fit for the variable
}
  \item{randomStruct}{
The random-effects structure to use  
}
  \item{fixedInteractions}{
Specific interaction terms to consider in the model, specified as a vector of strings with interacting terms separated by a ':'
}
  \item{fitInteractions}{
Whether to fit all two-way interactions between the fixed effects in the model. Default is FALSE
}
\item{alpha}{
The threshold P value used to determine the statistical significance of terms
}
  \item{verbose}{
Whether to report progress in detail
}
  \item{saveVars}{
Any variables in the original data frame to retain in the model data frame for later analysis  
}
  \item{optimizer}{
The GLMER optimizer to use. Options are 'bobyqa' (the default) and 'Nelder_Mead'
}
  \item{maxIters}{
The maximum number of iterations to allow by the optimizer (default is 10,000) 
}
}
\details{
The specified random-effects structures is fixed.

The model-selection routine starts with the most complex fixed-effects structure possible given the specified combination of explanatory variables and their interactions, and performs backward stepwise selection to obtain the minimum adequate model. Comparison of the fit of different models is based on likelihood-ratio tests against a specified threshold P value (alpha), which defaults to 0.05. Interaction terms are tested first, and then removed to test main effects. All main effects that are part of significant interaction terms are retained in the final model regardless of their significance as main effects. ML estimation is used during selection of terms, and then REML is used for fitting the final model.
}
\value{
model: the final minimum adequate model

data: the dataset used in fitting the models, i.e. a subset of the original data frame, containing only the variables fit in the model, variables specified to be saved, and with any rows containing NA values removed

stats: a table of statistics relating to each term considered

final.call: the call used to generate the final model
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

# Fit a model of log-transformed total abundance as a function of land use,
# human population density and distance to nearest road. Consider quadratic
# polynomials or simpler for the continuous effects. Do not consider
# interactions among fixed effects
m1 <- GLMERSelect(modelData = PREDICTSSites,responseVar = "LogAbund",
                                fitFamily = "gaussian",fixedFactors = "LandUse",
                                fixedTerms = list(logHPD.rs=2,logDistRd.rs=2),
                                randomStruct = "(1|SS)+(1|SSB)",verbose = TRUE)

# Fit a model of log-transformed total abundance as a function of land use,
# human population density and distance to nearest road. Consider quadratic
# polynomials or simpler for the continuous effects. Consider a single two-way
# interaction between land use and a quadratic effect of human population density
m2 <- GLMERSelect(modelData = PREDICTSSites,responseVar = "LogAbund",
                                fitFamily = "gaussian",fixedFactors = "LandUse",
                                fixedTerms = list(logHPD.rs=2,logDistRd.rs=2),
                                fixedInteractions=c("LandUse:poly(logHPD.rs,2)"),
                                randomStruct = "(1|SS)+(1|SSB)",verbose = TRUE)
                                
# Fit a model of log-transformed total abundance as a function of land use,
# human population density and distance to nearest road. Consider quadratic
# polynomials or simpler for the continuous effects. Consider all two-way
# interactions among fixed effects
m3 <- GLMERSelect(modelData = PREDICTSSites,responseVar = "LogAbund",
                                fitFamily = "gaussian",fixedFactors = "LandUse",
                                fixedTerms = list(logHPD.rs=2,logDistRd.rs=2),
                                fitInteractions=TRUE,
                                randomStruct = "(1|SS)+(1|SSB)",verbose = TRUE)
}