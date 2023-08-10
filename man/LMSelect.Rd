\name{LMSelect}
\alias{LMSelect}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Backward stepwise selection of a GLM}
\description{
Performs backward stepwise selection of terms in a generalized linear model. Tests interaction terms first, and then drops them to test main effects. Main effects that are part of interaction terms will be retained, regardless of their significance as main effects
}
\usage{
LMSelect=(modelData,responseVar,fitFamily,factors=
                       character(0),contEffects=list(),
                     interactions=character(0),
                     allInteractions=FALSE,
                     saveVars=character(0))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{modelData}{
A data frame containing the response variable, and all terms to be considered
}
  \item{responseVar}{
The response variable to fit in the model
}
  \item{fitFamily}{
The family to use for the generalized linear model
}
  \item{factors}{
The factors to consider in the model, specified as a vector of strings that correspond to the column names in modelData
}
  \item{contEffects}{
The continuous variables to consider in the model, specified as a list where the item names correspond to the column names in modelData and the values are integers specifying the maximum complexity of the polynomial term to fit for the variable
}
  \item{interactions}{
Specific interaction terms to consider in the model, specified as a vector of strings with interacting terms separated by a ':'
}
  \item{allInteractions}{
Whether to fit all two-way interactions between the fixed effects in the model. Default is FALSE
}
\item{alpha}{
The threshold P value used to determine the statistical significance of terms
}
  \item{saveVars}{
Any variables in the original data frame to retain in the model data frame for later analysis  
}
}
\details{

The model-selection routine starts with the most complex structure possible given the specified combination of explanatory variables and their interactions, and performs backward stepwise selection to obtain the minimum adequate model. Comparison of the fit of different models is based on likelihood-ratio tests, against a specified threshold P value (alpha), which defaults to 0.05. Interaction terms are tested first, and then removed to test main effects. All main effects that are part of significant interaction terms are retained in the final model regardless of their significance as main effects.
}
\value{
model: the final minimum adequate model

data: the dataset used in fitting the models, i.e. a subset of the original data frame, containing only the variables fit in the model, variables specified to be saved, and with any rows containing NA values removed

stats: a table of statistics relating to each term considered

final.call: the call used to generate the final model

family: the family of generalized linear model used - gaussian, poisson, binomial etc.
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
}