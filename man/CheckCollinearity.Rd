\name{CheckCollinearity}
\alias{CheckCollinearity}
\title{Check collinearity among variables}
\description{
Check collinearity among pairs of variables in a data frame
}
\usage{
  CheckCollinearity(data,allTerms)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
    A data frame
  }
  \item{allTerms}{
    The terms among which to check collinearity
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

# Check the collinearity among some of the explanatory variables in this dataset
# (land use, and rescaled variables describing human population density and 
# distance to nearest road)
CheckCollinearity(data = PREDICTSSites,allTerms = c("logHPD.rs","logDistRd.rs","LandUse"))
  }
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
    