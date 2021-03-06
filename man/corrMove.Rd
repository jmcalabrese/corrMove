\name{corrMove}
\alias{corrMove}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate the movement correlation indices
}
\description{
Selects the AICc best model for each partition and returns the corresponding point estimate and 95\% confidence intervals for each MCI conditional on the selected model.
}
\usage{
corrMove(data, prts)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
  A corrData object generated by as.corrData.
}
  \item{prts}{
  A vector of partition points, generated by findPrts.}
}
\details{
One estimate for each index, plus confidence intervals on each estimate, is produced for each partition in the data. These single values of each partition are repeated for each timestamp within a given partition for plotting purposes. Returns a value of 0 for both the point estimate and confidence interval for any MCI that cannot be produced by the selected model. For example, if the CU (correlated drift, uncorrelated diffusion) model is selected for a given partition, the diffusive correlation index point estimate and confidence interval limits will all be returned as 0, as there is no diffusive correlation in the CU model.
}
\value{
A dataframe consisting of the estimated values of each of the three MCIs and their corresponding confidence intervals, the selected model, and the partition number, for each timestamp in the dataset.
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Justin M. Calabrese
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
%\examples{
%##---- Should be DIRECTLY executable !! ----
%##-- ==>  Define data, use random,
%##--	or do  help(data=index)  for the standard data sets.

%## The function is currently defined as
%function (x)
%{
%  }
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
