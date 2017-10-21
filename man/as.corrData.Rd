\name{as.corrData}
\alias{as.corrData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Coerce a multi-individual movement dataset into a corrMove object.
}
\description{
Currently this only converts the timestamp field to POSIXct.
}
\usage{
as.corrMove(object, timeformat="", timezone="UTC")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
A multi-individual movement dataset, where each individual is observed at each timestamp.
}
\item{timeformat}{
Format arguement for strptime.
}
\item{timezone}{
Timezone arguement for strptime.
}
}
\details{
The data file must consist of a timestamp column with the name "timestamp", followed by a column of the x coordinates for the first individual, a column of the x coordinates of the second, individual, etc., then a column of the y coordinates of the first indivdiual, then a column of the y coordinates of the second individual. The timestamps can be in any standard date/time format provided the user supplies the appropriate formating string. Currently, the x and y coordinates must be in the project the user wishes to do the analysis in (e.g., UTM), as no reprojecting is done within the package.
}
%\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
%}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Chris H. Fleming and Justin M. Calabrese
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
