\name{findPrts}
\alias{findPrts}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Find partitions in a multi-individual movement dataset
}
\description{
Sequentially introduces partition points into a dataset by finding the partition point that results in the largest IC reduction. The algorithm stops either when no further IC reduction is possible, or when the minimum window length W is exceeded for the first time.
}
\usage{
findPrts(data, W, IC=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
    A corrData object from as.corrData.
}
\item{W}{
The user-specified minimum window width.
}
  \item{IC}{Information criterion. Current choices are AICc=1 (the default) or BIC=2.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A vector of loctions of the partition points in the dataset.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Justin M. Calabrese
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
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
\keyword{ Partition }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
