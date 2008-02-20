\name{RMA preprocess BufferedMatrix}
\alias{bg.correct.BufferedMatrix}
\alias{normalize.BufferedMatrix.quantiles}
\alias{median.polish.summarize}
\alias{BufferedMatrix.bg.correct.normalize.quantiles}
\alias{median.polish.summarize,BufferedMatrix-method}

\title{RMA preprocessing functions that work on BufferedMatrix objects}
\description{This group of functions can be used to apply the RMA
  background correction, Quantile normalization and Median polish
  summarization to data stored in a BufferedMatrix object.
}
\usage{
bg.correct.BufferedMatrix(x, copy=TRUE)
normalize.BufferedMatrix.quantiles(x,copy=TRUE)
BufferedMatrix.bg.correct.normalize.quantiles(x,copy=TRUE)
%%%median.polish.summarize(x,nProbeSets,ProbeNames)
}
\arguments{
\item{x}{a \code{\link[BufferedMatrix]{BufferedMatrix}} containing data
to be processed}
\item{copy}{should the \code{\link[BufferedMatrix]{BufferedMatrix}} be
copied or should the input object be changed on output}
%%%\item{nProbeSets}{Number of probesets stored in the
%%%\code{\link[BufferedMatrix]{BufferedMatrix}}}
%%%\item{ProbeNames}{Contains row names linking probes with probesets}
}
\value{In the case of \code{normalize.BufferedMatrix.quantiles} and
\code{bg.correct.BufferedMatrix} a
\code{\link[BufferedMatrix]{BufferedMatrix}} is returned. The function
\code{median.polish.summarize} returns a \code{\link{matrix}}.

The function \code{BufferedMatrix.bg.correct.normalize.quantiles} carries out both pre-processing steps with a single command.
}
\seealso{\code{\link[affy]{rma}}
}
\author{B. M. Bolstad \email{bmb@bmbolstad.com}}
\keyword{manip}
