\name{BufferedMatrix.read.celfiles}
\alias{BufferedMatrix.read.celfiles}
\title{Read CEL file data into PM or MM BufferedMatrix}
\description{
  Read CEL data into
  \code{\link[BufferedMatrix:BufferedMatrix-class]{BufferedMatrix}}
  objects. 
}
\usage{
BufferedMatrix.read.celfiles(..., filenames = character(0),celfile.path=NULL)

}
\arguments{
  \item{\dots}{ file names separated by comma.}
  \item{filenames}{ file names in a character vector.}
  \item{celfile.path}{path where CEL files are located}
}

\value{
  A \code{\link[BufferedMatrix:BufferedMatrix-class]{BufferedMatrix}}
  object containing the CEL file intensities. 
}

\author{Ben Bolstad \email{bmb@bmbolstad.com}}

\seealso{\code{\link[BufferedMatrix:BufferedMatrix-class]{BufferedMatrix}},
  \code{\link[BufferedMatrixMethods:BM_affyDataInput]{BufferedMatrix.read.probematrix}}}

