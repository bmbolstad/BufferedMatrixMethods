\name{BufferedMatrix.read.probematrix}
\alias{BufferedMatrix.read.probematrix}
\title{Read CEL file data into PM or MM BufferedMatrix}
\description{
  Read CEL data into \code{\link[BufferedMatrix:BufferedMatrix-class]{BufferedMatrix}} objects.
}
\usage{
BufferedMatrix.read.probematrix(..., filenames = character(0),celfile.path=NULL,rm.mask = FALSE, rm.outliers = FALSE, rm.extra = FALSE, verbose = FALSE,which="pm",cdfname = NULL)

}
\arguments{
  \item{\dots}{ file names separated by comma.}
  \item{filenames}{ file names in a character vector.}
  \item{celfile.path}{path where CEL files are located}
  \item{rm.mask}{ should the spots marked as 'MASKS' set to \code{NA} ? }
  \item{rm.outliers}{ should the spots marked as 'OUTLIERS' set to \code{NA} }
  \item{rm.extra}{ if \code{TRUE}, overrides what is in \code{rm.mask} and \code{rm.oultiers} }
  \item{verbose}{ verbosity flag }
  \item{which}{should be either "pm", "mm" or "both"}
  \item{cdfname}{Used to specify the name of an alternative cdf package. If set to
     \code{NULL}, the usual cdf package based on Affymetrix' mappings
     will be used.}
  
}

\value{
  A list of one or two \code{\link[BufferedMatrix:BufferedMatrix-class]{BufferedMatrix}} objects. Each \code{\link[BufferedMatrix:BufferedMatrix-class]{BufferedMatrix}} objects is either PM or MM data. No \code{\link[affy:AffyBatch-class]{AffyBatch}} is created.
}

\author{Ben Bolstad \email{bmb@bmbolstad.com}}

\seealso{\code{\link[affy:AffyBatch-class]{AffyBatch}}, \code{\link[affy:read.affybatch]{read.affybatch}}}

