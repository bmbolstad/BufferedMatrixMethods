\name{BufferedMatrix.justRMA}
\alias{BufferedMatrix.justRMA}
\title{Use BufferedMatrix objects to facilitate RMA computation with low memory overhead}
\description{
  Read CEL data into \code{\link[BufferedMatrix:BufferedMatrix-class]{BufferedMatrix}} objects.
}
\usage{
BufferedMatrix.justRMA(..., filenames=character(0),celfile.path=NULL,
                     phenoData=new("AnnotatedDataFrame"),
                     description=NULL,
                     notes="",
                     verbose=FALSE, background=TRUE, normalize=TRUE,
                     cdfname = NULL)
}
\arguments{
  \item{\dots}{ file names separated by comma.}
  \item{filenames}{ file names in a character vector.}
  \item{celfile.path}{path where CEL files are located}
  \item{phenoData}{ a
    \code{\link[Biobase:class.AnnotatedDataFrame]{AnnotatedDataFrame}}
    object} 
  \item{description}{ a \code{\link[Biobase:class.MIAME]{MIAME}} object }
  \item{notes}{ notes }
  \item{verbose}{ verbosity flag }
  \item{normalize}{logical value. If \code{TRUE} normalize data using
    quantile normalization}
  \item{background}{logical value. If \code{TRUE} background correct
    using RMA background correction}
  \item{cdfname}{Used to specify the name of an alternative cdf
    package. If set to 
     \code{NULL}, the usual cdf package based on Affymetrix' mappings
     will be used.} 


}

\value{
    An \code{ExpressionSet} object, containing expression values identical to
  what one would get from running \code{\link[affy:rma]{rma}} on an
    \code{\link[affy:AffyBatch-class]{AffyBatch}}. 

}

\author{Ben Bolstad \email{bmb@bmbolstad.com}}

\seealso{\code{\link[BufferedMatrix:BufferedMatrix-class]{BufferedMatrix}},
  \code{\link[BufferedMatrixMethods:BM_affyDataInput]{BufferedMatrix.read.probematrix}}}  

