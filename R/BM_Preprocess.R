###
###
### file: BM_Preprocess.R
###
### aim: Pre-procesing functions applied to BufferedMatrix objects
###
### Copyright (C) 2006 B. M. Bolstad
###
### Written by; B. M. Bolstad <bmb@bmbolstad.com>
###
### History: 
### June 27, 2006 - Initial version including support only for quantile normalization
### Aug 9, 2006 - add in the RMA background correction
### Aug 10, 2006 - add in median polish summarization



normalize.BufferedMatrix.quantiles <- function(x,copy=TRUE){

  if (!is(x,"BufferedMatrix")){
    stop("Need BufferedMatrix")
  }

  if (copy){
    x.copy <- duplicate(x)
  } else {
    x.copy <- x
  }
  


  x.copy@rawBufferedMatrix <- .Call("R_bm_quantile_normalize",x.copy@rawBufferedMatrix,PACKAGE="BufferedMatrixMethods")

  return (x.copy)
}




bg.correct.BufferedMatrix <- function(x, copy=TRUE){

  
  if (!is(x,"BufferedMatrix")){
    stop("Need BufferedMatrix")
  }
  
  if (copy){
    x.copy <- duplicate(x)
  } else {
    x.copy <- x
  }

  bg.dens <- function(x) {
    density(x, kernel = "epanechnikov", n = 2^14)
  }


  x.copy@rawBufferedMatrix <- .Call("R_bm_rma_bg_correct", x.copy@rawBufferedMatrix, body(bg.dens), new.env(), PACKAGE="BufferedMatrixMethods")

  return (x.copy)

 
}


setGeneric("median.polish.summarize", function(x,...) standardGeneric("median.polish.summarize"))

setMethod("median.polish.summarize", "BufferedMatrix", function(x,nProbeSets,ProbeNames){
####median.polish.summarize.BufferedMatrix <- function(x,nProbeSets,ProbeNames){

  if (length(ProbeNames) == dim(x)[1]){
    return(.Call("R_bm_summarize_medianpolish", x@rawBufferedMatrix, nProbeSets, ProbeNames, PACKAGE="BufferedMatrixMethods"))
  } else {
    stop("ProbeNames argument is of incorrect length")
  }
})
