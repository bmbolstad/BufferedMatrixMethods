###
### File: BM_justRMA.R
### 
### aim: implementation of RMA using BufferedMatrix objects
###
### History
### Feb 18, 2008 - Initial version (but really a port of justRMALite() that was in
###                AffyExtensions for several years previous)
###



BufferedMatrix.justRMA <- function(..., filenames=character(0),celfile.path=NULL,
                     phenoData=new("AnnotatedDataFrame"),
                     description=NULL,
                     notes="",
                     verbose=FALSE, background=TRUE, normalize=TRUE,
                     cdfname = NULL) {

  ### these have never been allowed to be TRUE for decent RMA performance
  ### although they have been lying around in the parameter list for the various versions		     
  rm.mask <- FALSE
  rm.outliers <- FALSE
  rm.extra <- FALSE		     

  require("affy")

  n <- length(filenames)


  pdata <- pData(phenoData)
  ##try to read sample names form phenoData. if not there use CEL filenames
  if(dim(pdata)[1]!=n){#if empty pdata filename are samplenames
    warning("Incompatible phenoData object. Created a new one.\n")

    samplenames <- gsub("^/?([^/]*/)*", "", unlist(filenames))
    pdata <- data.frame(sample=1:n,row.names=samplenames)
    phenoData <- new("AnnotatedDataFrame",
                     data=pdata,
                     varMetadata=data.frame(
                       labelDescription="arbitrary numbering",
                       row.names="sample"))
  }
  else samplenames <- rownames(pdata)

  if (is.null(description))
    {
      description <- new("MIAME")
      description@preprocessing$filenames <- filenames
      description@preprocessing$affyversion <- library(help=affy)$info[[2]][[2]][2]
    }


  if(is.null(cdfname)){
        headdetails <- read.celfile.header(filenames[[1]])
	cdfname <- headdetails[[1]]
  }
  tmp <- new("AffyBatch",
             cdfName=cdfname,
             annotation=cleancdfname(cdfname, addcdf=FALSE))
  pmIndex <- pmindex(tmp)
  probenames <- rep(names(pmIndex), unlist(lapply(pmIndex,length)))
  pmIndex <- unlist(pmIndex)

  tmp.buffmat <- BufferedMatrix.read.probematrix(..., filenames = filenames,celfile.path=celfile.path,
                             rm.mask = rm.mask, rm.outliers = rm.outliers, rm.extra = rm.extra, verbose = verbose,which="pm",cdfname = cdfname)	
			     
  ngenes <- length(geneNames(tmp))

  if (background & normalize){
  if (verbose){
       cat("Background Correcting and Normalizing\n")
    }
    tmp.buffmat <- BufferedMatrix.bg.correct.normalize.quantiles(tmp.buffmat,copy=FALSE)
  } else {
    ## Background correction
    if (background){
       if (verbose){
       	  cat("Background correcting\n")
       }	  
       tmp.buffmat <- bg.correct.BufferedMatrix(tmp.buffmat,copy=FALSE)
    }
    ## Normalization

    if (normalize){ 
       if (verbose){
       	  cat("Normalizing\n")
       }
       tmp.buffmat <- normalize.BufferedMatrix.quantiles(tmp.buffmat,copy=FALSE)
    }
  }
  ## Summarization

  set.buffer.dim(tmp.buffmat,cols=1,rows=5000)
  RowMode(tmp.buffmat)
  
   if (verbose){
     cat("Summarizing\n")
   }
  exprs <- median.polish.summarize(tmp.buffmat,ngenes,probenames)

  colnames(exprs) <- samplenames
  se.exprs <- array(NA, dim(exprs),
                    dimnames=list(rownames(exprs), colnames(exprs)))

  annotation <- annotation(tmp)
  notes(description) <- notes
  new("ExpressionSet",
      phenoData = phenoData,
      annotation = annotation,
      experimentData = description,
      exprs = exprs, se.exprs = se.exprs)  
}
