###
### File: BM_affyDataInput.R
### 
### aim: implement reading of Affymetrix CEL file data into BufferedMatrix object
###
### History
### Feb 18, 2008 - Initial version
###

checkValidFilenames <- function(filenames) {
    ## Returns TRUE if filenames is a character vector containing
    ## paths to files that exist (directories don't count).
    ## A suitable error message is printed via stop() if invalid
    ## file names are encountered.
    if (!is.character(filenames))
      stop(strwrap(paste("file names must be specified using a character",
                         "vector, not a", sQuote(typeof(filenames)))),
           call.=FALSE)

    if (length(filenames) == 0)
      stop("no file names provided")

    if (any(sapply(filenames, nchar) < 1))
      stop("empty file names are not allowed")

    finfo <- file.info(filenames)
    whBad <- sapply(finfo[["isdir"]], function(x) !identical(FALSE, x))
    if (any(whBad)) {
        msg <- paste("the following are not valid files:\n",
                     paste("  ", filenames[whBad], collapse="\n"))
        stop(msg, call.=FALSE)
    }
    TRUE
}

list.celfiles <-   function(...){
  files <- list.files(...)
  return(files[grep("\\.[cC][eE][lL]\\.gz$|\\.[cC][eE][lL]$", files)])
}

GetCelNames <- function(...,filenames=character(0), celfile.path=NULL){
  auxnames <- unlist(as.list(substitute(list(...)))[-1])

  if(!is.null(celfile.path)){
    auxnames <- file.path(celfile.path, auxnames)
    filenames <- file.path(celfile.path, filenames)
  }

  filenames <- c(filenames, auxnames)

  if(length(filenames)==0){
    if(is.null(celfile.path)) celfile.path <- getwd()
    filenames <- list.celfiles(celfile.path,full.names=TRUE)
  }
  if(length(filenames)==0) stop("No cel filennames specified and no cel files in specified directory:",celfile.path,"\n")
  
  filenames
}


BufferedMatrix.read.celfiles <- function(...,filenames=character(0),
                     celfile.path=NULL){

  require(affyio)			     

  filenames <- GetCelNames(..., filenames=filenames,celfile.path=celfile.path)

			         
  checkValidFilenames(filenames)

  samplenames <- gsub("^/?([^/]*/)*", "", unlist(filenames), extended=TRUE    )

  headdetails <- .Call("ReadHeader",as.character(filenames[[1]]), PACKAGE="affyio")

  for (filenm in filenames[-1]){
      temp.head <- .Call("ReadHeader",as.character(filenm[[1]]), PACKAGE="affyio")

      if (all(headdetails[[2]] != temp.head[[2]])){
        stop(paste("Dimension mismatch in CEL file",filenm[[1]],"\n"))
      }
      if (headdetails[[1]] != temp.head[[1]]){
        stop(paste("Chip type mismatch in CEL file",filenm[[1]], "Found:",temp.head[[1]],"Expected:",headdetails[[1]]))
      }
  }


  tmp.buffmat <- createBufferedMatrix(rows=headdetails[[2]][1]*headdetails[[2]][1],cols=1)

  tmp.buffmat[,1] <- read.celfile(filenames[[1]],intensity.means.only=TRUE )$INTENSITY$MEAN

  i <- 2
  for (filenm in filenames[-1]){
    AddColumn(tmp.buffmat)					 					 
    tmp.buffmat[,i] <- read.celfile(filenm[[1]],intensity.means.only=TRUE)$INTENSITY$MEAN
    i <- i+1
  }
  colnames(tmp.buffmat) <- samplenames
  tmp.buffmat 
}




BufferedMatrix.read.probematrix <- function(..., filenames = character(0),celfile.path=NULL,
                             rm.mask = FALSE, rm.outliers = FALSE, rm.extra = FALSE, verbose = FALSE,which="pm",
                             cdfname = NULL){
  
  which <- match.arg(which,c("pm","mm","both"))

  require(affy)
  require(affyio)			     
  

  filenames <- GetCelNames(..., filenames=filenames,celfile.path=celfile.path)

			         
  checkValidFilenames(filenames)
 
  samplenames <- gsub("^/?([^/]*/)*", "", unlist(filenames), extended=TRUE    )

  headdetails <- .Call("ReadHeader",as.character(filenames[[1]]), PACKAGE="affyio")

  for (filenm in filenames[-1]){
      temp.head <- .Call("ReadHeader",as.character(filenm[[1]]), PACKAGE="affyio")

      if (all(headdetails[[2]] != temp.head[[2]])){
        stop(paste("Dimension mismatch in CEL file",filenm[[1]],"\n"))
      }
      if (headdetails[[1]] != temp.head[[1]]){
        stop(paste("Chip type mismatch in CEL file",filenm[[1]], "Found:",temp.head[[1]],"Expected:",headdetails[[1]]))
      }
  }

  dim.intensity <- headdetails[[2]]
  ref.cdfName <- headdetails[[1]]
  ## Allow for usage of alternative cdfs
  if(is.null(cdfname)){	
      Data <- new("AffyBatch", cdfName = ref.cdfName, annotation = cleancdfname(ref.cdfName,addcdf = FALSE))
  } else {
      Data <- new("AffyBatch", cdfName = cdfname, annotation = cleancdfname(ref.cdfName, addcdf = FALSE))
  }
  cdfInfo <- as.list(getCdfInfo(Data))
  cdfInfo <- cdfInfo[order(names(cdfInfo))]
  
  
  tmp <- .Call("read_probeintensities", filenames[[1]],
            				rm.mask, rm.outliers, rm.extra, ref.cdfName,
  	            			dim.intensity, verbose, cdfInfo,which, PACKAGE="affyio")
 
  if (which == "pm"){
    tmp.buffmat <- createBufferedMatrix(length(tmp$pm),cols=1)
    tmp.buffmat[,1] <- tmp$pm 
  } else if (which == "mm"){
    tmp.buffmat <- createBufferedMatrix(length(tmp$mm),cols=1)
    tmp.buffmat[,1] <- tmp$mm 
  } else {
      tmp.buffmat <- createBufferedMatrix(length(tmp$pm),cols=1)
      tmp.buffmat[,1] <- tmp$pm 
      tmp.buffmat2 <- createBufferedMatrix(length(tmp$mm),cols=1)
      tmp.buffmat2[,1] <- tmp$mm 
  }

  i <- 2
  for (filenm in filenames[-1]){
      if (which == "pm"){
      	 AddColumn(tmp.buffmat)					 					 
    	 tmp.buffmat[,i] <-  .Call("read_probeintensities", filenm[[1]],
            				rm.mask, rm.outliers, rm.extra, ref.cdfName,
              			dim.intensity, verbose, cdfInfo,which, PACKAGE="affyio")$pm
      } else if (which == "mm"){
	 AddColumn(tmp.buffmat)					 					 
    	 tmp.buffmat[,i] <-  .Call("read_probeintensities", filenm[[1]],
            				rm.mask, rm.outliers, rm.extra, ref.cdfName,
              			dim.intensity, verbose, cdfInfo,which, PACKAGE="affyio")$mm
      } else {
         tmp <-  .Call("read_probeintensities", filenm[[1]],
            				rm.mask, rm.outliers, rm.extra, ref.cdfName,
              			dim.intensity, verbose, cdfInfo,which, PACKAGE="affyio")
      	 AddColumn(tmp.buffmat)					 					 
    	 tmp.buffmat[,i] <- tmp$pm
     	 AddColumn(tmp.buffmat2)					 					 
    	 tmp.buffmat2[,i] <- tmp$mm
      }
     i <- i+1
  }
  if (which == "pm" || which == "mm"){
    colnames(tmp.buffmat) <- samplenames	
    return(tmp.buffmat)
  } else {
    colnames(tmp.buffmat) <- samplenames 
    colnames(tmp.buffmat2) <- samplenames	
    return(list(pm=tmp.buffmat,mm=tmp.buffmat2))
  }
 
}
