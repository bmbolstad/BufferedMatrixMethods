/*****************************************************
 **
 ** file: preprocess_bm.c
 **
 ** Copyright (C) 2006    B. M. Bolstad <bmb@bmbolstad.com>
 **
 ** aim: Pre-process routines that cope with Buffered
 **      Matrix obejcts
 **
 **  History
 **  Jun 27, 2006 - Initial version
 **  Aug  9, 2006 - develop background correction implementation
 **  Aug 10, 2006 - median polish summarization
 **  Nov 13, 2006 - make max algorithm faster
 ** 
 *****************************************************/


#include "doubleBufferedMatrix.h"


#include <R.h>
#include <Rdefines.h>
#include <Rmath.h>
#include <Rinternals.h>

/**********************************************************
 **
 ** int sort_double(const void *a1,const void *a2)
 ** 
 ** a comparison function used when sorting doubles.
 **
 **********************************************************/

int sort_double(const double *a1,const double *a2){
  if (*a1 < *a2)
    return (-1);
  if (*a1 > *a2)
    return (1);
  return 0;
}




/***********************************************************
 **
 ** double find_max(double *x, int length)
 **
 ** this function returns the max of x
 **
 ************************************************************/

double find_max(double *x,int length){
  int i;
  /*  double *buffer = Calloc(length,double); */
  double max;

  /*
    for (i=0; i < length; i++){
    buffer[i] = x[i];
    }
  */

  /* qsort(buffer,length,sizeof(double),(int(*)(const void*, const void*))sort_double); */

  max = x[0];
  for (i=1; i < length; i++){
    if (x[i] > max){
      max = x[i];
    }
  }
  


  /* printf("max is %f \n", max); */
  /* Free(buffer); */
  return max;
}

/***************************************************************
 **
 ** double get_sd(double *MM, double MMmax, int rows, int cols, int column)
 **
 ** double *PM - pm matrix
 ** double PMmax - value of mode of PMs for column
 ** int rows,cols - dimensions of matrix
 ** int column - column (chip) of interest
 **
 ** estimate the sigma parameter given vector MM value of maximum of density
 ** of MM, dimensions of MM matrix and column of interest
 **
 **
 ***************************************************************/

double get_sd(double *MM, double MMmax, int rows){
  double sigma;
  double tmpsum = 0.0;
  int numtop=0;
  int i;


  for (i=0; i < rows; i++){
    if (MM[i] < MMmax){
      tmpsum = tmpsum + (MM[i] - MMmax)*(MM[i] - MMmax);
      numtop++;
    }
  }
  sigma = sqrt(tmpsum/(numtop -1))*sqrt(2.0)/0.85;
  return sigma;

}






/**************************************************************************************
 **
 ** double max_density(double *z,int rows,int cols,int column, SEXP fn,SEXP rho)
 **
 ** double *z - matrix of dimension rows*cols
 ** int cols - matrix dimension
 ** int rows - matrix dimension
 ** int column - column of interest
 ** SEXP fn - R function for estimation of density
 ** SEXP rho - an R environment to work within
 **
 *************************************************************************************/

double max_density(double *z,int rows, SEXP fn,SEXP rho){

  int i;
  SEXP x;
  SEXP results;

  double *dens_x;
  double *dens_y;
  double max_y,max_x;
  
  PROTECT(x = allocVector(REALSXP,rows));

  for (i=0; i< rows; i++){
    REAL(x)[i] = z[i];
  }

  /*  for(i=0; i < rows; i++)
      printf("%f\n",REAL(x)[i]); */
  
  defineVar(install("x"),x,rho);
  PROTECT(results = eval(fn,rho));
  
  
  dens_x = NUMERIC_POINTER(VECTOR_ELT(results,0));
  dens_y = NUMERIC_POINTER(VECTOR_ELT(results,1));

  max_y = find_max(dens_y,16384);
  
  i = 0;
  do {
    if (dens_y[i] == max_y)
      break;
    i++;

  } while(1);
  
  max_x = dens_x[i];
  

  /*  PROTECT(names = getAttrib(results,R_NamesSymbol)); */
  

  /* for (i = 0; i < length(results); i++)
  //  printf("%S \n",CHAR(STRING_ELT(names,i)));

  //printf("max_x: %f\n",max_x);
  */

  UNPROTECT(2);


  return max_x;

}










/**********************************************************************************
 **
 ** double Phi(double x)
 **
 ** Compute the standard normal distribution function
 **
 *********************************************************************************/

double Phi(double x){
   return pnorm5(x,0.0,1.0,1,0);
}

/***********************************************************************************
 **
 ** double phi(double x)
 **
 ** compute the standard normal density.
 **
 **
 **********************************************************************************/

double phi(double x){
  return dnorm4(x,0.0,1.0,0);
}


/***************************************************************
 **
 ** double  get_alpha2(double *PM,double PMmax, int rows,int cols,int column)
 **
 ** estimate the alpha parameter given vector PM value of maximum of density
 ** of PM, dimensions of MM matrix and column of interest using method proposed
 ** in affy2
 **
 **
 ***************************************************************/

double get_alpha2(double *PM, double PMmax, int length,SEXP fn,SEXP rho){
  double alpha;
  /* double tmpsum = 0.0;
     int numtop=0; */
  int i;

  for (i=0; i < length; i++){
    PM[i] = PM[i] - PMmax;
  }

  alpha = max_density(PM,length, fn,rho);

  alpha = 1.0/alpha;
  return alpha ;  
}

/********************************************************************************
 **
 ** void bg_parameters2(double *PM,double *MM, double *param, int rows, int cols, int column,SEXP fn,SEXP rho)
 **
 ** estimate the parameters for the background, they will be returned in *param
 ** param[0] is alpha, param[1] is mu, param[2] is sigma.
 **
 ** parameter estimates are same as those given by affy in bg.correct.rma (Version 1.1 release of affy)
 **
 *******************************************************************************/


void bg_parameters2(double *PM, double *param, int rows, SEXP fn,SEXP rho){
  int i = 0;
  double PMmax;
 
  double sd,alpha;
  int n_less=0,n_more=0;
  double *tmp_less = (double *)Calloc(rows,double);
  double *tmp_more = (double *)Calloc(rows,double);
 
 
  PMmax = max_density(PM,rows,fn,rho);
  
  for (i=0; i < rows; i++){
    if (PM[i] < PMmax){
      tmp_less[n_less] = PM[i];
      n_less++;
    }
    
  }  
  
  PMmax = max_density(tmp_less,n_less,fn,rho);
  sd = get_sd(PM,PMmax,rows)*0.85; 

  for (i=0; i < rows; i++){
    if (PM[i] > PMmax) {
      tmp_more[n_more] = PM[i];
      n_more++;
    }
  }

  /* the 0.85 is to fix up constant in above */
  alpha = get_alpha2(tmp_more,PMmax,n_more,fn,rho);

  param[0] = alpha;
  param[1] = PMmax;
  param[2] = sd;

  /*  Rprintf("%f %f %f\n",param[0],param[1],param[2]); */

  
  Free(tmp_less);
  Free(tmp_more);
}




/************************************************************************************
 **
 ** void bg_adjust(double *PM,double *MM, double *param, int rows, int cols, int column)
 **
 ** double *PM - PM matrix of dimension rows by cols
 ** double *MM - MM matrix of dimension rows by cols
 ** double *param - background model parameters
 ** int rows, cols - dimension of matrix
 ** int column - which column to adjust
 **
 ** note we will assume that param[0] is alpha, param[1] is mu, param[2] is sigma
 **
 ***********************************************************************************/

void bg_adjust(double *PM, double *param, int rows){
  int i;
  double a;
  
  for (i=0; i < rows; i++){
    a = PM[i] - param[1] - param[0]*param[2]*param[2]; 
    PM[i] = a + param[2] * phi(a/param[2])/Phi(a/param[2]);
  }
  
}







void bm_rma_bg_correct(doubleBufferedMatrix Matrix, SEXP fn,SEXP rho){

  int rows, cols;
  int i,j;
  double *params;
  double *datvec;

  rows = dbm_getRows(Matrix);
  cols = dbm_getCols(Matrix);
  
  params = Calloc(3,double);

  datvec = (double *)Calloc(rows,double);
  
  for (j = 0; j < cols; j++){
    dbm_getValueColumn(Matrix,&j,datvec,1);
    bg_parameters2(datvec, params,rows,fn,rho);
    bg_adjust(datvec,params,rows);
    dbm_setValueColumn(Matrix,&j,datvec,1); 
  }

  Free(params);
  Free(datvec);

}







SEXP R_bm_rma_bg_correct(SEXP R_BufferedMatrix, SEXP fn,SEXP rho){



  doubleBufferedMatrix Matrix;
  int current_mode;



  Matrix =  R_ExternalPtrAddr(R_BufferedMatrix);

  if (Matrix == NULL){
    return R_BufferedMatrix;
  }





  bm_rma_bg_correct(Matrix, fn, rho);



  return R_BufferedMatrix;

}











/*************************************************************
 **
 ** the dataitem record is used to keep track of data indicies 
 ** along with data value when sorting and unsorting in the 
 ** quantile algorithm.
 **
 ************************************************************/

typedef struct{
  double data;
  int rank;
} dataitem;
  


/*************************************************************
 **
 ** the dataitem_block record is used to keep track of data indicies 
 ** along with data value when sorting and unsorting in the 
 ** quantile algorithm in blocks
 **
 ************************************************************/

typedef struct{
  double data;
  int rank;
  int block;
} dataitem_block;







/***********************************************************
 **  
 ** int min(int x1, int x2)							    
 **
 ** returns the minimum of x1 and x2
 **		    
 **********************************************************/

static int min(int x1,int x2){
  if (x1 > x2)
    return x2;
  else
    return x1;
}

/**********************************************************
 **
 ** int sort_fn(const void *a1,const void *a2)
 **
 ** a comparison function for sorting objects of the dataitem type.
 **
 **
 **********************************************************/

static int sort_fn(const void *a1,const void *a2){
  dataitem *s1, *s2;
  s1 = (dataitem *)a1;
  s2 = (dataitem *)a2;
  
  if (s1->data < s2->data)
    return (-1);
  if (s1 ->data > s2->data)
    return (1);
  return 0;
}




/************************************************************
 **
 ** double *get_ranks(dataitem *x,int n)
 **
 ** get ranks in the same manner as R does. Assume that *x is
 ** already sorted
 **
 *************************************************************/

static void get_ranks(double *rank, dataitem *x,int n){
  int i,j,k;
   
  i = 0;

  while (i < n) {
    j = i;
    while ((j < n - 1) && (x[j].data  == x[j + 1].data))
      j++;
    if (i != j) {
      for (k = i; k <= j; k++)
	rank[k] = (i + j + 2) / 2.0;
    }
    else
      rank[i] = i + 1;
    i = j + 1;
  }
  /*return rank;*/
}



void bm_quantile_normalize(doubleBufferedMatrix Matrix){


  int i,j; /* indexing variables */

  int rows, cols; /* Dimensions of the matrix */
  int ind;
  dataitem **dimat;
  double *row_mean;
  double *datvec; 
  double *ranks;

  rows = dbm_getRows(Matrix);
  cols = dbm_getCols(Matrix);

  datvec = (double *)Calloc(rows,double);
  row_mean = (double *)Calloc(rows,double);
    
  for (i =0; i < rows; i++){
    row_mean[i] = 0.0;
    } 


  for (j = 0; j < cols; j++){
 
    dbm_getValueColumn(Matrix,&j,datvec,1);
 
    qsort(datvec,rows,sizeof(double),(int(*)(const void*, const void*))sort_double); 
 
    for (i =0; i < rows; i++){
      row_mean[i] += datvec[i]/((double)cols);
    }
  }
 
  ranks = (double *)Calloc(rows,double);
  /* now assign back distribution */
  dimat = (dataitem **)Calloc(1,dataitem *);
  dimat[0] = (dataitem *)Calloc(rows,dataitem);

  for (j = 0; j < cols; j++){ 
 
    dbm_getValueColumn(Matrix,&j,datvec,1); 
 
    for (i =0; i < rows; i++){
      dimat[0][i].data = datvec[i];
      dimat[0][i].rank = i;
    }
 
    qsort(dimat[0],rows,sizeof(dataitem),sort_fn);
 
    get_ranks(ranks,dimat[0],rows);
 
    for (i =0; i < rows; i++){
      ind = dimat[0][i].rank;
      dbm_setValue(Matrix,ind,j,row_mean[(int)floor(ranks[i])-1]);
    }
  }
  
  Free(ranks);
  Free(datvec);
  Free(dimat[0]);

  Free(dimat);
  Free(row_mean);
  
}



SEXP R_bm_quantile_normalize(SEXP R_BufferedMatrix){



  doubleBufferedMatrix Matrix;
  int current_mode;



  Matrix =  R_ExternalPtrAddr(R_BufferedMatrix);

  if (Matrix == NULL){
    return R_BufferedMatrix;
  }





  bm_quantile_normalize(Matrix);



  return R_BufferedMatrix;

}






/**************************************************************************
 **
 ** double median(double *x, int length)
 **
 ** double *x - vector
 ** int length - length of *x
 **
 ** returns the median of *x
 **
 *************************************************************************/

double  median(double *x, int length){
  int i;
  int half;
  double med;
  double *buffer = Calloc(length,double);

  for (i = 0; i < length; i++)
    buffer[i] = x[i];

  qsort(buffer,length,sizeof(double), (int(*)(const void*, const void*))sort_double);
  half = (length + 1)/2;
  if (length % 2 == 1){
    med = buffer[half - 1];
  } else {
    med = (buffer[half] + buffer[half-1])/2.0;
  }

  Free(buffer);
  return med;
}

/*******************************************************************************
 **
 ** double sum_abs(double *z, int rows, int cols)
 **
 ** double *z - matrix of doubles
 ** int rows - dimension of matrix
 ** int cols - dimension of matrix
 **
 ** returns the sum of the absolute values of elements of the matrix *z
 **
 ******************************************************************************/

double sum_abs(double *z, int rows, int cols){

  int i, j;
  double sum = 0.0;

  for (i=0; i < rows; i++)
    for (j=0; j < cols; j++)
      sum+=fabs(z[j*rows+i]);

  return sum;
}

/********************************************************************************
 **
 ** void get_row_median(double *z, double *rdelta, int rows, int cols)
 **
 ** double *z - matrix of dimension  rows*cols
 ** double *rdelta - on output will contain row medians (vector of length rows)
 ** int rows, cols - dimesion of matrix
 **
 ** get the row medians of a matrix
 **
 ********************************************************************************/

void get_row_median(double *z, double *rdelta, int rows, int cols){
  int i,j;
  double *buffer = Calloc(cols,double);

  for (i = 0; i < rows; i++){
    for (j = 0; j < cols; j++){
      buffer[j] = z[j*rows + i];
    }
    rdelta[i] = median(buffer,cols);
  }

  Free(buffer);
}

/********************************************************************************
 **
 ** void get_col_median(double *z, double *cdelta, int rows, int cols)
 **
 ** double *z - matrix of dimension  rows*cols
 ** double *cdelta - on output will contain col medians (vector of length cols)
 ** int rows, cols - dimesion of matrix
 **
 ** get the col medians of a matrix
 **
 ********************************************************************************/

void get_col_median(double *z, double *cdelta, int rows, int cols){

  int i, j;

  double *buffer = Calloc(rows,double);
  for (j = 0; j < cols; j++){
    for (i = 0; i < rows; i++){
      buffer[i] = z[j*rows + i];
    }
    cdelta[j] = median(buffer,rows);
  }

  Free(buffer);

}

/***********************************************************************************
 **
 ** void subtract_by_row(double *z, double *rdelta, int rows, int cols)
 **
 ** double *z - matrix of dimension rows by cols
 ** double *rdelta - vector of length rows
 ** int rows, cols dimensions of matrix
 **
 ** subtract the elements of *rdelta off each row of *z
 **
 ***********************************************************************************/

void subtract_by_row(double *z, double *rdelta, int rows, int cols){

  int i,j;

  for (i = 0; i < rows; i++){
    for (j = 0; j < cols; j++){
      z[j*rows +i]-= rdelta[i];
    }
  }
}


/***********************************************************************************
 **
 ** void subtract_by_col(double *z, double *cdelta, int rows, int cols)
 **
 ** double *z - matrix of dimension rows by cols
 ** double *cdelta - vector of length rows
 ** int rows, cols dimensions of matrix
 **
 ** subtract the elements of *cdelta off each col of *z
 **
 ***********************************************************************************/

void subtract_by_col(double *z, double *cdelta, int rows, int cols){

  int i,j;
  for (j = 0; j < cols; j++){
    for (i = 0; i < rows; i++){
      z[j*rows +i]-= cdelta[j];
    }
  }

}

/***********************************************************************************
 **
 ** void rmod(double *r, double *rdelta, int rows)
 **
 ** double *r - vector of length rows
 ** double *rdelta - vector of length rows
 ** int rows, cols dimensions of matrix
 **
 ** add elementwise *rdelta to *r
 **
 ***********************************************************************************/


void rmod(double *r, double *rdelta, int rows){
  int i;

  for (i = 0; i < rows; i++){
    r[i]= r[i] + rdelta[i];
  }
}

/***********************************************************************************
 **
 ** void cmod(double *c, double *cdelta, int cols)
 **
 ** double *c - vector of length rows
 ** double *cdelta - vector of length rows
 ** int cols length of vector
 **
 ** add elementwise *cdelta to *c
 **
 ***********************************************************************************/

void cmod(double *c, double *cdelta, int cols){
  int j;

  for (j = 0; j < cols; j++){
    c[j]= c[j] + cdelta[j];
  }
}


/*************************************************************************************
 **
 ** void median_polish(double *data, int rows, int cols, int *cur_rows, double *results, int nprobes)
 **
 ** double *data - a data matrix of dimension rows by cols (the entire PM matrix)
 ** int rows, cols - rows and columns dimensions of matrix
 ** int cur_rows - vector of length nprobes containg row indicies of *data matrix which apply for a
 **                particular probeset
 ** double *results - a vector of length cols already allocated. on output contains expression values
 ** int nprobes - number of probes in current probeset.
 **
 ** a function to do median polish.
 **
 *************************************************************************************/

void median_polish(doubleBufferedMatrix Matrix, int rows, int cols, int *cur_rows, double *results, int nprobes){

  int i,j,iter;
  int maxiter = 10;
  double eps=0.01;
  double oldsum = 0.0,newsum = 0.0;
  double t = 0.0;
  double delta;
  double *rdelta = Calloc(nprobes,double);
  double *cdelta = Calloc(cols,double);

  double *r = Calloc(nprobes,double);
  double *c = Calloc(cols,double);
  double *z = Calloc(nprobes*cols,double);




  dbm_getValueRow(Matrix,cur_rows,z,nprobes);

  for (j = 0; j < cols; j++){
    for (i =0; i < nprobes; i++){
      z[j*nprobes + i] = log(z[j*nprobes + i])/log(2.0);
    }
  }


  for (iter = 1; iter <= maxiter; iter++){
    get_row_median(z,rdelta,nprobes,cols);
    subtract_by_row(z,rdelta,nprobes,cols);
    rmod(r,rdelta,nprobes);
    delta = median(c,cols);
    for (j = 0; j < cols; j++){
      c[j] = c[j] - delta;
    }
    t = t + delta;
    get_col_median(z,cdelta,nprobes,cols);
    subtract_by_col(z,cdelta,nprobes,cols);
    cmod(c,cdelta,cols);
    delta = median(r,nprobes);
    for (i =0; i < nprobes; i ++){
      r[i] = r[i] - delta;
    }
    t = t+delta;
    newsum = sum_abs(z,nprobes,cols);
    if (newsum == 0.0 || fabs(1.0 - oldsum/newsum) < eps)
      break;
    oldsum = newsum;
  }

  for (j=0; j < cols; j++){
    results[j] =  t + c[j];
  }

  Free(rdelta);
  Free(cdelta);
  Free(r);
  Free(c);
  Free(z);
}





/************************************************************************************
 **
 **  void do_RMA(double *PM, char **ProbeNames, int *rows, int * cols)
 **
 ** double *PM - matrix of dimension rows by cols (probes by chips) should already be
 **              normalized and background corrected.
 ** char **ProbeNames - Probeset names, one for each probe.
 ** int *rows, *cols - dimensions of matrix
 **
 ** perform the multichip averaging. PM should be background corrected and normalized
 **
 ** assumed that Probes are sorted, by ProbeNames, so that we can just look at
 ** consecutive rows in PM matrix when doing the median polish
 **
 ** each item is then used to create a matrix that is median polished to give
 ** expression estimates.
 **
 ************************************************************************************/

void do_RMA_buffmat(doubleBufferedMatrix Matrix, char **ProbeNames, int *rows, int *cols, double *results, char **outNames, int nps){
  int j = 0;
  int i = 0;
  int k = 0;
  int size;
  char *first;
  int first_ind;
  int max_nrows = 1000;


  /* buffers of size 200 should be enough. */

  int *cur_rows=Calloc(max_nrows,int);
  int nprobes=0;

  double *cur_exprs = Calloc(*cols,double);

  /* double *OLDPM = NULL; */

  first = ProbeNames[0];
  first_ind = 0;
  i = 0;     /* indexes current probeset */
  j = 0;    /* indexes current row in PM matrix */
  k = 0;    /* indexes current probe in probeset */
  while ( j < *rows){
    if (strcmp(first,ProbeNames[j]) == 0){
      if (k >= max_nrows){
        max_nrows = 2*max_nrows;
        cur_rows = Realloc(cur_rows, max_nrows, int);
      }
      cur_rows[k] = j;
      k++;
      j++;

    } else {
      nprobes = k;
      median_polish(Matrix, *rows, *cols, cur_rows, cur_exprs, nprobes);
      for (k =0; k < *cols; k++){
        results[k*nps + i] = cur_exprs[k];
      }
      size = strlen(first);
      outNames[i] = Calloc(size+1,char);
      strcpy(outNames[i],first);
      i++;
      first = ProbeNames[j];
      k = 0;
    }
  }
  nprobes = k;
  median_polish(Matrix, *rows, *cols, cur_rows, cur_exprs, nprobes);
  for (k =0; k < *cols; k++){
    results[k*nps + i] = cur_exprs[k];
  }
  size = strlen(first);
  outNames[i] = Calloc(size+1,char);
  strcpy(outNames[i],first);


  Free(cur_exprs);
  Free(cur_rows);
}
















SEXP R_bm_summarize_medianpolish(SEXP R_BufferedMatrix, SEXP N_probes, SEXP ProbeNamesVec){



  doubleBufferedMatrix Matrix;
  int rows, cols;
  double *outexpr;
  char **outnames;
  char **ProbeNames;
  int i,nprobesets;

  SEXP outvec; /* ,outnamesvec; */
  SEXP dimnames,names;

  SEXP temp;


  Matrix =  R_ExternalPtrAddr(R_BufferedMatrix);

  if (Matrix == NULL){
    return R_BufferedMatrix;
  }


 
  rows = dbm_getRows(Matrix);
  cols = dbm_getCols(Matrix);

  nprobesets=INTEGER(N_probes)[0];

  ProbeNames = Calloc(rows,char *);

  for (i =0; i < rows; i++)
    ProbeNames[i] = CHAR(STRING_ELT(ProbeNamesVec,i));


  outnames = Calloc(nprobesets,char *);

  /* PROTECT(outvec = NEW_NUMERIC(nprobesets*cols)); */

  PROTECT(outvec = allocMatrix(REALSXP, nprobesets, cols));


  outexpr = NUMERIC_POINTER(outvec);

  /* printf("Calculating Expression\n"); */
  Rprintf("Calculating Expression\n");


  do_RMA_buffmat(Matrix, ProbeNames, &rows, &cols,outexpr,outnames,nprobesets);



  /* now lets put names on the matrix */

  PROTECT(dimnames = allocVector(VECSXP,2));
  PROTECT(names = allocVector(STRSXP,nprobesets));

  for ( i =0; i < nprobesets; i++){
    PROTECT(temp = mkChar(outnames[i]));
    SET_VECTOR_ELT(names,i,temp); /* was a direct mkChar prior to Sep 2, 2005*/
    UNPROTECT(1);
  }
  SET_VECTOR_ELT(dimnames,0,names);
  setAttrib(outvec, R_DimNamesSymbol, dimnames);
  UNPROTECT(2);
  for (i =0; i < nprobesets; i++)
    Free(outnames[i]);

  Free(outnames);
  Free(ProbeNames);
  UNPROTECT(1);
  return outvec;
  



}
