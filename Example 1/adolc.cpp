#include <cmath>
#include <cstdlib>
#include <cstdio>

#include <adolc/adolc.h>
#include <adolc/adolc_sparse.h>

#define tag 1

extern "C" {
  void myevalf(int n, double *x, double *f, int *flag);
	void myevalg(int n, double *x, double *g, int *flag);
	void myevalh(int n, double *x, int *hrow, int *hcol, double *hval, int *hnnz,
		       int lim, bool *lmem, int *flag);
	void myevalc(int n, double *x, int ind, double *c, int *flag);
	void myevaljac(int n, double *x, int ind, int *jcvar, double *jcval,
		       int *jcnnz, int lim, bool *lmem, int *flag);
	void myevalhc(int n, double *x, int ind, int *hcrow, int *hccol, double *hcval,
				 	 int *hcnnz, int lim, bool *lmem, int *flag);
	void myevalfc(int n, double *x, double *f, int m, double *c, int *flag);
	void myevalgjac(int n, double *x, double *g, int m, int *jcfun, int *jcvar,
					 double *jcval, int *jcnnz, int lim, bool *lmem, int *flag);
	void myevalgjacp(int n, double *x, double *g, int m, double *p, double *q,
					 char work, bool *gotj, int *flag);
	void myevalhl(int n, double *x, int m, double *lambda, double scalef,
					 double *scalec, int *hlrow, int *hlcol, double *hlval,
				   int *hlnnz, int lim, bool *lmem, int *flag);
	void myevalhlp(int n, double *x, int m, double *lambda, double scalef,
					 double *scalec, double *p, double *hp, bool *goth,
				   int *flag);
}



void myevalf(int n, double *x, double *f, int *flag) {

   *flag = 0;

   *f = 4 * pow( x[0],2 ) + 10 * pow( x[1], 2 );
}

/* ******************************************************************
   ****************************************************************** */

void myevalfad(int n, adouble *x, adouble *f, int *flag) {

	*flag = 0;

   *f = 4 * pow( x[0],2 ) + 10 * pow( x[1], 2 );

}


/* ******************************************************************
   ****************************************************************** */

void myevalg(int n, double *x, double *g, int *flag) {

	double f;
	adouble *xad = new adouble[n];
	adouble fad;
	int i;

	trace_on(tag);
	for (i=0;i<n;i++)
		xad[i] <<= x[i];

		myevalfad(n,xad,&fad,flag);

	 fad >>= f;
	 trace_off();

	 gradient(tag,n,x,g);

   *flag = 0;

}

/* ******************************************************************
   ****************************************************************** */

void myevalh(int n, double *x, int *hrow, int *hcol, double *hval, int *hnnz,
	     int lim, bool *lmem, int *flag) {

		double f;
		adouble *xad = new adouble[n];
		adouble fad;
		unsigned int    *rind  = NULL;
    unsigned int    *cind  = NULL;
    double *values = NULL;
		int i,nnz,options[2];

    options[0] = 0;          /*                               safe mode (default) */
    options[1] = 0;          /*                       indirect recovery (default) */

		trace_on(tag);
		for (i=0;i<n;i++)
			xad[i] <<= x[i];

 		myevalfad(n,xad,&fad,flag);

 	 fad >>= f;
 	 trace_off();

	 sparse_hess(tag, n, 0, x, &nnz, &rind, &cind, &values, options);

	 if( nnz > lim) {
		 *lmem = 1;
		 return;
	 }
   *lmem = 0;

	 for (i=0;i<nnz;i++) {
		 hval[i] = values[i];
		 hcol[i] = rind[i];
		 hrow[i] = cind[i];
	 }

	 *hnnz = nnz;

   *flag = 0;

	 free(rind);
	 free(cind);
	 free(values);
}

/* ******************************************************************
   ****************************************************************** */

void myevalc(int n, double *x, int ind, double *c, int *flag) {

   *flag = 0;

   if ( ind == 0 )
     *c = - pow( x[0],2 ) - pow( x[1],2 );

   else if ( ind == 1 )
     *c = pow( x[0],2 ) + pow( x[1],2 ) - 4.0;

   else
     *flag = -1;
}

/* ******************************************************************
   ****************************************************************** */

void myevalcad(int n, adouble *x, int ind, adouble *c, int *flag) {

   *flag = 0;

   if ( ind == 0 )
     *c = - pow( x[0],2 ) - pow( x[1],2 );

   else if ( ind == 1 )
     *c = pow( x[0],2 ) + pow( x[1],2 ) - 4.0;

   else
     *flag = -1;
}

/* ******************************************************************
   ****************************************************************** */

void myevaljac(int n, double *x, int ind, int *jcvar, double *jcval,
	       int *jcnnz, int lim, bool *lmem, int *flag) {

   *flag = 0;
   *lmem = 0;
   *jcnnz = 2;


     if( *jcnnz > lim ) {
       *lmem = 1;
       return;
     }
		 if ( ind == 0 ) {
     jcvar[0] = 0;
     jcval[0] = -2.0 * x[0];
     jcvar[1] = 1;
     jcval[1] = -2.0 * x[1];
   }
   else if ( ind == 1 ) {
     jcvar[0] = 0;
     jcval[0] = 2.0 * x[0];
     jcvar[1] = 1;
     jcval[1] = 2.0 * x[1];
   }
   else
     *flag = -1;
}

/* ******************************************************************
   ****************************************************************** */

void myevalhc(int n, double *x, int ind, int *hcrow, int *hccol, double *hcval,
	      int *hcnnz, int lim, bool *lmem, int *flag) {

	 double c;
	 adouble *xad = new adouble[n];
	 adouble cad;
	 unsigned int    *rind  = NULL;
	 unsigned int    *cind  = NULL;
	 double *values = NULL;
	 int i,nnz,options[2];

		options[0] = 0;          /*                               safe mode (default) */
		options[1] = 0;          /*                       indirect recovery (default) */

		trace_on(tag);
		for (i=0;i<n;i++)
		  xad[i] <<= x[i];

	  myevalcad(n,xad,ind,&cad,flag);

		cad >>= c;
		trace_off();

		sparse_hess(tag, n, 0, x, &nnz, &rind, &cind, &values, options);


		if( nnz > lim) {
			*lmem = 1;
			return;
		}
    *lmem = 0;

		for (i=0;i<nnz;i++) {
			hcval[i] = values[i];
			hccol[i] = cind[i];
			hcrow[i] = rind[i];
			 	}

		*hcnnz = nnz;

		*flag = 0;

		free(rind);
		free(cind);
		free(values);
}

/* *****************************************************************
   ***************************************************************** */

void myevalfc(int n, double *x, double *f, int m, double *c, int *flag) {

   *flag = -1;
}

/* *****************************************************************
   ***************************************************************** */

void myevalgjac(int n, double *x, double *g, int m, int *jcfun, int *jcvar,
		double *jcval, int *jcnnz, int lim, bool *lmem, int *flag) {

   *flag = -1;
}

/* *****************************************************************
   ***************************************************************** */

void myevalgjacp(int n, double *x, double *g, int m, double *p, double *q,
		 char work, bool *gotj, int *flag) {

   *flag = -1;
}

/* *****************************************************************
   ***************************************************************** */

void myevalhl(int n, double *x, int m, double *lambda, double scalef,
	      double *scalec, int *hlrow, int *hlcol, double *hlval,
	      int *hlnnz, int lim, bool *lmem, int *flag) {

   *flag = -1;
}

/* *****************************************************************
   ***************************************************************** */

void myevalhlp(int n, double *x, int m, double *lambda, double scalef,
	       double *scalec, double *p, double *hp, bool *goth,
	       int *flag) {

   *flag = -1;
}

/* ******************************************************************
   ****************************************************************** */
