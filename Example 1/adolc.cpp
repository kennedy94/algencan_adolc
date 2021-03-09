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
		double *scalec, double *p, double *hp, bool *goth, int *flag);
}

extern int ndim,nite;
extern double r,objrad;

/* ******************************************************************
****************************************************************** */

void myevalfad(int n, adouble *x, adouble *f, int *flag) {
	int i,j,k,ind1,ind2;
	adouble sum,aux,tmp;

	*flag = 0;

	*f = 0.0;
	for(i=0;i<nite;i++) {
		ind1 = ndim * i;

		for(j=i+1;j<nite;j++) {
			ind2 = ndim * j;
			sum = 0.0;
			for(k=0;k<ndim;k++) {
				sum = sum + pow( x[ind1+k] - x[ind2+k], 2 );
			}
			aux = 0.0;
			tmp = pow( 2.0 * r, 2 ) - sum;
			if (pow( 2.0 * r, 2 ) - sum > 0.0)
				aux = pow( 2.0 * r, 2 ) - sum;
			*f = *f + pow( aux, 2);
		}
	}
}

/* ******************************************************************
****************************************************************** */

void myevalcad(int n, adouble *x, int m, adouble *c, int *flag) {

	int i,k,ini;
	adouble sum;

	*flag = 0;

	for(i=0;i<m;i++) {
		ini = ndim * i;
		sum = 0.0;
		for(k=0;k<ndim;k++) {
			sum = sum + pow( x[ini+k], 2 );
		}

		c[i] = sum - pow( objrad - r, 2 );
	}

}

/* ******************************************************************
	****************************************************************** */

void myevalf(int n, double *x, double *f, int *flag) {

	*flag = -1;
}

/* ******************************************************************
	****************************************************************** */

void myevalg(int n, double *x, double *g, int *flag) {

	*flag = -1;
}

/* ******************************************************************
	****************************************************************** */

void myevalh(int n, double *x, int *hrow, int *hcol, double *hval, int *hnnz, int lim, bool *lmem, int *flag) {

	*flag = -1;
}

/* ******************************************************************
	****************************************************************** */

void myevalc(int n, double *x, int ind, double *c, int *flag) {

	*flag = -1;
}

/* ******************************************************************
	****************************************************************** */

void myevaljac(int n, double *x, int ind, int *jcvar, double *jcval, int *jcnnz, int lim, bool *lmem, int *flag) {
	*flag = -1;
}

/* ******************************************************************
	****************************************************************** */

void myevalhc(int n, double *x, int ind, int *hcrow, int *hccol, double *hcval, int *hcnnz, int lim, bool *lmem, int *flag) {

	*flag = -1;
}

/* *****************************************************************
	***************************************************************** */

void myevalfc(int n, double *x, double *f, int m, double *c, int *flag) {
	int i,j,k,ind1,ind2,ini;
	double sum,aux,tmp;

	*f = 0.0;
	for(i=0;i<nite;i++) {

		ind1 = ndim * i;

		for(j=i+1;j<nite;j++) {

			ind2 = ndim * j;
			sum = 0.0;
			for(k=0;k<ndim;k++) {
				sum = sum + pow( x[ind1+k] - x[ind2+k], 2 );
			}
			aux = 0.0;
			tmp = pow( 2.0 * r, 2 ) - sum;

			if (pow( 2.0 * r, 2 ) - sum > 0.0)
				aux = pow( 2.0 * r, 2 ) - sum;

			*f = *f + pow( aux, 2);
		}
	}

	for(i=0;i<m;i++) {
		ini = ndim * i;
		sum = 0.0;
		for(k=0;k<ndim;k++) {
			sum = sum + pow( x[ini+k], 2 );
		}
		c[i] = sum - pow( objrad - r, 2 );
	}
	*flag = 0;
}

/* *****************************************************************
	***************************************************************** */

void myevalgjac(int n, double *x, double *g, int m, int *jcfun, int *jcvar, double *jcval, int *jcnnz, int lim, bool *lmem, int *flag) {

	int i,nnz,options[4];
	double f, *c = new double[m], *values = NULL;
	adouble fad, *xad = new adouble[n], *cad = new adouble[m];
	unsigned int *rind  = NULL;
	unsigned int *cind  = NULL;

	trace_on(1);

	for (i=0;i<n;i++)
		xad[i] <<= x[i];

	myevalfad(n,xad,&fad,flag);

	fad >>= f;

	trace_off();

	gradient(1,n,x,g);

	options[0] = 0;          /* sparsity pattern by index domains (default) */
	options[1] = 0;          /*                         safe mode (default) */
	options[2] = 0;          /*              not required if options[0] = 0 */
	options[3] = 0;          /*                column compression (default) */

	trace_on(2);
	for(i=0;i<n;i++)
		xad[i] <<= x[i];

	myevalcad(n, xad, m, cad, flag);

	for(i=0;i<m;i++)
		cad[i] >>= c[i];
	trace_off();

	sparse_jac(2, m, n, 0, x, &nnz, &rind, &cind, &values, options);

	if( nnz > lim) {
		*lmem = 1;
		return;
	}
	*lmem = 0;

	for (i=0;i<nnz;i++) {
		jcval[i] = values[i];
		jcfun[i] = rind[i];
		jcvar[i] = cind[i];
	}

	*jcnnz = nnz;

	*flag = 0;

	free(rind);
	free(cind);
	free(values);
}

/* *****************************************************************
	***************************************************************** */

void myevalgjacp(int n, double *x, double *g, int m, double *p, double *q, char work, bool *gotj, int *flag) {
	*flag = -1;
}

/* *****************************************************************
	***************************************************************** */

void myevalhl(int n, double *x, int m, double *lambda, double scalef, double *scalec, int *hlrow, int *hlcol, double *hlval, int *hlnnz, int lim, bool *lmem, int *flag) {

	int i,nnz,options[2];
	double l, *values = NULL;
	adouble lad, *cad = new adouble[m], *xad = new adouble[n];
	unsigned int    *rind  = NULL;
	unsigned int    *cind  = NULL;

	options[0] = 0;          /*                               safe mode (default) */
	options[1] = 0;          /*                       indirect recovery (default) */

	trace_on(3);
	for (i=0;i<n;i++)
		xad[i] <<= x[i];

	myevalfad (n, xad, &lad, flag);

	lad = lad * scalef ;

	myevalcad(n, xad, m, cad, flag);

		for (i=0;i<m;i++)
			lad = lad + scalec[i] * lambda[i] * cad[i];

	lad >>= l;
	trace_off();

	sparse_hess(3, n, 0, x, &nnz, &rind, &cind, &values, options);

	if( nnz > lim) {
		*lmem = 1;
		return;
	}
	*lmem = 0;

	for (i=0;i<nnz;i++) {
		hlval[i] = values[i];
		hlcol[i] = rind[i];
		hlrow[i] = cind[i];
	}

	*hlnnz = nnz;

	*flag = 0;

	free(rind);
	free(cind);
	free(values);
}

/* *****************************************************************
***************************************************************** */

void myevalhlp(int n, double *x, int m, double *lambda, double scalef, double *scalec, double *p, double *hp, bool *goth, int *flag) {

	*flag = -1;
}

/* ******************************************************************
****************************************************************** */
