#include <cmath>
#include <cstdlib>
#include <cstdio>

#include <adolc/adolc.h>
#include <adolc/adolc_sparse.h>

#define tag 1

extern "C" {
	void myevalf(int n, double* x, double* f, int* flag);
	void myevalg(int n, double* x, double* g, int* flag);
	void myevalh(int n, double* x, int* hrow, int* hcol, double* hval, int* hnnz,
		int lim, bool* lmem, int* flag);
	void myevalc(int n, double* x, int ind, double* c, int* flag);
	void myevaljac(int n, double* x, int ind, int* jcvar, double* jcval,
		int* jcnnz, int lim, bool* lmem, int* flag);
	void myevalhc(int n, double* x, int ind, int* hcrow, int* hccol, double* hcval,
		int* hcnnz, int lim, bool* lmem, int* flag);
	void myevalfc(int n, double* x, double* f, int m, double* c, int* flag);
	void myevalgjac(int n, double* x, double* g, int m, int* jcfun, int* jcvar,
		double* jcval, int* jcnnz, int lim, bool* lmem, int* flag);
	void myevalgjacp(int n, double* x, double* g, int m, double* p, double* q,
		char work, bool* gotj, int* flag);
	void myevalhl(int n, double* x, int m, double* lambda, double scalef,
		double* scalec, int* hlrow, int* hlcol, double* hlval,
		int* hlnnz, int lim, bool* lmem, int* flag);
	void myevalhlp(int n, double* x, int m, double* lambda, double scalef,
		double* scalec, double* p, double* hp, bool* goth,
		int* flag);
}

extern int ndim, nite;
extern double r, objrad;

/* ******************************************************************
	 ****************************************************************** */

void myevalf(int n, double* x, double* f, int* flag) {
	int i, j, k, ind1, ind2;
	double sum;

	*flag = 0;

	*f = 0.0;
	for (i = 0;i < nite;i++) {
		ind1 = ndim * i;
		for (j = i + 1;j < nite;j++) {
			ind2 = ndim * j;
			sum = 0.0;
			for (k = 0;k < ndim;k++) {
				sum = sum + pow(x[ind1 + k] - x[ind2 + k], 2);
			}
			*f = *f + pow(std::max(0.0, pow(2.0 * r, 2) - sum), 2);
		}
	}
}

/* ******************************************************************
	 ****************************************************************** */

void myevalfad(int n, adouble* x, adouble* f, int* flag) {
	int i, j, k, ind1, ind2;
	adouble sum, aux;

	*flag = 0;

	*f = 0.0;
	for (i = 0;i < nite;i++) {
		ind1 = ndim * i;
		for (j = i + 1;j < nite;j++) {
			ind2 = ndim * j;
			sum = 0.0;
			for (k = 0;k < ndim;k++) {
				sum = sum + pow(x[ind1 + k] - x[ind2 + k], 2);
			}
			aux = 0.0;
			if (pow(2.0 * r, 2) - sum > 0.0)
				aux = pow(2.0 * r, 2) - sum;
			*f = *f + pow(aux, 2);
		}
	}
}

/* ******************************************************************
	 ****************************************************************** */

void myevalg(int n, double* x, double* g, int* flag) {

	double f;
	adouble* xad = new adouble[n];
	adouble fad;
	int i;

	trace_on(tag);
	for (i = 0;i < n;i++)
		xad[i] <<= x[i];

	myevalfad(n, xad, &fad, flag);

	fad >>= f;
	trace_off();

	gradient(tag, n, x, g);

	*flag = 0;
	delete[] xad;

}

/* ******************************************************************
	 ****************************************************************** */

void myevalh(int n, double* x, int* hrow, int* hcol, double* hval, int* hnnz,
	int lim, bool* lmem, int* flag) {

	double f;
	adouble* xad = new adouble[n];
	adouble fad;
	unsigned int* rind = NULL;
	unsigned int* cind = NULL;
	double* values = NULL;
	int i, nnz, options[2];

	options[0] = 0;          /*                               safe mode (default) */
	options[1] = 0;          /*                       indirect recovery (default) */

	trace_on(tag);
	for (i = 0;i < n;i++)
		xad[i] <<= x[i];

	myevalfad(n, xad, &fad, flag);

	fad >>= f;
	trace_off();

	sparse_hess(tag, n, 0, x, &nnz, &rind, &cind, &values, options);

	delete[] xad;

	if (nnz > lim) {
		*lmem = 1;
		return;
	}
	*lmem = 0;

	for (i = 0;i < nnz;i++) {
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

void myevalc(int n, double* x, int ind, double* c, int* flag) {

	int k, ini;
	double sum;

	*flag = 0;

	ini = ndim * ind;

	sum = 0.0;
	for (k = 0;k < ndim;k++) {
		sum = sum + pow(x[ini + k], 2);
	}
	*c = sum - pow(objrad - r, 2);
}

/* ******************************************************************
	 ****************************************************************** */

void myevalcad(int n, adouble* x, int ind, adouble* c, int* flag) {

	int k, ini;
	adouble sum;

	*flag = 0;

	ini = ndim * ind;

	sum = 0.0;
	for (k = 0;k < ndim;k++) {
		sum = sum + pow(x[ini + k], 2);
	}
	*c = sum - pow(objrad - r, 2);
}

/* ******************************************************************
	 ****************************************************************** */
void myevaljac(int n, double* x, int ind, int* jcvar, double* jcval,
	int* jcnnz, int lim, bool* lmem, int* flag) {

	int i, ini;

	*flag = 0;
	*lmem = 0;
	*jcnnz = ndim;

	if (*jcnnz > lim) {
		*lmem = 1;
		return;
	}

	ini = ndim * ind;

	for (i = 0;i < ndim;i++) {
		jcvar[i] = ini + i;
		jcval[i] = 2.0 * x[ini + i];
	}

}

/* ******************************************************************
	 ****************************************************************** */

void myevalhc(int n, double* x, int ind, int* hcrow, int* hccol, double* hcval,
	int* hcnnz, int lim, bool* lmem, int* flag) {

	double c;
	adouble* xad = new adouble[n];
	adouble cad;
	unsigned int* rind = NULL;
	unsigned int* cind = NULL;
	double* values = NULL;
	int i, nnz, options[2];

	options[0] = 0;          /*                               safe mode (default) */
	options[1] = 0;          /*                       indirect recovery (default) */

	trace_on(tag);
	for (i = 0;i < n;i++)
		xad[i] <<= x[i];

	myevalcad(n, xad, ind, &cad, flag);

	cad >>= c;
	trace_off();

	sparse_hess(tag, n, 0, x, &nnz, &rind, &cind, &values, options);

	delete[] xad;
	if (nnz > lim) {
		*lmem = 1;
		return;
	}
	*lmem = 0;

	for (i = 0;i < nnz;i++) {
		hcval[i] = values[i];
		hccol[i] = rind[i];
		hcrow[i] = cind[i];
	}

	*hcnnz = nnz;

	*flag = 0;

	free(rind);
	free(cind);
	free(values);
}

/* *****************************************************************
	 ***************************************************************** */

void myevalfc(int n, double* x, double* f, int m, double* c, int* flag) {

	*flag = -1;
}

/* *****************************************************************
	 ***************************************************************** */

void myevalgjac(int n, double* x, double* g, int m, int* jcfun, int* jcvar,
	double* jcval, int* jcnnz, int lim, bool* lmem, int* flag) {

	*flag = -1;
}

/* *****************************************************************
	 ***************************************************************** */

void myevalgjacp(int n, double* x, double* g, int m, double* p, double* q,
	char work, bool* gotj, int* flag) {

	*flag = -1;
}

/* *****************************************************************
	 ***************************************************************** */

void myevalhl(int n, double* x, int m, double* lambda, double scalef,
	double* scalec, int* hlrow, int* hlcol, double* hlval,
	int* hlnnz, int lim, bool* lmem, int* flag) {

	*flag = -1;
}

/* *****************************************************************
	 ***************************************************************** */

void myevalhlp(int n, double* x, int m, double* lambda, double scalef,
	double* scalec, double* p, double* hp, bool* goth,
	int* flag) {

	*flag = -1;
}
