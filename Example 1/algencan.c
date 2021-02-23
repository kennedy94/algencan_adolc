#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

void c_algencan(void *myevalf, void *myevalg, void *myevalh, void *myevalc,
	void *myevaljac, void *myevalhc, void *myevalfc, void *myevalgjac,
	void *myevalgjacp, void *myevalhl, void *myevalhlp, int jcnnzmax,
       	int hnnzmax,double *epsfeas, double *epsopt, double *efstin,
	double *eostin, double *efacc, double *eoacc, char *outputfnm,
	char *specfnm, int nvparam,char **vparam, int n, double *x,
	double *l, double *u, int m, double *lambda, _Bool *equatn,
	_Bool *linear, _Bool *coded, _Bool checkder, double *f,
	double *cnorm, double *snorm, double *nlpsupn,int *inform);

	void myevalf(int n, double *x, double *f, int *flag);
	void myevalg(int n, double *x, double *g, int *flag);
	void myevalh(int n, double *x, int *hrow, int *hcol, double *hval, int *hnnz,
		      int lim, _Bool *lmem, int *flag);
	void myevalc(int n, double *x, int ind, double *c, int *flag);
	void myevaljac(int n, double *x, int ind, int *jcvar, double *jcval,
		      int *jcnnz, int lim, _Bool *lmem, int *flag);
	void myevalhc(int n, double *x, int ind, int *hcrow, int *hccol, double *hcval,
					int *hcnnz, int lim, _Bool *lmem, int *flag);
	void myevalfc(int n, double *x, double *f, int m, double *c, int *flag);
	void myevalgjac(int n, double *x, double *g, int m, int *jcfun, int *jcvar,
							double *jcval, int *jcnnz, int lim, _Bool *lmem, int *flag);
	void myevalgjacp(int n, double *x, double *g, int m, double *p, double *q,
									char work, _Bool *gotj, int *flag);
	void myevalhl(int n, double *x, int m, double *lambda, double scalef,
										     double *scalec, int *hlrow, int *hlcol, double *hlval,
										     int *hlnnz, int lim, _Bool *lmem, int *flag);
	void myevalhlp(int n, double *x, int m, double *lambda, double scalef,
														      double *scalec, double *p, double *hp, _Bool *goth,
														      int *flag);

	int ndim,nite;
	double r,objrad;

/* ******************************************************************
   ****************************************************************** */

int main() {
  _Bool  checkder;
  int    hnnzmax,hnnzmax1,hnnzmax2,hnnzmax3,i,j,jcnnzmax,inform,m,n,nvparam,ncomp;
  double cnorm,efacc,efstin,eoacc,eostin,epsfeas,epsopt,f,nlpsupn,snorm;

  char   *specfnm, *outputfnm, **vparam;
  _Bool  coded[11],*equatn,*linear;
  double *l,*lambda,*u,*x;

  ndim   = 3;
  nite   = 100;
  r      = 1.0;
  objrad = 15.0;

  printf( "Density = %lf\n", nite * pow( ( r / objrad ), 3 ) );

  n = ndim * nite;
  m = nite;

  /* Memory allocation */
  x      = (double *) malloc(n * sizeof(double));
  l      = (double *) malloc(n * sizeof(double));
  u      = (double *) malloc(n * sizeof(double));
  lambda = (double *) malloc(m * sizeof(double));
  equatn = (_Bool  *) malloc(m * sizeof(_Bool ));
  linear = (_Bool  *) malloc(m * sizeof(_Bool ));

  if (     x == NULL ||      l == NULL ||      u == NULL ||
      lambda == NULL || equatn == NULL || linear == NULL ) {

    printf( "\nC ERROR IN MAIN PROGRAM: It was not possible to allocate memory.\n" );
    exit( 0 );

  }


  /* Initial point */
  srand (time(NULL));
  for(i = 0; i < n; i++) x[i] = -objrad + 2.0 * objrad * rand() / RAND_MAX;
  /* Lower and upper bounds */
  for(i=0;i<n;i++) {
    l[i] = - 1.0e20;
    u[i] =   1.0e20;
  }
  /* For each constraint i, set equatn[i] = 1. if it is an equality
     constraint of the form c_i(x) = 0, and set equatn[i] = 0 if it is
     an inequality constraint of the form c_i(x) <= 0. */
  for(i=0;i<m;i++) equatn[i] = 0;


  /* For each constraint i, set linear[i] = 1 if it is a linear
     constraint, otherwise set linear[i] = 0 */
  for(i=0;i<m;i++) linear[i] = 0;


  /* Lagrange multipliers approximation. */
  for( i = 0; i < m; i++ ) lambda[i] = 0.0;

  /* In this C interface evalf, evalg, evalh, evalc, evaljac and
     evalhc are present. evalfc, evalgjac, evalhl and evalhlp are
     not. */

  coded[0]  = 0; /* fsub     */
  coded[1]  = 0; /* gsub     */
  coded[2]  = 0; /* hsub     */
  coded[3]  = 0; /* csub     */
  coded[4]  = 0; /* jacsub   */
  coded[5]  = 0; /* hcsub    */
  coded[6]  = 1; /* fcsub    */
  coded[7]  = 1; /* gjacsub  */
  coded[8]  = 0; /* gjacpsub */
  coded[9]  = 1; /* hlsub    */
  coded[10] = 0; /* hlpsub   */

  /* Upper bounds on the number of sparse-matrices non-null
     elements */
  jcnnzmax = nite * ndim;
	hnnzmax1 = nite * ndim * ( nite * ndim + 1 ) / 2;
  hnnzmax2 = nite * ndim;
  hnnzmax3 = nite * ( ndim * ( ndim + 1 ) / 2 );
  hnnzmax  = hnnzmax1 + hnnzmax2 + hnnzmax3;

  /* Check derivatives? */
  checkder = 1;

  /* Parameters setting */
  epsfeas = 1.0e-08;
  epsopt  = 1.0e-08;
  efstin  = sqrt( epsfeas );
  eostin  = pow( epsopt, 1.5 );
  efacc   = sqrt( epsfeas );
  eoacc   = sqrt( epsopt );

  outputfnm = "algencan.out";
  specfnm   = "";

  nvparam = 0;


  /* Allocates VPARAM array */
  // vparam = ( char ** ) malloc( nvparam * sizeof( char * ) );

  /* Set algencan parameters */
  // vparam[0] = "TRUNCATED-NEWTON-LINE-SEARCH-INNER-SOLVER";
  // vparam[1] = "SKIP-ACCELERATION-STEP";
  // vparam[2] = "ITERATIONS-OUTPUT-DETAIL 11";

  // /* Optimize */
 c_algencan(&myevalf,&myevalg,&myevalh,&myevalc,&myevaljac,&myevalhc,&myevalfc,
      &myevalgjac,&myevalgjacp,&myevalhl,&myevalhlp,jcnnzmax,hnnzmax,
      &epsfeas,&epsopt,&efstin,&eostin,&efacc,&eoacc,outputfnm,specfnm,
      nvparam,vparam,n,x,l,u,m,lambda,equatn,linear,coded,checkder,
      &f,&cnorm,&snorm,&nlpsupn,&inform);

  /* Memory deallocation */
  free(x     );
  free(l     );
  free(u     );
  free(lambda);
  free(equatn);
  free(linear);
}
