/* based on golub.f -- translated by f2c (version 19970805).
	converted further by Einar Kjartansson, may 1999 	
	limit on maimum value of m removed by using calloc and free 
*/

#include <stdlib.h>
#include <math.h>

void  golubC(double *a, double *x, double *b, int m, int n)
{
    /* System generated locals */
    int a_dim1, a_offset, i__1, i__2, i__3;

    /* Local variables */
    double fact;
    int nmax;
    double dsum;
    int irow, i__, j, k;
    double dfact, *u;
    int kplus;
    double dsigma, dai, daj, dbi;
    int lim;
    double sum;


/* 	double precision version */
/* 	a(m,n) ; b(m) given wiht m >= n  solves for x(n) such that */
/* 	 || b - ax || = minium */
/* 	method of golub, from claerbout 1976 */

    /* Parameter adjustments */
    u = ( double *) calloc(m,sizeof(*u)) ;
    --b;
    --x;
    a_dim1 = m;
    a_offset = a_dim1 + 1;
    a -= a_offset;

    /* Function Body */
    nmax = n;
    if (m == n) {
	nmax = m - 1;
    }
    i__1 = nmax;
    for (k = 1; k <= i__1; ++k) {
	dsum = 0.;
	i__2 = m;
	for (i__ = k; i__ <= i__2; ++i__) {
	    daj = a[i__ + k * a_dim1];
	    dsum += daj * daj;
	}
	dai = a[k + k * a_dim1];
	dsigma = sqrt(dsum) ;
	if(dai < 0.0 ) dsigma = -dsigma ;

	dbi = sqrt(dai / dsigma + 1.);
	dfact = 1. / (dsigma * dbi);
	u[k - 1] = dbi;
	fact = dfact;
	kplus = k + 1;
	i__2 = m;
	for (i__ = kplus; i__ <= i__2; ++i__) {
	    u[i__ - 1] = fact * a[i__ + k * a_dim1];
	}

/* ......... I = UU' is a symmetric, orthogonal matrix which when appl
ied */
/* ......... to a(.,.) will annihilate the elements below the diagonal
 k */

	i__2 = n;
	for (j = k; j <= i__2; ++j) {
	    fact = 0.;
	    i__3 = m;
	    for (i__ = k; i__ <= i__3; ++i__) {
		fact += u[i__ - 1] * a[i__ + j * a_dim1];
	    }
	    i__3 = m;
	    for (i__ = k; i__ <= i__3; ++i__) {
		a[i__ + j * a_dim1] -= fact * u[i__ - 1];
	    }
	}
	fact = 0.;
	i__2 = m;
	for (i__ = k; i__ <= i__2; ++i__) {
	    fact += u[i__ - 1] * b[i__];
	}
	i__2 = m;
	for (i__ = k; i__ <= i__2; ++i__) {
	    b[i__] -= fact * u[i__ - 1];
	}
    }
/* ........Back substitute recursively */
    x[n] = b[n] / a[n + n * a_dim1];
    lim = n - 1;
    i__1 = lim;
    for (i__ = 1; i__ <= i__1; ++i__) {
	irow = n - i__;
	sum = 0.;
	i__2 = i__;
	for (j = 1; j <= i__2; ++j) {
	    sum += x[n - j + 1] * a[irow + (n - j + 1) * a_dim1];
	}
	x[irow] = (b[irow] - sum) / a[irow + irow * a_dim1];
    }
    free(u) ;
} /* golub_ */

