/* $Id: nnls.c,v 1.7 2000/11/07 16:29:30 tgkolda Exp $ */
/* $Source: /usr/local/cvsroot/appspack/apps/src/nnls.c,v $ */

/* Distributed with ASYNCHRONOUS PARALLEL PATTERN SEARCH (APPS) */

/* The routines in this file have been translated from Fortran to C by
   f2c. Additional modifications have been made to remove the
   dependencies on the f2c header file and library. The original
   Fortran 77 code accompanies the SIAM Publications printing of
   "Solving Least Squares Problems," by C. Lawson and R. Hanson and is
   freely available at www.netlib.org/lawson-hanson/all. */

/* nnls.F -- translated by f2c (version 19970805).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

/* The next line was removed after the f2c translation */
/* #include "f2c.h" */

/* The next lines were added after the f2c translation. Also swapped
   abs for nnls_abs and max for nnls_max to avoid confusion with some
   compilers. */
#include <stdio.h>
#include <math.h>
#define nnls_max(a,b) ((a) >= (b) ? (a) : (b))
#define nnls_abs(x) ((x) >= 0 ? (x) : -(x))
typedef int integer;
typedef double doublereal;

/* The following subroutine was added after the f2c translation */
double d_sign(double *a, double *b)
{
  double x;
  x = (*a >= 0 ? *a : - *a);
  return( *b >= 0 ? x : -x);
}

/* Table of constant values */

static integer c__1 = 1;
static integer c__0 = 0;
static integer c__2 = 2;

int h12_(integer *mode, integer *lpivot, integer *l1, const integer *m, doublereal *u, integer *iue, doublereal *up, doublereal *c__, integer *ice, const integer *icv, integer *ncv);
doublereal diff_(doublereal *x, doublereal *y);
int g1_(doublereal *a, doublereal *b, doublereal *cterm, doublereal *sterm, doublereal *sig);
int nnls_(doublereal *a, const integer *mda, const integer *m, const integer *n, doublereal *b, doublereal *x, doublereal *rnorm, doublereal *w, doublereal *zz, integer *index, integer *mode);


/*     SUBROUTINE NNLS  (A,MDA,M,N,B,X,RNORM,W,ZZ,INDEX,MODE) */

/*  Algorithm NNLS: NONNEGATIVE LEAST SQUARES */

/*  The original version of this code was developed by */
/*  Charles L. Lawson and Richard J. Hanson at Jet Propulsion Laboratory */
/*  1973 JUN 15, and published in the book */
/*  "SOLVING LEAST SQUARES PROBLEMS", Prentice-HalL, 1974. */
/*  Revised FEB 1995 to accompany reprinting of the book by SIAM. */

/*     GIVEN AN M BY N MATRIX, A, AND AN M-VECTOR, B,  COMPUTE AN */
/*     N-VECTOR, X, THAT SOLVES THE LEAST SQUARES PROBLEM */

/*                      A * X = B  SUBJECT TO X .GE. 0 */
/*     ------------------------------------------------------------------ */
/*                     Subroutine Arguments */

/*     A(),MDA,M,N     MDA IS THE FIRST DIMENSIONING PARAMETER FOR THE */
/*                     ARRAY, A().   ON ENTRY A() CONTAINS THE M BY N */
/*                     MATRIX, A.           ON EXIT A() CONTAINS */
/*                     THE PRODUCT MATRIX, Q*A , WHERE Q IS AN */
/*                     M BY M ORTHOGONAL MATRIX GENERATED IMPLICITLY BY */
/*                     THIS SUBROUTINE. */
/*     B()     ON ENTRY B() CONTAINS THE M-VECTOR, B.   ON EXIT B() CON- */
/*             TAINS Q*B. */
/*     X()     ON ENTRY X() NEED NOT BE INITIALIZED.  ON EXIT X() WILL */
/*             CONTAIN THE SOLUTION VECTOR. */
/*     RNORM   ON EXIT RNORM CONTAINS THE EUCLIDEAN NORM OF THE */
/*             RESIDUAL VECTOR. */
/*     W()     AN N-ARRAY OF WORKING SPACE.  ON EXIT W() WILL CONTAIN */
/*             THE DUAL SOLUTION VECTOR.   W WILL SATISFY W(I) = 0. */
/*             FOR ALL I IN SET P  AND W(I) .LE. 0. FOR ALL I IN SET Z */
/*     ZZ()     AN M-ARRAY OF WORKING SPACE. */
/*     INDEX()     AN INTEGER WORKING ARRAY OF LENGTH AT LEAST N. */
/*                 ON EXIT THE CONTENTS OF THIS ARRAY DEFINE THE SETS */
/*                 P AND Z AS FOLLOWS.. */

/*                 INDEX(1)   THRU INDEX(NSETP) = SET P. */
/*                 INDEX(IZ1) THRU INDEX(IZ2)   = SET Z. */
/*                 IZ1 = NSETP + 1 = NPP1 */
/*                 IZ2 = N */
/*     MODE    THIS IS A SUCCESS-FAILURE FLAG WITH THE FOLLOWING */
/*             MEANINGS. */
/*             1     THE SOLUTION HAS BEEN COMPUTED SUCCESSFULLY. */
/*             2     THE DIMENSIONS OF THE PROBLEM ARE BAD. */
/*                   EITHER M .LE. 0 OR N .LE. 0. */
/*             3    ITERATION COUNT EXCEEDED.  MORE THAN 3*N ITERATIONS. */

/*     ------------------------------------------------------------------ */
int nnls_(doublereal *a, const integer *mda, const integer *m, const integer *n, doublereal *b, doublereal *x, doublereal *rnorm, doublereal *w, doublereal *zz, integer *index, integer *mode)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;
    doublereal d__1, d__2;

    /* Builtin functions */
    /* The following lines were commented out after the f2c translation */
    /* double sqrt(); */
    /* integer s_wsfe(), do_fio(), e_wsfe(); */

    /* Local variables */
    static integer iter;
    static doublereal temp, wmax;
    static integer i__, j, l;
    static doublereal t, alpha, asave;
    static integer itmax, izmax, nsetp;
    static doublereal dummy, unorm, ztest, cc;
    static integer ii, jj, ip;
    static doublereal sm;
    static integer iz, jz;
    static doublereal up, ss;
    static integer rtnkey, iz1, iz2, npp1;

    /* Fortran I/O blocks */
    /* The following line was commented out after the f2c translation */
    /* static cilist io___22 = { 0, 6, 0, "(/a)", 0 }; */


/*     ------------------------------------------------------------------ 
*/
/*     integer INDEX(N) */
/*     double precision A(MDA,N), B(M), W(N), X(N), ZZ(M) */
/*     ------------------------------------------------------------------ 
*/
    /* Parameter adjustments */
    a_dim1 = *mda;
    a_offset = a_dim1 + 1;
    a -= a_offset;
    --b;
    --x;
    --w;
    --zz;
    --index;

    /* Function Body */
    *mode = 1;
    if (*m <= 0 || *n <= 0) {
	*mode = 2;
	return 0;
    }
    iter = 0;
    itmax = *n * 30;

/*                    INITIALIZE THE ARRAYS INDEX() AND X(). */

    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	x[i__] = 0.;
/* L20: */
	index[i__] = i__;
    }

    iz2 = *n;
    iz1 = 1;
    nsetp = 0;
    npp1 = 1;
/*                             ******  MAIN LOOP BEGINS HERE  ****** */
L30:
/*                  QUIT IF ALL COEFFICIENTS ARE ALREADY IN THE SOLUTION. 
*/
/*                        OR IF M COLS OF A HAVE BEEN TRIANGULARIZED. */

    if (iz1 > iz2 || nsetp >= *m) {
	goto L350;
    }

/*         COMPUTE COMPONENTS OF THE DUAL (NEGATIVE GRADIENT) VECTOR W(). 
*/

    i__1 = iz2;
    for (iz = iz1; iz <= i__1; ++iz) {
	j = index[iz];
	sm = 0.;
	i__2 = *m;
	for (l = npp1; l <= i__2; ++l) {
/* L40: */
	    sm += a[l + j * a_dim1] * b[l];
	}
	w[j] = sm;
/* L50: */
    }
/*                                   FIND LARGEST POSITIVE W(J). */
L60:
    wmax = 0.;
    i__1 = iz2;
    for (iz = iz1; iz <= i__1; ++iz) {
	j = index[iz];
	if (w[j] > wmax) {
	    wmax = w[j];
	    izmax = iz;
	}
/* L70: */
    }

/*             IF WMAX .LE. 0. GO TO TERMINATION. */
/*             THIS INDICATES SATISFACTION OF THE KUHN-TUCKER CONDITIONS. 
*/

    if (wmax <= 0.) {
	goto L350;
    }
    iz = izmax;
    j = index[iz];

/*     THE SIGN OF W(J) IS OK FOR J TO BE MOVED TO SET P. */
/*     BEGIN THE TRANSFORMATION AND CHECK NEW DIAGONAL ELEMENT TO AVOID */
/*     NEAR LINEAR DEPENDENCE. */

    asave = a[npp1 + j * a_dim1];
    i__1 = npp1 + 1;
    h12_(&c__1, &npp1, &i__1, m, &a[j * a_dim1 + 1], &c__1, &up, &dummy, &
	    c__1, &c__1, &c__0);
    unorm = 0.;
    if (nsetp != 0) {
	i__1 = nsetp;
	for (l = 1; l <= i__1; ++l) {
/* L90: */
/* Computing 2nd power */
	    d__1 = a[l + j * a_dim1];
	    unorm += d__1 * d__1;
	}
    }
    unorm = sqrt(unorm);
    d__2 = unorm + (d__1 = a[npp1 + j * a_dim1], nnls_abs(d__1)) * .01;
    if (diff_(&d__2, &unorm) > 0.) {

/*        COL J IS SUFFICIENTLY INDEPENDENT.  COPY B INTO ZZ, UPDATE Z
Z */
/*        AND SOLVE FOR ZTEST ( = PROPOSED NEW VALUE FOR X(J) ). */

	i__1 = *m;
	for (l = 1; l <= i__1; ++l) {
/* L120: */
	    zz[l] = b[l];
	}
	i__1 = npp1 + 1;
	h12_(&c__2, &npp1, &i__1, m, &a[j * a_dim1 + 1], &c__1, &up, &zz[1], &
		c__1, &c__1, &c__1);
	ztest = zz[npp1] / a[npp1 + j * a_dim1];

/*                                     SEE IF ZTEST IS POSITIVE */

	if (ztest > 0.) {
	    goto L140;
	}
    }

/*     REJECT J AS A CANDIDATE TO BE MOVED FROM SET Z TO SET P. */
/*     RESTORE A(NPP1,J), SET W(J)=0., AND LOOP BACK TO TEST DUAL */
/*     COEFFS AGAIN. */

    a[npp1 + j * a_dim1] = asave;
    w[j] = 0.;
    goto L60;

/*     THE INDEX  J=INDEX(IZ)  HAS BEEN SELECTED TO BE MOVED FROM */
/*     SET Z TO SET P.    UPDATE B,  UPDATE INDICES,  APPLY HOUSEHOLDER */
/*     TRANSFORMATIONS TO COLS IN NEW SET Z,  ZERO SUBDIAGONAL ELTS IN */
/*     COL J,  SET W(J)=0. */

L140:
    i__1 = *m;
    for (l = 1; l <= i__1; ++l) {
/* L150: */
	b[l] = zz[l];
    }

    index[iz] = index[iz1];
    index[iz1] = j;
    ++iz1;
    nsetp = npp1;
    ++npp1;

    if (iz1 <= iz2) {
	i__1 = iz2;
	for (jz = iz1; jz <= i__1; ++jz) {
	    jj = index[jz];
	    h12_(&c__2, &nsetp, &npp1, m, &a[j * a_dim1 + 1], &c__1, &up, &a[
		    jj * a_dim1 + 1], &c__1, mda, &c__1);
/* L160: */
	}
    }

    if (nsetp != *m) {
	i__1 = *m;
	for (l = npp1; l <= i__1; ++l) {
/* L180: */
	    a[l + j * a_dim1] = 0.;
	}
    }

    w[j] = 0.;
/*                                SOLVE THE TRIANGULAR SYSTEM. */
/*                                STORE THE SOLUTION TEMPORARILY IN ZZ(). 
*/
    rtnkey = 1;
    goto L400;
L200:

/*                       ******  SECONDARY LOOP BEGINS HERE ****** */

/*                          ITERATION COUNTER. */

L210:
    ++iter;
    if (iter > itmax) {
	*mode = 3;
	/* The following lines were replaced after the f2c translation */
	/* s_wsfe(&io___22); */
	/* do_fio(&c__1, " NNLS quitting on iteration count.", 34L); */
	/* e_wsfe(); */
	fprintf(stdout, "\n NNLS quitting on iteration count.\n");
	fflush(stdout);
	goto L350;
    }

/*                    SEE IF ALL NEW CONSTRAINED COEFFS ARE FEASIBLE. */
/*                                  IF NOT COMPUTE ALPHA. */

    alpha = 2.;
    i__1 = nsetp;
    for (ip = 1; ip <= i__1; ++ip) {
	l = index[ip];
	if (zz[ip] <= 0.) {
	    t = -x[l] / (zz[ip] - x[l]);
	    if (alpha > t) {
		alpha = t;
		jj = ip;
	    }
	}
/* L240: */
    }

/*          IF ALL NEW CONSTRAINED COEFFS ARE FEASIBLE THEN ALPHA WILL */
/*          STILL = 2.    IF SO EXIT FROM SECONDARY LOOP TO MAIN LOOP. */

    if (alpha == 2.) {
	goto L330;
    }

/*          OTHERWISE USE ALPHA WHICH WILL BE BETWEEN 0. AND 1. TO */
/*          INTERPOLATE BETWEEN THE OLD X AND THE NEW ZZ. */

    i__1 = nsetp;
    for (ip = 1; ip <= i__1; ++ip) {
	l = index[ip];
	x[l] += alpha * (zz[ip] - x[l]);
/* L250: */
    }

/*        MODIFY A AND B AND THE INDEX ARRAYS TO MOVE COEFFICIENT I */
/*        FROM SET P TO SET Z. */

    i__ = index[jj];
L260:
    x[i__] = 0.;

    if (jj != nsetp) {
	++jj;
	i__1 = nsetp;
	for (j = jj; j <= i__1; ++j) {
	    ii = index[j];
	    index[j - 1] = ii;
	    g1_(&a[j - 1 + ii * a_dim1], &a[j + ii * a_dim1], &cc, &ss, &a[j 
		    - 1 + ii * a_dim1]);
	    a[j + ii * a_dim1] = 0.;
	    i__2 = *n;
	    for (l = 1; l <= i__2; ++l) {
		if (l != ii) {

/*                 Apply procedure G2 (CC,SS,A(J-1,L),A(J,
L)) */

		    temp = a[j - 1 + l * a_dim1];
		    a[j - 1 + l * a_dim1] = cc * temp + ss * a[j + l * a_dim1]
			    ;
		    a[j + l * a_dim1] = -ss * temp + cc * a[j + l * a_dim1];
		}
/* L270: */
	    }

/*                 Apply procedure G2 (CC,SS,B(J-1),B(J)) */

	    temp = b[j - 1];
	    b[j - 1] = cc * temp + ss * b[j];
	    b[j] = -ss * temp + cc * b[j];
/* L280: */
	}
    }

    npp1 = nsetp;
    --nsetp;
    --iz1;
    index[iz1] = i__;

/*        SEE IF THE REMAINING COEFFS IN SET P ARE FEASIBLE.  THEY SHOULD 
*/
/*        BE BECAUSE OF THE WAY ALPHA WAS DETERMINED. */
/*        IF ANY ARE INFEASIBLE IT IS DUE TO ROUND-OFF ERROR.  ANY */
/*        THAT ARE NONPOSITIVE WILL BE SET TO ZERO */
/*        AND MOVED FROM SET P TO SET Z. */

    i__1 = nsetp;
    for (jj = 1; jj <= i__1; ++jj) {
	i__ = index[jj];
	if (x[i__] <= 0.) {
	    goto L260;
	}
/* L300: */
    }

/*         COPY B( ) INTO ZZ( ).  THEN SOLVE AGAIN AND LOOP BACK. */

    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L310: */
	zz[i__] = b[i__];
    }
    rtnkey = 2;
    goto L400;
L320:
    goto L210;
/*                      ******  END OF SECONDARY LOOP  ****** */

L330:
    i__1 = nsetp;
    for (ip = 1; ip <= i__1; ++ip) {
	i__ = index[ip];
/* L340: */
	x[i__] = zz[ip];
    }
/*        ALL NEW COEFFS ARE POSITIVE.  LOOP BACK TO BEGINNING. */
    goto L30;

/*                        ******  END OF MAIN LOOP  ****** */

/*                        COME TO HERE FOR TERMINATION. */
/*                     COMPUTE THE NORM OF THE FINAL RESIDUAL VECTOR. */

L350:
    sm = 0.;
    if (npp1 <= *m) {
	i__1 = *m;
	for (i__ = npp1; i__ <= i__1; ++i__) {
/* L360: */
/* Computing 2nd power */
	    d__1 = b[i__];
	    sm += d__1 * d__1;
	}
    } else {
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
/* L380: */
	    w[j] = 0.;
	}
    }
    *rnorm = sqrt(sm);
    return 0;

/*     THE FOLLOWING BLOCK OF CODE IS USED AS AN INTERNAL SUBROUTINE */
/*     TO SOLVE THE TRIANGULAR SYSTEM, PUTTING THE SOLUTION IN ZZ(). */

L400:
    i__1 = nsetp;
    for (l = 1; l <= i__1; ++l) {
	ip = nsetp + 1 - l;
	if (l != 1) {
	    i__2 = ip;
	    for (ii = 1; ii <= i__2; ++ii) {
		zz[ii] -= a[ii + jj * a_dim1] * zz[ip + 1];
/* L410: */
	    }
	}
	jj = index[ip];
	zz[ip] /= a[ip + jj * a_dim1];
/* L430: */
    }
    switch ((int)rtnkey) {
	case 1:  goto L200;
	case 2:  goto L320;
    }

    /* The next line was added after the f2c translation to keep
       compilers from complaining about a void return from a non-void
       function. */
    return 0;

} /* nnls_ */

int g1_(doublereal *a, doublereal *b, doublereal *cterm, doublereal *sterm, doublereal *sig)
{
    /* System generated locals */
    doublereal d__1;

    /* Builtin functions */
    /* The following line was commented out after the f2c translation */
    /* double sqrt(), d_sign(); */

    /* Local variables */
    static doublereal xr, yr;


/*     COMPUTE ORTHOGONAL ROTATION MATRIX.. */

/*  The original version of this code was developed by */
/*  Charles L. Lawson and Richard J. Hanson at Jet Propulsion Laboratory 
*/
/*  1973 JUN 12, and published in the book */
/*  "SOLVING LEAST SQUARES PROBLEMS", Prentice-HalL, 1974. */
/*  Revised FEB 1995 to accompany reprinting of the book by SIAM. */

/*     COMPUTE.. MATRIX   (C, S) SO THAT (C, S)(A) = (SQRT(A**2+B**2)) */
/*                        (-S,C)         (-S,C)(B)   (   0          ) */
/*     COMPUTE SIG = SQRT(A**2+B**2) */
/*        SIG IS COMPUTED LAST TO ALLOW FOR THE POSSIBILITY THAT */
/*        SIG MAY BE IN THE SAME LOCATION AS A OR B . */
/*     ------------------------------------------------------------------ 
*/
/*     ------------------------------------------------------------------ 
*/
    if (nnls_abs(*a) > nnls_abs(*b)) {
	xr = *b / *a;
/* Computing 2nd power */
	d__1 = xr;
	yr = sqrt(d__1 * d__1 + 1.);
	d__1 = 1. / yr;
	*cterm = d_sign(&d__1, a);
	*sterm = *cterm * xr;
	*sig = nnls_abs(*a) * yr;
	return 0;
    }
    if (*b != 0.) {
	xr = *a / *b;
/* Computing 2nd power */
	d__1 = xr;
	yr = sqrt(d__1 * d__1 + 1.);
	d__1 = 1. / yr;
	*sterm = d_sign(&d__1, b);
	*cterm = *sterm * xr;
	*sig = nnls_abs(*b) * yr;
	return 0;
    }
    *sig = 0.;
    *cterm = 0.;
    *sterm = 1.;
    return 0;
} /* g1_ */

/*     SUBROUTINE H12 (MODE,LPIVOT,L1,M,U,IUE,UP,C,ICE,ICV,NCV) */

/*  CONSTRUCTION AND/OR APPLICATION OF A SINGLE */
/*  HOUSEHOLDER TRANSFORMATION..     Q = I + U*(U**T)/B */

/*  The original version of this code was developed by */
/*  Charles L. Lawson and Richard J. Hanson at Jet Propulsion Laboratory */
/*  1973 JUN 12, and published in the book */
/*  "SOLVING LEAST SQUARES PROBLEMS", Prentice-HalL, 1974. */
/*  Revised FEB 1995 to accompany reprinting of the book by SIAM. */
/*     ------------------------------------------------------------------ */
/*                     Subroutine Arguments */

/*     MODE   = 1 OR 2   Selects Algorithm H1 to construct and apply a */
/*            Householder transformation, or Algorithm H2 to apply a */
/*            previously constructed transformation. */
/*     LPIVOT IS THE INDEX OF THE PIVOT ELEMENT. */
/*     L1,M   IF L1 .LE. M   THE TRANSFORMATION WILL BE CONSTRUCTED TO */
/*            ZERO ELEMENTS INDEXED FROM L1 THROUGH M.   IF L1 GT. M */
/*            THE SUBROUTINE DOES AN IDENTITY TRANSFORMATION. */
/*     U(),IUE,UP    On entry with MODE = 1, U() contains the pivot */
/*            vector.  IUE is the storage increment between elements. */
/*            On exit when MODE = 1, U() and UP contain quantities */
/*            defining the vector U of the Householder transformation. */
/*            on entry with MODE = 2, U() and UP should contain */
/*            quantities previously computed with MODE = 1.  These will */
/*            not be modified during the entry with MODE = 2. */
/*     C()    ON ENTRY with MODE = 1 or 2, C() CONTAINS A MATRIX WHICH */
/*            WILL BE REGARDED AS A SET OF VECTORS TO WHICH THE */
/*            HOUSEHOLDER TRANSFORMATION IS TO BE APPLIED. */
/*            ON EXIT C() CONTAINS THE SET OF TRANSFORMED VECTORS. */
/*     ICE    STORAGE INCREMENT BETWEEN ELEMENTS OF VECTORS IN C(). */
/*     ICV    STORAGE INCREMENT BETWEEN VECTORS IN C(). */
/*     NCV    NUMBER OF VECTORS IN C() TO BE TRANSFORMED. IF NCV .LE. 0 */
/*            NO OPERATIONS WILL BE DONE ON C(). */
/*     ------------------------------------------------------------------ */
int h12_(integer *mode, integer *lpivot, integer *l1, const integer *m, doublereal *u, integer *iue, doublereal *up, doublereal *c__, integer *ice, const integer *icv, integer *ncv)
{
    /* System generated locals */
    integer u_dim1, u_offset, i__1, i__2;
    doublereal d__1, d__2;

    /* Builtin functions */
    /* The following line was commented out after the f2c translation */
    /* double sqrt(); */

    /* Local variables */
    static integer incr;
    static doublereal b;
    static integer i__, j;
    static doublereal clinv;
    static integer i2, i3, i4;
    static doublereal cl, sm;

/*     ------------------------------------------------------------------ 
*/
/*     double precision U(IUE,M) */
/*     ------------------------------------------------------------------ 
*/
    /* Parameter adjustments */
    u_dim1 = *iue;
    u_offset = u_dim1 + 1;
    u -= u_offset;
    --c__;

    /* Function Body */
    if (0 >= *lpivot || *lpivot >= *l1 || *l1 > *m) {
	return 0;
    }
    cl = (d__1 = u[*lpivot * u_dim1 + 1], nnls_abs(d__1));
    if (*mode == 2) {
	goto L60;
    }
/*                            ****** CONSTRUCT THE TRANSFORMATION. ****** 
*/
    i__1 = *m;
    for (j = *l1; j <= i__1; ++j) {
/* L10: */
/* Computing MAX */
	d__2 = (d__1 = u[j * u_dim1 + 1], nnls_abs(d__1));
	cl = nnls_max(d__2,cl);
    }
    if (cl <= 0.) {
	goto L130;
    } else {
	goto L20;
    }
L20:
    clinv = 1. / cl;
/* Computing 2nd power */
    d__1 = u[*lpivot * u_dim1 + 1] * clinv;
    sm = d__1 * d__1;
    i__1 = *m;
    for (j = *l1; j <= i__1; ++j) {
/* L30: */
/* Computing 2nd power */
	d__1 = u[j * u_dim1 + 1] * clinv;
	sm += d__1 * d__1;
    }
    cl *= sqrt(sm);
    if (u[*lpivot * u_dim1 + 1] <= 0.) {
	goto L50;
    } else {
	goto L40;
    }
L40:
    cl = -cl;
L50:
    *up = u[*lpivot * u_dim1 + 1] - cl;
    u[*lpivot * u_dim1 + 1] = cl;
    goto L70;
/*            ****** APPLY THE TRANSFORMATION  I+U*(U**T)/B  TO C. ****** 
*/

L60:
    if (cl <= 0.) {
	goto L130;
    } else {
	goto L70;
    }
L70:
    if (*ncv <= 0) {
	return 0;
    }
    b = *up * u[*lpivot * u_dim1 + 1];
/*                       B  MUST BE NONPOSITIVE HERE.  IF B = 0., RETURN. 
*/

    if (b >= 0.) {
	goto L130;
    } else {
	goto L80;
    }
L80:
    b = 1. / b;
    i2 = 1 - *icv + *ice * (*lpivot - 1);
    incr = *ice * (*l1 - *lpivot);
    i__1 = *ncv;
    for (j = 1; j <= i__1; ++j) {
	i2 += *icv;
	i3 = i2 + incr;
	i4 = i3;
	sm = c__[i2] * *up;
	i__2 = *m;
	for (i__ = *l1; i__ <= i__2; ++i__) {
	    sm += c__[i3] * u[i__ * u_dim1 + 1];
/* L90: */
	    i3 += *ice;
	}
	if (sm != 0.) {
	    goto L100;
	} else {
	    goto L120;
	}
L100:
	sm *= b;
	c__[i2] += sm * *up;
	i__2 = *m;
	for (i__ = *l1; i__ <= i__2; ++i__) {
	    c__[i4] += sm * u[i__ * u_dim1 + 1];
/* L110: */
	    i4 += *ice;
	}
L120:
	;
    }
L130:
    return 0;
} /* h12_ */

doublereal diff_(doublereal *x, doublereal *y)
{
    /* System generated locals */
    doublereal ret_val;


/*  Function used in tests that depend on machine precision. */

/*  The original version of this code was developed by */
/*  Charles L. Lawson and Richard J. Hanson at Jet Propulsion Laboratory 
*/
/*  1973 JUN 7, and published in the book */
/*  "SOLVING LEAST SQUARES PROBLEMS", Prentice-HalL, 1974. */
/*  Revised FEB 1995 to accompany reprinting of the book by SIAM. */

    ret_val = *x - *y;
    return ret_val;
} /* diff_ */


/* The following subroutine was added after the f2c translation */
int nnls_c(double* a, const int* mda, const int* m, const int* n, double* b, 
	 double* x, double* rnorm, double* w, double* zz, int* index, 
	 int* mode)
{
  return (nnls_(a, mda, m, n, b, x, rnorm, w, zz, index, mode));
}

