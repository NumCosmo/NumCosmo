
#include <math.h>
#include <stdio.h>
#include <time.h>

#include "math/ncm_lapack.h"
#include "math/ncm_flapack.h"

//#include <R.h>
//#include <R_ext/Lapack.h>
/*
#include "mkl_lapack.h"
#include "mkl_blas.h"
*/

#ifdef _MKL_LAPACK_H_
    #define CTOF(x) x
#else
    #define CTOF(x) x ## _
#endif

/*#define MAX(A,B) ( (A) > (B) ? (A):(B))*/
/*#define MIN(A,B) ( (A) < (B) ? (A):(B))*/

#define LU   1
#define CHOL 2
#define SMW  3
#define PFCF 4

#define PRED 1
#define CORR 2

/*****************************************************************************/

/* Global Variables */
#define EPSTERM   5.0E-11
#define EPSROUND  1.0E-8
#define EPSINIT   1.0E-2
#define EPSIPM    1.0E-2
#define EPSPERT   1.0E-14

/*****************************************************************************/

void PrintMatrix( char* name, double* vals, int* rows, int* cols )
{
    int i, j;
    printf("%8s = [\n", name);

    for (i=0;i<(*rows);i++)
    {
        for (j=0;j<(*cols);j++)
        {
            printf("%15.10e ", vals[i+j*(*rows)]);
        }
        printf(";\n");
    }
    printf("];\n");
}

/*****************************************************************************/

double VectorAbsSum( double* v, int *n )
{
    int one=1;
    return CTOF(dasum)( n, v, &one );
    /* return dasum( n, v, &one ); */
}

/*****************************************************************************/

void VectorVectorCopy( double* lhs, double* rhs, int* n )
{
    int one=1;
    CTOF(dcopy)( n, rhs, &one, lhs, &one );
    /* dcopy( n, rhs, &one, lhs, &one ); */
}

/*****************************************************************************/

void VectorVectorDivide( double* top, double* bot, double* res, int* n )
{
    int i;
    for (i=0;i<(*n);i++) res[i] = top[i]/bot[i];
}
/*****************************************************************************/

void VectorVectorMult( double* alpha, double* x, double* y, int* n )
{
    int one=1;
    CTOF(daxpy)( n, alpha, x, &one, y, &one );
    /* daxpy( n, alpha, x, &one, y, &one ); */
}

/*****************************************************************************/

void VectorVectorMinus( double* x, double* y, double* res, int* n )
{
    double alpha = -1.0;
    VectorVectorCopy( res, x, n );
    VectorVectorMult( &alpha, y, res, n );
}

/*****************************************************************************/

double VectorVectorDot( double* x, double* y, int* n )
{
    int one=1;
    return CTOF(ddot)( n, x, &one, y, &one );
    /* return ddot( n, x, &one, y, &one ); */
}

/*****************************************************************************/

void MatrixMatrixDiagSolve( double* D, double* A, double* res, int* rows, int*
cols )
{
    int i,j;
    int n = (*rows);
    int m = (*cols);
    double temp;
    for (i=0;i<n;i++)
    {
        temp = D[i];
        for (j=0;j<m;j++) res[i + j*n] = A[i + j*n] / temp;
    }
}

/*****************************************************************************/

void MatrixVectorMult( double* alpha, double* A, int trans, double* x, double*
beta, double* b, int* rows, int* cols )
{
    int one=1;
    if (trans)
    {
        CTOF(dgemv)("T", rows, cols, alpha, A, rows, x, &one, beta, b, &one );
        /* dgemv('T', *rows, *cols, *alpha, A, *rows, x, one, *beta, b, one ); */
    }
    else
    {
        CTOF(dgemv)("N", rows, cols, alpha, A, rows, x, &one, beta, b, &one );
        /* dgemv('N', *rows, *cols, *alpha, A, *rows, x, one, *beta, b, one ); */
    }
}

/*****************************************************************************/

void MatrixConstantPlusDiag( double* A, double c, int* n )
{
    int i;
    for (i=0;i<(*n);i++) A[i+i*(*n)] += c;
}

/*****************************************************************************/

void MatrixCholFactorize( double* A, int* n, int* info )
{
    CTOF(dpotrf)( "L", n, A, n, info );
    /* dpotrf( 'L', *n, A, *n, info ); */
}

/*****************************************************************************/

void MatrixCholSolve( double* A, int* n, double* rhs, int *nrhs, int* info )
{
    CTOF(dpotrs)("L", n, nrhs, A, n, rhs, n, info );
    /* dpotrs('L', *n, *nrhs, A, *n, rhs, *n, info ); */
}

/*****************************************************************************/

void MatrixLUFactorize( double* A, int* n, int* ipiv, int* info )
{
    CTOF(dgetrf)( n, n, A, n, ipiv, info );
    /* dgetrf( *n, *n, A, *n, ipiv, info ); */
}

/*****************************************************************************/

void MatrixLUSolve( double* A, int* n, int* ipiv, double* rhs, int *nrhs )
{
    int info = 0;
    /*int one  = 1;*/
    CTOF(dgetrs)("N", n, nrhs, A, n, ipiv, rhs, n, &info );
    /* dgetrs('N', *n, *nrhs, A, *n, ipiv, rhs, *n, &info ); */
}

/*****************************************************************************/

void MatrixMatrixPlusDiag( double* A, double* D, int* n )
{
    int i;
    for (i=0;i<(*n);i++) A[i+i*(*n)] += D[i];
}

/*****************************************************************************/

void MatrixMatrixCopy( double* lhs, double* rhs, int* rows, int* cols )
{
    int i;
    int len = (*rows)*(*cols);
    for (i=0;i<len;i++) lhs[i] = rhs[i];
}

/*****************************************************************************/

void MatrixMatrixMinus( double* x, double* y, double* res, int* rows, int* cols)
{
    int i;
    int len = (*rows)*(*cols);
    for (i=0;i<len;i++) res[i] = x[i] - y[i];
}

/*****************************************************************************/

void MatrixMatrixPlus( double* x, double* y, double* res, int* rows, int* cols)
{
    int i;
    int len = (*rows)*(*cols);
    for (i=0;i<len;i++) res[i] = x[i] + y[i];
}

/*****************************************************************************/

void MatrixMatrixPlusEquals( double* x, double* res, int* rows, int* cols)
{
    int i;
    int len = (*rows)*(*cols);
    for (i=0;i<len;i++) res[i] += x[i] ;
}

/*****************************************************************************/

void MatrixConstantSet( double* A, double c, int* rows, int* cols)
{
    int i;
    int len = (*rows)*(*cols);
    for (i=0;i<len;i++) A[i] = c;
}

/*****************************************************************************/

void MatrixMatrixMult( double *alpha, double* A, int transA, double* B,
    int transB, double* beta, double* C, int* rA, int *cA, int *rB, int *cB,
    int* rC, int* cC  )
{
    if (transA)
    {
        if (transB)
        {
            CTOF(dgemm)( "T", "T", rC, cC, cB, alpha, A, rA, B, rB, beta, C, rC );
            /* dgemm( 'T', 'T', *rC, *cC, *cB, *alpha, A, *rA, B, *rB, *beta, C, *rC ); */
        }
        else
        {
            CTOF(dgemm)( "T", "N", rC, cC, rB, alpha, A, rA, B, rB, beta, C, rC );
            /* dgemm( 'T', 'N', *rC, *cC, *rB, *alpha, A, *rA, B, *rB, *beta, C,*rC ); */
        }
    }
    else
    {
        if (transB)
        {
            CTOF(dgemm)( "N", "T", rC, cC, cA, alpha, A, rA, B, rB, beta, C, rC );
            /* dgemm( 'N', 'T', *rC, *cC, *cA, *alpha, A, *rA, B, *rB, *beta, C, *rC ); */
        }
        else
        {
            CTOF(dgemm)( "N", "N", rC, cC, cA, alpha, A, rA, B, rB, beta, C, rC );
            /* dgemm( 'N', 'N', *rC, *cC, *cA, *alpha, A, *rA, B, *rB, *beta, C, *rC ); */
        }
    }
}

/*****************************************************************************/

void LRQPHeader()
{
    printf("ITER  PRIM            DUAL            COMP            GAP           TERM\n");
}

/*****************************************************************************/

void LRQPInitPoint( int *n, int *m, int *p, double *Q, double *c, double *A,
    double *b, double *u, double *alpha, double* beta, double *xi, double *zeta,
    double *w, double *temp )
{
    int i;
    double pone =  1.0;
    double mone = -1.0;
    double zero =  0.0;
    for (i=0;i<(*n);i++) alpha[i] = MIN(EPSINIT,u[i]*EPSINIT);
    for (i=0;i<(*p);i++) beta[i]  = 0.0;
    MatrixVectorMult( &pone, Q, 1, alpha, &zero, w, n, m );
    if ( (*n)!=(*m) ) MatrixVectorMult( &mone, Q, 0, w, &zero, temp, n, m );
    else              for (i=0;i<(*n);i++) temp[i] += -w[i];
    VectorVectorMult( &mone, c, temp, n );
    for (i=0;i<(*n);i++)
    {
        xi[i]   = MAX(EPSINIT,temp[i]);
        zeta[i] = MAX(EPSINIT,xi[i]-temp[i]);
    }
}

/*****************************************************************************
*
*  MATLAB CODE FOR QpIpmCalcStats
*
*    w               = Q'*alpha;
*    UminusAlpha     = (u - alpha);
*    XiOnUminusAlpha = xi./UminusAlpha;
*    ZetaOnAlpha     = zeta./alpha;
*
*    if (n==m)
*        r1 = - c - w   - A*beta - xi + zeta;
*        quad = alpha'*w;
*    else
*        r1 = - c - Q*w   - A*beta - xi + zeta;
*        quad = w'*w;
*    end;
*
*    r2      = b - A'*alpha;
*    comp    = alpha'*zeta + UminusAlpha'*xi;
*    cTalpha = c'*alpha;
*    gap     = abs(quad + cTalpha + u'*xi + b'*beta);
*    term    = comp / (abs( 0.5*quad + cTalpha ) + 1.0);
*    t       = comp*(((1.0-mult+eps)/(10+mult))^2)/(2*n);
*    D       = spdiags( XiOnUminusAlpha + ZetaOnAlpha + 1.0e-14, 0, n, n );
*
******************************************************************************/

void LRQPCalcStats( int *n, int *m, int *p, double *Q, double *c, double *A,
    double *b, double *u, double *alpha, double* beta, double *xi, 
    double *zeta, double *dalpha, double* dbeta, double *dxi, double *dzeta,
    double *UminusAlpha, double *XiOnUminusAlpha, double *ZetaOnAlpha, 
    double* w, double *r1, double *r2, double *D, double *prim, double *dual, 
    double *comp, double *gap, double *term, double * mult, double *t)
{
    int i;
    int one = 1;
    double quad;
    double cTalpha;
    double temp;
    double pone =  1.0;
    double mone = -1.0;
    double zero =  0.0;
    MatrixVectorMult( &pone, Q, 1, alpha, &zero, w, n, m );
    VectorVectorMinus( u, alpha, UminusAlpha, n );
    VectorVectorDivide( xi, UminusAlpha, XiOnUminusAlpha, n );
    VectorVectorDivide( zeta, alpha, ZetaOnAlpha, n );
    if ( (*n)==(*m) )
    {   
        if (*p) MatrixVectorMult( &mone, A, 0, beta, &zero, r1, n, p );
        else    MatrixConstantSet( r1, 0.0, n, &one );
        VectorVectorMult( &mone, w,  r1, n );
        VectorVectorMult( &mone, c,  r1, n );
        VectorVectorMult( &mone, xi, r1, n );
        VectorVectorMult( &pone, zeta, r1, n );
        quad = VectorVectorDot( alpha, w, n );
    }
    else
    {
        if (*p) MatrixVectorMult( &mone, A, 0, beta, &zero, r1, n, p );
        else    MatrixConstantSet( r1, 0.0, n, &one );
        MatrixVectorMult( &mone, Q, 0, w, &pone, r1, n, m );
        VectorVectorMult( &mone, c,  r1, n );
        VectorVectorMult( &mone, xi, r1, n );
        VectorVectorMult( &pone, zeta, r1, n );
        quad = VectorVectorDot( w, w, m );
    }
    if (*p)
    {
        VectorVectorCopy( r2, b, p );
        MatrixVectorMult( &mone, A, 1, alpha, &pone, r2, n, p );
        *dual = VectorAbsSum( r2, p );
    }
    *prim   = VectorAbsSum( r1, n );
    *comp   = VectorVectorDot( alpha, zeta, n ) + VectorVectorDot( UminusAlpha, xi, n );
    cTalpha = VectorVectorDot( c, alpha, n );
    if (*p) *gap = fabs( quad + cTalpha + VectorVectorDot( u, xi, n ) + VectorVectorDot( b, beta, p ) );
    else    *gap = fabs( quad + cTalpha + VectorVectorDot( u, xi, n ) );
    *term   = *comp / ( fabs( 0.5*quad + cTalpha) + 1.0);
    temp    = (1.0 - *mult + EPSIPM)/(10.0 + *mult);
    *t      = *comp*(temp*temp)/(2*(*n));
    for (i=0;i<(*n);i++) D[i] = XiOnUminusAlpha[i] + ZetaOnAlpha[i] + EPSPERT;
}

/*****************************************************************************/

void LRQPDisplay( int i, double *prim, double *dual, double *comp, double *gap,
                  double *term )
{
    printf("%3d %15.7e %15.7e %15.7e %15.7e %15.7e \n", i, *prim, *dual, *comp, *gap, *term );
}

/*****************************************************************************/

void LRQPSummary( int i, int niter, int method, int n, int m, double *prim, 
    double *dual, double *comp, double *gap, double *term )
{
    if (i==niter)
    {
        printf("LowRankQP FAILED TO CONVERGE\n");
        if (n==m)
        {
            if (method==CHOL) printf("    Try increasing niter, using method=LU, or rescaling problem.\n");
            else              printf("    Try increasing niter, or rescaling problem.\n");
        }
        else
        {
            if (method==SMW)  printf("    Try increasing niter, using method=PFCF, using method=LU, or rescaling problem.\n");
            else              printf("    Try increasing niter, or rescaling problem.\n");
        }
    }
    else
    {
        printf("LowRankQP CONVERGED IN %d ITERATIONS\n\n", i+1 );
        printf("    Primal Feasibility    = %15.7e\n", *prim);
        printf("    Dual Feasibility      = %15.7e\n", *dual);
        printf("    Complementarity Value = %15.7e\n", *comp);
        printf("    Duality Gap           = %15.7e\n", *gap);
        printf("    Termination Condition = %15.7e\n", *term);
    }
}

/****************************************************************************/

void PfcfSolve( int *n, double *p, double *beta, double *r, int trans )
{
    int j;
    double sigma;
    int nm1 = (*n)-1;
    if (trans)
    {
        sigma = r[nm1]*p[nm1];
        for (j=(nm1-1);j>0;j--)
        {
            r[j]  -= sigma*beta[j];
            sigma += r[j]*p[j];
        }
        r[0] -= sigma*beta[0];
    }
    else
    {
        sigma = r[0]*beta[0];
        for (j=1;j<nm1;j++)
        {
            r[j]  -= sigma*p[j];
            sigma += r[j]*beta[j];
        }
        r[nm1] -= sigma*p[nm1];
    }
}

/****************************************************************************/

void PfcfFactorize( int *n, int *m, double *Q, double *D, double *P,
    double *Beta, double *Lambda, double *LambdaTemp, double *t  )
{
    int i, j/*, k*/;
    double infinity = 1.0e12;
    double epsilon  = 1.0e-12;
    double temp;

    /*int one = 1;*/

    int ind;
    /*int indi;*/
    /*int indj;*/

    VectorVectorCopy( LambdaTemp, D, n );
    MatrixMatrixCopy( P, Q, n, m );
    for (i=0;i<(*m);i++)
    {
        for (j=0;j<i;j++) PfcfSolve( n, P+j*(*n), Beta+j*(*n), P+i*(*n), 0 );
        t[0] = 1.0;
        for (j=0;j<(*n);j++)
        {
            ind  = j+i*(*n);
            temp = P[ind];   
            if ( fabs(t[j]) < infinity )
            {
                if ( fabs(LambdaTemp[j]) > epsilon )
                {
                    t[j+1]    = t[j] + ( temp*temp ) / LambdaTemp[j];
                    Lambda[j] = LambdaTemp[j]*t[j+1] / t[j];
                    Beta[ind] = temp / ( LambdaTemp[j] * t[j+1]);
                }
                else
                {
                    if ( fabs(temp) > epsilon )
                    {
                        t[j+1]    = infinity;
                        Lambda[j] = temp*temp / t[j];
                        Beta[ind] = 1.0 / temp;
                    }
                    else
                    {
                        t[j+1]    = t[j];
                        Lambda[j] = 0.0;
                        Beta[ind] = 0.0;
                    }
                }
            }
            else
            {
                t[j+1]    = infinity;
                Lambda[j] = LambdaTemp[j];
                Beta[ind] = 0.0;
            }
        }
        VectorVectorCopy( LambdaTemp, Lambda, n );
    }
}

/*****************************************************************************/

void LRQPFactorize( int *n, int *m, int *method, double *Q, double *D,
    double *M, int* pivN, double *buffNxM, double *P, double *Beta, 
    double *Lambda, double *LambdaTemp, double *T )
{
    int    info =  0;
    double pone =  1.0;
    /*double mone = -1.0;*/
    double zero =  0.0;

    if ((*method==LU)||(*method==CHOL))
    {
        if ((*n)==(*m)) MatrixMatrixCopy( M, Q, n, n );
        else            MatrixMatrixMult( &pone, Q, 0, Q, 1, &zero, M, n, m, n, m, n, n );
        MatrixMatrixPlusDiag( M, D, n );
        if (*method==LU)        MatrixLUFactorize( M, n, pivN, &info );
        else if (*method==CHOL) MatrixCholFactorize( M, n, &info );
    }
    else
    {
        if (*method==SMW)
        {
            MatrixMatrixDiagSolve( D, Q, buffNxM, n, m );
            MatrixMatrixMult( &pone, Q, 1, buffNxM, 0, &zero, M, n, m, n, m, m, m );
            MatrixConstantPlusDiag( M, pone, m );
            MatrixCholFactorize( M, m, &info );
        }
        else if (*method==PFCF)
        {
            PfcfFactorize( n, m, Q, D, P, Beta, Lambda, LambdaTemp, T );
        }
    }
}

/******************************************************************************
*
*  MATLAB CODE FOR IpmSolve
*
*    if (strcmp(info.method,'FULL'))
*        sol = info.M \ rhs;
*    elseif (strcmp(info.method,'CHOL'))
*        z   = info.M'\rhs;
*        sol = info.M\z;
*    elseif (strcmp(info.method,'SMW'))
*        z   = D \ rhs;
*        t   = info.M'\(Q'*z);
*        t   = info.M \ t;
*        sol = z - D \ (Q*t);
*    elseif (strcmp(info.method,'PFCF'))
*        PfcfSolve()
*    end;
*
******************************************************************************/

void LRQPSolve( int *n, int *m, int *nrhs, int *method, double *Q, double *D,
    double *rhs, double *sol, double *M, int* pivN, double *buffMxP, double *P,
    double *Beta, double *Lambda )
{
    int i, j;
    int info    = 0;
    int one     = 1;
    double pone =  1.0;
    double mone = -1.0;
    double zero =  0.0;

    
    /*long int start;*/
    /*long int finish;*/

    MatrixMatrixCopy( sol, rhs, n, nrhs );
    if (*method==LU)
    {
        MatrixLUSolve( M, n, pivN, sol, nrhs );
    }
    else if (*method==CHOL)
    {
        MatrixCholSolve( M, n, sol, nrhs, &info );
    }
    else if (*method==SMW)
    {
        MatrixMatrixDiagSolve( D, sol, sol, n, nrhs );
        MatrixMatrixMult( &pone, Q, 1, sol, 0, &zero, buffMxP, n, m, n, nrhs, m, nrhs );
        MatrixCholSolve( M, m, buffMxP, nrhs, &info );
        MatrixMatrixMult( &mone, Q, 0, buffMxP, 0, &zero, sol, n, m, m, nrhs, n,nrhs );
        MatrixMatrixPlus( sol, rhs, sol, n, nrhs );
        MatrixMatrixDiagSolve( D, sol, sol, n, nrhs );
    }
    else if (*method==PFCF)
    {
        for (i=0;i<(*nrhs);i++)
        {
            for (j=0;j<(*m);j++) PfcfSolve( n, P+(j*(*n)), Beta+(j*(*n)), sol+(i*(*n)), 0);
            MatrixMatrixDiagSolve( Lambda, sol+(i*(*n)), sol+(i*(*n)), n, &one );
            for (j=((*m)-1);j>=0;j--) PfcfSolve( n, P+(j*(*n)), Beta+(j*(*n)),sol+(i*(*n)), 1);
        }
    }
}

/******************************************************************************
*
*  MATLAB CODE FOR IpmSolve
*
*    R      = QpIpmSolve( PREDinfo, Q, D, A );
*    r3     = - zeta;    OR    r3 = (t - Dalpha.*Dzeta)./alpha - zeta;
*    r4     = - xi;      OR    r4 = (t + Dalpha.*Dxi )./UminusAlpha - xi;
*    r5     = r1 + r3 - r4
*    r      = QpIpmSolve( info, Q, D, r5 );
*    Dbeta  = (A'*R) \ (A'*r - r2);
*    Dalpha = r - R*Dbeta;
*    Dzeta  = r3 - ZetaOnAlpha.*Dalpha;
*    Dxi    = r4 + XiOnUminusAlpha.*Dalpha;
*
******************************************************************************/

void LRQPCalcDx( int *n, int *m, int *p, int *method, double *Q, double *c,
    double *A, double *b, double * u, double *alpha, double* beta, double *xi,
    double *zeta, double *dalpha, double* dbeta, double *dxi, double *dzeta,
    double *UminusAlpha, double *ZetaOnAlpha, double *XiOnUminusAlpha,
    double *buffMxP, double *buffMx1, double* buffPxP, double *buffPx1,
    int *pivN, double *R, double *r, double *r1, double* r2, double *r3,
    double *r4, double* r5, double *D, double *M, double *t, double *P,
    double *Beta, double *Lambda, double *LambdaTemp, double *T, int predcorr)
{

    int i/*, j*/;
    int    info = 0;
    int    one  = 1;
    double pone =  1.0;
    double mone = -1.0;
    double zero =  0.0;
    if ((predcorr==PRED)&&(*p))
    {
        LRQPSolve( n, m, p, method, Q, D, A, R, M, pivN, buffMxP, P, Beta, Lambda );
    }
    for (i=0;i<(*n);i++)
    {
        r3[i] = - zeta[i];
        r4[i] = - xi[i];
    }
    if (predcorr==CORR)
    {
        for (i=0;i<(*n);i++) r3[i] += ( *t - (dalpha[i] * dzeta[i]) )/alpha[i];
        for (i=0;i<(*n);i++) r4[i] += ( *t + (dalpha[i] * dxi[i]) )/UminusAlpha[i];
    }
    for (i=0;i<(*n);i++) r5[i] = r1[i] + r3[i] - r4[i];

    if (*p)
    {
        LRQPSolve( n, m, &one, method, Q, D, r5, r, M, pivN, buffMx1, P, Beta, Lambda );
        VectorVectorCopy( buffPx1, r2, p );
        MatrixVectorMult( &pone, A, 1, r, &mone, buffPx1, n, p );
        MatrixMatrixMult( &pone, A, 1, R, 0, &zero, buffPxP, n, p, n, p, p, p );
        MatrixCholFactorize( buffPxP, p, &info );
        MatrixCholSolve( buffPxP, p, buffPx1, &one, &info);
        VectorVectorCopy( dbeta, buffPx1, p );
        VectorVectorCopy( dalpha, r , n);
        MatrixVectorMult( &mone, R, 0, dbeta, &pone, dalpha, n, p );
    }
    else
    {
        LRQPSolve( n, m, &one, method, Q, D, r5, dalpha, M, pivN, buffMx1, P, Beta, Lambda );
    }
    for (i=0;i<(*n);i++) dzeta[i] = r3[i] - (ZetaOnAlpha[i] * dalpha[i]);
    for (i=0;i<(*n);i++) dxi[i]   = r4[i] + (XiOnUminusAlpha[i] * dalpha[i]);

}

/*****************************************************************************
*
*  MATLAB CODE FOR QpIpmStep
*
*    temp = [ 1.0 (-alpha./Dalpha)' (UminusAlpha./Dalpha)' (-xi./Dxi)'(-zeta./Dzeta)' ];
*    mult = 0.99*min( temp( find(temp>0.0) ) );
*    alpha = alpha + mult*Dalpha;
*    beta  = beta  + mult*Dbeta;
*    xi    = xi    + mult*Dxi;
*    zeta  = zeta  + mult*Dzeta;
*
*****************************************************************************/

void LRQPStep( int *n, int *p, double *alpha, double* beta, double *xi,
    double *zeta, double *dalpha, double* dbeta, double *dxi, double *dzeta,
    double *UminusAlpha, double *mult)
{
    int i;
    *mult= 1.0;
    for (i=0;i<(*n);i++)
    {
        if (dalpha[i]<0.0) *mult = MIN(*mult,(-alpha[i]/dalpha[i]));
        if (dalpha[i]>0.0) *mult = MIN(*mult,(UminusAlpha[i]/dalpha[i]));
        if (dxi[i]<0.0)    *mult = MIN(*mult,(-xi[i]/dxi[i]));
        if (dzeta[i]<0.0)  *mult = MIN(*mult,(-zeta[i]/dzeta[i]));
    }
    *mult *= 0.99;
    VectorVectorMult( mult, dalpha, alpha, n );
    if (*p) VectorVectorMult( mult, dbeta,  beta,  p );
    VectorVectorMult( mult, dxi,    xi,    n );
    VectorVectorMult( mult, dzeta,  zeta,  n );
}

/******************************************************************************/

void LowRankQP( int *n, int *m, int *p, int* method, int* verbose, int* niter, 
    double *Q, double *c, double *A, double *b, double *u, double *alpha,
    double* beta, double *xi, double *zeta)
{
    int i;
    
    /*long int start;*/
    /*long int finish;*/

    /* Iteration Display variables */
    double mult = 0.0;
    double prim = 0.0;
    double dual = 0.0;
    double comp = 0.0;
    double gap  = 0.0;
    double term = 0.0;
    double t    = 0.0;

    /* Step direction vectors */
    double *dalpha = (double *) calloc( (*n), sizeof(double) );
    double *dbeta  = NULL;
    double *dxi    = (double *) calloc( (*n), sizeof(double) );
    double *dzeta  = (double *) calloc( (*n), sizeof(double) );

    /* Some repeatedly occuring vectors */
    double *UminusAlpha     = (double *) calloc( *n, sizeof(double) );
    double *XiOnUminusAlpha = (double *) calloc( *n, sizeof(double) );
    double *ZetaOnAlpha     = (double *) calloc( *n, sizeof(double) );

    /* Some vectors used during calculations */
    double *w  = (double *) calloc( *m, sizeof(double) );
    double *r1 = (double *) calloc( *n, sizeof(double) );
    double *r2 = NULL;
    double *r3 = (double *) calloc( *n, sizeof(double) );
    double *r4 = (double *) calloc( *n, sizeof(double) );
    double *r5 = (double *) calloc( *n, sizeof(double) );
    double *D  = (double *) calloc( *n, sizeof(double) );
    double *r  = (double *) calloc( *n, sizeof(double) );
    double *R  = NULL;

    /* Various Buffers */
    double *buffMxP = NULL;
    double *buffPxP = NULL;
    double *buffPx1 = NULL;

    double *M = NULL;
    int    *pivN = NULL;

    double *buffNxM = NULL;
    double *buffMx1 = NULL;

    double *P = NULL;
    double *Beta = NULL;
    double *Lambda = NULL;
    double *LambdaTemp = NULL;
    double *T = NULL;

    /* Vectors to be created if p!=0 */
    if (*p)
    {
        dbeta   = (double *) calloc( (*p), sizeof(double) );
        r2      = (double *) calloc( (*p), sizeof(double) );
        R       = (double *) calloc( (*n)*(*p), sizeof(double) );
        buffMxP = (double *) calloc( (*m)*(*p), sizeof(double) );
        buffPxP = (double *) calloc( (*p)*(*p), sizeof(double) );
        buffPx1 = (double *) calloc( (*p), sizeof(double) );
    }

    if ((*method==LU)||(*method==CHOL))
    {
        M    = (double *) calloc( (*n)*(*n), sizeof(double) );
        pivN = (int *) calloc( *n, sizeof(int) );
    }
    else
    {
        if (*method==SMW)
        {
            buffNxM  = (double *) calloc( (*n)*(*m), sizeof(double) );
            M        = (double *) calloc( (*m)*(*m), sizeof(double) );
            buffMx1  = (double *) calloc( (*m), sizeof(double) );
        }
        else if (*method==PFCF)
        {
            P          = (double *) calloc( (*n)*(*m), sizeof(double) );
            Beta       = (double *) calloc( (*n)*(*m), sizeof(double) );
            Lambda     = (double *) calloc( (*n), sizeof(double) );
            LambdaTemp = (double *) calloc( (*n), sizeof(double) );
            T          = (double *) calloc(1+(*n), sizeof(double) );
        }
    }

    /* Main Loop */
    if ( *verbose ) LRQPHeader();
    LRQPInitPoint( n, m, p, Q, c, A, b, u, alpha, beta, xi, zeta, w, r1 );

    for (i=0;i<(*niter);i++)
    {
        /* start = clock(); */
        LRQPCalcStats( n, m, p, Q, c, A, b, u, alpha, beta, xi, zeta, dalpha,
            dbeta, dxi, dzeta, UminusAlpha, XiOnUminusAlpha, ZetaOnAlpha, w, r1,
            r2, D, &prim, &dual, &comp, &gap, &term, &mult, &t );
        /* finish = clock();
        TimeIpmCalcStats += (finish - start)/CLK_TCK; */

        if ( *verbose ) LRQPDisplay( i+1, &prim, &dual, &comp, &gap, &term  );
        if ( term < EPSTERM ) break;

        LRQPFactorize( n, m, method, Q, D, M, pivN, buffNxM, P, Beta, Lambda,
            LambdaTemp, T );

        LRQPCalcDx( n, m, p, method, Q, c, A, b, u, alpha, beta, xi, zeta,
            dalpha, dbeta, dxi, dzeta, UminusAlpha, ZetaOnAlpha, 
            XiOnUminusAlpha, buffMxP, buffMx1, buffPxP, buffPx1, pivN, R, r, r1,
            r2, r3, r4, r5, D, M, &t, P, Beta, Lambda, LambdaTemp, T, PRED);

        LRQPCalcDx( n, m, p, method, Q, c, A, b, u, alpha, beta, xi, zeta,
            dalpha, dbeta, dxi, dzeta, UminusAlpha, ZetaOnAlpha,
            XiOnUminusAlpha, buffMxP, buffMx1, buffPxP, buffPx1, pivN, R, r, r1,
            r2, r3, r4, r5, D, M, &t, P, Beta, Lambda, LambdaTemp, T, CORR);

        LRQPStep( n, p, alpha, beta, xi, zeta, dalpha, dbeta, dxi, dzeta,
            UminusAlpha, &mult );
    }

    if ( *verbose ) LRQPSummary( i, *niter, *method, *n, *m, &prim, &dual, &comp, &gap, &term );

    /* Free Memory */
    free(dalpha);      free(dxi);             free(dzeta);
    free(UminusAlpha); free(XiOnUminusAlpha); free(ZetaOnAlpha);
    free(r1); free(r3); free(r4); free(r5);
    free(D);  free(w);  free(r);  
        
    if (*p)
    {
        free( dbeta );
        free( r2 );
        free( R );
        free( buffMxP );
        free( buffPxP );
        free( buffPx1 );
    }

    if ((*method==LU)||(*method==CHOL))
    {
        free( M );
        free( pivN );
    }
    else
    {
        if (*method==SMW)
        {
            free( buffMx1 );
            free( buffNxM );
            free( M  );
        }
        else if (*method==PFCF)
        {
            free( P );
            free( Beta );
            free( Lambda );
            free( LambdaTemp );
            free( T );
        }
    }
}
