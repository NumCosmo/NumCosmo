#define dtrsv dtrsv_
#define dpotrf dpotrf_
#define dpotrs dpotrs_
#define dpotri dpotri_
#define dtrtri dtrtri_
#define dtrmm dtrmm_
#define dtrmv dtrmv_
#define dgeqrf dgeqrf_
#define dormqr dormqr_
#define dsyev dsyev_
#define dgesvd dgesvd_
#define dsymv dsymv_
#define dgemv dgemv_
#define dgemm dgemm_
#define dsyrk dsyrk_
#define dsyr2k dsyr2k_
#define daxpy daxpy_
#define dtrsm dtrsm_
#define dsymm dsymm_
#define dsyr dsyr_
#define ddot ddot_

void dtrsv(const char *uplo, const char *trans, const char *diag, const int  *n,
           const double *a, const int *lda, double *x, const int *incx);
void dpotrf( char* uplo, int * n, double* a, int * lda, int * info );
void dpotri( char* uplo, int * n, double* a, int * lda, int * info );
void dgemv(const char *trans, const int *m, const int *n, const double *alpha,
           const double *a, const int *lda, const double *x, const int *incx,
           const double *beta, double *y, const int *incy);
void dsyrk(const char *uplo, const char *trans, const int *n, const int *k,
           const double *alpha, const double *a, const int *lda, const double *beta,
           double *c, const int *ldc);
void dsyr2k(const char *uplo, const char *trans, const int *n, const int *k,
            const double *alpha, const double *a, const int *lda, const double *b, const int *ldb,
            const double *beta, double *c, const int *ldc);
void dgesvd( char* jobu, char* jobvt, int * m, int * n, double* a, int * lda, double* s, double* u, int * ldu, double* vt, int * ldvt, double* work, int * lwork, int * info );
void dgemm(const char *transa, const char *transb, const int *m, const int *n, const int *k,
           const double *alpha, const double *a, const int *lda, const double *b, const int *ldb,
           const double *beta, double *c, const int *ldc);
void dtrtri( char* uplo, char* diag, int * n, double* a, int * lda, int * info );
void dtrmm(const char *side, const char *uplo, const char *transa, const char *diag,
           const int *m, const int *n, const double *alpha, const double *a, const int *lda,
           double *b, const int *ldb);
void dtrmv(const char *uplo, const char *transa, const char *diag, const int *n,
           const double *a, const int *lda, double *b, const int *incx);
void dgeqrf( int * m, int * n, double* a, int * lda, double* tau, double* work, int * lwork, int * info );
void dormqr( char* side, char* trans, int * m, int * n, int * k, double* a, int * lda, double* tau, double* c, int * ldc, double* work, int * lwork, int * info );
void dsyev( char* jobz, char* uplo, int * n, double* a, int * lda, double* w, double* work, int * lwork, int * info );
void dsymv(const char *uplo, const int *n, const double *alpha, const double *a, const int *lda,
           const double *x, const int *incx, const double *beta, double *y, const int *incy);
void daxpy(const int *n, const double *alpha, const double *x, const int *incx, double *y, const int *incy);
void dtrsm(const char *side, const char *uplo, const char *transa, const char *diag,
           const int *m, const int *n, const double *alpha, const double *a, const int *lda,
           double *b, const int *ldb);
void dsyr(const char *uplo, const int *n, const double *alpha, const double *x, const int *incx,
         double *a, const int *lda);
void dsymm(const char *side, const char *uplo, const int *m, const int *n,
           const double *alpha, const double *a, const int *lda, const double *b, const int *ldb,
           const double *beta, double *c, const int *ldc);
double ddot(int* N,double *DX, int* INCX,double *DY,int* INCY);
void dpotrs(char* UPLO,int * N,int * NRHS,double* A,int* LDA,double* B,int* LDB,int* INFO );
           
         

