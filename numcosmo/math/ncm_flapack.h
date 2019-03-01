
#ifndef NUMCOSMO_GIR_SCAN

extern double dasum_ (int *, double *, int *);
extern int idamax_ (const int *n, const double *dx, const int *incx);
extern void dscal_ (const int *n, const double *da, const double *dx, const int *incx);
extern void daxpy_ (const int *n, const double *alpha, const double *dx, const int *incx, double *dy, const int *incy);
extern void dcopy_ (const int *n, const double *dx, const int *incx, double *dy, const int *incy);
extern double ddot_ (const int *n, const double *dx, const int *incx, const double *dy, const int *incy);
extern void dgemv_ (const char *trans, const int *m, const int *n, const double *alpha, const double *a, const int *lda, const double *x, const int *incx, const double *beta, double *y, const int *incy);
extern void dgemm_ (const char *transa, const char *transb, const int *m, const int *n, const int *k, const double *alpha, const double *a, const int *lda, const double *b, const int *ldb, const double *beta, double *c, const int *ldc);
extern void dgetrf_ (int *m, int *n, double *a, int *lda, int *ipiv, int *info);
extern void dgetrs_ (char *trans, int *n, int *nrhs, double *a, int *lda, int *ipiv, double *b, int *ldb, int *info);

extern void dptsv_ (gint *n, gint *nrhs, gdouble *d, gdouble *e, gdouble *b, gint *ldb, gint *info);
extern void dpotrf_ (const gchar *uplo, const gint *n, gdouble *a, const gint *lda, gint *info);
extern void dpotri_ (const gchar *uplo, const gint *n, gdouble *a, const gint *lda, gint *info);
extern void dpotrs_ (const gchar *uplo, const gint *n, const gint *nrhs, gdouble *a, const gint *lda, gdouble *b, const gint *ldb, gint *info);
extern void dposv_ (const gchar *uplo, const gint *n, const gint *nrhs, gdouble *a, const gint *lda, gdouble *b, const gint *ldb, gint *info);

extern void dsytrf_ (const gchar *uplo, gint *n, gdouble *a, gint *lda, gint *ipiv, gdouble *work, gint *lwork, gint *info);
extern void dsytrs_ (const gchar *uplo, const gint *n, const gint *nrhs, gdouble *a, const gint *lda, const gint *ipiv, gdouble *b, const gint *ldb, gint *info);
extern void dsytri_ (const gchar *uplo, const gint *n, gdouble *a, const gint *lda, const gint *ipiv, gdouble *work, gint *info);
extern void dsysvx_ (const gchar *fact, gchar *uplo, const gint* n, const gint *nrhs, gdouble *a, const gint *lda, gdouble *af, const gint *ldaf, gint *ipiv, gdouble *b, const gint *ldb, gdouble *x, const gint *ldx, gdouble *rcond, gdouble *ferr, gdouble *berr, gdouble *work, gint *lwork, gint *iwork, gint *info);
extern void dsysvxx_ (const gchar *fact, gchar *uplo, const gint* n, const gint *nrhs, gdouble *a, const gint *lda, gdouble *af, const gint *ldaf, gint *ipiv, gchar *equed, gdouble *s, gdouble *b, const gint *ldb, gdouble *x, const gint *ldx, gdouble *rcond, gdouble *rpvgrw, gdouble *berr, const gint *n_err_bnds, gdouble *err_bnds_norm, gdouble *err_bnds_comp, const gint *nparams, gdouble *params, gdouble *work, gint *iwork, gint *info);

extern void dsyevr_ (const gchar *jobz, const gchar *range, const gchar *uplo, const gint *n, gdouble *a, const gint *lda, const gdouble *vl, const gdouble *vu, const gint *il, const gint *iu, const gdouble *abstol, gint *m, gdouble *w, double *z, const gint *ldz, gint *isuppz, gdouble *work, const gint *lwork, gint *iwork, const gint *liwork, gint *info);
extern void dsyevd_ (const gchar *jobz, const gchar *uplo, const gint *n, gdouble *a, const gint *lda, gdouble *w, gdouble *work, const gint *lwork, gint *iwork, const gint *liwork, gint *info);

extern void dgeev_ (const gchar *jobvl, const gchar *jobvr, gint *n, gdouble *a, gint *lda, gdouble *wr, gdouble *wi, gdouble *vl, gint *ldvl, gdouble *vr, gint *ldvr, gdouble *work, gint *lwork, gint *info);
extern void dgeevx_ (const gchar *balanc, const gchar *jobvl, const gchar *jobvr, const gchar *sense, const gint *n, gdouble *a, const gint *lda, gdouble *wr, gdouble *wi, gdouble *vl, const gint *ldvl, gdouble *vr, const gint *ldvr, gint *ilo, gint *ihi, gdouble *scale, gdouble *abnrm, gdouble *rconde, gdouble *rcondv, gdouble *work, const gint *lwork, gint *iwork, gint *info);

extern void dggglm_ (const gint *n, const gint *M, const gint *P, gdouble *X, const gint *lda, gdouble *L, const gint *ldb, gdouble *d, gdouble *p, gdouble *y, gdouble *work, const gint *lwork, gint *info);

extern void dgeqrf_ (const gint *m, const gint *n, gdouble *a, const gint *lda, gdouble *tau, gdouble *work, const gint *lwork, gint *info);
extern void dgerqf_ (const gint *m, const gint *n, gdouble *a, const gint *lda, gdouble *tau, gdouble *work, const gint *lwork, gint *info);
extern void dgeqlf_ (const gint *m, const gint *n, gdouble *a, const gint *lda, gdouble *tau, gdouble *work, const gint *lwork, gint *info);
extern void dgelqf_ (const gint *m, const gint *n, gdouble *a, const gint *lda, gdouble *tau, gdouble *work, const gint *lwork, gint *info);

extern void dtrsv_ (const gchar *uplo, gchar *t, gchar *d, gint *n, gdouble *a, gint *lda, gdouble *x, gint *incx);

#endif /* NUMCOSMO_GIR_SCAN */
