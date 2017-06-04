

typedef long int ptrdiff_t;
typedef long unsigned int size_t;
typedef int wchar_t;
typedef struct {
  long long __max_align_ll __attribute__((__aligned__(__alignof__(long long))));
  long double __max_align_ld __attribute__((__aligned__(__alignof__(long double))));
} max_align_t;

typedef double (*nlopt_func)(unsigned n, const double *x,
        double *gradient,
        void *func_data);

typedef void (*nlopt_mfunc)(unsigned m, double *result,
       unsigned n, const double *x,
        double *gradient,
        void *func_data);



typedef void (*nlopt_precond)(unsigned n, const double *x, const double *v,
         double *vpre, void *data);

typedef enum {
     NLOPT_GN_DIRECT = 0,
     NLOPT_GN_DIRECT_L,
     NLOPT_GN_DIRECT_L_RAND,
     NLOPT_GN_DIRECT_NOSCAL,
     NLOPT_GN_DIRECT_L_NOSCAL,
     NLOPT_GN_DIRECT_L_RAND_NOSCAL,

     NLOPT_GN_ORIG_DIRECT,
     NLOPT_GN_ORIG_DIRECT_L,

     NLOPT_GD_STOGO,
     NLOPT_GD_STOGO_RAND,

     NLOPT_LD_LBFGS_NOCEDAL,

     NLOPT_LD_LBFGS,

     NLOPT_LN_PRAXIS,

     NLOPT_LD_VAR1,
     NLOPT_LD_VAR2,

     NLOPT_LD_TNEWTON,
     NLOPT_LD_TNEWTON_RESTART,
     NLOPT_LD_TNEWTON_PRECOND,
     NLOPT_LD_TNEWTON_PRECOND_RESTART,

     NLOPT_GN_CRS2_LM,

     NLOPT_GN_MLSL,
     NLOPT_GD_MLSL,
     NLOPT_GN_MLSL_LDS,
     NLOPT_GD_MLSL_LDS,

     NLOPT_LD_MMA,

     NLOPT_LN_COBYLA,

     NLOPT_LN_NEWUOA,
     NLOPT_LN_NEWUOA_BOUND,

     NLOPT_LN_NELDERMEAD,
     NLOPT_LN_SBPLX,

     NLOPT_LN_AUGLAG,
     NLOPT_LD_AUGLAG,
     NLOPT_LN_AUGLAG_EQ,
     NLOPT_LD_AUGLAG_EQ,

     NLOPT_LN_BOBYQA,

     NLOPT_GN_ISRES,



     NLOPT_AUGLAG,
     NLOPT_AUGLAG_EQ,
     NLOPT_G_MLSL,
     NLOPT_G_MLSL_LDS,

     NLOPT_LD_SLSQP,

     NLOPT_LD_CCSAQ,

     NLOPT_GN_ESCH,

     NLOPT_NUM_ALGORITHMS
} NcmFitNloptAlgorithm;

extern const char * NcmFitNloptAlgorithm_name(nlopt_algorithm a);

typedef enum {
     NLOPT_FAILURE = -1,
     NLOPT_INVALID_ARGS = -2,
     NLOPT_OUT_OF_MEMORY = -3,
     NLOPT_ROUNDOFF_LIMITED = -4,
     NLOPT_FORCED_STOP = -5,
     NLOPT_SUCCESS = 1,
     NLOPT_STOPVAL_REACHED = 2,
     NLOPT_FTOL_REACHED = 3,
     NLOPT_XTOL_REACHED = 4,
     NLOPT_MAXEVAL_REACHED = 5,
     NLOPT_MAXTIME_REACHED = 6
} NcmFitNloptResult;



extern void nlopt_srand(unsigned long seed);
extern void nlopt_srand_time(void);

extern void nlopt_version(int *major, int *minor, int *bugfix);
struct nlopt_opt_s;
typedef struct nlopt_opt_s *nlopt_opt;




extern nlopt_opt nlopt_create(NcmFitNloptAlgorithm algorithm, unsigned n);
extern void nlopt_destroy(nlopt_opt opt);
extern nlopt_opt nlopt_copy(const nlopt_opt opt);

extern NcmFitNloptResult nlopt_optimize(nlopt_opt opt, double *x,
      double *opt_f);

extern NcmFitNloptResult nlopt_set_min_objective(nlopt_opt opt, nlopt_func f,
        void *f_data);
extern NcmFitNloptResult nlopt_set_max_objective(nlopt_opt opt, nlopt_func f,
        void *f_data);

extern NcmFitNloptResult nlopt_set_precond_min_objective(nlopt_opt opt, nlopt_func f, nlopt_precond pre, void *f_data);
extern NcmFitNloptResult nlopt_set_precond_max_objective(nlopt_opt opt, nlopt_func f, nlopt_precond pre, void *f_data);

extern NcmFitNloptAlgorithm nlopt_get_algorithm(const nlopt_opt opt);
extern unsigned nlopt_get_dimension(const nlopt_opt opt);



extern NcmFitNloptResult nlopt_set_lower_bounds(nlopt_opt opt,
       const double *lb);
extern NcmFitNloptResult nlopt_set_lower_bounds1(nlopt_opt opt, double lb);
extern NcmFitNloptResult nlopt_get_lower_bounds(const nlopt_opt opt,
       double *lb);
extern NcmFitNloptResult nlopt_set_upper_bounds(nlopt_opt opt,
       const double *ub);
extern NcmFitNloptResult nlopt_set_upper_bounds1(nlopt_opt opt, double ub);
extern NcmFitNloptResult nlopt_get_upper_bounds(const nlopt_opt opt,
       double *ub);

extern NcmFitNloptResult nlopt_remove_inequality_constraints(nlopt_opt opt);
extern NcmFitNloptResult nlopt_add_inequality_constraint(nlopt_opt opt,
         nlopt_func fc,
         void *fc_data,
         double tol);
extern NcmFitNloptResult nlopt_add_precond_inequality_constraint(
     nlopt_opt opt, nlopt_func fc, nlopt_precond pre, void *fc_data,
     double tol);
extern NcmFitNloptResult nlopt_add_inequality_mconstraint(nlopt_opt opt,
           unsigned m,
           nlopt_mfunc fc,
           void *fc_data,
           const double *tol);

extern NcmFitNloptResult nlopt_remove_equality_constraints(nlopt_opt opt);
extern NcmFitNloptResult nlopt_add_equality_constraint(nlopt_opt opt,
       nlopt_func h,
       void *h_data,
       double tol);
extern NcmFitNloptResult nlopt_add_precond_equality_constraint(
     nlopt_opt opt, nlopt_func h, nlopt_precond pre, void *h_data,
     double tol);
extern NcmFitNloptResult nlopt_add_equality_mconstraint(nlopt_opt opt,
         unsigned m,
         nlopt_mfunc h,
         void *h_data,
         const double *tol);



extern NcmFitNloptResult nlopt_set_stopval(nlopt_opt opt, double stopval);
extern double nlopt_get_stopval(const nlopt_opt opt);

extern NcmFitNloptResult nlopt_set_ftol_rel(nlopt_opt opt, double tol);
extern double nlopt_get_ftol_rel(const nlopt_opt opt);
extern NcmFitNloptResult nlopt_set_ftol_abs(nlopt_opt opt, double tol);
extern double nlopt_get_ftol_abs(const nlopt_opt opt);

extern NcmFitNloptResult nlopt_set_xtol_rel(nlopt_opt opt, double tol);
extern double nlopt_get_xtol_rel(const nlopt_opt opt);
extern NcmFitNloptResult nlopt_set_xtol_abs1(nlopt_opt opt, double tol);
extern NcmFitNloptResult nlopt_set_xtol_abs(nlopt_opt opt, const double *tol);
extern NcmFitNloptResult nlopt_get_xtol_abs(const nlopt_opt opt,
          double *tol);

extern NcmFitNloptResult nlopt_set_maxeval(nlopt_opt opt, int maxeval);
extern int nlopt_get_maxeval(const nlopt_opt opt);

extern NcmFitNloptResult nlopt_set_maxtime(nlopt_opt opt, double maxtime);
extern double nlopt_get_maxtime(const nlopt_opt opt);

extern NcmFitNloptResult nlopt_force_stop(nlopt_opt opt);
extern NcmFitNloptResult nlopt_set_force_stop(nlopt_opt opt, int val);
extern int nlopt_get_force_stop(const nlopt_opt opt);



extern NcmFitNloptResult nlopt_set_local_optimizer(nlopt_opt opt,
          const nlopt_opt local_opt);

extern NcmFitNloptResult nlopt_set_population(nlopt_opt opt, unsigned pop);
extern unsigned nlopt_get_population(const nlopt_opt opt);

extern NcmFitNloptResult nlopt_set_vector_storage(nlopt_opt opt, unsigned dim);
extern unsigned nlopt_get_vector_storage(const nlopt_opt opt);

extern NcmFitNloptResult nlopt_set_default_initial_step(nlopt_opt opt,
        const double *x);
extern NcmFitNloptResult nlopt_set_initial_step(nlopt_opt opt,
       const double *dx);
extern NcmFitNloptResult nlopt_set_initial_step1(nlopt_opt opt, double dx);
extern NcmFitNloptResult nlopt_get_initial_step(const nlopt_opt opt,
       const double *x, double *dx);





typedef void* (*nlopt_munge)(void *p);
extern void nlopt_set_munge(nlopt_opt opt,
      nlopt_munge munge_on_destroy,
      nlopt_munge munge_on_copy);
typedef void* (*nlopt_munge2)(void *p, void *data);
extern void nlopt_munge_data(nlopt_opt opt,
                                    nlopt_munge2 munge, void *data);
typedef double (*nlopt_func_old)(int n, const double *x,
     double *gradient,
     void *func_data);

extern NcmFitNloptResult nlopt_minimize(
     NcmFitNloptAlgorithm algorithm,
     int n, nlopt_func_old f, void *f_data,
     const double *lb, const double *ub,
     double *x,
     double *minf,
     double minf_max, double ftol_rel, double ftol_abs,
     double xtol_rel, const double *xtol_abs,
     int maxeval, double maxtime) __attribute__((deprecated));

extern NcmFitNloptResult nlopt_minimize_constrained(
     NcmFitNloptAlgorithm algorithm,
     int n, nlopt_func_old f, void *f_data,
     int m, nlopt_func_old fc, void *fc_data, ptrdiff_t fc_datum_size,
     const double *lb, const double *ub,
     double *x,
     double *minf,
     double minf_max, double ftol_rel, double ftol_abs,
     double xtol_rel, const double *xtol_abs,
     int maxeval, double maxtime) __attribute__((deprecated));

extern NcmFitNloptResult nlopt_minimize_econstrained(
     NcmFitNloptAlgorithm algorithm,
     int n, nlopt_func_old f, void *f_data,
     int m, nlopt_func_old fc, void *fc_data, ptrdiff_t fc_datum_size,
     int p, nlopt_func_old h, void *h_data, ptrdiff_t h_datum_size,
     const double *lb, const double *ub,
     double *x,
     double *minf,
     double minf_max, double ftol_rel, double ftol_abs,
     double xtol_rel, const double *xtol_abs,
     double htol_rel, double htol_abs,
     int maxeval, double maxtime) __attribute__((deprecated));

extern void nlopt_get_local_search_algorithm(NcmFitNloptAlgorithm *deriv,
          NcmFitNloptAlgorithm *nonderiv,
          int *maxeval) __attribute__((deprecated));
extern void nlopt_set_local_search_algorithm(NcmFitNloptAlgorithm deriv,
          NcmFitNloptAlgorithm nonderiv,
          int maxeval) __attribute__((deprecated));

extern int nlopt_get_stochastic_population(void) __attribute__((deprecated));
extern void nlopt_set_stochastic_population(int pop) __attribute__((deprecated));
