#include <stdbool.h>

#ifndef PLENSLIKE
#define PLENSLIKE

#define min(X,Y) ((X) < (Y) ? (X) : (Y))
#define max(X,Y) ((X) > (Y) ? (X) : (Y))
#define min3(X,Y,Z) ( min( X, min(Y, Z) ) )

#ifndef perr
#define perr(...) do { \
    printf ("@ %s (%d): ", __FILE__, __LINE__); \
    printf (__VA_ARGS__); \
} while (0)
#endif


typedef struct {
  int     ntrm;
  int     lmax;
  int     **s12L;
  double ***w12L;
} qest;

typedef struct {
  int     nbins;
  int     lmax;
  int    *bin_lmins;
  int    *bin_lmaxs;
  double *bin_vals;
  double *mat_sigma;
  double *mat_sigma_inv;
  double *clpp_fid;
  double *cltt_fid;
  double *bl_fid;
  double *fl;
  double *vl_inv;
  double *al_inv;
} plenslike_dat_mono;

typedef struct {
  int     nbins;
  int     lmax;
  int     lmaxt;
  int     lmax1;
  int     lmax2;
  int     lmax3;
  int     lmax4;
  double  s4hat;
  double  s4std;
  int    *bin_lmins;
  int    *bin_lmaxs;
  double *bin_vals;
  double *mat_sigma;
  double *mat_sigma_inv;
  double *clpp_fid;
  double *vl_inv;
  double *rl_inv;
  double *sl_fid;
  double *cltt_fid;
  double *bl1n1_fid;
  double *bl2n1_fid;
  double *bl3n1_fid;
  double *bl4n1_fid;
  double *fl1;
  double *fl2;
  double *fl3;
  double *fl4;
  qest   *qe12;
  qest   *qe34;
} plenslike_dat_quad;

typedef struct {
  int     nbins;
  int     lmaxcmb;
  int     lmaxphi;
  int     nqe;
  int     nx;

  int    *bin_lmins;
  int    *bin_lmaxs;
  double *bin_vals;
  double *mat_sigma;
  double *mat_sigma_inv;

  double *cltt_fid;
  double *clee_fid;
  double *clte_fid;

  double *clpp_fid;
  double *qlpp_fid;
  double *vlpp_inv;
  double *qlpp_inv;

  qest  **qes;
  int    *qefs;
  int    *qe12;
  int    *qe34;
} plenslike_dat_qecl;

typedef struct {
  int     nbins;
  int     lmaxcmb;
  int     lmaxphi;
  int     n1lqbins;
  int     n1lpbins;
  int     nqe;
  int     nx;

  double  s4fid;
  double  s4std;

  int    *bin_lmins;
  int    *bin_lmaxs;
  double *bin_vals;
  double *mat_sigma;
  double *mat_sigma_inv;

  double *cltt_fid;
  double *clee_fid;
  double *clte_fid;

  double *clpp_fid;
  double *qlpp_fid;
  double *qlss_fid;
  double *n1pp_fid;
  double *vlpp_inv;
  double *qlpp_inv;

  double *n1lqs;
  double *n1lps;
  double *mat_n1;

  qest  **qes;
  int    *qefs;
  int    *qe12;
  int    *qe34;
} plenslike_dat_full;

// qest.c
void free_qe( qest *qe );

void fill_qe_resp( int lmax, double *resp, qest *qee, qest *qes, 
		   double *f1, int f1lmax, double *f2, int f2lmax );
void fill_qe_covn( int lmax, double *covn, qest *q12, qest *q34, 
		   double *c1, int c1lmax, double *c2, int c2lmax, 
		   bool switch_34, bool mult_n1_l1l2L );

void init_qe_plm( qest *qe, int lmax, double *cltt);
void init_qe_slm( qest *qe, int lmax );
void init_qe_mlm( qest *qe, int lmax, double *cltt, double *fbl1, double *fbl2 );
void init_qe_nlm( qest *qe, int lmax, double *fbl1, double *fbl2 );

void init_qe_plm_tt( qest *qe, int lmax, double *cltt);
void init_qe_plm_ee( qest *qe, int lmax, double *clee);
void init_qe_plm_te( qest *qe, int lmax, double *clte);
void init_qe_plm_tb( qest *qe, int lmax, double *clte);
void init_qe_plm_eb( qest *qe, int lmax, double *clee);

void switch_qe( qest *qe );
void init_qe_plm_et( qest *qe, int lmax, double *clte);
void init_qe_plm_bt( qest *qe, int lmax, double *clte);
void init_qe_plm_be( qest *qe, int lmax, double *clee);

void init_qls( qest **qls, int lmaxcmb, double *cltt, double *clee, double *clte );
void free_qls( qest **qls );

// wignerd.c
void init_gauss_legendre_quadrature( int n, double *x, double *w );
void wignerd_cf_from_cl( int s1, int s2, int nfunc, int ntheta, int lmax, const double *cos_theta, 
			 double *out_cf, const double *in_cl );
void wignerd_cl_from_cf( int s1, int s2, int nfunc, int ntheta, int lmax, const double *cos_theta, 
			 const double *integration_weights, double *out_cl, const double *in_cf );

// plenslike_dat_mono.c
void   load_plenslike_dat_mono( plenslike_dat_mono *dat, char *tfname );
void   free_plenslike_dat_mono( plenslike_dat_mono *dat );

void   fill_plenslike_mono_bins( plenslike_dat_mono *dat, double *bins, double *clpp );
void   fill_qe_plm_resp_plm_mono( int lmax, double *resp, double *cltt_fid, double *bl_fid, double *fl, double *cltt, double *bl);

double calc_plenslike_mono( plenslike_dat_mono *dat, double *clpp );
double calc_plenslike_mono_renorm( plenslike_dat_mono *dat, double *clpp, double *cltt, double *bl );

// plenslike_dat_quad.c
void   load_plenslike_dat_quad( plenslike_dat_quad *dat, char *tfname );
void   free_plenslike_dat_quad( plenslike_dat_quad *dat );

void   fill_plenslike_quad_bins( plenslike_dat_quad *dat, double *bins, double *clpp );

void   fill_quad_resp_pp_cltt( int lmax, double *rl, plenslike_dat_quad *dat, double *cltt );
void   fill_quad_resp_ss_bfid( int lmax, double *rl, plenslike_dat_quad *dat );
void   fill_quad_resp_qq_bfid( int lmax, double *rl, plenslike_dat_quad *dat, qest *qq );

double calc_plenslike_quad( plenslike_dat_quad *dat, double *clpp );
double calc_plenslike_quad_renorm_cltt( plenslike_dat_quad *dat, double *clpp, double *cltt );

// plenslike_dat_qecl.c
int  load_plenslike_dat_qecl( plenslike_dat_qecl *dat, char *tfname );
void free_plenslike_dat_qecl( plenslike_dat_qecl *dat );
void fill_plenslike_qecl_bins( plenslike_dat_qecl *dat, double *bins, double *clpp );

void fill_qecl_resp_pp_qls( int lmax, double *rl, plenslike_dat_qecl *dat, qest **qls );
void full_qecl_resp_pp_cls( int lmax, double *rl, plenslike_dat_qecl *dat, double *cltt, double *clee, double *clte);

double calc_plenslike_qecl( plenslike_dat_qecl *dat, double *clpp );
double calc_plenslike_qecl_renorm( plenslike_dat_qecl *dat, double *clpp, double *cltt, double *clee, double *clte );

// plenslike_dat_full.c
void fill_full_n1( int lmax, double *n1, plenslike_dat_full *dat, double *clpp );
int  load_plenslike_dat_full( plenslike_dat_full *dat, char *tfname );
void free_plenslike_dat_full( plenslike_dat_full *dat );
void fill_plenslike_full_bins( plenslike_dat_full *dat, double *bins, double *clpp );

void fill_full_resp_pp_qls( int lmax, double *rl, plenslike_dat_full *dat, qest **qls );
void full_full_resp_pp_cls( int lmax, double *rl, plenslike_dat_full *dat, double *cltt, double *clee, double *clte);
void fill_full_n1( int lmax, double *n1, plenslike_dat_full *dat, double *clpp );

double calc_plenslike_full( plenslike_dat_full *dat, double *clpp );
double calc_plenslike_full_renorm( plenslike_dat_full *dat, double *clpp, double *cltt, double *clee, double *clte );

double calc_plenslike_full_ren1( plenslike_dat_full *dat, double *clpp );
double calc_plenslike_full_renorm_ren1( plenslike_dat_full *dat, double *clpp, double *cltt, double *clee, double *clte );

#endif
