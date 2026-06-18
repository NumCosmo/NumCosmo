/*************************************************************************/
/* Copyright (C) 1996-1997 Marlis Hochbruck                              */
/* All rights reserved.                                                  */
/*                                                                       */
/* C code written by Mathias Froehlich.                                  */
/*                                                                       */
/* This code is part of a copyrighted package. For details, see the file */
/* "COPYING" in the top-level directory.                                 */
/*************************************************************************/
#define GAMMA_UP(i)    schur_par[i]
#define ALPHA_LEF(i)   schur_par[(i)+N]
#define BETA_LEF(i)    schur_par[(i)+2*N]
#define ALPHA_UP(i)    schur_par[(i)+3*N]
#define BETA_UP(i)     schur_par[(i)+4*N]
#define ALPHA(i)       ALPHA_LEF((i)-1)
#define BETA(i)        GAMMA_UP((i)-1)
#define CR_GAMMA_UP(i) BETA_UP(i)
#define CR_GAMMA(i)    ALPHA_UP(i)

extern int update_x(fdouble *work,finteger lwork,polynom *q,polynom *q_up,
      finteger *ldq,fdouble *b,fdouble *x,fdouble *d__,finteger *ldd,
      finteger *d_perm,fdouble *d_tau,finteger n,finteger k);

extern fdouble get_singularity(fdouble *work,finteger lwork,fdouble *m,
      finteger *ldm,finteger k);

extern void update_d(fdouble *d__,finteger *ldd,fdouble *q__,finteger *ldq,
      fdouble *e_up__,finteger *lde,finteger n,finteger k);

extern void update_d_schur(fdouble *d__,finteger *ldd,fdouble *e_up,
      finteger *lde,fdouble *p,finteger *ldp,finteger k);

extern int update_y(fdouble *work,finteger lwork,fdouble *b,fdouble *y,
      fdouble *e_up__,finteger *lde,finteger N,finteger n,finteger k);

extern int update_Dy_schur(fdouble *work,finteger iwork,fdouble *b,
      fdouble *y,fdouble *e_up__,fdouble *d,finteger *ldd,finteger N,
      finteger n,finteger k,finteger *ldq);

extern int update_y_schur(fdouble *work,finteger iwork,fdouble *b,fdouble *y,
      fdouble *e_up__,finteger N,finteger n,finteger k,finteger *ldq);

extern int update_x_from_y(polynom *q_up,fdouble *y,fdouble *x,
      finteger n,finteger k,finteger *ldq);

extern int update_x_from_Dy_schur(fdouble *work,finteger lwork,fdouble *L__,
      fdouble *y,fdouble *x,finteger N,finteger *reg_index,
      finteger reg_pos,finteger *ldq);

extern int update_x_from_y_schur(fdouble *work,finteger lwork,
      fdouble *schur_par,fdouble *y,fdouble *x,finteger N,finteger *reg_index,
      finteger reg_pos,finteger *ldq);

extern void get_M(fdouble *M__,finteger *ldM,fdouble *e_lef,
      fdouble *e_up,fdouble *p_lef,fdouble *p_up,finteger k);

extern void get_colreg_pair(polynom *res,polynom *res_up,toeplitz *h,
      polynom *q,polynom *q_up,fdouble *e,fdouble *e_up,fdouble *p,
      fdouble *p_up,fdouble *v,fdouble *v_up,finteger n);

extern void get_colreg_pair_schurlev(polynom *res,polynom *res_up,toeplitz *h,
      polynom *q,polynom *q_up,fdouble *e,fdouble *e_up,fdouble *eres,
      fdouble *eres_up,fdouble *p,fdouble *p_up,fdouble *pres,fdouble *pres_up,
      finteger *ldsmall,fdouble *v,fdouble *v_up,finteger n);

extern void get_inner_q_up(polynom *res_up,polynom *q_lef,polynom *q_up,
      fdouble *v_up,finteger k);

extern void get_inner_x_up_schur(fdouble *eres_up,fdouble *e_lef,
      fdouble *e_up,fdouble *pres_up,fdouble *p_lef,fdouble *p_up,
      fdouble *v_up,finteger N,finteger n,finteger k);

extern void get_inner_q_up_schurlev(polynom *res_up,polynom *q_lef,
      polynom *q_up,fdouble *eres_up,fdouble *e_lef,fdouble *e_up,
      fdouble *pres_up,fdouble *p_lef,fdouble *p_up,
      fdouble *v_up,finteger N,finteger n,finteger k);

extern void get_inner_q(polynom *res,polynom *q_lef,polynom *q,fdouble *v,
      finteger k);

extern void get_inner_x_schur(fdouble *eres,fdouble *e_lef,fdouble *e,
      fdouble *pres,fdouble *p_lef,fdouble *p,
      fdouble *v,finteger N,finteger n,finteger k);

extern void get_inner_q_schurlev(polynom *res,polynom *q_lef,polynom *q,
      fdouble *eres,fdouble *e_lef,fdouble *e,
      fdouble *pres,fdouble *p_lef,fdouble *p,
      fdouble *v,finteger N,finteger n,finteger k);

extern void get_new_coefs(fdouble *e__,fdouble *e_up__,fdouble *p__,
      fdouble *p_up__,fdouble *p_lef,fdouble *e_lef,finteger *ldsmall,
      toeplitz *h,polynom *q,polynom *q_up,polynom *q_lef,fdouble *v,
      fdouble *v_up,finteger k,finteger n);

extern void get_new_e_pup_coefs(fdouble *e__,fdouble *p_up__,fdouble *p_lef,
      fdouble *e_lef,finteger *ldsmall,fdouble *v,fdouble *v_up,finteger k,
      finteger n);

extern void get_anti_reg_coefs(fdouble *e__,fdouble *e_up__,
      fdouble *p__,fdouble *p_up__,fdouble *p_lef,fdouble *e_lef,
      finteger *ldsmall,toeplitz *h,polynom *q,polynom *q_up,polynom *q_lef,
      fdouble *v,fdouble *v_up,finteger k,finteger n);

extern int get_anti_regular_up(fdouble *work,finteger lwork,polynom *q_lef,
      polynom *q_up,polynom *q,fdouble *e_up__,fdouble *e_lef,
      fdouble *p_up__,fdouble *d__,finteger *d_perm,fdouble *d_tau,
      fdouble *p_lef,finteger *ldsmall,finteger *ldq,finteger n,finteger k);

extern int get_anti_regular_up_schurlev(fdouble *work,finteger lwork,
      fdouble *theta_,
      polynom *q_lef,polynom *q_up,polynom *q,fdouble *e_up__,
      fdouble *e_lef,fdouble *p_up__,fdouble *p_lef,fdouble *d__,
      finteger *d_perm,fdouble *d_tau,finteger *ldsmall,finteger *ldq,
      finteger N,finteger n,finteger k);

extern int get_anti_regular_lef(fdouble *work,finteger lwork,polynom *q,
      polynom *q_lef,fdouble *e__,fdouble *e_lef,fdouble *p__,
      fdouble *d__,finteger *d_perm,fdouble *d_tau,finteger *ldsmall,
      finteger *ldq,finteger n,finteger k);

extern int get_anti_regular_lef_schurlev(fdouble *work,finteger lwork,
      fdouble *theta,
      polynom *q,polynom *q_lef,fdouble *e__,fdouble *e_lef,fdouble *p__,
      fdouble *p_lef,fdouble *d__,finteger *d_perm,fdouble *d_tau,
      finteger *ldsmall,finteger *ldq,finteger N,finteger n,finteger k);

extern void get_regular_q(polynom *q,polynom *q_up,polynom *q_lef,
      fdouble *e_lef,fdouble *e_up);

extern void get_regular_q_schurlev(polynom *q,polynom *q_up,polynom *q_lef,
      fdouble *e,fdouble *e_lef,fdouble *e_up,fdouble *p,
      fdouble *p_lef,fdouble *p_up,finteger N,finteger n,finteger k);

extern fdouble get_e_lef_k(polynom* q_lef, polynom* q_up, fdouble* e_lef,
      fdouble* e_up__, finteger *ldsmall, finteger k);

extern fdouble get_p_lef_km1(polynom* q_lef, polynom* q_up,
      fdouble* p_lef,fdouble* p_up__, finteger *ldsmall, finteger k);

extern int get_anti_regular_pair(fdouble *work,finteger lwork,fdouble *M,
      finteger *ldM,polynom *q_lef,polynom *q_up,finteger *ldq,fdouble *e_up,
      fdouble *e_lef,fdouble *p_up,fdouble *p_lef,finteger n,
      finteger k);

extern int get_anti_regular_pair_schurlev(fdouble *work,finteger lwork,
      fdouble *M,finteger *ldM,polynom *q_lef,polynom *q_up,finteger *ldq,
      fdouble *e_up,fdouble *e_lef,fdouble *p_up,fdouble *p_lef,
      fdouble *schur_par,finteger N,finteger n,finteger k);

extern void get_alpha_beta_up_from_theta(fdouble *schur_par,fdouble *theta,
      finteger N,finteger n,finteger k);

extern void get_alpha_beta_lef_from_theta(fdouble *schur_par,fdouble *theta,
      fdouble *v,finteger N,finteger n,finteger k);

extern int decomp_matrix(fdouble *lawork,finteger *k,fdouble *mat,
      finteger *ldmat,finteger *perm,fdouble *tau);
