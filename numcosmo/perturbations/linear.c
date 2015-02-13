/***************************************************************************
 *            linear.c
 *
 *  Sat Oct 25 21:02:36 2008
 *  Copyright  2008  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2012 <sandro@lapsandro>
 * numcosmo is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * numcosmo is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * SECTION:linear
 * @title: NcLinearPert
 * @short_description: Linear perturbations.
 *
 * FIXME
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "perturbations/linear.h"
#include "math/grid_one.h"
#include "math/ncm_sf_sbessel_int.h"
#include "math/integral.h"
#include "math/ncm_spline_cubic_notaknot.h"
#include "math/memory_pool.h"

guint _itheta_table[3]   = {NC_PERT_THETA0, NC_PERT_THETA1, NC_PERT_THETA2};
guint _itheta_p_table[3] = {NC_PERT_THETA_P0, NC_PERT_THETA_P1, NC_PERT_THETA_P2};
guint _nc_default_los_init[] = {2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 14, 16, 18, 20, 0};
guint _nc_default_los_step[] = {5, 50, 10, 100, 20, 200, 25, 300, 50, 0};

/**
 * nc_pert_linear_create_los_table: (skip)
 * @lmax_los: FIXME
 * @los_ini: FIXME
 * @los_step: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
*/
GArray *
nc_pert_linear_create_los_table (guint lmax_los, guint *los_ini, guint *los_step)
{
  GArray *los_table = g_array_sized_new (FALSE, FALSE, sizeof (guint), 1000);
  guint last_i, i = 0;
  while (los_ini[i] != 0) i++;
  last_i = los_ini[i - 1];
  /* lmax_los_size = i; */

  g_assert (lmax_los > 30);

  g_array_append_vals (los_table, los_ini, i);

  if (last_i < lmax_los)
  {
    guint j;
    i = 0;
    while (last_i < lmax_los)
    {
      if (los_step[i + 1] == 0)
      {
        /* lmax_los_size += (lmax_los - last_i) / los_step[i]; */
        for (j = last_i + los_step[i]; j <= lmax_los; j += los_step[i])
          g_array_append_val (los_table, j);
        break;
      }
      else
      {
        guint stest = GSL_MIN(los_step[i + 1], lmax_los);
        /* lmax_los_size += (stest - last_i) / los_step[i]; */
        for (j = last_i + los_step[i]; j <= stest; j += los_step[i])
          g_array_append_val (los_table, j);

        last_i = stest;
      }
      i += 2;
    }
  }
  return los_table;
}

/**
 * nc_pert_linear_workspace_new: (skip)
 * @pert: a #NcLinearPert
 *
 * FIXME
 *
 * Returns: FIXME
*/
NcLinearPertWorkSpace *
nc_pert_linear_workspace_new (NcLinearPert *pert)
{
  NcLinearPertWorkSpace *pws = g_slice_new (NcLinearPertWorkSpace);
  pws->pert = pert;
  pws->tight_coupling = FALSE;
  pws->tight_coupling_end = FALSE;
  pws->k = 0.0;

  return pws;
}

/**
 * nc_pert_linear_new: (skip)
 * @cosmo: a #NcHICosmo.
 * @recomb: FIXME
 * @lmax: FIXME
 * @tc_reltol: FIXME
 * @reltol: FIXME
 * @tc_abstol: FIXME
 * @abstol: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
NcLinearPert *
nc_pert_linear_new (NcHICosmo *cosmo, NcRecomb *recomb, guint lmax, gdouble tc_reltol, gdouble reltol, gdouble tc_abstol, gdouble abstol)
{
  NcLinearPert *pert = g_slice_new (NcLinearPert);

  g_assert (lmax > 2);

  pert->cosmo = cosmo;
  pert->recomb = nc_recomb_ref (recomb);
  pert->a = nc_scale_factor_new (NC_TIME_CONFORMAL, NC_PERTURBATION_START_X - 1.0);

  nc_recomb_prepare_if_needed (pert->recomb, cosmo);
  nc_scale_factor_prepare_if_needed (pert->a, cosmo);

  pert->lambdai = NC_PERTURBATIONS_X2LAMBDA (NC_PERTURBATION_START_X);
  pert->lambdaf = NC_PERTURBATIONS_X2LAMBDA (1.0);

  nc_recomb_v_tau_lambda_features (pert->recomb, cosmo, 2.0 * M_LN10, 
                                   &pert->lambda_rec, 
                                   &pert->lambda_rec_10m2_max[0], 
                                   &pert->lambda_rec_10m2_max[1]);
  pert->lambda_opt_cutoff = nc_recomb_tau_cutoff (pert->recomb, cosmo);
  
  pert->eta0 = nc_scale_factor_t_z (pert->a, 0.0);

  pert->lmax = lmax;

  pert->reltol = reltol;
  pert->abstol = abstol;
  pert->tc_reltol = tc_reltol;
  pert->tc_abstol = tc_abstol;

  pert->sys_size = (NC_PERT_THETA_P2 + 1) + 2 * (lmax + 1 - 3);

  pert->solver = g_slice_new (NcLinearPertOdeSolver);
  pert->solver = cvodes_solver;
	//pert->solver = ncm_gsl_odeiv2_solver;

  pert->solver->data = pert->solver->create (pert);

  pert->pws = nc_pert_linear_workspace_new (pert);

  printf ("# lmax %d | eta0 %g\n", lmax, pert->eta0);

  return pert;
}

/**
 * nc_pert_linear_free:
 * @pert: a #NcLinearPert.
 *
 * FIXME
 *
*/
void
nc_pert_linear_free (NcLinearPert *pert)
{
  NCM_UNUSED (pert);
  g_assert_not_reached ();
}

/**
 * nc_pert_linear_clear:
 * @pert: a #NcLinearPert.
 *
 * FIXME
 *
*/
void
nc_pert_linear_clear (NcLinearPert **pert)
{
  g_clear_object (pert);
}

/**
 * nc_pert_linear_splines_free:
 * @pspline: a #NcLinearPertSplines.
 *
 * FIXME
 *
*/
void 
nc_pert_linear_splines_free (NcLinearPertSplines *pspline)
{
  NCM_UNUSED (pspline);
  g_assert_not_reached ();
}

/**
 * nc_pert_linear_splines_clear:
 * @pspline: a #NcLinearPertSplines.
 *
 * FIXME
 *
*/
void 
nc_pert_linear_splines_clear (NcLinearPertSplines **pspline)
{
  nc_pert_linear_splines_free (*pspline);
  *pspline = NULL;
}

typedef struct _local_los_int_data
{
  gint l;
  NcmSpline *S0;
  NcmSpline *S1;
  NcmSpline *S2;
  gint count;
  gdouble k;
} local_los_int_data;

static gdouble local_los_int (gdouble deta, gpointer params);

typedef struct _Nc_k_integrand_data
{
  glong l;
  NcLinearPertSplines *pspline;
  glong count;
} Nc_k_integrand_data;

static gdouble
_Nc_k_integrand (gdouble k, gpointer params)
{
  Nc_k_integrand_data *Nc_data = (Nc_k_integrand_data *)params;
  NcLinearPertSplines *pspline = Nc_data->pspline;
  const glong l = Nc_data->l;
  glong nppo = 10;
  gdouble xf = 5250.0;
  glong np = ceil(xf / (2.0 * M_PI) * nppo);
  NcmGridSection x_sec[2] = {{NCM_GRID_NODES_BOTH,  0, np, 0.0, xf}, {0}};
  NcmSFSphericalBesselIntSpline *int_xnjl_spline = ncm_sf_sbessel_jl_xj_integrate_spline_cached_new (l, x_sec, TRUE);
  gdouble theta_l_k;
//printf ("# Uhu![%ld] %.15g\n", l, k);
  ncm_sf_sbessel_jl_xj_integrate_spline_goto (int_xnjl_spline, l);
  nc_pert_linear_spline_set_source_at (pspline, k);
  theta_l_k = ncm_sf_sbessel_jl_xj_integral_spline (int_xnjl_spline, pspline->Sg[0], pspline->Sg[1], pspline->Sg[2], k);
  Nc_data->count++;
  return theta_l_k * theta_l_k / k;
}

/**
 * nc_pert_linear_los_integrate:
 * @pspline: a #NcLinearPertSplines
 * @l: FIXME
 * @k: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
*/
gdouble
nc_pert_linear_los_integrate (NcLinearPertSplines *pspline, glong l, gdouble k)
{
  static gsl_integration_workspace **w = NULL;
  gsl_function F;
  local_los_int_data ldata;
  gdouble res_int = 0.0, err;
  gdouble res_rules = 0.0;
//  gint i;

  NcmSFSphericalBesselIntSpline *int_xnjl_spline;

  {
    glong nppo = 10;
    gdouble xf = 5250.0;
    glong np = ceil(xf / (2.0 * M_PI) * nppo);
    NcmGridSection x_sec[2] = {{NCM_GRID_NODES_BOTH,  0, np, 0.0, xf}, {0}};
    int_xnjl_spline = ncm_sf_sbessel_jl_xj_integrate_spline_cached_new (l, x_sec, TRUE);
    ncm_sf_sbessel_jl_xj_integrate_spline_goto (int_xnjl_spline, l);
    //printf ("# Rules x_sec [%.15g %.15g] / %ld\n", 0.0, xf, np);
  }
  nc_pert_linear_spline_set_source_at (pspline, k);

  return ncm_sf_sbessel_jl_xj_integral_spline (int_xnjl_spline, pspline->Sg[0], pspline->Sg[1], pspline->Sg[2], k);

  if (w == NULL)
    w = ncm_integral_get_workspace();

  F.params = &ldata;
  F.function = &local_los_int;

  ldata.k = k;
  ldata.l = l;
  ldata.S0 = pspline->Sg[0];
  ldata.S1 = pspline->Sg[1];
  ldata.S2 = pspline->Sg[2];
  ldata.count = 0;
/*
  for (i = 0; i < ldata.S0->len - 1; i++)
  {
    const gdouble aa = ncm_mpsf_sbessel_integrate (int_jlspline, pspline->Sg[0], l, ki, i, 0);
    const gdouble bb = ncm_mpsf_sbessel_integrate (int_jlspline, pspline->Sg[1], l, ki, i, 1);
    const gdouble cc = ncm_mpsf_sbessel_integrate (int_jlspline, pspline->Sg[2], l, ki, i, 2);
    res_rules += aa + bb + cc;
  }
*/
  if (FALSE)
  {
    gdouble deta_opt = ncm_vector_get (ldata.S0->xv, ncm_vector_len(ldata.S0->xv) - 1);
    gsl_integration_qag (&F, 0.0, deta_opt, 0.0, 1e-13, NCM_INTEGRAL_PARTITION, 1, *w, &res_int, &err);fflush (stdout);
    printf ("# RES % .15e % .15e %.4e %.4e %08d\n",
      res_int, res_rules, fabs((res_rules - res_int) / res_int), fabs(err / res_int), ldata.count);fflush (stdout);
  }

//  printf ("%.15g %.15g\n", res_rules, ncm_sf_sbessel_jl_xj_integral_spline (int_xnjl_spline, pspline->Sg[0], pspline->Sg[1], pspline->Sg[2], k));

  return res_rules;
}

/**
 * nc_pert_transfer_function_new: (skip)
 * @pert: a #NcLinearPert
 * @k0: FIXME
 * @k1: FIXME
 * @np: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
*/
NcLinearPertTF *
nc_pert_transfer_function_new (NcLinearPert *pert, gdouble k0, gdouble k1, gulong np)
{
  NcLinearPertTF *perttf = g_slice_new (NcLinearPertTF);
  NcmVector *logkv = ncm_vector_new (np);
  NcmVector *logPhiv = ncm_vector_new (np);

  perttf->logk0 = log (k0 * ncm_c_hubble_radius ());
  perttf->logk1 = log (k1 * ncm_c_hubble_radius ());
  perttf->np = np;
  perttf->logPhi_logk = ncm_spline_cubic_notaknot_new_full (logkv, logPhiv, FALSE);
  perttf->pert = pert;

  ncm_vector_free (logkv);
  ncm_vector_free (logPhiv);
  
  return perttf;
}

/**
 * nc_pert_transfer_function_prepare:
 * @perttf: a #NcLinearPertTF
 *
 * FIXME
 *
*/
void
nc_pert_transfer_function_prepare (NcLinearPertTF *perttf)
{
  guint i;

  for (i = 0; i < perttf->np; i++)
  {
    gdouble logk = perttf->logk0 + (perttf->logk1 - perttf->logk0) / (perttf->np - 1.0) * i;
    gdouble logPhi;
    perttf->pert->pws->k = exp (logk);
    perttf->pert->solver->init (perttf->pert);
    perttf->pert->solver->evol (perttf->pert, perttf->pert->lambdaf);
    //perttf->pert->solver->evol (perttf->pert, perttf->pert->g_opt_cutoff);

    logPhi = log (perttf->pert->solver->get (perttf->pert, NC_PERT_PHI));
		ncm_vector_set (perttf->logPhi_logk->xv, i, logk);
		ncm_vector_set (perttf->logPhi_logk->yv, i, logPhi);
    printf ("# Calculated mode % 20.15g\n", perttf->pert->pws->k);
    fflush (stdout);
  }

  ncm_spline_prepare (perttf->logPhi_logk);
}

/**
 * nc_pert_transfer_function_get:
 * @perttf: a #NcLinearPertTF
 * @kh: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
*/
gdouble
nc_pert_transfer_function_get (NcLinearPertTF *perttf, gdouble kh)
{
  const gdouble logk = log (kh * ncm_c_hubble_radius ());
  const gdouble logPhi = ncm_spline_eval (perttf->logPhi_logk, logk);
  return exp (logPhi);
}

static void
_nc_pert_linear_prepare_deta_grid (NcLinearPert *pert, gdouble *deta_grid, glong ni, glong nc, glong nt)
{
  const gdouble x_opt_cutoff    = NC_PERTURBATIONS_LAMBDA2X (pert->lambda_opt_cutoff);
  const gdouble x_rec_10m2_max0 = NC_PERTURBATIONS_LAMBDA2X (pert->lambda_rec_10m2_max[0]);
  const gdouble x_rec_10m2_max1 = NC_PERTURBATIONS_LAMBDA2X (pert->lambda_rec_10m2_max[1]);
  const glong n = nc + nt + ni;
  gdouble last_eta;
  gint i;

  for (i = 0; i < ni; i++)
  {
    gdouble x_i = x_opt_cutoff + (x_rec_10m2_max1 - x_opt_cutoff) / (ni - 1.0) * i;
    gdouble deta_i = pert->eta0 - nc_scale_factor_t_x (pert->a, x_i);
    deta_grid[n - 1 - i] = deta_i;
  }

  for (i = 0; i < nc; i++)
  {
    gdouble x_i = x_rec_10m2_max1 + (x_rec_10m2_max0 - x_rec_10m2_max1) / (nc * 1.0) * (i + 1.0);
    gdouble deta_i = pert->eta0 - nc_scale_factor_t_x (pert->a, x_i);
    deta_grid[n - 1 - (i + ni)] = deta_i;
  }
  last_eta = deta_grid[n - nc - ni];

  for (i = 0; i < nt; i++)
  {
    gdouble deta_i = pert->eta0 - (last_eta + (pert->eta0 - last_eta) / (nt * 1.0) * (i + 1.0));
    deta_grid[n - 1 - (i + nc + ni)] = deta_i;
  }
}

static void
_nc_pert_linear_prepare_k_grid (NcLinearPert *pert, gdouble *k_grid, glong ni, glong nt, gdouble k0, gdouble k1)
{
  const gdouble log_k0 = log(k0);
  gint i;

  NCM_UNUSED (pert);
  
  for (i = 0; i < ni; i++)
    k_grid[i] = k0 * exp ((M_LN10 - log_k0) / (ni - 1.0) * i);
  for (i = 0; i < nt; i++)
    k_grid[ni + i] = 10.0 + (k1 - 10.0) / (nt * 1.0) * (i + 1.0);
}

/**
 * nc_pert_linear_splines_new: (skip)
 * @pert: a #NcLinearPert
 * @types: a #NcLinearPertSplineTypes
 * @n_deta: FIXME
 * @n_evol: FIXME
 * @k0: FIXME
 * @k1: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
*/
NcLinearPertSplines *
nc_pert_linear_splines_new (NcLinearPert *pert, NcLinearPertSplineTypes types, gulong n_deta, gulong n_evol, gdouble k0, gdouble k1)
{
  NcLinearPertSplines *pspline = g_slice_new (NcLinearPertSplines);
  gint i;
  pspline->pert = pert;
  pspline->n_deta = n_deta;
  pspline->n_evol = n_evol;
  pspline->types = types;
  pspline->k0 = k0;
  pspline->k1 = k1;

  pspline->ga = ncm_vector_new (n_deta);
  pspline->ka = ncm_vector_new (n_evol);

	for (i = 0; i < 3; i++)
  {
    guint j;

		pspline->Sg_data[i] = ncm_vector_new (n_deta);

    pspline->Sg[i] = ncm_spline_cubic_notaknot_new_full (pspline->ga, pspline->Sg_data[i], FALSE);

    pspline->Sk_data[i] = ncm_matrix_new (n_deta, n_evol);

		pspline->Sk[i] = g_slice_alloc (sizeof (NcmSpline *) * n_deta);
    for (j = 0; j < n_deta; j++)
    {
      NcmVector *Sk_row_ij = ncm_matrix_get_row (pspline->Sk_data[i], j);
      pspline->Sk[i][j] = ncm_spline_cubic_notaknot_new_full (pspline->ka, Sk_row_ij, FALSE);
      ncm_vector_free (Sk_row_ij);
    }
  }

  for (i = 0; i < NC_PERTURBATION_BASE_SIZE; i++)
  {
    pspline->sdata[i] = ncm_vector_new (n_evol);
    pspline->s[i] = ncm_spline_cubic_notaknot_new_full (pspline->ka, pspline->sdata[i], FALSE);
  }

  pspline->Nc = NULL;

  g_assert (n_deta >= 60);
  g_assert (n_evol >= 220);

  _nc_pert_linear_prepare_deta_grid (pert, ncm_vector_ptr (pspline->ga, 0), 20, n_deta - 40, 20);
  _nc_pert_linear_prepare_k_grid (pert, ncm_vector_ptr (pspline->ka, 0), 20, n_evol - 20, k0, k1);

  return pspline;
}

/**
 * nc_pert_linear_prepare_splines:
 * @pspline: a #NcLinearPertSplines
 *
 * FIXME
 *
*/
void
nc_pert_linear_prepare_splines (NcLinearPertSplines *pspline)
{
  NcLinearPert *pert = pspline->pert;
  guint i;
  printf ("# Evoling [%lu] modes in: [%.5e %.5e]\n#\n", pspline->n_evol, pspline->k0, pspline->k1);

  for (i = 0; i < pspline->n_evol; i++)
  {
    guint j;
    pert->pws->k = ncm_vector_get (pspline->ka, i);
    pert->solver->init (pert);
    pert->solver->evol (pert, pert->lambda_opt_cutoff);

    if (pspline->types & NC_LINEAR_PERTURBATIONS_SPLINE_SOURCES)
    {
      for (j = 0; j < pspline->n_deta; j++)
      {
        gdouble S0, S1, S2;
        gint jj = pspline->n_deta - 1 - j;
        gdouble detaj = ncm_vector_get (pspline->ga, jj);
        gdouble xj = nc_scale_factor_z_t (pert->a, pert->eta0 - detaj) + 1.0;
        gdouble lambdaj = -log(xj);

        pert->solver->evol (pert, lambdaj);
        pert->solver->get_sources (pert, &S0, &S1, &S2);
				ncm_matrix_set (pspline->Sk_data[0], jj, i, S0);
				ncm_matrix_set (pspline->Sk_data[1], jj, i, S1);
				ncm_matrix_set (pspline->Sk_data[2], jj, i, S2);
      }
    }
    else
      pert->solver->evol (pert, pert->lambdaf);

    for (j = 0; j < NC_PERTURBATION_BASE_SIZE; j++)
      if (pspline->types & (1 << j))
        ncm_vector_set (pspline->sdata[j], i, pert->solver->get (pert, j));
    printf ("# integrated mode % 20.15g\n", pert->pws->k);fflush(stdout);
  }

  if (pspline->types & NC_LINEAR_PERTURBATIONS_SPLINE_SOURCES)
  {
    for (i = 0; i < pspline->n_deta; i++)
    {
      ncm_spline_prepare (pspline->Sk[0][i]);
      ncm_spline_prepare (pspline->Sk[1][i]);
      ncm_spline_prepare (pspline->Sk[2][i]);
    }
  }

  for (i = 0; i < NC_PERTURBATION_BASE_SIZE; i++)
    if (pspline->types & (1 << i))
      ncm_spline_prepare (pspline->s[i]);

  printf ("# All splines ready.\n");
}

/**
 * nc_pert_linear_spline_set_source_at:
 * @pspline: a #NcLinearPertSplines
 * @k: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
*/
gboolean
nc_pert_linear_spline_set_source_at (NcLinearPertSplines *pspline, gdouble k)
{
  guint i;
  //printf ("# Setting source at %.15g [%.15g %.15g]\n", k, pspline->k0, pspline->k1);fflush(stdout);
  g_assert (pspline->types & NC_LINEAR_PERTURBATIONS_SPLINE_SOURCES);
  g_assert ((k/pspline->k0 - 1.0 >= -1e-13) && (k/pspline->k1 - 1.0 <= 1e-13));

  for (i = 0; i < ncm_vector_len (pspline->ga); i++)
  {
		ncm_vector_set (pspline->Sg_data[0], i, ncm_spline_eval (pspline->Sk[0][i], k));
		ncm_vector_set (pspline->Sg_data[1], i, ncm_spline_eval (pspline->Sk[1][i], k));
		ncm_vector_set (pspline->Sg_data[2], i, ncm_spline_eval (pspline->Sk[2][i], k));
  }

  ncm_spline_prepare (pspline->Sg[0]);
  ncm_spline_prepare (pspline->Sg[1]);
  ncm_spline_prepare (pspline->Sg[2]);

  return TRUE;
}

/**
 * nc_pert_linear_calc_Nc_spline: (skip)
 * @pspline: a #NcLinearPertSplines
 * @pw_spline: a #NcmSpline
 * @los_table: FIXME
 * @n_interp: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
*/
gboolean
nc_pert_linear_calc_Nc_spline (NcLinearPertSplines *pspline, NcmSpline *pw_spline, GArray *los_table, gulong n_interp)
{
  guint i;
  g_assert (pspline->types & NC_LINEAR_PERTURBATIONS_SPLINE_SOURCES);
  NcmVector *la, *Nca;

  NCM_UNUSED (pw_spline);
  
  if (pspline->Nc == NULL)
  {
    la = ncm_vector_new (200);
    Nca = ncm_vector_new (200);
    pspline->Nc = ncm_spline_cubic_notaknot_new_full (la, Nca, FALSE);
    ncm_vector_free (la);
    ncm_vector_free (Nca);
  }
  else
  {
    la = pspline->Nc->xv;
    Nca = pspline->Nc->yv;
		g_assert (ncm_vector_len (la) == los_table->len);
		g_assert (ncm_vector_len (Nca) == los_table->len);
  }

  for (i = 0; i < los_table->len; i++)
  {
    guint j;
    //gint l = g_array_index (los_table, guint, los_table->len - 1 - i);
    gint l = g_array_index (los_table, guint, i);
    gdouble tot_res;
    Nc_k_integrand_data Nc_data = {l, pspline, 0};
    gsl_function F = {&_Nc_k_integrand, &Nc_data};

    printf ("# Integrating Nc[%d]\n", l);fflush (stdout);

    if (TRUE)
    {
      gsl_integration_workspace **w = ncm_integral_get_workspace();
      gdouble error;
      gsl_integration_qag (&F,
                           pspline->k0,
                           pspline->k1, 0.0, 1e-11, NCM_INTEGRAL_PARTITION, 6, *w, ncm_vector_ptr (Nca, i), &error);
      ncm_memory_pool_return (w);
      ncm_vector_set (la, i, l);
      printf ("# Integrate with %ld evals\n", Nc_data.count);
      //printf ("%d %.15g\n", l, g_array_index (Nca, gdouble, i));
      //continue;
    }

#ifdef HAVE_GSL_GLF
		{
      gsl_integration_glfixed_table *glt = gsl_integration_glfixed_table_alloc (9);
      gdouble logki = log (pspline->k0);
      tot_res = 0.0;

      for (j = 1; j < n_interp; j++)
      {
        gdouble piece = 0.0;
        gint m;
        const gdouble logkf = log (pspline->k0) + log (pspline->k1 / pspline->k0) / (n_interp - 1.0) * j;
        const gdouble logkfki_2 = (logkf + logki) / 2.0;
        const gdouble logkf_ki_2 = (logkf - logki) / 2.0;
        {
          const gdouble thetal0 = nc_pert_linear_los_integrate (pspline, l, exp (logkfki_2));
          piece += thetal0 * thetal0 * glt->w[0];
        }
        for (m = 1; m < 5; m++)
        {
          const gdouble logkp = logkfki_2 + logkf_ki_2 * glt->x[m];
          const gdouble logkm = logkfki_2 - logkf_ki_2 * glt->x[m];
          const gdouble thetalp = nc_pert_linear_los_integrate (pspline, l, exp (logkp));
          const gdouble thetalm = nc_pert_linear_los_integrate (pspline, l, exp (logkm));
          piece += (thetalp * thetalp + thetalm * thetalm) * glt->w[m];
        }
        piece *= logkf_ki_2;
        logki = logkf;
        tot_res += piece;
      }
      printf ("%d %.15g %.15g %.5e\n", l, tot_res, ncm_vector_get (Nca, i), fabs((tot_res - ncm_vector_get (Nca, i)) / ncm_vector_get (Nca, i)));
      ncm_vector_set (la, i, l);
      ncm_vector_set (Nca, i, tot_res);
      printf ("\n\n");
    }
#endif /* HAVE_GSL_GLF */
    //printf ("\n\n");
  }
  ncm_spline_prepare (pspline->Nc);

  return TRUE;
}

static gdouble
local_los_int (gdouble deta, gpointer params)
{
  local_los_int_data *ldata = (local_los_int_data *)params;
  const glong l = ldata->l;
  gdouble kdeta = ldata->k * deta;
  gdouble kdeta2 = kdeta * kdeta;
  gdouble jl, jlp1, djl, d2jl;

  jl = ncm_sf_sbessel (l, kdeta);
  jlp1 = ncm_sf_sbessel (l + 1, kdeta);

  djl = (kdeta != 0) ? (l * jl / kdeta - jlp1) : (0.0);
  d2jl = (kdeta != 0) ? (((l * (l-1.0) - kdeta2) * jl + 2.0 * kdeta * jlp1) / kdeta2) : (0.0);

  ldata->count++;
  return
    (
      ncm_spline_eval (ldata->S0, deta) * jl +
      ncm_spline_eval (ldata->S1, deta) * djl +
      ncm_spline_eval (ldata->S2, deta) * d2jl
      );
}
