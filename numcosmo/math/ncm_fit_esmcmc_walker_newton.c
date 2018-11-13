/***************************************************************************
 *            ncm_fit_esmcmc_walker_aps.c
 *
 *  Sat October 27 13:08:13 2018
 *  Copyright  2018  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * ncm_fit_esmcmc_walker_aps.c
 * Copyright (C) 2018 Sandro Dias Pinto Vitenti <sandro@isoftware.com.br>
 *
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
 * SECTION:ncm_fit_esmcmc_walker_aps
 * @title: NcmFitESMCMCWalkerAPS
 * @short_description: Ensemble sampler Markov Chain Monte Carlo walker - aps move.
 *
 * Implementing aps move walker for #NcmFitESMCMC (affine invariant).
 * 
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_fit_esmcmc_walker.h"
#include "math/ncm_fit_esmcmc_walker_aps.h"

#include "math/ncm_fit_esmcmc.h"
#include "math/ncm_stats_dist_nd_gauss.h"

enum
{
  PROP_0,
};

struct _NcmFitESMCMCWalkerAPSPrivate
{
  guint size;
  guint size_2;
  guint nparams;
  NcmVector *norm_box;
  gchar *desc;
  NcmStatsDistNdKDEGauss dndg0;
  NcmStatsDistNdKDEGauss dndg1;
};

G_DEFINE_TYPE_WITH_CODE (NcmFitESMCMCWalkerAPS, ncm_fit_esmcmc_walker_aps, NCM_TYPE_FIT_ESMCMC_WALKER, G_ADD_PRIVATE (NcmFitESMCMCWalkerAPS));

static void
ncm_fit_esmcmc_walker_aps_init (NcmFitESMCMCWalkerAPS *aps)
{
  NcmFitESMCMCWalkerAPSPrivate * const self = aps->priv = G_TYPE_INSTANCE_GET_PRIVATE (aps, NCM_TYPE_FIT_ESMCMC_WALKER_APS, NcmFitESMCMCWalkerAPSPrivate);

  self->size     = 0;
  self->size_2   = 0;
  self->nparams  = 0;
  self->norm_box = NULL;
  self->desc     = "APS-Move";
  self->dndg0    = NULL:
  self->dndg1    = NULL:
}

static void
_ncm_fit_esmcmc_walker_aps_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  /*NcmFitESMCMCWalkerAPS *aps = NCM_FIT_ESMCMC_WALKER_APS (object);*/
  g_return_if_fail (NCM_IS_FIT_ESMCMC_WALKER_APS (object));

  switch (prop_id)
  {
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_ncm_fit_esmcmc_walker_aps_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  /*NcmFitESMCMCWalkerAPS *aps = NCM_FIT_ESMCMC_WALKER_APS (object);*/
  g_return_if_fail (NCM_IS_FIT_ESMCMC_WALKER_APS (object));

  switch (prop_id)
  {
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_ncm_fit_esmcmc_walker_aps_dispose (GObject *object)
{
  NcmFitESMCMCWalkerAPS *aps = NCM_FIT_ESMCMC_WALKER_APS (object);
  NcmFitESMCMCWalkerAPSPrivate * const self = aps->priv;
  
  ncm_vector_clear (&self->norm_box);
  ncm_stats_dist_nd_kde_gauss_clear (self->dndg0);
  ncm_stats_dist_nd_kde_gauss_clear (self->dndg1);

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_fit_esmcmc_walker_aps_parent_class)->dispose (object);
}

static void
_ncm_fit_esmcmc_walker_aps_finalize (GObject *object)
{
  /*NcmFitESMCMCWalkerAPS *aps = NCM_FIT_ESMCMC_WALKER_APS (object);*/
  /*NcmFitESMCMCWalkerAPSPrivate * const self = aps->priv;*/
  
  /* Chain up : end */
  G_OBJECT_CLASS (ncm_fit_esmcmc_walker_aps_parent_class)->finalize (object);
}

static void _ncm_fit_esmcmc_walker_aps_set_size (NcmFitESMCMCWalker *walker, guint size);
static guint _ncm_fit_esmcmc_walker_aps_get_size (NcmFitESMCMCWalker *walker);
static void _ncm_fit_esmcmc_walker_aps_set_nparams (NcmFitESMCMCWalker *walker, guint nparams);
static guint _ncm_fit_esmcmc_walker_aps_get_nparams (NcmFitESMCMCWalker *walker);
static void _ncm_fit_esmcmc_walker_aps_setup (NcmFitESMCMCWalker *walker, GPtrArray *theta, GPtrArray *m2lnL, guint ki, guint kf, NcmRNG *rng);
static void _ncm_fit_esmcmc_walker_aps_step (NcmFitESMCMCWalker *walker, GPtrArray *theta, GPtrArray *m2lnL, NcmVector *thetastar, guint k);
static gdouble _ncm_fit_esmcmc_walker_aps_prob (NcmFitESMCMCWalker *walker, GPtrArray *theta, GPtrArray *m2lnL, NcmVector *thetastar, guint k, const gdouble m2lnL_cur, const gdouble m2lnL_star);
static gdouble _ncm_fit_esmcmc_walker_aps_prob_norm (NcmFitESMCMCWalker *walker, GPtrArray *theta, GPtrArray *m2lnL, NcmVector *thetastar, guint k);
static void _ncm_fit_esmcmc_walker_aps_clean (NcmFitESMCMCWalker *walker, guint ki, guint kf);
static const gchar *_ncm_fit_esmcmc_walker_aps_desc (NcmFitESMCMCWalker *walker);

static void
ncm_fit_esmcmc_walker_aps_class_init (NcmFitESMCMCWalkerAPSClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);
  NcmFitESMCMCWalkerClass *walker_class = NCM_FIT_ESMCMC_WALKER_CLASS (klass);

  object_class->set_property = _ncm_fit_esmcmc_walker_aps_set_property;
  object_class->get_property = _ncm_fit_esmcmc_walker_aps_get_property;
  object_class->dispose      = _ncm_fit_esmcmc_walker_aps_dispose;
  object_class->finalize     = _ncm_fit_esmcmc_walker_aps_finalize;

/*  
   g_object_class_install_property (object_class,
                                   PROP_GCONST,
                                   g_param_spec_double ("G",
                                                        NULL,
                                                        "`Gravitation' constant G",
                                                        1.0e-15, G_MAXDOUBLE, _G,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
*/  
  walker_class->set_size    = &_ncm_fit_esmcmc_walker_aps_set_size;
  walker_class->get_size    = &_ncm_fit_esmcmc_walker_aps_get_size;
  walker_class->set_nparams = &_ncm_fit_esmcmc_walker_aps_set_nparams;
  walker_class->get_nparams = &_ncm_fit_esmcmc_walker_aps_get_nparams;
  walker_class->setup       = &_ncm_fit_esmcmc_walker_aps_setup;
  walker_class->step        = &_ncm_fit_esmcmc_walker_aps_step;
  walker_class->prob        = &_ncm_fit_esmcmc_walker_aps_prob;
  walker_class->prob_norm   = &_ncm_fit_esmcmc_walker_aps_prob_norm;
  walker_class->clean       = &_ncm_fit_esmcmc_walker_aps_clean;
  walker_class->desc        = &_ncm_fit_esmcmc_walker_aps_desc;
}

static void 
_ncm_fit_esmcmc_walker_aps_set_sys (NcmFitESMCMCWalker *walker, guint size, guint nparams)
{
  NcmFitESMCMCWalkerAPS *aps = NCM_FIT_ESMCMC_WALKER_APS (walker);
  NcmFitESMCMCWalkerAPSPrivate * const self = aps->priv;
  
  g_assert_cmpuint (size, >, 0);
  g_assert_cmpuint (nparams, >, 0);
  
  if (self->size != size || self->nparams != nparams)
  {
    guint i;
    
    ncm_stats_dist_nd_kde_gauss_clear (&self->dndg0);
    ncm_stats_dist_nd_kde_gauss_clear (&self->dndg1);
    ncm_vector_clear (&self->norm_box);
    
    g_assert (size % 2 == 0);
    self->size    = size;
    self->size_2  = size / 2;
    self->nparams = nparams;
    
    self->norm_box = ncm_vector_new (self->size);
    self->dndg0    = ncm_stats_dist_nd_kde_gauss_new (self->nparams);
    self->dndg1    = ncm_stats_dist_nd_kde_gauss_new (self->nparams);
  }
}

static void 
_ncm_fit_esmcmc_walker_aps_set_size (NcmFitESMCMCWalker *walker, guint size)
{
  NcmFitESMCMCWalkerAPS *aps = NCM_FIT_ESMCMC_WALKER_APS (walker);
  NcmFitESMCMCWalkerAPSPrivate * const self = aps->priv;

  g_assert_cmpuint (size, >, 0);
  
  if (self->nparams != 0)
    _ncm_fit_esmcmc_walker_aps_set_sys (walker, size, self->nparams);
  else
    self->size = size;  
}

static guint 
_ncm_fit_esmcmc_walker_aps_get_size (NcmFitESMCMCWalker *walker)
{
  NcmFitESMCMCWalkerAPS *aps = NCM_FIT_ESMCMC_WALKER_APS (walker);
  NcmFitESMCMCWalkerAPSPrivate * const self = aps->priv;

  return self->size;
}

static void 
_ncm_fit_esmcmc_walker_aps_set_nparams (NcmFitESMCMCWalker *walker, guint nparams)
{
  NcmFitESMCMCWalkerAPS *aps = NCM_FIT_ESMCMC_WALKER_APS (walker);
  NcmFitESMCMCWalkerAPSPrivate * const self = aps->priv;

  g_assert_cmpuint (nparams, >, 0);
  
  if (self->size != 0)
    _ncm_fit_esmcmc_walker_aps_set_sys (walker, self->size, nparams);
  else
    self->nparams = nparams;
}

static guint 
_ncm_fit_esmcmc_walker_aps_get_nparams (NcmFitESMCMCWalker *walker)
{
  NcmFitESMCMCWalkerAPS *aps = NCM_FIT_ESMCMC_WALKER_APS (walker);
  NcmFitESMCMCWalkerAPSPrivate * const self = aps->priv;

  return self->nparams;
}

static void 
_ncm_fit_esmcmc_walker_aps_setup (NcmFitESMCMCWalker *walker, GPtrArray *theta, GPtrArray *m2lnL, guint ki, guint kf, NcmRNG *rng)
{
  NcmFitESMCMCWalkerAPS *aps = NCM_FIT_ESMCMC_WALKER_APS (walker);
  NcmFitESMCMCWalkerAPSPrivate * const self = aps->priv;
  /*const gdouble nf = sqrt (1.0 * self->nparams);*/
  guint k;

  /*printf ("# Setup: %u %u %u %u\n", ki, kf, self->size, self->size_2);*/
    
  for (k = ki; k < kf; k++)
  {
    guint i;
    for (i = 0; i < self->nparams; i++)
    {
      ncm_matrix_set (self->vel_norm, k, i, ncm_rng_ugaussian_gen (rng) / (_NF));
    }
  }
}

static gdouble 
_ncm_fit_esmcmc_walker_aps_theta_to_x (NcmFitESMCMCWalkerAPS *aps, guint i, const gdouble theta_i)
{
  NcmFitESMCMCWalkerAPSPrivate * const self = aps->priv;
  
  const gdouble lb      = ncm_matrix_get (self->box, i, 0);
  const gdouble ub      = ncm_matrix_get (self->box, i, 1);
  const gdouble bsize   = ub - lb;
  const gdouble tb      = theta_i - lb;
  const gdouble maxatan = - 0.5 * log (GSL_DBL_EPSILON * 0.5);

  if (tb == 0.0)             
    return -maxatan;
  else if (tb == bsize)        
    return maxatan;
  else
    return atanh (2.0 * tb / bsize - 1.0);
}

static gdouble 
_ncm_fit_esmcmc_walker_aps_x_to_theta (NcmFitESMCMCWalkerAPS *aps, guint i, const gdouble x_i)
{
  NcmFitESMCMCWalkerAPSPrivate * const self = aps->priv;
  
  const gdouble lb     = ncm_matrix_get (self->box, i, 0);
  const gdouble ub     = ncm_matrix_get (self->box, i, 1);
  const gdouble bsize  = ub - lb;
  const gdouble middle = 0.5 * (ub + lb);

  return middle + bsize * tanh (x_i) * 0.5;
}

static gdouble 
_ncm_fit_esmcmc_walker_aps_norm (NcmFitESMCMCWalkerAPS *aps, guint i, const gdouble thetastar_i, const gdouble theta_i)
{
  NcmFitESMCMCWalkerAPSPrivate * const self = aps->priv;
  const gdouble lb    = ncm_matrix_get (self->box, i, 0);
  const gdouble ub    = ncm_matrix_get (self->box, i, 1);
  const gdouble tb    = theta_i - lb;
  const gdouble bt    = ub - theta_i;
  const gdouble tsb   = thetastar_i - lb;
  const gdouble bts   = ub - thetastar_i;
  const gdouble norm  = tsb * bts / (tb * bt);
  
  return log (norm);
}

#define _POS(i) ((i) * 2) 
#define _VEL(i) ((i) * 2 + 1) 

static gint
_ncm_fit_esmcmc_walker_aps_f (realtype x, N_Vector y, N_Vector ydot, gpointer f_data)
{
  NcmFitESMCMCWalkerAPSPrivate * const self = (NcmFitESMCMCWalkerAPSPrivate *) f_data;
  guint i, j;

  for (j = 0; j < self->nparams; j++)
  {
    NV_Ith_S (ydot, _POS (j)) = NV_Ith_S (y, _VEL (j));
    NV_Ith_S (ydot, _VEL (j)) = 0.0;
  }
  
  
  for (i = 0; i < self->size_2; i++)
  {
    const gdouble mass_i = ncm_vector_get (self->mass, i);
    gdouble norm_i = 0.0;

    if (mass_i == 0.0)
      continue;
    
    for (j = 0; j < self->nparams; j++)
    {
      const gdouble psi_j       = NV_Ith_S (y, _POS (j));
      const gdouble Xij_m_psi_i = ncm_matrix_get (self->Xij, i, j) - psi_j;

      norm_i += gsl_pow_2 (Xij_m_psi_i);
    }

    /*norm_i = pow (norm_i + 1.0 * 1.0e-2, 0.5 * self->nparams);*/
    norm_i = pow (norm_i + 1.0 * _EPSILON * _EPSILON, 1.5);

    for (j = 0; j < self->nparams; j++)
    {
      const gdouble psi_j       = NV_Ith_S (y, _POS (j));
      const gdouble Xij_m_psi_i = ncm_matrix_get (self->Xij, i, j) - psi_j;

      NV_Ith_S (ydot, _VEL (j)) += mass_i * self->G * Xij_m_psi_i / (norm_i);
    }    
  }
  
  return 0;
}

static void
_ncm_fit_esmcmc_walker_aps_step (NcmFitESMCMCWalker *walker, GPtrArray *theta, GPtrArray *m2lnL, NcmVector *thetastar, guint k)
{
  NcmFitESMCMCWalkerAPS *aps = NCM_FIT_ESMCMC_WALKER_APS (walker);
  NcmFitESMCMCWalkerAPSPrivate * const self = aps->priv;
  NcmVector *theta_k      = g_ptr_array_index (theta, k);
  const gdouble m2lnL_k   = ncm_vector_get (g_ptr_array_index (m2lnL, k), 0);
  const guint subensemble = (k < self->size_2) ? self->size_2 : 0;
  gdouble XtoY_YtoX = 0.0;
  gdouble tf        = 0.0;
  gdouble m2lnL_min = 1.0e300;
  guint i, l;
  gint flag;

  ncm_vector_set (self->norm_box, k, 0.0);
    
  ncm_stats_vec_reset (self->walker_stats, TRUE);

  for (i = 0; i < self->size_2; i++)
  {
    guint j               = i + subensemble;
    NcmVector *theta_j    = g_ptr_array_index (theta, j);
    const gdouble m2lnL_j = ncm_vector_get (g_ptr_array_index (m2lnL, j), 0);

    m2lnL_min = MIN (m2lnL_min, m2lnL_j);
    
    ncm_vector_set (self->mass, i, m2lnL_j);

    /*printf ("# setting mass %3u (%3u) % 22.15g | % 22.15g (% 22.15g)\n", i, j, -0.5 * m2lnL_j, -0.5 * m2lnL_k, -0.5 * (m2lnL_j - m2lnL_k));*/
    
    for (l = 0; l < self->nparams; l++)
    {
      const gdouble theta_j_l = ncm_vector_get (theta_j, l);
      const gdouble x_l       = g_array_index (self->use_box, gboolean, l) ? _ncm_fit_esmcmc_walker_aps_theta_to_x (aps, l, theta_j_l) : theta_j_l;

      ncm_stats_vec_set (self->walker_stats, l, x_l);
    }
    ncm_stats_vec_update (self->walker_stats);
  }

  for (i = 0; i < self->size_2; i++)
  {
    const gdouble A = -0.5 * (ncm_vector_get (self->mass, i) - m2lnL_min);
    ncm_vector_set (self->mass, i, exp (A));
    
    /*printf ("# mass %4u % 22.15g % 22.15g\n", i, ncm_vector_get (self->mass, i), A);*/
  }

  for (l = 0; l < self->nparams; l++)
  {
    const gdouble sd_l = ncm_stats_vec_get_sd (self->walker_stats, l);
    
    NV_Ith_S (self->y, _POS (l)) = ncm_vector_get (theta_k, l) / sd_l;
    NV_Ith_S (self->y, _VEL (l)) = ncm_matrix_get (self->vel_norm, k, l);

    XtoY_YtoX += gsl_pow_2 (NV_Ith_S (self->y, _VEL (l)));

    /*printf ("POS %u % 22.15g : % 22.15g / % 22.15g\n", l, NV_Ith_S (self->y, _POS (l)), ncm_vector_get (theta_k, l), sd_l);*/
    
    for (i = 0; i < self->size_2; i++)
    {
      guint j                 = i + subensemble;
      NcmVector *theta_j      = g_ptr_array_index (theta, j);
      const gdouble theta_j_l = ncm_vector_get (theta_j, l);
      const gdouble x_l       = g_array_index (self->use_box, gboolean, l) ? _ncm_fit_esmcmc_walker_aps_theta_to_x (aps, l, theta_j_l) : theta_j_l;

      ncm_matrix_set (self->Xij, i, l, x_l / sd_l);
    }
  }
/*
  ncm_matrix_log_vals (self->Xij, "    Xij: ", "% 22.15g");

  printf ("INI POS:");
  for (i = 0; i < self->nparams; i++)
    printf (" % 22.15g", NV_Ith_S (self->y, _POS (i)));
  printf ("\n");

  printf ("INI VEL:");
  for (i = 0; i < self->nparams; i++)
    printf (" % 22.15g", NV_Ith_S (self->y, _VEL (i)));
  printf ("\n");
*/  
  if (!self->cvode_init)
  {
    flag = CVodeInit (self->cvode, &_ncm_fit_esmcmc_walker_aps_f, 0.0, self->y);
    NCM_CVODE_CHECK (&flag, "CVodeInit", 1, );
    
    self->cvode_init = TRUE;
  }
  else
  {
    flag = CVodeReInit (self->cvode, 0.0, self->y);
    NCM_CVODE_CHECK (&flag, "CVodeReInit", 1, );
  }

  flag = CVodeSStolerances (self->cvode, 1.0e-13, 0.0);
  NCM_CVODE_CHECK (&flag, "CVodeSStolerances", 1, );
  
  flag = CVodeSetMaxNumSteps (self->cvode, 100000000);
  NCM_CVODE_CHECK (&flag, "CVodeSetMaxNumSteps", 1, );
  
  flag = CVodeSetUserData (self->cvode, self);
  NCM_CVODE_CHECK (&flag, "CVodeSetUserData", 1, );

  flag = CVodeSetStopTime (self->cvode, MAX_TIME);
  NCM_CVODE_CHECK (&flag, "CVodeSetStopTime", 1, );

#if HAVE_SUNDIALS_MAJOR == 2
  flag = CVDense (self->cvode, NCM_HOAA_VAR_SYS_SIZE);
  NCM_CVODE_CHECK (&flag, "CVDense", 1, );

  /*flag = CVDlsSetDenseJacFn (self->cvode, J);*/
  /*NCM_CVODE_CHECK (&flag, "CVDlsSetDenseJacFn", 1, );*/
#elif HAVE_SUNDIALS_MAJOR == 3
  flag = CVDlsSetLinearSolver (self->cvode, self->LS, self->A);
  NCM_CVODE_CHECK (&flag, "CVDlsSetLinearSolver", 1, );

  /*flag = CVDlsSetJacFn (self->cvode, J);*/
  /*NCM_CVODE_CHECK (&flag, "CVDlsSetJacFn", 1, );*/
#endif
  
  flag = CVode (self->cvode, MAX_TIME, self->y, &tf, CV_NORMAL);
  NCM_CVODE_CHECK (&flag, "ncm_ode_spline_prepare[CVode]", 1, );
/*
  printf ("END POS:");
  for (i = 0; i < self->nparams; i++)
    printf (" % 22.15g", NV_Ith_S (self->y, _POS (i)));
  printf ("\n");

  printf ("END VEL:");
  for (i = 0; i < self->nparams; i++)
    printf (" % 22.15g", NV_Ith_S (self->y, _VEL (i)));
  printf ("\n");
*/  
  for (i = 0; i < self->nparams; i++)
  {
    const gdouble sd_l = ncm_stats_vec_get_sd (self->walker_stats, i);
    const gdouble x_i  = NV_Ith_S (self->y, _POS (i)) * sd_l;
    gdouble thetastar_i;

    XtoY_YtoX -= gsl_pow_2 (NV_Ith_S (self->y, _VEL (i)));
    
    if (g_array_index (self->use_box, gboolean, i))
    {
      const gdouble theta_k_i = ncm_vector_get (theta_k, i);
      thetastar_i = _ncm_fit_esmcmc_walker_aps_x_to_theta (aps, i, x_i);

      ncm_vector_addto (self->norm_box, k, 
                        _ncm_fit_esmcmc_walker_aps_norm (aps, i, thetastar_i, theta_k_i));
    }
    else
      thetastar_i = x_i;

    ncm_vector_set (thetastar, i, thetastar_i);
  }

  ncm_vector_addto (self->norm_box, k, 0.5 * XtoY_YtoX);

  printf ("XtoY_YtoX: % 22.15g\n", 0.5 * XtoY_YtoX);

  ncm_vector_log_vals (theta_k,   "    THETA: ", "% 22.15g", TRUE);
  ncm_vector_log_vals (thetastar, "THETASTAR: ", "% 22.15g", TRUE);
}

static gdouble 
_ncm_fit_esmcmc_walker_aps_prob (NcmFitESMCMCWalker *walker, GPtrArray *theta, GPtrArray *m2lnL, NcmVector *thetastar, guint k, const gdouble m2lnL_cur, const gdouble m2lnL_star)
{
  NcmFitESMCMCWalkerAPS *aps = NCM_FIT_ESMCMC_WALKER_APS (walker);
  NcmFitESMCMCWalkerAPSPrivate * const self = aps->priv;

  return exp ((m2lnL_cur - m2lnL_star) * 0.5 + ncm_vector_get (self->norm_box, k));
}

static gdouble 
_ncm_fit_esmcmc_walker_aps_prob_norm (NcmFitESMCMCWalker *walker, GPtrArray *theta, GPtrArray *m2lnL, NcmVector *thetastar, guint k)
{
  NcmFitESMCMCWalkerAPS *aps = NCM_FIT_ESMCMC_WALKER_APS (walker);
  NcmFitESMCMCWalkerAPSPrivate * const self = aps->priv;

  return ncm_vector_get (self->norm_box, k);
}

static void 
_ncm_fit_esmcmc_walker_aps_clean (NcmFitESMCMCWalker *walker, guint ki, guint kf)
{
  /* Nothing to do. */
}

const gchar *
_ncm_fit_esmcmc_walker_aps_desc (NcmFitESMCMCWalker *walker)
{
  NcmFitESMCMCWalkerAPS *aps = NCM_FIT_ESMCMC_WALKER_APS (walker);  
  NcmFitESMCMCWalkerAPSPrivate * const self = aps->priv;
  
  return self->desc;
}

/**
 * ncm_fit_esmcmc_walker_aps_new:
 * @nwalkers: number of walkers
 * @nparams: number of parameters
 * 
 * Creates a new #NcmFitESMCMCWalkerAPS to be used
 * with @nwalkers.
 * 
 * Returns: (transfer full): a new #NcmFitESMCMCWalkerAPS.
 */
NcmFitESMCMCWalkerAPS *
ncm_fit_esmcmc_walker_aps_new (guint nwalkers, guint nparams)
{
  NcmFitESMCMCWalkerAPS *aps = g_object_new (NCM_TYPE_FIT_ESMCMC_WALKER_APS,
                                                     "size", nwalkers,
                                                     "nparams", nparams,
                                                     NULL);

  return aps;
}

/**
 * ncm_fit_esmcmc_walker_aps_set_G:
 * @aps: a #NcmFitESMCMCWalkerAPS
 * @G: new `gravitation' constant $G > 10^{-15}$
 * 
 * Sets the value of the `gravitation' constant $G$.
 * 
 */
void
ncm_fit_esmcmc_walker_aps_set_G (NcmFitESMCMCWalkerAPS *aps, const gdouble G)
{
  NcmFitESMCMCWalkerAPSPrivate * const self = aps->priv;
  
  g_assert_cmpfloat (G, >, 1.0e-15);
  self->G = G;
}

/**
 * ncm_fit_esmcmc_walker_aps_get_G:
 * @aps: a #NcmFitESMCMCWalkerAPS
 * 
 * Gets the value of the `gravitation' constant $G$.
 * 
 * Returns: current value of $G$.
 */
gdouble
ncm_fit_esmcmc_walker_aps_get_G (NcmFitESMCMCWalkerAPS *aps)
{
  NcmFitESMCMCWalkerAPSPrivate * const self = aps->priv;
  
  return self->G;
}

/**
 * ncm_fit_esmcmc_walker_aps_set_box:
 * @aps: a #NcmFitESMCMCWalkerAPS
 * @n: parameter index
 * @lb: lower bound
 * @ub: upper bound
 * 
 * Sets box sampling for the @n-th parameter using @lb as lower bound
 * and @ub as upper bound.
 * 
 */
void 
ncm_fit_esmcmc_walker_aps_set_box (NcmFitESMCMCWalkerAPS *aps, guint n, const gdouble lb, const gdouble ub)
{
  NcmFitESMCMCWalkerAPSPrivate * const self = aps->priv;
  
  g_assert_cmpuint (n, <, self->nparams);
  g_assert_cmpfloat (lb, <, ub);

  g_array_index (self->use_box, gboolean, n) = TRUE;
  ncm_matrix_set (self->box, n, 0, lb);
  ncm_matrix_set (self->box, n, 1, ub);
}

/**
 * ncm_fit_esmcmc_walker_aps_set_box_mset:
 * @aps: a #NcmFitESMCMCWalkerAPS
 * @mset: a #NcmMSet
 * 
 * Sets box sampling for the parameters using bounds from
 * @mset.
 * 
 */
void 
ncm_fit_esmcmc_walker_aps_set_box_mset (NcmFitESMCMCWalkerAPS *aps, NcmMSet *mset)
{
  NcmFitESMCMCWalkerAPSPrivate * const self = aps->priv;
  const guint fparams_len = ncm_mset_fparams_len (mset);
  guint i;
  
  g_assert_cmpuint (self->nparams, ==, fparams_len);
  for (i = 0; i < fparams_len; i++)
  {
    ncm_fit_esmcmc_walker_aps_set_box (aps, i, 
                                          ncm_mset_fparam_get_lower_bound (mset, i),
                                          ncm_mset_fparam_get_upper_bound (mset, i));
  }  
}
