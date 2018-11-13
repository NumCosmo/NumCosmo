/***************************************************************************
 *            ncm_fit_esmcmc_walker_newton.c
 *
 *  Sat October 27 13:08:13 2018
 *  Copyright  2018  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * ncm_fit_esmcmc_walker_newton.c
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
 * SECTION:ncm_fit_esmcmc_walker_newton
 * @title: NcmFitESMCMCWalkerNewton
 * @short_description: Ensemble sampler Markov Chain Monte Carlo walker - newton move.
 *
 * Implementing newton move walker for #NcmFitESMCMC (affine invariant).
 * 
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_fit_esmcmc_walker.h"
#include "math/ncm_fit_esmcmc_walker_newton.h"

#include "math/ncm_fit_esmcmc.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <cvodes/cvodes.h>
#if HAVE_SUNDIALS_MAJOR == 2
#include <cvodes/cvodes_dense.h>
#include <cvodes/cvodes_band.h>
#elif HAVE_SUNDIALS_MAJOR == 3
#include <cvodes/cvodes_direct.h> 
#endif
#include <nvector/nvector_serial.h> 

#if HAVE_SUNDIALS_MAJOR == 3
#include <sundials/sundials_matrix.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <sunlinsol/sunlinsol_dense.h>
#define SUN_DENSE_ACCESS SM_ELEMENT_D
#endif 

#include <sundials/sundials_types.h> 

#include <gsl/gsl_linalg.h>
#endif /* NUMCOSMO_GIR_SCAN */

enum
{
  PROP_0,
  PROP_GCONST,
};

struct _NcmFitESMCMCWalkerNewtonPrivate
{
  guint size;
  guint size_2;
  guint nparams;
  gdouble G;
  NcmMatrix *box;
  NcmVector *norm_box;
  GArray *use_box;
  GArray *indices;
  GArray *numbers;
  gchar *desc;
  NcmVector *mass;
  NcmMatrix *Xij;
  NcmMatrix *vel_norm;
  NcmStatsVec *walker_stats;
  N_Vector y;
  gpointer cvode;
  gboolean cvode_init;
#if HAVE_SUNDIALS_MAJOR == 3
  SUNMatrix A;
  SUNLinearSolver LS;
#endif
};

G_DEFINE_TYPE_WITH_CODE (NcmFitESMCMCWalkerNewton, ncm_fit_esmcmc_walker_newton, NCM_TYPE_FIT_ESMCMC_WALKER, G_ADD_PRIVATE (NcmFitESMCMCWalkerNewton));

#define _G (1.0e2)
#define _NF (1.0)
#define _EPSILON (1.0e-2)
#define MAX_TIME (1.0)

static void
ncm_fit_esmcmc_walker_newton_init (NcmFitESMCMCWalkerNewton *newton)
{
  NcmFitESMCMCWalkerNewtonPrivate * const self = newton->priv = G_TYPE_INSTANCE_GET_PRIVATE (newton, NCM_TYPE_FIT_ESMCMC_WALKER_NEWTON, NcmFitESMCMCWalkerNewtonPrivate);

  self->size         = 0;
  self->size_2       = 0;
  self->nparams      = 0;
  self->G            = 0.0;
  self->box          = NULL;
  self->norm_box     = NULL;
  self->use_box      = g_array_new (TRUE, TRUE, sizeof (gboolean));
  self->indices      = g_array_new (TRUE, TRUE, sizeof (guint));
  self->numbers      = g_array_new (TRUE, TRUE, sizeof (guint));
  self->desc         = "Newton-Move";
  self->mass         = NULL;
  self->Xij          = NULL;
  self->vel_norm     = NULL;
  self->walker_stats = NULL;
  self->y            = NULL;
  self->cvode        = NULL;
  self->cvode_init   = FALSE;

#if HAVE_SUNDIALS_MAJOR == 3
  self->A            = NULL;
  self->LS           = NULL;
#endif  
}

static void
_ncm_fit_esmcmc_walker_newton_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmFitESMCMCWalkerNewton *newton = NCM_FIT_ESMCMC_WALKER_NEWTON (object);
  g_return_if_fail (NCM_IS_FIT_ESMCMC_WALKER_NEWTON (object));

  switch (prop_id)
  {
    case PROP_GCONST:
      ncm_fit_esmcmc_walker_newton_set_G (newton, g_value_get_double (value));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_ncm_fit_esmcmc_walker_newton_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmFitESMCMCWalkerNewton *newton = NCM_FIT_ESMCMC_WALKER_NEWTON (object);
  g_return_if_fail (NCM_IS_FIT_ESMCMC_WALKER_NEWTON (object));

  switch (prop_id)
  {
    case PROP_GCONST:
      g_value_set_double (value, ncm_fit_esmcmc_walker_newton_get_G (newton));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_ncm_fit_esmcmc_walker_newton_dispose (GObject *object)
{
  NcmFitESMCMCWalkerNewton *newton = NCM_FIT_ESMCMC_WALKER_NEWTON (object);
  NcmFitESMCMCWalkerNewtonPrivate * const self = newton->priv;
  
  ncm_stats_vec_clear (&self->walker_stats);

  ncm_vector_clear (&self->norm_box);
  ncm_vector_clear (&self->mass);
  ncm_matrix_clear (&self->vel_norm);
  ncm_matrix_clear (&self->Xij);
  ncm_matrix_clear (&self->box);

  g_clear_pointer (&self->use_box, g_array_unref);
  g_clear_pointer (&self->indices, g_array_unref);
  g_clear_pointer (&self->numbers, g_array_unref);

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_fit_esmcmc_walker_newton_parent_class)->dispose (object);
}

static void
_ncm_fit_esmcmc_walker_newton_finalize (GObject *object)
{
  NcmFitESMCMCWalkerNewton *newton = NCM_FIT_ESMCMC_WALKER_NEWTON (object);
  NcmFitESMCMCWalkerNewtonPrivate * const self = newton->priv;

  if (self->cvode != NULL)
  {
    CVodeFree (&self->cvode);
    self->cvode = NULL;
  }
  if (self->y != NULL)
  {
    N_VDestroy (self->y);
    self->y = NULL;
  }

#if HAVE_SUNDIALS_MAJOR == 3
  if (self->A != NULL)
  {
    SUNMatDestroy (self->A);
    self->A = NULL;
  }

  if (self->LS != NULL)
  {
    SUNLinSolFree (self->LS);
    self->LS = NULL;
  }
#endif
  
  /* Chain up : end */
  G_OBJECT_CLASS (ncm_fit_esmcmc_walker_newton_parent_class)->finalize (object);
}

static void _ncm_fit_esmcmc_walker_newton_set_size (NcmFitESMCMCWalker *walker, guint size);
static guint _ncm_fit_esmcmc_walker_newton_get_size (NcmFitESMCMCWalker *walker);
static void _ncm_fit_esmcmc_walker_newton_set_nparams (NcmFitESMCMCWalker *walker, guint nparams);
static guint _ncm_fit_esmcmc_walker_newton_get_nparams (NcmFitESMCMCWalker *walker);
static void _ncm_fit_esmcmc_walker_newton_setup (NcmFitESMCMCWalker *walker, GPtrArray *theta, GPtrArray *m2lnL, guint ki, guint kf, NcmRNG *rng);
static void _ncm_fit_esmcmc_walker_newton_step (NcmFitESMCMCWalker *walker, GPtrArray *theta, GPtrArray *m2lnL, NcmVector *thetastar, guint k);
static gdouble _ncm_fit_esmcmc_walker_newton_prob (NcmFitESMCMCWalker *walker, GPtrArray *theta, GPtrArray *m2lnL, NcmVector *thetastar, guint k, const gdouble m2lnL_cur, const gdouble m2lnL_star);
static gdouble _ncm_fit_esmcmc_walker_newton_prob_norm (NcmFitESMCMCWalker *walker, GPtrArray *theta, GPtrArray *m2lnL, NcmVector *thetastar, guint k);
static void _ncm_fit_esmcmc_walker_newton_clean (NcmFitESMCMCWalker *walker, guint ki, guint kf);
static const gchar *_ncm_fit_esmcmc_walker_newton_desc (NcmFitESMCMCWalker *walker);

static void
ncm_fit_esmcmc_walker_newton_class_init (NcmFitESMCMCWalkerNewtonClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);
  NcmFitESMCMCWalkerClass *walker_class = NCM_FIT_ESMCMC_WALKER_CLASS (klass);

  object_class->set_property = _ncm_fit_esmcmc_walker_newton_set_property;
  object_class->get_property = _ncm_fit_esmcmc_walker_newton_get_property;
  object_class->dispose      = _ncm_fit_esmcmc_walker_newton_dispose;
  object_class->finalize     = _ncm_fit_esmcmc_walker_newton_finalize;

  g_object_class_install_property (object_class,
                                   PROP_GCONST,
                                   g_param_spec_double ("G",
                                                        NULL,
                                                        "`Gravitation' constant G",
                                                        1.0e-15, G_MAXDOUBLE, _G,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  
  walker_class->set_size    = &_ncm_fit_esmcmc_walker_newton_set_size;
  walker_class->get_size    = &_ncm_fit_esmcmc_walker_newton_get_size;
  walker_class->set_nparams = &_ncm_fit_esmcmc_walker_newton_set_nparams;
  walker_class->get_nparams = &_ncm_fit_esmcmc_walker_newton_get_nparams;
  walker_class->setup       = &_ncm_fit_esmcmc_walker_newton_setup;
  walker_class->step        = &_ncm_fit_esmcmc_walker_newton_step;
  walker_class->prob        = &_ncm_fit_esmcmc_walker_newton_prob;
  walker_class->prob_norm   = &_ncm_fit_esmcmc_walker_newton_prob_norm;
  walker_class->clean       = &_ncm_fit_esmcmc_walker_newton_clean;
  walker_class->desc        = &_ncm_fit_esmcmc_walker_newton_desc;
}

static void 
_ncm_fit_esmcmc_walker_newton_set_sys (NcmFitESMCMCWalker *walker, guint size, guint nparams)
{
  NcmFitESMCMCWalkerNewton *newton = NCM_FIT_ESMCMC_WALKER_NEWTON (walker);
  NcmFitESMCMCWalkerNewtonPrivate * const self = newton->priv;
  
  g_assert_cmpuint (size, >, 0);
  g_assert_cmpuint (nparams, >, 0);
  
  if (self->size != size || self->nparams != nparams)
  {
    guint i;
    
    ncm_stats_vec_clear (&self->walker_stats);

    ncm_vector_clear (&self->norm_box);
    ncm_vector_clear (&self->mass);
    
    ncm_matrix_clear (&self->box);
    ncm_matrix_clear (&self->Xij);
    ncm_matrix_clear (&self->vel_norm);

    if (self->cvode != NULL)
    {
      CVodeFree (&self->cvode);
      self->cvode = NULL;
    }

    if (self->y != NULL)
    {
      N_VDestroy (self->y);
      self->y = NULL;
    }
    
    g_assert (size % 2 == 0);
    self->size    = size;
    self->size_2  = size / 2;
    self->nparams = nparams;
    
    self->walker_stats = ncm_stats_vec_new (nparams, NCM_STATS_VEC_VAR, FALSE);
    self->norm_box     = ncm_vector_new (self->size);
    self->mass         = ncm_vector_new (self->size_2);
    self->box          = ncm_matrix_new (self->nparams, 2);
    self->vel_norm     = ncm_matrix_new (self->size,   self->nparams);
    self->Xij          = ncm_matrix_new (self->size_2, self->nparams);

    g_array_set_size (self->use_box, self->nparams);
    g_array_set_size (self->indices, self->size * self->nparams);
    g_array_set_size (self->numbers, self->size);

    for (i = 0; i < size; i++)
      g_array_index (self->numbers, guint, i) = i;

    self->y          = N_VNew_Serial (self->nparams * 2);
    self->cvode      = CVodeCreate (CV_ADAMS, CV_FUNCTIONAL);//CVodeCreate (CV_BDF, CV_NEWTON);
    self->cvode_init = FALSE;

#if HAVE_SUNDIALS_MAJOR == 3
    if (self->A != NULL)
    {
      SUNMatDestroy (self->A);
      self->A = NULL;
    }

    if (self->LS != NULL)
    {
      SUNLinSolFree (self->LS);
      self->LS = NULL;
    }
        
    self->A  = SUNDenseMatrix (self->nparams * 2, self->nparams * 2);
    self->LS = SUNDenseLinearSolver (self->y, self->A);

    NCM_CVODE_CHECK ((gpointer)self->A, "SUNDenseMatrix", 0, );
    NCM_CVODE_CHECK ((gpointer)self->LS, "SUNDenseLinearSolver", 0, );
#endif
  }
}

static void 
_ncm_fit_esmcmc_walker_newton_set_size (NcmFitESMCMCWalker *walker, guint size)
{
  NcmFitESMCMCWalkerNewton *newton = NCM_FIT_ESMCMC_WALKER_NEWTON (walker);
  NcmFitESMCMCWalkerNewtonPrivate * const self = newton->priv;

  g_assert_cmpuint (size, >, 0);
  
  if (self->nparams != 0)
    _ncm_fit_esmcmc_walker_newton_set_sys (walker, size, self->nparams);
  else
    self->size = size;  
}

static guint 
_ncm_fit_esmcmc_walker_newton_get_size (NcmFitESMCMCWalker *walker)
{
  NcmFitESMCMCWalkerNewton *newton = NCM_FIT_ESMCMC_WALKER_NEWTON (walker);
  NcmFitESMCMCWalkerNewtonPrivate * const self = newton->priv;

  return self->size;
}

static void 
_ncm_fit_esmcmc_walker_newton_set_nparams (NcmFitESMCMCWalker *walker, guint nparams)
{
  NcmFitESMCMCWalkerNewton *newton = NCM_FIT_ESMCMC_WALKER_NEWTON (walker);
  NcmFitESMCMCWalkerNewtonPrivate * const self = newton->priv;

  g_assert_cmpuint (nparams, >, 0);
  
  if (self->size != 0)
    _ncm_fit_esmcmc_walker_newton_set_sys (walker, self->size, nparams);
  else
    self->nparams = nparams;
}

static guint 
_ncm_fit_esmcmc_walker_newton_get_nparams (NcmFitESMCMCWalker *walker)
{
  NcmFitESMCMCWalkerNewton *newton = NCM_FIT_ESMCMC_WALKER_NEWTON (walker);
  NcmFitESMCMCWalkerNewtonPrivate * const self = newton->priv;

  return self->nparams;
}

static void 
_ncm_fit_esmcmc_walker_newton_setup (NcmFitESMCMCWalker *walker, GPtrArray *theta, GPtrArray *m2lnL, guint ki, guint kf, NcmRNG *rng)
{
  NcmFitESMCMCWalkerNewton *newton = NCM_FIT_ESMCMC_WALKER_NEWTON (walker);
  NcmFitESMCMCWalkerNewtonPrivate * const self = newton->priv;
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
_ncm_fit_esmcmc_walker_newton_theta_to_x (NcmFitESMCMCWalkerNewton *newton, guint i, const gdouble theta_i)
{
  NcmFitESMCMCWalkerNewtonPrivate * const self = newton->priv;
  
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
_ncm_fit_esmcmc_walker_newton_x_to_theta (NcmFitESMCMCWalkerNewton *newton, guint i, const gdouble x_i)
{
  NcmFitESMCMCWalkerNewtonPrivate * const self = newton->priv;
  
  const gdouble lb     = ncm_matrix_get (self->box, i, 0);
  const gdouble ub     = ncm_matrix_get (self->box, i, 1);
  const gdouble bsize  = ub - lb;
  const gdouble middle = 0.5 * (ub + lb);

  return middle + bsize * tanh (x_i) * 0.5;
}

static gdouble 
_ncm_fit_esmcmc_walker_newton_norm (NcmFitESMCMCWalkerNewton *newton, guint i, const gdouble thetastar_i, const gdouble theta_i)
{
  NcmFitESMCMCWalkerNewtonPrivate * const self = newton->priv;
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
_ncm_fit_esmcmc_walker_newton_f (realtype x, N_Vector y, N_Vector ydot, gpointer f_data)
{
  NcmFitESMCMCWalkerNewtonPrivate * const self = (NcmFitESMCMCWalkerNewtonPrivate *) f_data;
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
_ncm_fit_esmcmc_walker_newton_step (NcmFitESMCMCWalker *walker, GPtrArray *theta, GPtrArray *m2lnL, NcmVector *thetastar, guint k)
{
  NcmFitESMCMCWalkerNewton *newton = NCM_FIT_ESMCMC_WALKER_NEWTON (walker);
  NcmFitESMCMCWalkerNewtonPrivate * const self = newton->priv;
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
      const gdouble x_l       = g_array_index (self->use_box, gboolean, l) ? _ncm_fit_esmcmc_walker_newton_theta_to_x (newton, l, theta_j_l) : theta_j_l;

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
      const gdouble x_l       = g_array_index (self->use_box, gboolean, l) ? _ncm_fit_esmcmc_walker_newton_theta_to_x (newton, l, theta_j_l) : theta_j_l;

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
    flag = CVodeInit (self->cvode, &_ncm_fit_esmcmc_walker_newton_f, 0.0, self->y);
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
      thetastar_i = _ncm_fit_esmcmc_walker_newton_x_to_theta (newton, i, x_i);

      ncm_vector_addto (self->norm_box, k, 
                        _ncm_fit_esmcmc_walker_newton_norm (newton, i, thetastar_i, theta_k_i));
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
_ncm_fit_esmcmc_walker_newton_prob (NcmFitESMCMCWalker *walker, GPtrArray *theta, GPtrArray *m2lnL, NcmVector *thetastar, guint k, const gdouble m2lnL_cur, const gdouble m2lnL_star)
{
  NcmFitESMCMCWalkerNewton *newton = NCM_FIT_ESMCMC_WALKER_NEWTON (walker);
  NcmFitESMCMCWalkerNewtonPrivate * const self = newton->priv;

  return exp ((m2lnL_cur - m2lnL_star) * 0.5 + ncm_vector_get (self->norm_box, k));
}

static gdouble 
_ncm_fit_esmcmc_walker_newton_prob_norm (NcmFitESMCMCWalker *walker, GPtrArray *theta, GPtrArray *m2lnL, NcmVector *thetastar, guint k)
{
  NcmFitESMCMCWalkerNewton *newton = NCM_FIT_ESMCMC_WALKER_NEWTON (walker);
  NcmFitESMCMCWalkerNewtonPrivate * const self = newton->priv;

  return ncm_vector_get (self->norm_box, k);
}

static void 
_ncm_fit_esmcmc_walker_newton_clean (NcmFitESMCMCWalker *walker, guint ki, guint kf)
{
  /* Nothing to do. */
}

const gchar *
_ncm_fit_esmcmc_walker_newton_desc (NcmFitESMCMCWalker *walker)
{
  NcmFitESMCMCWalkerNewton *newton = NCM_FIT_ESMCMC_WALKER_NEWTON (walker);  
  NcmFitESMCMCWalkerNewtonPrivate * const self = newton->priv;
  
  return self->desc;
}

/**
 * ncm_fit_esmcmc_walker_newton_new:
 * @nwalkers: number of walkers
 * @nparams: number of parameters
 * 
 * Creates a new #NcmFitESMCMCWalkerNewton to be used
 * with @nwalkers.
 * 
 * Returns: (transfer full): a new #NcmFitESMCMCWalkerNewton.
 */
NcmFitESMCMCWalkerNewton *
ncm_fit_esmcmc_walker_newton_new (guint nwalkers, guint nparams)
{
  NcmFitESMCMCWalkerNewton *newton = g_object_new (NCM_TYPE_FIT_ESMCMC_WALKER_NEWTON,
                                                     "size", nwalkers,
                                                     "nparams", nparams,
                                                     NULL);

  return newton;
}

/**
 * ncm_fit_esmcmc_walker_newton_set_G:
 * @newton: a #NcmFitESMCMCWalkerNewton
 * @G: new `gravitation' constant $G > 10^{-15}$
 * 
 * Sets the value of the `gravitation' constant $G$.
 * 
 */
void
ncm_fit_esmcmc_walker_newton_set_G (NcmFitESMCMCWalkerNewton *newton, const gdouble G)
{
  NcmFitESMCMCWalkerNewtonPrivate * const self = newton->priv;
  
  g_assert_cmpfloat (G, >, 1.0e-15);
  self->G = G;
}

/**
 * ncm_fit_esmcmc_walker_newton_get_G:
 * @newton: a #NcmFitESMCMCWalkerNewton
 * 
 * Gets the value of the `gravitation' constant $G$.
 * 
 * Returns: current value of $G$.
 */
gdouble
ncm_fit_esmcmc_walker_newton_get_G (NcmFitESMCMCWalkerNewton *newton)
{
  NcmFitESMCMCWalkerNewtonPrivate * const self = newton->priv;
  
  return self->G;
}

/**
 * ncm_fit_esmcmc_walker_newton_set_box:
 * @newton: a #NcmFitESMCMCWalkerNewton
 * @n: parameter index
 * @lb: lower bound
 * @ub: upper bound
 * 
 * Sets box sampling for the @n-th parameter using @lb as lower bound
 * and @ub as upper bound.
 * 
 */
void 
ncm_fit_esmcmc_walker_newton_set_box (NcmFitESMCMCWalkerNewton *newton, guint n, const gdouble lb, const gdouble ub)
{
  NcmFitESMCMCWalkerNewtonPrivate * const self = newton->priv;
  
  g_assert_cmpuint (n, <, self->nparams);
  g_assert_cmpfloat (lb, <, ub);

  g_array_index (self->use_box, gboolean, n) = TRUE;
  ncm_matrix_set (self->box, n, 0, lb);
  ncm_matrix_set (self->box, n, 1, ub);
}

/**
 * ncm_fit_esmcmc_walker_newton_set_box_mset:
 * @newton: a #NcmFitESMCMCWalkerNewton
 * @mset: a #NcmMSet
 * 
 * Sets box sampling for the parameters using bounds from
 * @mset.
 * 
 */
void 
ncm_fit_esmcmc_walker_newton_set_box_mset (NcmFitESMCMCWalkerNewton *newton, NcmMSet *mset)
{
  NcmFitESMCMCWalkerNewtonPrivate * const self = newton->priv;
  const guint fparams_len = ncm_mset_fparams_len (mset);
  guint i;
  
  g_assert_cmpuint (self->nparams, ==, fparams_len);
  for (i = 0; i < fparams_len; i++)
  {
    ncm_fit_esmcmc_walker_newton_set_box (newton, i, 
                                          ncm_mset_fparam_get_lower_bound (mset, i),
                                          ncm_mset_fparam_get_upper_bound (mset, i));
  }  
}
