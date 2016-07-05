/***************************************************************************
 *            nc_powspec_mnl_halofit.c
 *
 *  Thu March 17 14:57:40 2016
 *  Copyright  2016  Cyrille Doux
 *  <cdoux@apc.in2p3.fr>
 ****************************************************************************/
/*
 * nc_powspec_mnl_halofit.c
 * Copyright (C) 2016 Cyrille Doux <cdoux@apc.in2p3.fr>
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
 * SECTION:nc_powspec_mnl_halofit
 * @title: NcPowspecMNLHaloFit
 * @short_description: Class for linear matter power spectrum from a transfer function.
 *
 * Provides a linear matter power spectrum using a transfer function #NcTransferFunc.
 *
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "nc_hicosmo.h"
#include "nc_powspec_mnl_halofit.h"

#include "model/nc_hicosmo_de.h"

#include "math/ncm_spline_cubic_notaknot.h"
#include "math/integral.h"
#include "math/memory_pool.h"

#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_sf_exp.h>

struct _NcPowspecMNLHaloFitPrivate
{
  gdouble z;
  gdouble ksigma;
  gdouble an;
  gdouble bn; 
  gdouble cn;
  gdouble gamman;
  gdouble alphan;
  gdouble betan; 
  gdouble nun;
  gdouble f1; 
  gdouble f2; 
  gdouble f3;
};

enum
{
   PROP_0,
   PROP_PSML,
   PROP_ZMAXNL,
   PROP_PREC,
 };

G_DEFINE_TYPE (NcPowspecMNLHaloFit, nc_powspec_mnl_halofit, NC_TYPE_POWSPEC_MNL);

static void
nc_powspec_mnl_halofit_init (NcPowspecMNLHaloFit *pshf)
{
  pshf->psml   = NULL;
  pshf->zmaxnl = 0.0;
  pshf->prec   = 0.0;
  pshf->Rsigma = NULL;
  pshf->neff   = NULL;
  pshf->Cur    = NULL;
  pshf->priv   = G_TYPE_INSTANCE_GET_PRIVATE (pshf, NC_TYPE_POWSPEC_MNL_HALOFIT, NcPowspecMNLHaloFitPrivate);

  pshf->priv->z = HUGE_VAL;
}

static void
_nc_powspec_mnl_halofit_set_property (GObject* object, guint prop_id, const GValue* value, GParamSpec* pspec)
{
  NcPowspecMNLHaloFit* pshf = NC_POWSPEC_MNL_HALOFIT (object);
  g_return_if_fail (NC_IS_POWSPEC_MNL_HALOFIT (object));

  switch (prop_id)
  {
    case PROP_PSML:
      pshf->psml = g_value_dup_object (value);
      break;
    case PROP_ZMAXNL:
      pshf->zmaxnl = g_value_get_double (value);
      break;
    case PROP_PREC:
      pshf->prec = g_value_get_double (value);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_powspec_mnl_halofit_get_property (GObject* object, guint prop_id, GValue* value, GParamSpec* pspec)
{
  NcPowspecMNLHaloFit* pshf = NC_POWSPEC_MNL_HALOFIT (object);
  g_return_if_fail (NC_IS_POWSPEC_MNL_HALOFIT (object));

  switch (prop_id)
  {
    case PROP_PSML:
      g_value_set_object (value, pshf->psml);
      break;
    case PROP_ZMAXNL:
      g_value_set_double (value, pshf->zmaxnl);
      break;
    case PROP_PREC:
      g_value_set_double (value, pshf->prec);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_powspec_mnl_halofit_constructed (GObject* object)
{
  /* Chain up : start */
  G_OBJECT_CLASS (nc_powspec_mnl_halofit_parent_class)->constructed (object);
  {
    NcPowspecMNLHaloFit* pshf = NC_POWSPEC_MNL_HALOFIT (object);
    NcmVector *zv             = ncm_vector_new (NC_POWSPEC_MNL_HALOFIT_NKNOTS);
    guint i;

    if (NCM_POWSPEC (pshf)->zf < pshf->zmaxnl)
      ncm_powspec_require_zf (NCM_POWSPEC (pshf), pshf->zmaxnl);

    if (NCM_POWSPEC (pshf->psml)->zf < pshf->zmaxnl)
      ncm_powspec_require_zf (NCM_POWSPEC (pshf->psml), pshf->zmaxnl);

    pshf->Rsigma = ncm_spline_cubic_notaknot_new ();
    ncm_spline_set_len (pshf->Rsigma, NC_POWSPEC_MNL_HALOFIT_NKNOTS);

    pshf->neff = ncm_spline_cubic_notaknot_new ();
    ncm_spline_set_len (pshf->neff, NC_POWSPEC_MNL_HALOFIT_NKNOTS);

    pshf->Cur = ncm_spline_cubic_notaknot_new ();
    ncm_spline_set_len (pshf->Cur, NC_POWSPEC_MNL_HALOFIT_NKNOTS);

    for (i = 0; i < NC_POWSPEC_MNL_HALOFIT_NKNOTS; i++)
    {
      ncm_vector_set (zv, i, pshf->zmaxnl / (NC_POWSPEC_MNL_HALOFIT_NKNOTS - 1.0) * i);
    }

    ncm_spline_set_xv (pshf->Rsigma, zv, FALSE);
    ncm_spline_set_xv (pshf->neff, zv, FALSE);
    ncm_spline_set_xv (pshf->Cur, zv, FALSE);

    ncm_vector_free (zv);

    ncm_spline_acc (pshf->Rsigma, TRUE);
    ncm_spline_acc (pshf->neff, TRUE);
    ncm_spline_acc (pshf->Cur, TRUE);
  }
}

static void
_nc_powspec_mnl_halofit_dispose (GObject* object)
{
  NcPowspecMNLHaloFit* pshf = NC_POWSPEC_MNL_HALOFIT (object);

  nc_powspec_ml_clear (&pshf->psml);
  ncm_spline_clear (&pshf->Rsigma);
  ncm_spline_clear (&pshf->neff);
  ncm_spline_clear (&pshf->Cur);

  /* Chain up : end */
  G_OBJECT_CLASS (nc_powspec_mnl_halofit_parent_class)->dispose (object);
}

static void
_nc_powspec_mnl_halofit_finalize (GObject* object)
{

  /* Chain up : end */
  G_OBJECT_CLASS (nc_powspec_mnl_halofit_parent_class)->finalize (object);
}

static void _nc_powspec_mnl_halofit_prepare (NcmPowspec* powspec, NcmModel *model);
static gdouble _nc_powspec_mnl_halofit_eval (NcmPowspec* powspec, NcmModel *model, const gdouble z, const gdouble k);
static void _nc_powspec_mnl_halofit_eval_vec (NcmPowspec* powspec, NcmModel *model, const gdouble z, NcmVector* k, NcmVector* Pk);

static void
nc_powspec_mnl_halofit_class_init (NcPowspecMNLHaloFitClass* klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);
  NcmPowspecClass* powspec_class = NCM_POWSPEC_CLASS (klass);

  object_class->set_property = &_nc_powspec_mnl_halofit_set_property;
  object_class->get_property = &_nc_powspec_mnl_halofit_get_property;

  object_class->constructed  = &_nc_powspec_mnl_halofit_constructed;
  object_class->dispose      = &_nc_powspec_mnl_halofit_dispose;
  object_class->finalize     = &_nc_powspec_mnl_halofit_finalize;

  g_type_class_add_private (klass, sizeof (NcPowspecMNLHaloFitPrivate));

  g_object_class_install_property (object_class,
                                   PROP_PSML,
                                   g_param_spec_object ("power-spec",
                                                        NULL,
                                                        "Linear power spectrum.",
                                                        NC_TYPE_POWSPEC_ML,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_ZMAXNL,
                                   g_param_spec_double ("zmaxnl",
                                                        NULL,
                                                        "Max redshift for halofit correction",
                                                        0.0, 10000.0, 10.0,
                                                        G_PARAM_READWRITE |G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_PREC,
                                   g_param_spec_double ("prec",
                                                        NULL,
                                                        "Precision for halofit computations",
                                                        1e-15, 1.0, 1e-3,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));


  powspec_class->prepare  = &_nc_powspec_mnl_halofit_prepare;
  powspec_class->eval     = &_nc_powspec_mnl_halofit_eval;
  powspec_class->eval_vec = &_nc_powspec_mnl_halofit_eval_vec;
}

typedef struct _int_var_moment_params
{
  const gint n;
  const gdouble R;
  const gdouble z;
  NcPowspecML* ps;
  NcHICosmo* cosmo;

} int_var_moment_params;

///////////////////// INTEGRATION OVER K*R ////////////////////////////
static gdouble
_nc_powspec_mnl_halofit_var_moment_integrand (gdouble kR, gpointer params)
{
  int_var_moment_params* ts = (int_var_moment_params*)params;
  const gdouble k = kR / ts->R;
  const gdouble matter_P = ncm_powspec_eval (NCM_POWSPEC (ts->ps), NCM_MODEL (ts->cosmo), ts->z, k);
  const gdouble kR2 = kR * kR;
  const gdouble W2 = exp (-kR2);
  // printf ("%g %g %g %g %g %i %g\n", k, ts->R, ts->z, matter_P, kR2, ts->n, W2);
  return matter_P * gsl_pow_int (kR2, ts->n + 1) * W2;
}

static gdouble
_nc_powspec_mnl_halofit_var_moment (NcPowspecML* ps, NcHICosmo* cosmo, const gdouble R, const gdouble z, const gint n)
{
  gdouble result, error;
  gsl_function F;
  int_var_moment_params ivmps = { n, R, z, ps, cosmo };

  gsl_integration_workspace *w = gsl_integration_workspace_alloc (NCM_INTEGRAL_PARTITION);

  F.function = &_nc_powspec_mnl_halofit_var_moment_integrand;
  F.params = &ivmps;

  gsl_integration_qagiu (&F, NCM_DEFAULT_PRECISION, 0.0, NCM_DEFAULT_PRECISION, NCM_INTEGRAL_PARTITION, w, &result, &error);

  gsl_integration_workspace_free (w);

  return result / (gsl_pow_3 (R) * ncm_c_2_pi_2 ());
}

typedef struct _var_params
{
  NcPowspecML* ps;
  NcHICosmo* cosmo;
  const gdouble z;
} var_params;

static gdouble
_nc_powspec_mnl_halofit_varm1 (gdouble lnR, gpointer params)
{
  var_params* vps = (var_params*)params;
  const gdouble R = gsl_sf_exp (lnR);
  const gdouble sigma2_0 = _nc_powspec_mnl_halofit_var_moment (vps->ps, vps->cosmo, R, vps->z, 0);
  return sigma2_0 - 1.0;
}

static gdouble
_nc_powspec_mnl_halofit_varm1_deriv (gdouble lnR, gpointer params)
{
  var_params* vps = (var_params*)params;
  const gdouble R = gsl_sf_exp (lnR);
  const gdouble sigma2_1 = _nc_powspec_mnl_halofit_var_moment (vps->ps, vps->cosmo, R, vps->z, 1);

  return -2.0 * sigma2_1;
}

static void
_nc_powspec_mnl_halofit_varm1_fdf (gdouble lnR, gpointer params, gdouble* varm1, gdouble* dvarm1)
{
  var_params* vps = (var_params*)params;
  const gdouble R = gsl_sf_exp (lnR);
  const gdouble sigma2_0 = _nc_powspec_mnl_halofit_var_moment (vps->ps, vps->cosmo, R, vps->z, 0);
  const gdouble sigma2_1 = _nc_powspec_mnl_halofit_var_moment (vps->ps, vps->cosmo, R, vps->z, 1);

  *varm1 = sigma2_0 - 1.0;
  *dvarm1 = -2.0 * sigma2_1;
}

typedef struct _root_params
{
  NcPowspecML* ps;
  NcHICosmo* cosmo;
} root_params;

static gdouble
_nc_powspec_mnl_halofit_linear_scale (NcPowspecMNLHaloFit* pshf, NcHICosmo* cosmo, const gdouble z)
{
  /* Using gsl_root_fsolver_newton (with derivatives) */

  gint status;
  gint iter = 0, max_iter = 20;
  const gsl_root_fdfsolver_type* T;
  gsl_root_fdfsolver* s;

  gdouble lnR0 = 0.0;
  gdouble lnR  = (-z / 2.0 < NC_POWSPEC_MNL_HALOFIT_LOGRMIN) ? NC_POWSPEC_MNL_HALOFIT_LOGRMIN : -z / 2.0 + NCM_DEFAULT_PRECISION;
  gsl_function_fdf FDF;

  var_params vps = { pshf->psml, cosmo, z };

  FDF.f = &_nc_powspec_mnl_halofit_varm1;
  FDF.df = &_nc_powspec_mnl_halofit_varm1_deriv;
  FDF.fdf = &_nc_powspec_mnl_halofit_varm1_fdf;
  FDF.params = &vps;

  T = gsl_root_fdfsolver_newton;
  s = gsl_root_fdfsolver_alloc (T);
  gsl_root_fdfsolver_set (s, &FDF, lnR);

  do
  {
    iter++;
    status = gsl_root_fdfsolver_iterate (s);

    lnR0 = lnR;
    lnR = gsl_root_fdfsolver_root (s);
    status = gsl_root_test_delta (exp (lnR), exp (lnR0), 0.0, pshf->prec); // Compare R vs R0 !

  } while (status == GSL_CONTINUE && iter < max_iter);

  gsl_root_fdfsolver_free (s);

  return exp (lnR);
}

static void 
_nc_powspec_mnl_halofit_prepare_nl (NcPowspecMNLHaloFit* pshf, NcmModel* model)
{
  NcHICosmo* cosmo = NC_HICOSMO (model);
  guint i;
  gdouble ztemp, Rtemp, sigma2_0, sigma2_1, sigma2_2, d1, d2;

  NcmVector *zv      = ncm_spline_get_xv (pshf->Rsigma);
  NcmVector *Rsigmav = ncm_spline_get_yv (pshf->Rsigma);
  NcmVector *neffv   = ncm_spline_get_yv (pshf->neff);
  NcmVector *Curv    = ncm_spline_get_yv (pshf->Cur);

  for (i = 0; i < NC_POWSPEC_MNL_HALOFIT_NKNOTS; i++)
  {
    ztemp = ncm_vector_get (zv, i);
    Rtemp = _nc_powspec_mnl_halofit_linear_scale (pshf, cosmo, ztemp);

printf ("# z = % 20.15g, R = % 20.15g | % 20.15g % 20.15g % 20.15g\n", ztemp, Rtemp, 
        _nc_powspec_mnl_halofit_var_moment (pshf->psml, cosmo, Rtemp, ztemp, 0),
        _nc_powspec_mnl_halofit_var_moment (pshf->psml, cosmo, Rtemp, ztemp, 1),
        _nc_powspec_mnl_halofit_var_moment (pshf->psml, cosmo, Rtemp, ztemp, 2));
    
    ncm_vector_set (Rsigmav, i, Rtemp);

    sigma2_0 = _nc_powspec_mnl_halofit_var_moment (pshf->psml, cosmo, Rtemp, ztemp, 0); // = 1.0 by construction
    sigma2_1 = _nc_powspec_mnl_halofit_var_moment (pshf->psml, cosmo, Rtemp, ztemp, 1);
    sigma2_2 = _nc_powspec_mnl_halofit_var_moment (pshf->psml, cosmo, Rtemp, ztemp, 2);

    d1 = -2.0 * sigma2_1 / sigma2_0;
    d2 = -d1 * d1 + 4.0 * (-sigma2_1 + sigma2_2) / sigma2_0;

    ncm_vector_set (neffv, i, -3.0 - d1);
    ncm_vector_set (Curv, i, -d2);
  }

  ncm_spline_prepare (pshf->Rsigma);
  ncm_spline_prepare (pshf->neff);
  ncm_spline_prepare (pshf->Cur);

  ncm_vector_free (zv);
  ncm_vector_free (Rsigmav);
  ncm_vector_free (neffv);
  ncm_vector_free (Curv);
}

static void 
_nc_powspec_mnl_halofit_prepare (NcmPowspec* powspec, NcmModel* model)
{
  NcPowspecMNLHaloFit* pshf = NC_POWSPEC_MNL_HALOFIT (powspec);

  g_assert (NC_IS_HICOSMO (model));

  ncm_powspec_prepare (NCM_POWSPEC (pshf->psml), model);

  _nc_powspec_mnl_halofit_prepare_nl (pshf, model);
}

static void
_nc_powspec_mnl_halofit_preeval (NcPowspecMNLHaloFit *pshf, NcHICosmo *cosmo, const gdouble z)
{
  const gdouble E2     = nc_hicosmo_E2 (cosmo, z);
  const gdouble Rsigma = ncm_spline_eval (pshf->Rsigma, z);

  const gdouble neff = ncm_spline_eval (pshf->neff, z);
  const gdouble Cur  = ncm_spline_eval (pshf->Cur, z);

  const gdouble neff2 = neff * neff;
  const gdouble neff3 = neff2 * neff;
  const gdouble neff4 = neff2 * neff2;

  const gdouble Omega_de_onepp = NC_IS_HICOSMO_DE (cosmo) ? nc_hicosmo_de_E2Omega_de_onepw (NC_HICOSMO_DE (cosmo), z) / E2 : 0.0;
  const gdouble Omega_m        = nc_hicosmo_Omega_m0 (cosmo) * gsl_pow_3 (1.0 + z) / E2;

  pshf->priv->z      = z;
  
  pshf->priv->ksigma = 1.0 / Rsigma;
  
  pshf->priv->an     = ncm_util_exp10 (1.5222 + 2.8553 * neff + 2.3706 * neff2 + 0.9903 * neff3 + 0.2250 * neff4 - 0.6038 * Cur + 0.1749 * Omega_de_onepp);
  pshf->priv->bn     = ncm_util_exp10 (-0.5642 + 0.5864 * neff + 0.5716 * neff2 - 1.5474 * Cur + 0.2279 * Omega_de_onepp);
  pshf->priv->cn     = ncm_util_exp10 (0.3698 + 2.0404 * neff + 0.8161 * neff2 + 0.5869 * Cur);
  pshf->priv->gamman = 0.1971 - 0.0843 * neff + 0.8460 * Cur;
  pshf->priv->alphan = fabs (6.0835 + 1.3373 * neff - 0.1959 * neff2 - 5.5274 * Cur);
  pshf->priv->betan  = 2.0379 - 0.7354 * neff + 0.3157 * neff2 + 1.2490 * neff3 + 0.3980 * neff4 - 0.1682 * Cur; // + fnu*(1.081 + 0.395*pow(rneff,2)
  pshf->priv->nun    = ncm_util_exp10 (5.2105 + 3.6902 * neff);

  pshf->priv->f1     = pow (Omega_m, NC_POWSPEC_MNL_HALOFIT_F1POW);
  pshf->priv->f2     = pow (Omega_m, NC_POWSPEC_MNL_HALOFIT_F2POW);
  pshf->priv->f3     = pow (Omega_m, NC_POWSPEC_MNL_HALOFIT_F3POW);
}

static gdouble
_nc_powspec_mnl_halofit_eval (NcmPowspec* powspec, NcmModel* model, const gdouble z, const gdouble k)
{
  NcHICosmo *cosmo          = NC_HICOSMO (model);
  NcPowspecMNLHaloFit *pshf = NC_POWSPEC_MNL_HALOFIT (powspec);
  const gdouble Pklin       = ncm_powspec_eval (NCM_POWSPEC (pshf->psml), model, z, k);

  if (z > pshf->zmaxnl)
  {
    return Pklin;
  }
  else
  {
    if (z != pshf->priv->z)
      _nc_powspec_mnl_halofit_preeval (pshf, cosmo, z);

    const gdouble k3      = gsl_pow_3 (k);
    const gdouble k3o2pi2 = k3 / ncm_c_2_pi_2 ();

    const gdouble Delta_lin = k3o2pi2 * Pklin;

    const gdouble y = k / pshf->priv->ksigma;

    const gdouble P_Q = Pklin * (pow (1.0 + Delta_lin, pshf->priv->betan) / (1.0 + pshf->priv->alphan * Delta_lin)) * exp (-y / 4.0 - y * y / 8.0);

    const gdouble Delta_Hprime = pshf->priv->an * pow (y, 3.0 * pshf->priv->f1) / (1.0 + pshf->priv->bn * pow (y, pshf->priv->f2) + pow (pshf->priv->cn * pshf->priv->f3 * y, 3.0 - pshf->priv->gamman));
    const gdouble Delta_H      = Delta_Hprime / (1.0 + pshf->priv->nun / (y * y));

    const gdouble P_H = Delta_H / k3o2pi2;

    return P_Q + P_H;
  }
}

static void 
_nc_powspec_mnl_halofit_eval_vec (NcmPowspec* powspec, NcmModel* model, const gdouble z, NcmVector* k, NcmVector* Pk)
{
  NcHICosmo* cosmo = NC_HICOSMO (model);
  // NcHIPrim* prim = NC_HIPRIM (ncm_model_peek_submodel_by_mid (model, nc_hiprim_id ()));
  NcPowspecMNLHaloFit* pshf = NC_POWSPEC_MNL_HALOFIT (powspec);

  ncm_powspec_eval_vec (NCM_POWSPEC (pshf->psml), model, z, k, Pk);

  if (z > pshf->zmaxnl)
  {
    return;
  }
  else
  {
    const guint len = ncm_vector_len (k);
    guint i;

    if (z != pshf->priv->z)
      _nc_powspec_mnl_halofit_preeval (pshf, cosmo, z);
    
    for (i = 0; i < len; i++)
    {
      const gdouble ki       = ncm_vector_get (k, i);
      const gdouble ki3      = gsl_pow_3 (ki);
      const gdouble ki3o2pi2 = ki3 / ncm_c_2_pi_2 ();

      const gdouble Pklin    = ncm_vector_get (Pk, i);

      const gdouble Delta_lin = ki3o2pi2 * Pklin;

      const gdouble y = ki / pshf->priv->ksigma;

      const gdouble P_Q = Pklin * (pow (1.0 + Delta_lin, pshf->priv->betan) / (1.0 + pshf->priv->alphan * Delta_lin)) * exp (-y / 4.0 - y * y / 8.0);

      const gdouble Delta_Hprime = pshf->priv->an * pow (y, 3.0 * pshf->priv->f1) / (1.0 + pshf->priv->bn * pow (y, pshf->priv->f2) + pow (pshf->priv->cn * pshf->priv->f3 * y, 3.0 - pshf->priv->gamman));
      const gdouble Delta_H      = Delta_Hprime / (1.0 + pshf->priv->nun / (y * y));

      const gdouble P_H = Delta_H / ki3o2pi2;

      ncm_vector_set (Pk, i, P_Q + P_H);
    }
  }
}

/**
 * nc_powspec_mnl_halofit_new:
 * @psml: a #NcPowspecML
 * @zmaxnl: a gdouble
 * @prec: a gdouble
 *
 * Creates a new #NcPowspecMNLHaloFit from the transfer
 * function @tf.
 *
 * Returns: (transfer full): the newly created #NcPowspecMNLHaloFit.
 */
NcPowspecMNLHaloFit*
nc_powspec_mnl_halofit_new (NcPowspecML* psml, gdouble zmaxnl, gdouble prec)
{
  NcPowspecMNLHaloFit *pshf = g_object_new (NC_TYPE_POWSPEC_MNL_HALOFIT,
                                            "power-spec", psml, 
                                            "zmaxnl", zmaxnl, 
                                            "prec", prec,
                                            NULL);

  return pshf;
}
