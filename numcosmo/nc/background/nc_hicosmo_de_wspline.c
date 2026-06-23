/***************************************************************************
 *            nc_hicosmo_de_wspline.c
 *
 *  Mon Oct 11 16:22:12 2021
 *  Copyright  2021  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 *  Copyright  2021  Sanderson Carlos Ribeiro
 *  <sander23.ribeiro@uel.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2021 <vitenti@uel.br>
 * Copyright (C) Sanderson Carlos Ribeiro 2021 <sander23.ribeiro@uel.br>
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
 * NcHICosmoDEWSpline:
 *
 * Dark Energy -- spline equation of state.
 *
 * Dark Energy equation of state: $w(\alpha)$ approximated by a cubic spline.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "nc/background/nc_hicosmo_de_wspline.h"
#include "nc_enum_types.h"
#include "ncm/spline/ncm_spline_cubic_notaknot.h"
#include "ncm/spline/ncm_spline_gsl.h"
#include "ncm/integration/ncm_integrate.h"
#include "ncm/core/ncm_memory_pool.h"
#include "ncm/model/ncm_mset_func_list.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_integration.h>
#endif /* NUMCOSMO_GIR_SCAN */

struct _NcHICosmoDEWSplinePrivate
{
  guint nknots;
  guint size;
  gdouble z_1;
  gdouble z_f;
  gdouble alpha_f;
  gdouble w_f;
  gdouble int_f;
  NcHICosmoSplineKnots knots;
  NcmSpline *w_alpha;
};

G_DEFINE_TYPE_WITH_PRIVATE (NcHICosmoDEWSpline, nc_hicosmo_de_wspline, NC_TYPE_HICOSMO_DE)

enum
{
  PROP_0,
  PROP_Z_1,
  PROP_Z_F,
  PROP_KNOTS,
  PROP_SIZE,
};

static void
nc_hicosmo_de_wspline_init (NcHICosmoDEWSpline *wspline)
{
  NcHICosmoDEWSplinePrivate * const self = wspline->priv = nc_hicosmo_de_wspline_get_instance_private (wspline);

  self->nknots  = 0;
  self->size    = 0;
  self->z_1     = 0.0;
  self->z_f     = 0.0;
  self->alpha_f = 0.0;
  self->w_f     = 0.0;
  self->int_f   = 0.0;
  self->w_alpha = NULL;
}

static void
_nc_hicosmo_de_wspline_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcHICosmoDEWSpline *wspline            = NC_HICOSMO_DE_WSPLINE (object);
  NcHICosmoDEWSplinePrivate * const self = wspline->priv;

  g_return_if_fail (NC_IS_HICOSMO_DE_WSPLINE (object));

  switch (prop_id)
  {
    case PROP_Z_1:
      g_value_set_double (value, self->z_1);
      break;
    case PROP_Z_F:
      g_value_set_double (value, self->z_f);
      break;
    case PROP_KNOTS:
      g_value_set_enum (value, self->knots);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_nc_hicosmo_de_wspline_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcHICosmoDEWSpline *wspline            = NC_HICOSMO_DE_WSPLINE (object);
  NcHICosmoDEWSplinePrivate * const self = wspline->priv;

  g_return_if_fail (NC_IS_HICOSMO_DE_WSPLINE (object));

  switch (prop_id)
  {
    case PROP_Z_1:
      self->z_1 = g_value_get_double (value);
      break;
    case PROP_Z_F:
      self->z_f = g_value_get_double (value);
      break;
    case PROP_KNOTS:
      self->knots = g_value_get_enum (value);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_nc_hicosmo_de_wspline_constructed (GObject *object)
{
  /* Chain up : start */
  G_OBJECT_CLASS (nc_hicosmo_de_wspline_parent_class)->constructed (object);
  {
    NcHICosmoDEWSpline *wspline            = NC_HICOSMO_DE_WSPLINE (object);
    NcHICosmoDEWSplinePrivate * const self = wspline->priv;
    NcmModel *model                        = NCM_MODEL (wspline);
    NcmVector *orig_vec                    = ncm_model_orig_params_peek_vector (model);
    NcmModelClass *model_class             = NCM_MODEL_GET_CLASS (model);
    guint wz_size                          = ncm_model_vparam_len (model, NC_HICOSMO_DE_WSPLINE_W);
    const gdouble alphaf                   = log1p (self->z_f);
    NcmVector *alphav, *wv;
    guint i, wvi;

    self->nknots  = wz_size;
    self->size    = model_class->sparam_len + self->nknots;
    self->alpha_f = alphaf * (1.0 - GSL_DBL_EPSILON);

    g_assert_cmpuint (wz_size, >, 2);

    wvi = ncm_model_vparam_index (model, NC_HICOSMO_DE_WSPLINE_W, 0);

    alphav = ncm_vector_new (wz_size);
    wv     = ncm_vector_get_subvector (orig_vec, wvi, wz_size);

    {
      /* Knot placement in alpha = ln(1+z) over [0, alpha_f]. The default
       * Chebyshev grid clusters knots toward both endpoints (good for
       * interpolation conditioning), which puts the densest knots at high z
       * where expansion data are sparsest. The uniform alternative spreads the
       * resolution evenly, trading endpoint conditioning for a less degenerate
       * high-z boundary.
       */
      for (i = 0; i < wz_size; i++)
      {
        if (self->knots == NC_HICOSMO_SPLINE_KNOTS_UNIFORM)
        {
          const gdouble alpha = alphaf * i / (wz_size - 1.0);

          ncm_vector_set (alphav, i, alpha);
        }
        else
        {
          const gdouble alpha = 0.5 * alphaf + 0.5 * alphaf * cos (M_PI * (i * 1.0) / (wz_size - 1.0));

          ncm_vector_set (alphav, wz_size - 1 - i, alpha);
        }
      }
    }

    {
      NcmSpline *s;

      if (self->nknots >= 6)
        s = NCM_SPLINE (ncm_spline_cubic_notaknot_new ());
      else
        s = NCM_SPLINE (ncm_spline_gsl_new (gsl_interp_polynomial));

      self->w_alpha = ncm_spline_new (s, alphav, wv, TRUE);

      ncm_spline_free (s);
      ncm_vector_free (alphav);
      ncm_vector_free (wv);
    }
  }
}

static void
_nc_hicosmo_de_wspline_dispose (GObject *object)
{
  NcHICosmoDEWSpline *wspline            = NC_HICOSMO_DE_WSPLINE (object);
  NcHICosmoDEWSplinePrivate * const self = wspline->priv;

  ncm_spline_clear (&self->w_alpha);

  /* Chain up : end */
  G_OBJECT_CLASS (nc_hicosmo_de_wspline_parent_class)->dispose (object);
}

static void
_nc_hicosmo_de_wspline_finalize (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (nc_hicosmo_de_wspline_parent_class)->finalize (object);
}

static gdouble _nc_hicosmo_de_wspline_E2Omega_de (NcHICosmoDE *cosmo_de, gdouble z);
static gdouble _nc_hicosmo_de_wspline_dE2Omega_de_dz (NcHICosmoDE *cosmo_de, gdouble z);
static gdouble _nc_hicosmo_de_wspline_d2E2Omega_de_dz2 (NcHICosmoDE *cosmo_de, gdouble z);
static gdouble _nc_hicosmo_de_wspline_w_de (NcHICosmoDE *cosmo_de, gdouble z);

static void
nc_hicosmo_de_wspline_class_init (NcHICosmoDEWSplineClass *klass)
{
  GObjectClass *object_class     = G_OBJECT_CLASS (klass);
  NcHICosmoDEClass *parent_class = NC_HICOSMO_DE_CLASS (klass);
  NcmModelClass *model_class     = NCM_MODEL_CLASS (klass);

  model_class->set_property = &_nc_hicosmo_de_wspline_set_property;
  model_class->get_property = &_nc_hicosmo_de_wspline_get_property;
  object_class->constructed = &_nc_hicosmo_de_wspline_constructed;
  object_class->dispose     = &_nc_hicosmo_de_wspline_dispose;
  object_class->finalize    = &_nc_hicosmo_de_wspline_finalize;

  ncm_model_class_set_name_nick (model_class, "WSpline - DE EOS reconstruction w(z)", "WSpline");
  ncm_model_class_add_params (model_class, 0,
                              NC_HICOSMO_DE_WSPLINE_VPARAM_LEN - NC_HICOSMO_DE_VPARAM_LEN,
                              PROP_SIZE);

  g_object_class_install_property (object_class,
                                   PROP_Z_1,
                                   g_param_spec_double ("z1",
                                                        NULL,
                                                        "second redshift knot",
                                                        1.0e-3, 0.5, 1.0e-2,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  g_object_class_install_property (object_class,
                                   PROP_Z_F,
                                   g_param_spec_double ("zf",
                                                        NULL,
                                                        "final redshift",
                                                        1.0, 1.0e10, 2.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  g_object_class_install_property (object_class,
                                   PROP_KNOTS,
                                   g_param_spec_enum ("knots",
                                                      NULL,
                                                      "knot placement in alpha=ln(1+z)",
                                                      NC_TYPE_HICOSMO_SPLINE_KNOTS,
                                                      NC_HICOSMO_SPLINE_KNOTS_CHEBYSHEV,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /* Set w_0 param info */
  ncm_model_class_set_vparam (model_class, NC_HICOSMO_DE_WSPLINE_W, 6, "w", "w",
                              -5.0, 1.0, 1.0e-2,
                              NC_HICOSMO_DEFAULT_PARAMS_ABSTOL, NC_HICOSMO_DE_WSPLINE_DEFAULT_W0,
                              NCM_PARAM_TYPE_FREE);

  /* Check for errors in parameters initialization */
  ncm_model_class_check_params_info (model_class);

  nc_hicosmo_de_set_E2Omega_de_impl (parent_class, &_nc_hicosmo_de_wspline_E2Omega_de);
  nc_hicosmo_de_set_dE2Omega_de_dz_impl (parent_class, &_nc_hicosmo_de_wspline_dE2Omega_de_dz);
  nc_hicosmo_de_set_d2E2Omega_de_dz2_impl (parent_class, &_nc_hicosmo_de_wspline_d2E2Omega_de_dz2);
  nc_hicosmo_de_set_w_de_impl (parent_class, &_nc_hicosmo_de_wspline_w_de);
}

#define VECTOR (NCM_MODEL (cosmo_de))
#define OMEGA_X (ncm_model_orig_param_get (VECTOR, NC_HICOSMO_DE_OMEGA_X))

static void
_nc_hicosmo_de_wspline_prepare (NcHICosmoDEWSpline *wspline)
{
  NcHICosmoDEWSplinePrivate * const self = wspline->priv;

  if (!ncm_model_lstate_is_update (NCM_MODEL (wspline), 0))
  {
    ncm_spline_prepare (self->w_alpha);

    self->w_f   = ncm_spline_eval (self->w_alpha, self->alpha_f);
    self->int_f = ncm_spline_eval_integ (self->w_alpha, 0.0, self->alpha_f);

    ncm_model_lstate_set_update (NCM_MODEL (wspline), 0);
  }
  else
  {
    return;
  }
}

static gdouble
_nc_hicosmo_de_wspline_E2Omega_de (NcHICosmoDE *cosmo_de, gdouble z)
{
  NcHICosmoDEWSpline *wspline            = NC_HICOSMO_DE_WSPLINE (cosmo_de);
  NcHICosmoDEWSplinePrivate * const self = wspline->priv;

  _nc_hicosmo_de_wspline_prepare (wspline);

  if (z < self->z_f)
  {
    const gdouble alpha = log1p (z);

    return OMEGA_X * exp (3.0 * (alpha + ncm_spline_eval_integ (self->w_alpha, 0.0, alpha)));
  }
  else
  {
    const gdouble alpha = log1p (z);

    return OMEGA_X * exp (3.0 * (alpha + self->int_f + self->w_f * (alpha - self->alpha_f)));
  }
}

static gdouble
_nc_hicosmo_de_wspline_dE2Omega_de_dz (NcHICosmoDE *cosmo_de, gdouble z)
{
  NcHICosmoDEWSpline *wspline            = NC_HICOSMO_DE_WSPLINE (cosmo_de);
  NcHICosmoDEWSplinePrivate * const self = wspline->priv;

  _nc_hicosmo_de_wspline_prepare (wspline);

  if (z < self->z_f)
  {
    const gdouble alpha             = log1p (z);
    const gdouble w                 = ncm_spline_eval (self->w_alpha, alpha);
    const gdouble int_w             = ncm_spline_eval_integ (self->w_alpha, 0.0, alpha);
    const gdouble exp_malpha_OmegaX = OMEGA_X * exp (2.0 * alpha + 3.0 * int_w);

    return 3.0 * (1.0 + w) * exp_malpha_OmegaX;
  }
  else
  {
    const gdouble alpha = log1p (z);

    return OMEGA_X * 3.0 * (1.0 + self->w_f) * exp (2.0 * alpha + 3.0 * (self->int_f + self->w_f * (alpha - self->alpha_f)));
  }
}

static gdouble
_nc_hicosmo_de_wspline_d2E2Omega_de_dz2 (NcHICosmoDE *cosmo_de, gdouble z)
{
  NcHICosmoDEWSpline *wspline            = NC_HICOSMO_DE_WSPLINE (cosmo_de);
  NcHICosmoDEWSplinePrivate * const self = wspline->priv;

  _nc_hicosmo_de_wspline_prepare (wspline);

  if (z < self->z_f)
  {
    const gdouble alpha              = log1p (z);
    const gdouble w                  = ncm_spline_eval (self->w_alpha, alpha);
    const gdouble dw                 = ncm_spline_eval_deriv (self->w_alpha, alpha);
    const gdouble int_w              = ncm_spline_eval_integ (self->w_alpha, 0.0, alpha);
    const gdouble exp_m2alpha_OmegaX = OMEGA_X * exp (1.0 * alpha + 3.0 * int_w);

    return 3.0 * (2.0 + w * (5.0 + 3.0 * w) + dw) * exp_m2alpha_OmegaX;
  }
  else
  {
    const gdouble alpha              = log1p (z);
    const gdouble exp_m2alpha_OmegaX = OMEGA_X * exp (1.0 * alpha + 3.0 * (self->int_f + self->w_f * (alpha - self->alpha_f)));

    return 3.0 * (2.0 + self->w_f * (5.0 + 3.0 * self->w_f)) * exp_m2alpha_OmegaX;
  }
}

static gdouble
_nc_hicosmo_de_wspline_w_de (NcHICosmoDE *cosmo_de, gdouble z)
{
  NcHICosmoDEWSpline *wspline            = NC_HICOSMO_DE_WSPLINE (cosmo_de);
  NcHICosmoDEWSplinePrivate * const self = wspline->priv;

  _nc_hicosmo_de_wspline_prepare (wspline);

  if (z < self->z_f)
  {
    const gdouble alpha = z < self->z_f ? log1p (z) : log1p (self->z_f);
    const gdouble w     = ncm_spline_eval (self->w_alpha, alpha);

    return w;
  }
  else
  {
    return self->w_f;
  }
}

/**
 * nc_hicosmo_de_wspline_new:
 * @nknots: number of knots in the $w(z)$ spline
 * @z_f: redshift $z_f$ of the last knot.
 *
 * This function instantiates a new object of type #NcHICosmoDEWSpline.
 *
 * Returns: A new #NcHICosmoDEWSpline
 */
NcHICosmoDEWSpline *
nc_hicosmo_de_wspline_new (gsize nknots, const gdouble z_f)
{
  NcHICosmoDEWSpline *wspline = g_object_new (NC_TYPE_HICOSMO_DE_WSPLINE,
                                              "zf", z_f,
                                              "w-length", nknots,
                                              NULL);

  return wspline;
}

/**
 * nc_hicosmo_de_wspline_get_alpha:
 * @wspline: a #NcHICosmoDEWSpline
 *
 * Gets vector with the current knots in $\alpha$.
 *
 * Returns: (transfer none): A new #NcmVector
 */
NcmVector *
nc_hicosmo_de_wspline_get_alpha (NcHICosmoDEWSpline *wspline)
{
  NcHICosmoDEWSplinePrivate * const self = wspline->priv;

  return ncm_spline_peek_xv (self->w_alpha);
}

/**
 * nc_hicosmo_de_wspline_lp_norm:
 * @wspline: a #NcHICosmoDEWSpline
 * @ctype: a #NcmSplineCurvatureType
 * @p: the norm order $p > 0$
 *
 * Computes the domain-normalized $L_p$ norm of the $w(\alpha)$ curvature density
 * (selected by @ctype) over the active range $[0, \alpha_f]$. The case @p = 2
 * yields the root-mean-square curvature; large @p approaches the maximum (see
 * ncm_spline_curvature_lp_norm()).
 *
 * Returns: the curvature $L_p$ norm.
 */
gdouble
nc_hicosmo_de_wspline_lp_norm (NcHICosmoDEWSpline *wspline, NcmSplineCurvatureType ctype, const gdouble p)
{
  NcHICosmoDEWSplinePrivate * const self = wspline->priv;

  _nc_hicosmo_de_wspline_prepare (wspline);

  return ncm_spline_curvature_lp_norm (self->w_alpha, ctype, p, 0.0, self->alpha_f);
}

/**
 * nc_hicosmo_de_wspline_weighted_lp_norm:
 * @wspline: a #NcHICosmoDEWSpline
 * @ctype: a #NcmSplineCurvatureType
 * @p: the norm order $p > 0$
 * @weight: a #NcmSpline holding the non-negative weight density $W(\alpha)$
 *
 * Computes the weight-normalized $L_p$ norm of the $w(\alpha)$ curvature density
 * (selected by @ctype) over the active range $[0, \alpha_f]$, see
 * ncm_spline_curvature_weighted_lp_norm(). The @weight abscissa is
 * $\alpha = \ln(1 + z)$, matching the internal $w(\alpha)$ spline; a large
 * $W(\alpha)$ penalizes curvature near that $\alpha$, a small one tolerates it.
 *
 * Returns: the weighted curvature $L_p$ norm.
 */
gdouble
nc_hicosmo_de_wspline_weighted_lp_norm (NcHICosmoDEWSpline *wspline, NcmSplineCurvatureType ctype, const gdouble p, NcmSpline *weight)
{
  NcHICosmoDEWSplinePrivate * const self = wspline->priv;

  _nc_hicosmo_de_wspline_prepare (wspline);
  ncm_spline_prepare (weight);

  return ncm_spline_curvature_weighted_lp_norm (self->w_alpha, ctype, p, weight, 0.0, self->alpha_f);
}

/**
 * nc_hicosmo_de_wspline_mean_kappa:
 * @wspline: a #NcHICosmoDEWSpline
 *
 * Gets the mean value of $w(z)$ curvature, i.e. the root-mean-square geometric
 * curvature of $w(\alpha)$ over $[0, \alpha_f]$ (the $p = 2$,
 * #NCM_SPLINE_CURVATURE_GEOMETRIC case of nc_hicosmo_de_wspline_lp_norm()).
 *
 * Returns: The mean value of $w(z)$ curvature
 */
gdouble
nc_hicosmo_de_wspline_mean_kappa (NcHICosmoDEWSpline *wspline)
{
  return nc_hicosmo_de_wspline_lp_norm (wspline, NCM_SPLINE_CURVATURE_GEOMETRIC, 2.0);
}

static void
_nc_hicosmo_de_wspline_mean_kappa (NcmMSetFuncList *flist, NcmMSet *mset, const gdouble *x, gdouble *f)
{
  NcHICosmoDEWSpline *cosmo = NC_HICOSMO_DE_WSPLINE (ncm_mset_peek (mset, nc_hicosmo_id ()));

  g_assert (NC_IS_HICOSMO_DE_WSPLINE (cosmo));
  f[0] = nc_hicosmo_de_wspline_mean_kappa (cosmo);
}

static void
_nc_hicosmo_de_wspline_lp_kappa (NcmMSetFuncList *flist, NcmMSet *mset, const gdouble *x, gdouble *f)
{
  NcHICosmoDEWSpline *cosmo = NC_HICOSMO_DE_WSPLINE (ncm_mset_peek (mset, nc_hicosmo_id ()));

  g_assert (NC_IS_HICOSMO_DE_WSPLINE (cosmo));
  f[0] = nc_hicosmo_de_wspline_lp_norm (cosmo, NCM_SPLINE_CURVATURE_GEOMETRIC, x[0]);
}

static void
_nc_hicosmo_de_wspline_lp_w2 (NcmMSetFuncList *flist, NcmMSet *mset, const gdouble *x, gdouble *f)
{
  NcHICosmoDEWSpline *cosmo = NC_HICOSMO_DE_WSPLINE (ncm_mset_peek (mset, nc_hicosmo_id ()));

  g_assert (NC_IS_HICOSMO_DE_WSPLINE (cosmo));
  f[0] = nc_hicosmo_de_wspline_lp_norm (cosmo, NCM_SPLINE_CURVATURE_D2, x[0]);
}

static void
_nc_hicosmo_de_wspline_wlp_kappa (NcmMSetFuncList *flist, NcmMSet *mset, const gdouble *x, gdouble *f)
{
  NcHICosmoDEWSpline *cosmo = NC_HICOSMO_DE_WSPLINE (ncm_mset_peek (mset, nc_hicosmo_id ()));
  NcmSpline *weight         = NCM_SPLINE (ncm_mset_func_list_peek_obj (flist));

  g_assert (NC_IS_HICOSMO_DE_WSPLINE (cosmo));
  f[0] = nc_hicosmo_de_wspline_weighted_lp_norm (cosmo, NCM_SPLINE_CURVATURE_GEOMETRIC, x[0], weight);
}

static void
_nc_hicosmo_de_wspline_wlp_w2 (NcmMSetFuncList *flist, NcmMSet *mset, const gdouble *x, gdouble *f)
{
  NcHICosmoDEWSpline *cosmo = NC_HICOSMO_DE_WSPLINE (ncm_mset_peek (mset, nc_hicosmo_id ()));
  NcmSpline *weight         = NCM_SPLINE (ncm_mset_func_list_peek_obj (flist));

  g_assert (NC_IS_HICOSMO_DE_WSPLINE (cosmo));
  f[0] = nc_hicosmo_de_wspline_weighted_lp_norm (cosmo, NCM_SPLINE_CURVATURE_D2, x[0], weight);
}

void
_nc_hicosmo_de_wspline_register_functions (void)
{
  ncm_mset_func_list_register ("mean_kappa", "\\bar{\\kappa}", "NcHICosmoDEWSpline", "Mean w geometric curvature (L2)", G_TYPE_NONE, _nc_hicosmo_de_wspline_mean_kappa, 0, 1);
  ncm_mset_func_list_register ("lp_kappa", "\\Vert\\kappa\\Vert_p", "NcHICosmoDEWSpline", "Lp norm of w geometric curvature, x[0]=p", G_TYPE_NONE, _nc_hicosmo_de_wspline_lp_kappa, 1, 1);
  ncm_mset_func_list_register ("lp_w2", "\\Vert w''\\Vert_p", "NcHICosmoDEWSpline", "Lp norm of w'' (second derivative), x[0]=p", G_TYPE_NONE, _nc_hicosmo_de_wspline_lp_w2, 1, 1);
  ncm_mset_func_list_register ("wlp_kappa", "\\Vert\\kappa\\Vert_{p,W}", "NcHICosmoDEWSpline", "Weighted Lp norm of w geometric curvature, x[0]=p, obj=weight spline W(alpha)", NCM_TYPE_SPLINE, _nc_hicosmo_de_wspline_wlp_kappa, 1, 1);
  ncm_mset_func_list_register ("wlp_w2", "\\Vert w''\\Vert_{p,W}", "NcHICosmoDEWSpline", "Weighted Lp norm of w'' (second derivative), x[0]=p, obj=weight spline W(alpha)", NCM_TYPE_SPLINE, _nc_hicosmo_de_wspline_wlp_w2, 1, 1);
}

