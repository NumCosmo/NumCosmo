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
 * SECTION:nc_hicosmo_de_wspline
 * @title: NcHICosmoDEWSpline
 * @short_description: Dark Energy -- spline equation of state
 *
 * Dark Energy equation of state: $w(\alpha)$ approximated by a cubic spline.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "model/nc_hicosmo_de_wspline.h"
#include "math/ncm_spline_cubic_notaknot.h"
#include "math/ncm_spline_gsl.h"

struct _NcHICosmoDEWSplinePrivate
{
  guint nknots;
  guint size;
  gdouble z_1;
  gdouble z_f;
  gdouble alpha_f;
  gdouble w_f;
  gdouble int_f;
  NcmSpline *w_alpha;
};

G_DEFINE_TYPE_WITH_PRIVATE (NcHICosmoDEWSpline, nc_hicosmo_de_wspline, NC_TYPE_HICOSMO_DE)

enum
{
  PROP_0,
  PROP_Z_1,
  PROP_Z_F,
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
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
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
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
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
    const gdouble alpha1                   = log1p (self->z_1);
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

    ncm_vector_set (alphav, 0, 0.0);
    ncm_vector_set (alphav, 1, alpha1);

    {
      const gdouble dalpha = (log (alphaf) - log (alpha1)) / (wz_size - 2);

      for (i = 0; i < wz_size - 2; i++)
      {
        const gdouble alpha = alpha1 * exp (dalpha * (i + 1.0));

        ncm_vector_set (alphav, 2 + i, alpha);
      }
    }

    {
      NcmSpline *s;

      if (self->nknots >= 6)
        s = ncm_spline_cubic_notaknot_new ();
      else
        s = ncm_spline_gsl_new (gsl_interp_polynomial);

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
static gdouble _nc_hicosmo_de_wspline_dw_de (NcHICosmoDE *cosmo_de, gdouble z);
static gdouble _nc_hicosmo_de_wspline_ln_rho_rho0 (NcHICosmoDE *cosmo_de, gdouble z);

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

  ncm_model_class_set_name_nick (model_class, "XCDM - Constant EOS", "XCDM");
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
  nc_hicosmo_de_set_dw_de_impl (parent_class, &_nc_hicosmo_de_wspline_dw_de);
  nc_hicosmo_de_set_ln_rho_rho0_impl (parent_class, &_nc_hicosmo_de_wspline_ln_rho_rho0);
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

static gdouble
_nc_hicosmo_de_wspline_dw_de (NcHICosmoDE *cosmo_de, gdouble z)
{
  NcHICosmoDEWSpline *wspline            = NC_HICOSMO_DE_WSPLINE (cosmo_de);
  NcHICosmoDEWSplinePrivate * const self = wspline->priv;

  _nc_hicosmo_de_wspline_prepare (wspline);

  if (z < self->z_f)
  {
    const gdouble alpha = z < self->z_f ? log1p (z) : log1p (self->z_f);
    const gdouble dw    = ncm_spline_eval_deriv (self->w_alpha, alpha);

    return dw;
  }
  else
  {
    return 0.0;
  }
}

static gdouble
_nc_hicosmo_de_wspline_ln_rho_rho0 (NcHICosmoDE *cosmo_de, gdouble z)
{
  NcHICosmoDEWSpline *wspline            = NC_HICOSMO_DE_WSPLINE (cosmo_de);
  NcHICosmoDEWSplinePrivate * const self = wspline->priv;

  _nc_hicosmo_de_wspline_prepare (wspline);

  if (z < self->z_f)
  {
    const gdouble alpha = log1p (z);

    return 3.0 * (alpha + ncm_spline_eval_integ (self->w_alpha, 0.0, alpha));
  }
  else
  {
    const gdouble alpha = log1p (z);

    return 3.0 * (alpha + self->int_f + self->w_f * (alpha - self->alpha_f));
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

