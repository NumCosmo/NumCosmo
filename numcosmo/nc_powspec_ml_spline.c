/***************************************************************************
 *            nc_powspec_ml_spline.c
 *
 *  Sun Jun 18 10:26:38 2017
 *  Copyright  2017  Mariana Penna Lima
 *  <pennalima@gmail.com>
 ****************************************************************************/
/*
 * nc_powspec_ml_spline.c
 * Copyright (C) 2017 Mariana Penna Lima <pennalima@gmail.com>
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
 * NcPowspecMLSpline:
 *
 * Class for linear matter power spectrum from a 1d-spline.
 *
 * Provides a linear matter power spectrum as a function of mode $k$ at a redshift $z_0$
 * computing a spline, whose knots and the respective $P(k, z = z_0)$ values are
 * provided in a file. The time evolution of the power spectrum is computed using the
 * growth function.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "nc_powspec_ml_spline.h"
#include "nc_hiprim.h"

struct _NcPowspecMLSpline
{
  /*< private > */
  NcPowspecML parent_instance;
  NcmSerialize *ser;
  NcmSpline *Pk;
  NcGrowthFunc *gf;
};

enum
{
  PROP_0,
  PROP_SPLINE,
  PROP_SIZE
};

G_DEFINE_TYPE (NcPowspecMLSpline, nc_powspec_ml_spline, NC_TYPE_POWSPEC_ML)

static void
nc_powspec_ml_spline_init (NcPowspecMLSpline *ps_fs)
{
  ps_fs->ser = ncm_serialize_new (0);
  ps_fs->Pk  = NULL;
  ps_fs->gf  = nc_growth_func_new ();
}

static void
_nc_powspec_ml_spline_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcPowspecMLSpline *ps_fs = NC_POWSPEC_ML_SPLINE (object);

  g_return_if_fail (NC_IS_POWSPEC_ML_SPLINE (object));

  switch (prop_id)
  {
    case PROP_SPLINE:
      nc_powspec_ml_spline_set_spline (ps_fs, g_value_get_object (value));
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_nc_powspec_ml_spline_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcPowspecMLSpline *ps_fs = NC_POWSPEC_ML_SPLINE (object);

  g_return_if_fail (NC_IS_POWSPEC_ML_SPLINE (object));

  switch (prop_id)
  {
    case PROP_SPLINE:
      g_value_set_object (value, nc_powspec_ml_spline_peek_spline (ps_fs));
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_nc_powspec_ml_spline_constructed (GObject *object)
{
  /* Chain up : start */
  G_OBJECT_CLASS (nc_powspec_ml_spline_parent_class)->constructed (object);
}

static void
_nc_powspec_ml_spline_dispose (GObject *object)
{
  NcPowspecMLSpline *ps_fs = NC_POWSPEC_ML_SPLINE (object);

  ncm_spline_clear (&ps_fs->Pk);
  ncm_serialize_clear (&ps_fs->ser);
  nc_growth_func_clear (&ps_fs->gf);

  /* Chain up : end */
  G_OBJECT_CLASS (nc_powspec_ml_spline_parent_class)->dispose (object);
}

static void
_nc_powspec_ml_spline_finalize (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (nc_powspec_ml_spline_parent_class)->finalize (object);
}

static void _nc_powspec_ml_spline_prepare (NcmPowspec *powspec, NcmModel *model);
static gdouble _nc_powspec_ml_spline_eval (NcmPowspec *powspec, NcmModel *model, const gdouble z, const gdouble k);
static void _nc_powspec_ml_spline_eval_vec (NcmPowspec *powspec, NcmModel *model, const gdouble z, NcmVector *k, NcmVector *Pk);
static void _nc_powspec_ml_spline_get_nknots (NcmPowspec *powspec, guint *Nz, guint *Nk);

static void
nc_powspec_ml_spline_class_init (NcPowspecMLSplineClass *klass)
{
  GObjectClass *object_class     = G_OBJECT_CLASS (klass);
  NcmPowspecClass *powspec_class = NCM_POWSPEC_CLASS (klass);

  object_class->set_property = &_nc_powspec_ml_spline_set_property;
  object_class->get_property = &_nc_powspec_ml_spline_get_property;

  object_class->constructed = &_nc_powspec_ml_spline_constructed;

  object_class->dispose  = &_nc_powspec_ml_spline_dispose;
  object_class->finalize = &_nc_powspec_ml_spline_finalize;

  /**
   * NcPowspecMLSpline:filename:
   *
   * The file name containing $P(k, z=z_0)$ in a serialized #NcmSpline format.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_SPLINE,
                                   g_param_spec_object ("spline",
                                                        NULL,
                                                        "Spline representing Pk at z=0",
                                                        NCM_TYPE_SPLINE,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));


  powspec_class->prepare    = &_nc_powspec_ml_spline_prepare;
  powspec_class->eval       = &_nc_powspec_ml_spline_eval;
  powspec_class->eval_vec   = &_nc_powspec_ml_spline_eval_vec;
  powspec_class->get_nknots = &_nc_powspec_ml_spline_get_nknots;
}

static void
_nc_powspec_ml_spline_prepare (NcmPowspec *powspec, NcmModel *model)
{
  NcHICosmo *cosmo         = NC_HICOSMO (model);
  NcPowspecMLSpline *ps_fs = NC_POWSPEC_ML_SPLINE (powspec);

  g_assert (NC_IS_HICOSMO (cosmo));
  g_assert (ncm_model_peek_submodel_by_mid (model, nc_hiprim_id ()) != NULL);

  nc_growth_func_prepare_if_needed (ps_fs->gf, cosmo);
  ncm_spline_prepare (ps_fs->Pk);
}

static gdouble
_nc_powspec_ml_spline_eval (NcmPowspec *powspec, NcmModel *model, const gdouble z, const gdouble k)
{
  NcHICosmo *cosmo         = NC_HICOSMO (model);
  NcPowspecMLSpline *ps_fs = NC_POWSPEC_ML_SPLINE (powspec);
  const gdouble growth     = nc_growth_func_eval (ps_fs->gf, cosmo, z);
  const gdouble gf2        = gsl_pow_2 (growth);

  return ncm_spline_eval (ps_fs->Pk, k) * gf2;
}

static void
_nc_powspec_ml_spline_eval_vec (NcmPowspec *powspec, NcmModel *model, const gdouble z, NcmVector *k, NcmVector *Pk)
{
  const guint len = ncm_vector_len (k);
  guint i;

  for (i = 0; i < len; i++)
  {
    const gdouble ki   = ncm_vector_get (k, i);
    const gdouble Pk_z = _nc_powspec_ml_spline_eval (powspec, model, z, ki);

    ncm_vector_set (Pk, i, Pk_z);
  }
}

static void
_nc_powspec_ml_spline_get_nknots (NcmPowspec *powspec, guint *Nz, guint *Nk)
{
  NcPowspecMLSpline *ps_fs = NC_POWSPEC_ML_SPLINE (powspec);
  NcmPowspec *powspec_base = NCM_POWSPEC (ps_fs);
  const gdouble zi         = ncm_powspec_get_zi (powspec_base);
  const gdouble zf         = ncm_powspec_get_zf (powspec_base);

  Nz[0] = GSL_MAX ((zf - zi) * 50, 10);
  Nk[0] = ncm_spline_get_len (ps_fs->Pk);
}

/**
 * nc_powspec_ml_spline_new:
 * @Pk: a #NcmSpline
 *
 * Creates a new #NcPowspecMLSpline from a spline @Pk of the power spectrum, $P(k,
 * z=z_0)$ at a given redshift $z_0$.
 *
 * Returns: (transfer full): the newly created #NcPowspecMLSpline.
 */
NcPowspecMLSpline *
nc_powspec_ml_spline_new (NcmSpline *Pk)
{
  NcPowspecMLSpline *ps_fs = g_object_new (NC_TYPE_POWSPEC_ML_SPLINE,
                                           "spline", Pk,
                                           NULL);

  return ps_fs;
}

/**
 * nc_powspec_ml_spline_set_spline:
 * @ps_fs: a #NcPowspecMLSpline
 * @Pk: a #NcmSpline
 *
 * Attributes a spline of the power spectrum at a redshift $z_0$, P(k, z_0), from @Pk.
 *
 */
void
nc_powspec_ml_spline_set_spline (NcPowspecMLSpline *ps_fs, NcmSpline *Pk)
{
  ncm_spline_clear (&ps_fs->Pk);
  ps_fs->Pk = ncm_spline_ref (Pk);
}

/**
 * nc_powspec_ml_spline_peek_spline:
 * @ps_fs: a #NcPowspecMLSpline
 *
 * Returns: (transfer none): the spline of the power spectrum at a redshift $z_0$, P(k,
 * z_0).
 */
NcmSpline *
nc_powspec_ml_spline_peek_spline (NcPowspecMLSpline *ps_fs)
{
  return ps_fs->Pk;
}

