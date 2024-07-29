/***************************************************************************
 *            nc_powspec_ml_fix_spline.c
 *
 *  Sun Jun 18 10:26:38 2017
 *  Copyright  2017  Mariana Penna Lima
 *  <pennalima@gmail.com>
 ****************************************************************************/
/*
 * nc_powspec_ml_fix_spline.c
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
 * SECTION:nc_powspec_ml_fix_spline
 * @title: NcPowspecMLFixSpline
 * @short_description: Class for linear matter power spectrum from a fixed spline.
 * @stability: Stable
 * @include: numcosmo/nc_powspec_ml_fix_spline.h
 *
 * Provides a linear matter power spectrum as a function of mode $k$ at a redshift $z_0$ computing
 * a spline, whose knots and the respective $P(k, z = z_0)$ values are provided in a file.
 *
 *
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "nc_powspec_ml_fix_spline.h"
#include "nc_hiprim.h"

enum
{
  PROP_0,
  PROP_FILENAME,
  PROP_SIZE
};

G_DEFINE_TYPE (NcPowspecMLFixSpline, nc_powspec_ml_fix_spline, NC_TYPE_POWSPEC_ML);

static void
nc_powspec_ml_fix_spline_init (NcPowspecMLFixSpline *ps_fs)
{
  ps_fs->ser      = ncm_serialize_new (0);
  ps_fs->Pk       = NULL;
  ps_fs->gf       = nc_growth_func_new ();
  ps_fs->filename = NULL;
}

static void
_nc_powspec_ml_fix_spline_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcPowspecMLFixSpline *ps_fs = NC_POWSPEC_ML_FIX_SPLINE (object);
  
  g_return_if_fail (NC_IS_POWSPEC_ML_FIX_SPLINE (object));
  
  switch (prop_id)
  {
    case PROP_FILENAME:
      nc_powspec_ml_fix_spline_set_file (ps_fs, g_value_get_string (value));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_powspec_ml_fix_spline_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcPowspecMLFixSpline *ps_fs = NC_POWSPEC_ML_FIX_SPLINE (object);
  
  g_return_if_fail (NC_IS_POWSPEC_ML_FIX_SPLINE (object));
  
  switch (prop_id)
  {
    case PROP_FILENAME:
      g_value_set_string (value, ps_fs->filename);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_powspec_ml_fix_spline_constructed (GObject *object)
{
  /* Chain up : start */
  G_OBJECT_CLASS (nc_powspec_ml_fix_spline_parent_class)->constructed (object);
}

static void
_nc_powspec_ml_fix_spline_dispose (GObject *object)
{
  NcPowspecMLFixSpline *ps_fs = NC_POWSPEC_ML_FIX_SPLINE (object);
  
  g_clear_pointer (&ps_fs->filename, g_free);
  ncm_spline_clear (&ps_fs->Pk);
  ncm_serialize_clear (&ps_fs->ser);
  nc_growth_func_clear (&ps_fs->gf);
  
  /* Chain up : end */
  G_OBJECT_CLASS (nc_powspec_ml_fix_spline_parent_class)->dispose (object);
}

static void
_nc_powspec_ml_fix_spline_finalize (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (nc_powspec_ml_fix_spline_parent_class)->finalize (object);
}

static void _nc_powspec_ml_fix_spline_prepare (NcmPowspec *powspec, NcmModel *model);
static gdouble _nc_powspec_ml_fix_spline_eval (NcmPowspec *powspec, NcmModel *model, const gdouble z, const gdouble k);
static void _nc_powspec_ml_fix_spline_eval_vec (NcmPowspec *powspec, NcmModel *model, const gdouble z, NcmVector *k, NcmVector *Pk);
static void _nc_powspec_ml_fix_spline_get_nknots (NcmPowspec *powspec, guint *Nz, guint *Nk);

static void
nc_powspec_ml_fix_spline_class_init (NcPowspecMLFixSplineClass *klass)
{
  GObjectClass *object_class     = G_OBJECT_CLASS (klass);
  NcmPowspecClass *powspec_class = NCM_POWSPEC_CLASS (klass);
  
  object_class->set_property = &_nc_powspec_ml_fix_spline_set_property;
  object_class->get_property = &_nc_powspec_ml_fix_spline_get_property;
  
  object_class->constructed = &_nc_powspec_ml_fix_spline_constructed;
  
  object_class->dispose  = &_nc_powspec_ml_fix_spline_dispose;
  object_class->finalize = &_nc_powspec_ml_fix_spline_finalize;
  
  /**
   * NcPowspecMLFixSpline:filename:
   *
   * The file name containing $P(k, z=z_0)$ in a serialized #NcmSpline format.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_FILENAME,
                                   g_param_spec_string ("filename",
                                                        NULL,
                                                        "Filename",
                                                        NULL,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  
  
  powspec_class->prepare    = &_nc_powspec_ml_fix_spline_prepare;
  powspec_class->eval       = &_nc_powspec_ml_fix_spline_eval;
  powspec_class->eval_vec   = &_nc_powspec_ml_fix_spline_eval_vec;
  powspec_class->get_nknots = &_nc_powspec_ml_fix_spline_get_nknots;
}

static void
_nc_powspec_ml_fix_spline_prepare (NcmPowspec *powspec, NcmModel *model)
{
  NcHICosmo *cosmo            = NC_HICOSMO (model);
  NcPowspecMLFixSpline *ps_fs = NC_POWSPEC_ML_FIX_SPLINE (powspec);
  
  g_assert (NC_IS_HICOSMO (cosmo));
  g_assert (ncm_model_peek_submodel_by_mid (model, nc_hiprim_id ()) != NULL);
  
  nc_growth_func_prepare_if_needed (ps_fs->gf, cosmo);
  ncm_spline_prepare (ps_fs->Pk);
}

static gdouble
_nc_powspec_ml_fix_spline_eval (NcmPowspec *powspec, NcmModel *model, const gdouble z, const gdouble k)
{
  NcHICosmo *cosmo            = NC_HICOSMO (model);
  NcHIPrim *prim              = NC_HIPRIM (ncm_model_peek_submodel_by_mid (model, nc_hiprim_id ()));
  NcPowspecMLFixSpline *ps_fs = NC_POWSPEC_ML_FIX_SPLINE (powspec);
  const gdouble growth        = nc_growth_func_eval (ps_fs->gf, cosmo, z);
  const gdouble gf2           = gsl_pow_2 (growth);
  
  NCM_UNUSED (prim);
  
  return ncm_spline_eval (ps_fs->Pk, k) * gf2;
}

static void
_nc_powspec_ml_fix_spline_eval_vec (NcmPowspec *powspec, NcmModel *model, const gdouble z, NcmVector *k, NcmVector *Pk)
{
  const guint len = ncm_vector_len (k);
  guint i;
  
  for (i = 0; i < len; i++)
  {
    const gdouble ki   = ncm_vector_get (k, i);
    const gdouble Pk_z = _nc_powspec_ml_fix_spline_eval (powspec, model, z, ki);
    
    ncm_vector_set (Pk, i, Pk_z);
  }
}

static void
_nc_powspec_ml_fix_spline_get_nknots (NcmPowspec *powspec, guint *Nz, guint *Nk)
{
  Nz[0] = 20;
  Nk[0] = 1000;
}

/**
 * nc_powspec_ml_fix_spline_new:
 * @filename: a file containing a serialized #NcmSpline
 *
 * Creates a new #NcPowspecMLFixSpline from a spline of the power spectrum,
 * $P(k, z=z_0)$ at a given redshift $z_0$, which is present in the serialized object @filename.
 *
 * Returns: (transfer full): the newly created #NcPowspecMLFixSpline.
 */
NcPowspecMLFixSpline *
nc_powspec_ml_fix_spline_new (const gchar *filename)
{
  NcPowspecMLFixSpline *ps_fs = g_object_new (NC_TYPE_POWSPEC_ML_FIX_SPLINE,
                                              "filename", filename,
                                              NULL);
  
  return ps_fs;
}

/**
 * nc_powspec_ml_fix_spline_set_file:
 * @ps_fs: a #NcPowspecMLFixSpline
 * @filename: a file containing a serialized #NcmSpline
 *
 * Attributes a spline of the power spectrum at a redshift $z_0$, P(k, z_0), from @filename.
 *
 */
void
nc_powspec_ml_fix_spline_set_file (NcPowspecMLFixSpline *ps_fs, const gchar *filename)
{
  g_assert (filename != NULL);
  {
    GObject *obj = ncm_serialize_from_file (ps_fs->ser, filename);
    
    g_assert (NCM_IS_SPLINE (obj));
    
    g_clear_pointer (&ps_fs->filename, g_free);
    ncm_spline_clear (&ps_fs->Pk);
    
    ps_fs->filename = g_strdup (filename);
    ps_fs->Pk       = NCM_SPLINE (obj);
  }
}

