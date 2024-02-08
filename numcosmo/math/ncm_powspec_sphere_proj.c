/***************************************************************************
 *            ncm_powspec_sphere_proj.c
 *
 *  Mon April 01 11:07:28 2019
 *  Copyright  2019  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * ncm_powspec_sphere_proj.c
 * Copyright (C) 2019 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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
 * SECTION:ncm_powspec_sphere_proj
 * @title: NcmPowspecSphereProj
 * @short_description: Class to compute spherical projection of power spectra
 *
 * The #NcmPowspecSphereProj class computes the spherical projection of power spectra.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_powspec_sphere_proj.h"
#include "math/ncm_spline_cubic_notaknot.h"
#include "math/ncm_spline2d_bicubic.h"
#include "math/ncm_fftlog_sbessel_jljm.h"
#include "math/ncm_c.h"
#include "ncm_enum_types.h"

typedef struct _NcmPowspecSphereProjPrivate
{
  NcmPowspec *ps;
  GArray *w_array;
  GPtrArray *fftlog;
  gdouble lnr0;
  gdouble lnk0;
  gdouble Lk;
  gdouble xi_i;
  gdouble xi_f;
  gdouble k_pivot;
  gdouble Pk_pivot;
  gboolean calibrated;
  gdouble reltol;
  gdouble reltol_z;
  guint ell_min;
  guint ell_max;
  GPtrArray *lnk_array;
  GPtrArray *k2Pk_array;
  GPtrArray *Cell_spline;
  NcmModelCtrl *ctrl;
  gboolean constructed;
} NcmPowspecSphereProjPrivate;

enum
{
  PROP_0,
  PROP_TYPE,
  PROP_XI_I,
  PROP_XI_F,
  PROP_K_PIVOT,
  PROP_RELTOL,
  PROP_RELTOL_Z,
  PROP_POWERSPECTRUM,
  PROP_ELL_MIN,
  PROP_ELL_MAX,
  PROP_SIZE,
};

struct _NcmPowspecSphereProj
{
  GObject parent_instance;
};

G_DEFINE_TYPE_WITH_PRIVATE (NcmPowspecSphereProj, ncm_powspec_sphere_proj, G_TYPE_OBJECT)

static void
ncm_powspec_sphere_proj_init (NcmPowspecSphereProj *psp)
{
  NcmPowspecSphereProjPrivate * const self = ncm_powspec_sphere_proj_get_instance_private (psp);

  self->ps          = NULL;
  self->w_array     = g_array_new (FALSE, FALSE, sizeof (gdouble));
  self->lnr0        = 0.0;
  self->lnk0        = 0.0;
  self->Lk          = 0.0;
  self->xi_i        = 0.0;
  self->xi_f        = 0.0;
  self->k_pivot     = 0.0;
  self->Pk_pivot    = 0.0;
  self->reltol      = 0.0;
  self->reltol_z    = 0.0;
  self->lnk_array   = g_ptr_array_new ();
  self->k2Pk_array  = g_ptr_array_new ();
  self->fftlog      = g_ptr_array_new ();
  self->Cell_spline = g_ptr_array_new ();
  self->calibrated  = FALSE;
  self->ctrl        = ncm_model_ctrl_new (NULL);
  self->constructed = FALSE;

  {
    gint i;

    for (i = 0; i < 10; i++)
    {
      gdouble w = exp (log (0.90) / 9.0 * i);

      g_array_append_val (self->w_array, w);
    }
  }

  g_ptr_array_set_free_func (self->lnk_array,   (GDestroyNotify) ncm_vector_free);
  g_ptr_array_set_free_func (self->k2Pk_array,  (GDestroyNotify) ncm_vector_free);
  g_ptr_array_set_free_func (self->fftlog,      (GDestroyNotify) g_ptr_array_unref);
  g_ptr_array_set_free_func (self->Cell_spline, (GDestroyNotify) ncm_spline2d_free);
}

static void
_ncm_powspec_sphere_proj_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmPowspecSphereProj *psp                = NCM_POWSPEC_SPHERE_PROJ (object);
  NcmPowspecSphereProjPrivate * const self = ncm_powspec_sphere_proj_get_instance_private (psp);

  g_return_if_fail (NCM_IS_POWSPEC_SPHERE_PROJ (object));

  switch (prop_id)
  {
    case PROP_XI_I:
      ncm_powspec_sphere_proj_set_xi_i (psp, g_value_get_double (value));
      break;
    case PROP_XI_F:
      ncm_powspec_sphere_proj_set_xi_f (psp, g_value_get_double (value));
      break;
    case PROP_K_PIVOT:
      ncm_powspec_sphere_proj_set_k_pivot (psp, g_value_get_double (value));
      break;
    case PROP_RELTOL:
      self->reltol = g_value_get_double (value);
      break;
    case PROP_RELTOL_Z:
      self->reltol_z = g_value_get_double (value);
      break;
    case PROP_POWERSPECTRUM:
      self->ps = g_value_dup_object (value);
      break;
    case PROP_ELL_MIN:
      ncm_powspec_sphere_proj_set_ell_min (psp, g_value_get_uint (value));
      break;
    case PROP_ELL_MAX:
      ncm_powspec_sphere_proj_set_ell_max (psp, g_value_get_uint (value));
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_ncm_powspec_sphere_proj_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmPowspecSphereProj *psp                = NCM_POWSPEC_SPHERE_PROJ (object);
  NcmPowspecSphereProjPrivate * const self = ncm_powspec_sphere_proj_get_instance_private (psp);

  g_return_if_fail (NCM_IS_POWSPEC_SPHERE_PROJ (object));

  switch (prop_id)
  {
    case PROP_XI_I:
      g_value_set_double (value, self->xi_i);
      break;
    case PROP_XI_F:
      g_value_set_double (value, self->xi_f);
      break;
    case PROP_RELTOL:
      g_value_set_double (value, self->reltol);
      break;
    case PROP_K_PIVOT:
      g_value_set_double (value, self->k_pivot);
      break;
    case PROP_RELTOL_Z:
      g_value_set_double (value, self->reltol_z);
      break;
    case PROP_POWERSPECTRUM:
      g_value_set_object (value, self->ps);
      break;
    case PROP_ELL_MIN:
      g_value_set_uint (value, ncm_powspec_sphere_proj_get_ell_min (psp));
      break;
    case PROP_ELL_MAX:
      g_value_set_uint (value, ncm_powspec_sphere_proj_get_ell_max (psp));
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_ncm_powspec_sphere_proj_adjust_fftlog_array (NcmPowspecSphereProj *psp, guint ell_min, guint ell_max, gboolean reset)
{
  NcmPowspecSphereProjPrivate * const self = ncm_powspec_sphere_proj_get_instance_private (psp);
  guint i, w_i;

  g_assert_cmpint (ell_min, <=, ell_max);

  if (reset)
    g_ptr_array_set_size (self->fftlog, 0);

  if (self->fftlog->len < self->w_array->len)
  {
    for (i = self->fftlog->len; i < self->w_array->len; i++)
    {
      GPtrArray *fftlog_w = g_ptr_array_new ();

      g_ptr_array_set_free_func (fftlog_w, (GDestroyNotify) ncm_fftlog_free);
      g_ptr_array_add (self->fftlog, fftlog_w);
    }
  }
  else if (self->fftlog->len > self->w_array->len)
  {
    g_ptr_array_set_size (self->fftlog, self->w_array->len);
  }

  for (w_i = 0; w_i < self->w_array->len; w_i++)
  {
    const gdouble w     = g_array_index (self->w_array, gdouble, w_i);
    GPtrArray *fftlog_w = g_ptr_array_index (self->fftlog, w_i);

    if (fftlog_w->len == 0)
    {
      guint nell = ell_max - ell_min;

      self->ell_min = ell_min;
      self->ell_max = ell_max;

      for (i = 0; i <= nell; i++)
      {
        const guint ell            = i + self->ell_min;
        NcmFftlogSBesselJLJM *jljm = ncm_fftlog_sbessel_jljm_new (ell, 0, log (w), self->lnr0, self->lnk0, self->Lk, 100);

        ncm_fftlog_use_smooth_padding (NCM_FFTLOG (jljm), TRUE);
        ncm_fftlog_set_smooth_padding_scale (NCM_FFTLOG (jljm), -180.0);
        /*ncm_fftlog_set_padding (NCM_FFTLOG (jljm), 0.0); */
        ncm_fftlog_sbessel_jljm_set_best_lnr0 (jljm);
        g_ptr_array_add (fftlog_w, jljm);
      }
    }
    else
    {
      g_assert_cmpint (self->ell_max - self->ell_min + 1, ==, fftlog_w->len);

      if (ell_min < self->ell_min)
      {
        const guint nell = self->ell_min - ell_min;

        for (i = 0; i < nell; i++)
        {
          const gint ell             = self->ell_min - i - 1;
          NcmFftlogSBesselJLJM *jljm = ncm_fftlog_sbessel_jljm_new (ell, 0, log (w), self->lnr0, self->lnk0, self->Lk, 100);

          ncm_fftlog_use_smooth_padding (NCM_FFTLOG (jljm), TRUE);
          ncm_fftlog_set_smooth_padding_scale (NCM_FFTLOG (jljm), -180.0);

          /*ncm_fftlog_set_padding (NCM_FFTLOG (jljm), 0.0); */

          ncm_fftlog_sbessel_jljm_set_best_lnr0 (jljm);
          g_ptr_array_insert (fftlog_w, 0, jljm);
        }
      }
      else
      {
        const gint nell = ell_min - self->ell_min;

        g_ptr_array_remove_range (fftlog_w, 0, nell);
      }

      self->ell_min = ell_min;

      if (ell_max > self->ell_max)
      {
        const guint nell = ell_max - self->ell_max;

        for (i = 0; i < nell; i++)
        {
          const gint ell             = self->ell_max + i + 1;
          NcmFftlogSBesselJLJM *jljm = ncm_fftlog_sbessel_jljm_new (ell, 0, log (w), self->lnr0, self->lnk0, self->Lk, 100);

          ncm_fftlog_use_smooth_padding (NCM_FFTLOG (jljm), TRUE);
          ncm_fftlog_set_smooth_padding_scale (NCM_FFTLOG (jljm), -180.0);

          /*ncm_fftlog_set_padding (NCM_FFTLOG (jljm), 0.0); */

          ncm_fftlog_sbessel_jljm_set_best_lnr0 (jljm);
          g_ptr_array_add (fftlog_w, jljm);
        }
      }
      else
      {
        const gint nell = self->ell_max - ell_max;

        g_ptr_array_remove_range (fftlog_w, ell_max + 1 - self->ell_min, nell);
      }
    }
  }
}

static void
_ncm_powspec_sphere_proj_constructed (GObject *object)
{
  /* Chain up : start */
  G_OBJECT_CLASS (ncm_powspec_sphere_proj_parent_class)->constructed (object);
  {
    NcmPowspecSphereProj *psp                = NCM_POWSPEC_SPHERE_PROJ (object);
    NcmPowspecSphereProjPrivate * const self = ncm_powspec_sphere_proj_get_instance_private (psp);
    const gdouble lnk_min                    = log (ncm_powspec_get_kmin (self->ps));
    const gdouble lnk_max                    = log (ncm_powspec_get_kmax (self->ps));

    self->lnk0 = 0.5 * (lnk_max + lnk_min);
    self->Lk   = (lnk_max - lnk_min);

    _ncm_powspec_sphere_proj_adjust_fftlog_array (psp, self->ell_min, self->ell_max, TRUE);

    self->constructed = TRUE;
  }
}

static void
_ncm_powspec_sphere_proj_dispose (GObject *object)
{
  NcmPowspecSphereProj *psp                = NCM_POWSPEC_SPHERE_PROJ (object);
  NcmPowspecSphereProjPrivate * const self = ncm_powspec_sphere_proj_get_instance_private (psp);

  ncm_powspec_clear (&self->ps);

  ncm_model_ctrl_clear (&self->ctrl);

  g_clear_pointer (&self->w_array,     g_array_unref);
  g_clear_pointer (&self->lnk_array,   g_ptr_array_unref);
  g_clear_pointer (&self->k2Pk_array,  g_ptr_array_unref);
  g_clear_pointer (&self->fftlog,      g_ptr_array_unref);
  g_clear_pointer (&self->Cell_spline, g_ptr_array_unref);

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_powspec_sphere_proj_parent_class)->dispose (object);
}

static void
_ncm_powspec_sphere_proj_finalize (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (ncm_powspec_sphere_proj_parent_class)->finalize (object);
}

static void
ncm_powspec_sphere_proj_class_init (NcmPowspecSphereProjClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);

  object_class->set_property = &_ncm_powspec_sphere_proj_set_property;
  object_class->get_property = &_ncm_powspec_sphere_proj_get_property;
  object_class->constructed  = &_ncm_powspec_sphere_proj_constructed;
  object_class->dispose      = &_ncm_powspec_sphere_proj_dispose;
  object_class->finalize     = &_ncm_powspec_sphere_proj_finalize;

  g_object_class_install_property (object_class,
                                   PROP_XI_I,
                                   g_param_spec_double ("xi-i",
                                                        NULL,
                                                        "Output initial scale",
                                                        0.0, G_MAXDOUBLE, 1.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_XI_F,
                                   g_param_spec_double ("xi-f",
                                                        NULL,
                                                        "Output final scale",
                                                        0.0, G_MAXDOUBLE, 1.0e4,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_K_PIVOT,
                                   g_param_spec_double ("k-pivot",
                                                        NULL,
                                                        "Pivot k for growth computation",
                                                        0.0, G_MAXDOUBLE, 1.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_RELTOL,
                                   g_param_spec_double ("reltol",
                                                        NULL,
                                                        "Relative tolerance for calibration",
                                                        GSL_DBL_EPSILON, 1.0, 1.0e-5,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_RELTOL_Z,
                                   g_param_spec_double ("reltol-z",
                                                        NULL,
                                                        "Relative tolerance for calibration in the redshift direction",
                                                        GSL_DBL_EPSILON, 1.0, 1.0e-6,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_POWERSPECTRUM,
                                   g_param_spec_object ("powerspectrum",
                                                        NULL,
                                                        "NcmPowspec object",
                                                        NCM_TYPE_POWSPEC,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_ELL_MIN,
                                   g_param_spec_uint ("ell-min",
                                                      NULL,
                                                      "Minimum ell",
                                                      0, G_MAXUINT, 0,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_ELL_MAX,
                                   g_param_spec_uint ("ell-max",
                                                      NULL,
                                                      "Maximum ell",
                                                      1, G_MAXUINT, 1,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
}

/**
 * ncm_powspec_sphere_proj_new:
 * @ps: a #NcmPowspec
 * @ell_min: minimum $\ell$
 * @ell_max: maximum $\ell$
 *
 * Creates a new #NcmPowspecSphereProj from the power spectrum @ps.
 *
 * Returns: (transfer full): the newly created #NcmPowspecSphereProj.
 */
NcmPowspecSphereProj *
ncm_powspec_sphere_proj_new (NcmPowspec *ps, guint ell_min, guint ell_max)
{
  NcmPowspecSphereProj *psp = g_object_new (NCM_TYPE_POWSPEC_SPHERE_PROJ,
                                            "powerspectrum", ps,
                                            "ell-min",       ell_min,
                                            "ell-max",       ell_max,
                                            NULL);

  return psp;
}

/**
 * ncm_powspec_sphere_proj_ref:
 * @psp: a #NcmPowspecSphereProj
 *
 * Increases the reference count of @psp by one.
 *
 * Returns: (transfer full): @psp
 */
NcmPowspecSphereProj *
ncm_powspec_sphere_proj_ref (NcmPowspecSphereProj *psp)
{
  return g_object_ref (psp);
}

/**
 * ncm_powspec_sphere_proj_free:
 * @psp: a #NcmPowspecSphereProj
 *
 * Decreases the reference count of @psp by one.
 *
 */
void
ncm_powspec_sphere_proj_free (NcmPowspecSphereProj *psp)
{
  g_object_unref (psp);
}

/**
 * ncm_powspec_sphere_proj_clear:
 * @psp: a #NcmPowspecSphereProj
 *
 * If @psp is different from NULL, decreases the reference count of
 * @psp by one and sets @fftlog to NULL.
 *
 */
void
ncm_powspec_sphere_proj_clear (NcmPowspecSphereProj **psp)
{
  g_clear_object (psp);
}

typedef struct _NcmPowspecSphereProjArg
{
  NcmPowspecSphereProj *psp;
  NcmModel *model;
  gdouble z;
} NcmPowspecSphereProjArg;

static gdouble
_ncm_powspec_sphere_proj_k2Pk (gdouble k, gpointer userdata)
{
  NcmPowspecSphereProjArg *arg             = (NcmPowspecSphereProjArg *) userdata;
  NcmPowspecSphereProjPrivate * const self = ncm_powspec_sphere_proj_get_instance_private (arg->psp);
  const gdouble two_pi                     = 2.0 / ncm_c_pi ();
  const gdouble k2                         = k * k;
  const gdouble Pk                         = ncm_powspec_eval (self->ps, arg->model, arg->z, k);
  const gdouble f                          = Pk * k2 * two_pi;

  return f;
}

static gdouble
_ncm_powspec_sphere_proj_growth2 (NcmPowspecSphereProj *psp, NcmModel *model, const gdouble z)
{
  NcmPowspecSphereProjPrivate * const self = ncm_powspec_sphere_proj_get_instance_private (psp);

  return ncm_powspec_eval (self->ps, model, z, self->k_pivot) / self->Pk_pivot;
}

/**
 * ncm_powspec_sphere_proj_prepare:
 * @psp: a #NcmPowspecSphereProj
 * @model: a #NcmModel
 *
 * Prepares the object applying the filter to the power spectrum.
 *
 */
void
ncm_powspec_sphere_proj_prepare (NcmPowspecSphereProj *psp, NcmModel *model)
{
  NcmPowspecSphereProjPrivate * const self = ncm_powspec_sphere_proj_get_instance_private (psp);
  NcmPowspecSphereProjArg arg;
  gsl_function F;
  guint w_i;

  F.function = &_ncm_powspec_sphere_proj_k2Pk;
  F.params   = &arg;

  arg.psp   = psp;
  arg.model = model;
  arg.z     = 0.0;

  ncm_powspec_prepare_if_needed (self->ps, model);
  self->Pk_pivot = ncm_powspec_eval (self->ps, model, 0.0, self->k_pivot);

  {
    GPtrArray *fftlog_w0      = g_ptr_array_index (self->fftlog, 0);
    NcmFftlog *fftlog_w0_last = g_ptr_array_index (fftlog_w0, fftlog_w0->len - 1);
    const gdouble lnk_min     = log (ncm_powspec_get_kmin (self->ps));
    const gdouble lnk_max     = log (ncm_powspec_get_kmax (self->ps));

    self->lnk0 = 0.5 * (lnk_max + lnk_min);
    self->Lk   = (lnk_max - lnk_min);

    if ((self->lnk0 != ncm_fftlog_get_lnk0 (fftlog_w0_last)) || (self->Lk != ncm_fftlog_get_length (fftlog_w0_last)))
    {
      NcmFftlogSBesselJLJM *jljm0 = NCM_FFTLOG_SBESSEL_JLJM (fftlog_w0_last);

      ncm_fftlog_set_lnk0 (fftlog_w0_last, self->lnk0);
      ncm_fftlog_set_length (fftlog_w0_last, self->Lk);
      ncm_fftlog_sbessel_jljm_set_best_lnr0 (jljm0);

      self->calibrated = FALSE;
    }
  }

  if (!self->calibrated)
  {
    for (w_i = 0; w_i < self->w_array->len; w_i++)
    {
      const gdouble w          = g_array_index (self->w_array, gdouble, w_i);
      GPtrArray *fftlog_w      = g_ptr_array_index (self->fftlog, w_i);
      NcmFftlog *fftlog_w_last = g_ptr_array_index (fftlog_w, fftlog_w->len - 1);

      ncm_fftlog_set_eval_r_min (fftlog_w_last, self->xi_i / w);
      ncm_fftlog_set_eval_r_max (fftlog_w_last, self->xi_f * w);
      ncm_fftlog_use_eval_interval (fftlog_w_last, TRUE);

      ncm_fftlog_calibrate_size_gsl (fftlog_w_last, &F, self->reltol);
      {
        NcmVector *lnk;
        NcmVector *k2Pk;
        guint n = ncm_fftlog_get_size (fftlog_w_last);
        guint i;

        if (w_i < self->lnk_array->len)
        {
          ncm_vector_free (g_ptr_array_index (self->lnk_array, w_i));
          ncm_vector_free (g_ptr_array_index (self->k2Pk_array, w_i));

          g_ptr_array_index (self->lnk_array, w_i)  = lnk  = ncm_vector_new (n);
          g_ptr_array_index (self->k2Pk_array, w_i) = k2Pk = ncm_vector_new (n);
        }
        else
        {
          g_assert_cmpint (w_i, ==, self->lnk_array->len);

          lnk  = ncm_vector_new (n);
          k2Pk = ncm_vector_new (n);

          g_ptr_array_add (self->lnk_array, lnk);
          g_ptr_array_add (self->k2Pk_array, k2Pk);
        }

        ncm_fftlog_get_lnk_vector (fftlog_w_last, lnk);

        /*ncm_vector_log_vals (lnk, "lnk: ", "% 22.15g", TRUE);*/

        for (i = 0; i < fftlog_w->len - 1; i++)
        {
          NcmFftlog *fftlog_w_i          = g_ptr_array_index (fftlog_w, i);
          NcmFftlogSBesselJLJM *jljm_w_i = NCM_FFTLOG_SBESSEL_JLJM (fftlog_w_i);

          ncm_fftlog_set_size (fftlog_w_i, n);
          ncm_fftlog_sbessel_jljm_set_best_lnr0 (jljm_w_i);
          ncm_fftlog_set_eval_r_min (fftlog_w_i, self->xi_i / w);
          ncm_fftlog_set_eval_r_max (fftlog_w_i, self->xi_f * w);
          ncm_fftlog_use_eval_interval (fftlog_w_i, TRUE);

          g_assert_cmpint (ncm_fftlog_get_size (fftlog_w_i), ==, n);
        }
      }
    }

    self->calibrated = TRUE;
  }

  for (w_i = 0; w_i < self->w_array->len; w_i++)
  {
    GPtrArray *fftlog_w = g_ptr_array_index (self->fftlog, w_i);
    NcmVector *lnk      = g_ptr_array_index (self->lnk_array, w_i);
    NcmVector *k2Pk     = g_ptr_array_index (self->k2Pk_array, w_i);
    const guint n       = ncm_vector_len (lnk);
    guint i;

    for (i = 0; i < n; i++)
    {
      const gdouble lnk_i = ncm_vector_get (lnk, i);
      const gdouble k_i   = exp (lnk_i);

      ncm_vector_set (k2Pk, i, GSL_FN_EVAL (&F, k_i));
    }

    for (i = 0; i < fftlog_w->len; i++)
    {
      NcmFftlog *fftlog_w_i = g_ptr_array_index (fftlog_w, i);

      ncm_fftlog_eval_by_vector (fftlog_w_i, k2Pk);
      ncm_fftlog_prepare_splines (fftlog_w_i);
    }
  }

  {
    const guint nn = (self->ell_max - self->ell_min + 1);
    NcmVector *w_v = ncm_vector_new (self->w_array->len);
    guint i;

    for (i = 0; i < self->w_array->len; i++)
      ncm_vector_set (w_v, i, g_array_index (self->w_array, gdouble, self->w_array->len - 1 - i));

    g_ptr_array_set_size (self->Cell_spline, 0);

    for (i = 0; i < nn; i++)
    {
      GPtrArray *fftlog_w_last   = g_ptr_array_index (self->fftlog, self->w_array->len - 1);
      NcmFftlog *fftlog_w_last_i = g_ptr_array_index (fftlog_w_last, i);

      NcmSpline *Cell_s   = ncm_fftlog_peek_spline_Gr (fftlog_w_last_i, 0);
      NcmVector *lnxi_v   = ncm_spline_get_xv (Cell_s);
      const gint xi_v_len = ncm_vector_len (lnxi_v);

      NcmMatrix *Cell_m = ncm_matrix_new (self->w_array->len, xi_v_len);

      for (w_i = 0; w_i < self->w_array->len; w_i++)
      {
        GPtrArray *fftlog_w   = g_ptr_array_index (self->fftlog, w_i);
        NcmFftlog *fftlog_w_i = g_ptr_array_index (fftlog_w, i);
        NcmSpline *Cell_w_s   = ncm_fftlog_peek_spline_Gr (fftlog_w_i, 0);
        gint j;

        for (j = 0; j < xi_v_len; j++)
        {
          const gdouble lnxi_j = ncm_vector_get (lnxi_v, j);

          ncm_matrix_set (Cell_m, self->w_array->len - 1 - w_i, j, ncm_spline_eval (Cell_w_s, lnxi_j));
        }
      }

      {
        NcmSpline2d *s2d = ncm_spline2d_bicubic_notaknot_new ();

        ncm_spline2d_set (s2d, lnxi_v, w_v, Cell_m, TRUE);
        g_ptr_array_add (self->Cell_spline, s2d);
      }

      ncm_vector_free (lnxi_v);
    }

    ncm_vector_free (w_v);
  }
}

/**
 * ncm_powspec_sphere_proj_prepare_if_needed:
 * @psp: a #NcmPowspecSphereProj
 * @model: a #NcmModel
 *
 * Prepares (if necessary) the object applying the filter to the power spectrum.
 *
 */
void
ncm_powspec_sphere_proj_prepare_if_needed (NcmPowspecSphereProj *psp, NcmModel *model)
{
  NcmPowspecSphereProjPrivate * const self = ncm_powspec_sphere_proj_get_instance_private (psp);
  gboolean model_up                        = ncm_model_ctrl_update (self->ctrl, model);

  if (model_up)
    ncm_powspec_sphere_proj_prepare (psp, model);
}

/**
 * ncm_powspec_sphere_proj_set_xi_i:
 * @psp: a #NcmPowspecSphereProj
 * @xi_i: the output initial scale $\xi_i$
 *
 * Sets the initial scale $\xi$ to @xi_i.
 *
 */
void
ncm_powspec_sphere_proj_set_xi_i (NcmPowspecSphereProj *psp, gdouble xi_i)
{
  NcmPowspecSphereProjPrivate * const self = ncm_powspec_sphere_proj_get_instance_private (psp);

  if (self->xi_i != xi_i)
  {
    self->xi_i = xi_i;
    ncm_model_ctrl_force_update (self->ctrl);
    self->calibrated = FALSE;
  }
}

/**
 * ncm_powspec_sphere_proj_set_xi_f:
 * @psp: a #NcmPowspecSphereProj
 * @xi_f: the output final scale $\xi_f$
 *
 * Sets the final scale $\xi$ to @xi_f.
 *
 */
void
ncm_powspec_sphere_proj_set_xi_f (NcmPowspecSphereProj *psp, gdouble xi_f)
{
  NcmPowspecSphereProjPrivate * const self = ncm_powspec_sphere_proj_get_instance_private (psp);

  if (self->xi_f != xi_f)
  {
    self->xi_f = xi_f;
    ncm_model_ctrl_force_update (self->ctrl);
    self->calibrated = FALSE;
  }
}

/**
 * ncm_powspec_sphere_proj_set_k_pivot:
 * @psp: a #NcmPowspecSphereProj
 * @k_pivot: the pivot k
 *
 * Sets the pivot k to @k_pivot.
 *
 */
void
ncm_powspec_sphere_proj_set_k_pivot (NcmPowspecSphereProj *psp, gdouble k_pivot)
{
  NcmPowspecSphereProjPrivate * const self = ncm_powspec_sphere_proj_get_instance_private (psp);

  if (self->k_pivot != k_pivot)
  {
    self->k_pivot = k_pivot;
    ncm_model_ctrl_force_update (self->ctrl);
    self->calibrated = FALSE;
  }
}

/**
 * ncm_powspec_sphere_proj_get_r_min:
 * @psp: a #NcmPowspecSphereProj
 *
 * Gets the minimum $r$.
 *
 * Returns: the minimum $r$.
 */
gdouble
ncm_powspec_sphere_proj_get_r_min (NcmPowspecSphereProj *psp)
{
  NcmPowspecSphereProjPrivate * const self = ncm_powspec_sphere_proj_get_instance_private (psp);

  return exp (self->lnr0 - self->Lk * 0.5);
}

/**
 * ncm_powspec_sphere_proj_get_r_max:
 * @psp: a #NcmPowspecSphereProj
 *
 * Gets the maximum $r$.
 *
 * Returns: the maximum $r$.
 */
gdouble
ncm_powspec_sphere_proj_get_r_max (NcmPowspecSphereProj *psp)
{
  NcmPowspecSphereProjPrivate * const self = ncm_powspec_sphere_proj_get_instance_private (psp);

  return exp (self->lnr0 + self->Lk * 0.5);
}

/**
 * ncm_powspec_sphere_proj_set_ell_min:
 * @psp: a #NcmPowspecSphereProj
 * @ell_min: minimum $\ell$
 *
 * Sets the minimum $\ell$ to @ell_min.
 */
void
ncm_powspec_sphere_proj_set_ell_min (NcmPowspecSphereProj *psp, const guint ell_min)
{
  NcmPowspecSphereProjPrivate * const self = ncm_powspec_sphere_proj_get_instance_private (psp);

  if (self->constructed && (ell_min != self->ell_min))
  {
    _ncm_powspec_sphere_proj_adjust_fftlog_array (psp, ell_min, self->ell_max, FALSE);
    ncm_model_ctrl_force_update (self->ctrl);
    self->calibrated = FALSE;
  }
  else
  {
    self->ell_min = ell_min;
  }
}

/**
 * ncm_powspec_sphere_proj_set_ell_max:
 * @psp: a #NcmPowspecSphereProj
 * @ell_max: maximum $\ell$
 *
 * Sets the maximum $\ell$ to @ell_max.
 */
void
ncm_powspec_sphere_proj_set_ell_max (NcmPowspecSphereProj *psp, const guint ell_max)
{
  NcmPowspecSphereProjPrivate * const self = ncm_powspec_sphere_proj_get_instance_private (psp);

  if (self->constructed && (ell_max != self->ell_max))
  {
    _ncm_powspec_sphere_proj_adjust_fftlog_array (psp, self->ell_min, ell_max, FALSE);
    ncm_model_ctrl_force_update (self->ctrl);
    self->calibrated = FALSE;
  }
  else
  {
    self->ell_max = ell_max;
  }
}

/**
 * ncm_powspec_sphere_proj_get_ell_min:
 * @psp: a #NcmPowspecSphereProj
 *
 * Returns: the minimum $\ell$.
 */
guint
ncm_powspec_sphere_proj_get_ell_min (NcmPowspecSphereProj *psp)
{
  NcmPowspecSphereProjPrivate * const self = ncm_powspec_sphere_proj_get_instance_private (psp);

  return self->ell_min;
}

/**
 * ncm_powspec_sphere_proj_get_ell_max:
 * @psp: a #NcmPowspecSphereProj
 *
 * Returns: the maximum $\ell$.
 */
guint
ncm_powspec_sphere_proj_get_ell_max (NcmPowspecSphereProj *psp)
{
  NcmPowspecSphereProjPrivate * const self = ncm_powspec_sphere_proj_get_instance_private (psp);

  return self->ell_max;
}

/**
 * ncm_powspec_sphere_proj_get_w:
 * @psp: a #NcmPowspecSphereProj
 * @w_i: w index
 *
 * Returns: the value of the @w_i-th $w$.
 */
gdouble
ncm_powspec_sphere_proj_get_w (NcmPowspecSphereProj *psp, const guint w_i)
{
  NcmPowspecSphereProjPrivate * const self = ncm_powspec_sphere_proj_get_instance_private (psp);

  g_assert_cmpuint (w_i, <, self->w_array->len);

  return g_array_index (self->w_array, gdouble, w_i);
}

/**
 * ncm_powspec_sphere_proj_get_ell:
 * @psp: a #NcmPowspecSphereProj
 * @w_i: w index
 * @ell: the value of $\ell$
 * @lnxi: (out) (transfer full) (array) (element-type gdouble): array with the $\xi$ values
 * @Cell: (out) (transfer full) (array) (element-type gdouble): array with the $\xi$ values
 *
 * Returns the results for $\ell$.
 *
 */
void
ncm_powspec_sphere_proj_get_ell (NcmPowspecSphereProj *psp, const guint w_i, const gint ell, GArray **lnxi, GArray **Cell)
{
  NcmPowspecSphereProjPrivate * const self = ncm_powspec_sphere_proj_get_instance_private (psp);

  g_assert_cmpuint (w_i, <, self->w_array->len);
  g_assert_cmpint (ell, <=, self->ell_max);
  g_assert_cmpint (ell, >=, self->ell_min);

  {
    GPtrArray *fftlog_w = g_ptr_array_index (self->fftlog, w_i);
    NcmFftlog *fftlog   = g_ptr_array_index (fftlog_w, ell - self->ell_min);
    NcmSpline *Cell_s   = ncm_fftlog_peek_spline_Gr (fftlog, 0);
    NcmVector *lnxi_v   = ncm_spline_get_xv (Cell_s);
    NcmVector *Cell_v   = ncm_spline_get_yv (Cell_s);

    *lnxi = ncm_vector_dup_array (lnxi_v);
    *Cell = ncm_vector_dup_array (Cell_v);

    ncm_vector_free (lnxi_v);
    ncm_vector_free (Cell_v);
  }
}

/**
 * ncm_powspec_sphere_proj_peek_ell_spline:
 * @psp: a #NcmPowspecSphereProj
 * @w_i: w index
 * @ell: the value of $\ell$
 *
 * Returns: (transfer none): the resulting spline for $\ell$.
 */
NcmSpline *
ncm_powspec_sphere_proj_peek_ell_spline (NcmPowspecSphereProj *psp, const guint w_i, const gint ell)
{
  NcmPowspecSphereProjPrivate * const self = ncm_powspec_sphere_proj_get_instance_private (psp);
  GPtrArray *fftlog_w                      = g_ptr_array_index (self->fftlog, w_i);

  g_assert_cmpuint (w_i, <, self->w_array->len);
  g_assert_cmpint (ell, <=, self->ell_max);
  g_assert_cmpint (ell, >=, self->ell_min);

  {
    NcmFftlog *fftlog = g_ptr_array_index (fftlog_w, ell - self->ell_min);

    return ncm_fftlog_peek_spline_Gr (fftlog, 0);
  }
}

/**
 * ncm_powspec_sphere_proj_eval_Cell_xi1_xi2:
 * @psp: a #NcmPowspecSphereProj
 * @model: a #NcmModel
 * @ell: the value of $\ell$
 * @z1: the value of $z_1$
 * @z2: the value of $z_2$
 * @xi1: the value of $\xi_1$
 * @xi2: the value of $\xi_2$
 *
 * Returns: the value of $C_\ell(\xi_1, \xi_2)$.
 */
gdouble
ncm_powspec_sphere_proj_eval_Cell_xi1_xi2 (NcmPowspecSphereProj *psp, NcmModel *model, const gint ell, const gdouble z1, const gdouble z2, const gdouble xi1, const gdouble xi2)
{
  NcmPowspecSphereProjPrivate * const self = ncm_powspec_sphere_proj_get_instance_private (psp);
  NcmSpline2d *s2d                         = g_ptr_array_index (self->Cell_spline, ell - self->ell_min);
  const gdouble lnxi                       = 0.5 * log (xi1 * xi2);
  const gdouble w                          = (xi1 > xi2) ? sqrt (xi2 / xi1) : sqrt (xi1 / xi2);

  /*printf ("% 22.15g % 22.15g % 22.15g % 22.15g | % 22.15g\n", xi1, xi2, exp (lnxi), w, g_array_index (self->w_array, gdouble, self->w_array->len - 1)); */
  if (w < g_array_index (self->w_array, gdouble, self->w_array->len - 1))
    return 0.0;
  else
    return ncm_spline2d_eval (s2d, lnxi, w) * sqrt (_ncm_powspec_sphere_proj_growth2 (psp, model, z1) * _ncm_powspec_sphere_proj_growth2 (psp, model, z2));
}

