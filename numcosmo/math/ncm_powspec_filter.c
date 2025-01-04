/***************************************************************************
 *            ncm_powspec_filter.c
 *
 *  Fri June 17 10:12:06 2016
 *  Copyright  2016  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/* excerpt from: */

/***************************************************************************
 *            nc_window_gaussian.c
 *
 *  Mon Jun 28 15:09:13 2010
 *  Copyright  2010  Mariana Penna Lima
 *  <pennalima@gmail.com>
 ****************************************************************************/
/***************************************************************************
 *            nc_window_tophat.c
 *
 *  Mon Jun 28 15:09:13 2010
 *  Copyright  2010  Mariana Penna Lima
 *  <pennalima@gmail.com>
 ****************************************************************************/
/*
 * ncm_powspec_filter.c
 * Copyright (C) 2016 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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
 * NcmPowspecFilter:
 *
 * Class to compute filtered power spectrum.
 *
 * This class computes the filtered power spectrum, $\sigma^2(k, r)$, and its derivatives with respect to $\ln r$
 * (#ncm_powspec_filter_eval_dnvar_dlnrn()) using the FFTLog approach (see #NcmFftlog),
 * \begin{equation}\label{eq:variance}
 * \sigma^2(r, z) = \frac{1}{2\pi^2} \int_0^\infty k^2 \ P(k, z) \vert W(k,r) \vert^2 \ \mathrm{d}k,
 * \end{equation}
 * where $P(k, z)$ is the power spectrum at mode $k$ and redshift $z$ and $W(k, r)$ is the filter (or window function).
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_powspec_filter.h"
#include "math/ncm_spline_cubic_notaknot.h"
#include "math/ncm_spline2d_bicubic.h"
#include "math/ncm_fftlog_tophatwin2.h"
#include "math/ncm_fftlog_gausswin2.h"
#include "math/ncm_c.h"
#include "ncm_enum_types.h"

enum
{
  PROP_0,
  PROP_TYPE,
  PROP_LNR0,
  PROP_ZI,
  PROP_ZF,
  PROP_RELTOL,
  PROP_RELTOL_Z,
  PROP_MAX_K_KNOTS,
  PROP_MAX_Z_KNOTS,
  PROP_POWERSPECTRUM,
  PROP_SIZE,
};

struct _NcmPowspecFilter
{
  /*< private >*/
  GObject parent_instance;
  NcmPowspec *ps;
  NcmPowspecFilterType type;
  NcmFftlog *fftlog;
  gdouble lnr0;
  gdouble lnk0;
  gdouble Lk;
  gdouble zi;
  gdouble zf;
  gboolean calibrated;
  gdouble reltol;
  gdouble reltol_z;
  guint max_k_knots;
  guint max_z_knots;
  NcmSpline2d *var;
  NcmSpline2d *dvar;
  NcmModelCtrl *ctrl;
  gboolean constructed;
};

G_DEFINE_TYPE (NcmPowspecFilter, ncm_powspec_filter, G_TYPE_OBJECT)

static void
ncm_powspec_filter_init (NcmPowspecFilter *psf)
{
  psf->ps          = NULL;
  psf->lnr0        = 0.0;
  psf->lnk0        = 0.0;
  psf->Lk          = 0.0;
  psf->zi          = 0.0;
  psf->zf          = 0.0;
  psf->reltol      = 0.0;
  psf->reltol_z    = 0.0;
  psf->max_k_knots = 0;
  psf->max_z_knots = 0;
  psf->type        = NCM_POWSPEC_FILTER_TYPE_LEN;
  psf->fftlog      = NULL;
  psf->calibrated  = FALSE;
  psf->var         = ncm_spline2d_bicubic_notaknot_new ();
  psf->dvar        = ncm_spline2d_bicubic_notaknot_new ();
  psf->ctrl        = ncm_model_ctrl_new (NULL);
  psf->constructed = FALSE;
}

static void
_ncm_powspec_filter_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmPowspecFilter *psf = NCM_POWSPEC_FILTER (object);

  g_return_if_fail (NCM_IS_POWSPEC_FILTER (object));

  switch (prop_id)
  {
    case PROP_TYPE:
      ncm_powspec_filter_set_type (psf, g_value_get_enum (value));
      break;
    case PROP_LNR0:
      ncm_powspec_filter_set_lnr0 (psf, g_value_get_double (value));
      break;
    case PROP_ZI:
      ncm_powspec_filter_set_zi (psf, g_value_get_double (value));
      break;
    case PROP_ZF:
      ncm_powspec_filter_set_zf (psf, g_value_get_double (value));
      break;
    case PROP_RELTOL:
      psf->reltol = g_value_get_double (value);
      break;
    case PROP_RELTOL_Z:
      psf->reltol_z = g_value_get_double (value);
      break;
    case PROP_MAX_K_KNOTS:
      psf->max_k_knots = g_value_get_uint (value);
      break;
    case PROP_MAX_Z_KNOTS:
      psf->max_z_knots = g_value_get_uint (value);
      break;
    case PROP_POWERSPECTRUM:
      psf->ps = g_value_dup_object (value);
      psf->zi = ncm_powspec_get_zi (psf->ps);
      psf->zf = ncm_powspec_get_zf (psf->ps);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_ncm_powspec_filter_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmPowspecFilter *psf = NCM_POWSPEC_FILTER (object);

  g_return_if_fail (NCM_IS_POWSPEC_FILTER (object));

  switch (prop_id)
  {
    case PROP_TYPE:
      g_value_set_enum (value, psf->type);
      break;
    case PROP_LNR0:
      g_value_set_double (value, psf->lnr0);
      break;
    case PROP_ZI:
      g_value_set_double (value, psf->zi);
      break;
    case PROP_ZF:
      g_value_set_double (value, psf->zf);
      break;
    case PROP_RELTOL:
      g_value_set_double (value, psf->reltol);
      break;
    case PROP_RELTOL_Z:
      g_value_set_double (value, psf->reltol_z);
      break;
    case PROP_MAX_K_KNOTS:
      g_value_set_uint (value, psf->max_k_knots);
      break;
    case PROP_MAX_Z_KNOTS:
      g_value_set_uint (value, psf->max_z_knots);
      break;
    case PROP_POWERSPECTRUM:
      g_value_set_object (value, psf->ps);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_ncm_powspec_filter_constructed (GObject *object)
{
  /* Chain up : start */
  G_OBJECT_CLASS (ncm_powspec_filter_parent_class)->constructed (object);
  {
    NcmPowspecFilter *psf     = NCM_POWSPEC_FILTER (object);
    NcmPowspecFilterType type = psf->type;

    psf->constructed = TRUE;
    psf->type        = NCM_POWSPEC_FILTER_TYPE_LEN;

    ncm_powspec_filter_set_type (psf, type);
  }
}

static void
_ncm_powspec_filter_dispose (GObject *object)
{
  NcmPowspecFilter *psf = NCM_POWSPEC_FILTER (object);

  ncm_powspec_clear (&psf->ps);
  ncm_fftlog_clear (&psf->fftlog);

  ncm_spline2d_clear (&psf->var);
  ncm_spline2d_clear (&psf->dvar);

  ncm_model_ctrl_clear (&psf->ctrl);

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_powspec_filter_parent_class)->dispose (object);
}

static void
_ncm_powspec_filter_finalize (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (ncm_powspec_filter_parent_class)->finalize (object);
}

static void
ncm_powspec_filter_class_init (NcmPowspecFilterClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);

  object_class->set_property = &_ncm_powspec_filter_set_property;
  object_class->get_property = &_ncm_powspec_filter_get_property;
  object_class->constructed  = &_ncm_powspec_filter_constructed;
  object_class->dispose      = &_ncm_powspec_filter_dispose;
  object_class->finalize     = &_ncm_powspec_filter_finalize;

  /**
   * NcmPowspecFilter:lnr0:
   *
   * The output center value for $\ln(r)$.
   */
  g_object_class_install_property (object_class,
                                   PROP_LNR0,
                                   g_param_spec_double ("lnr0",
                                                        NULL,
                                                        "Output center value",
                                                        -G_MAXDOUBLE, G_MAXDOUBLE, 0.0,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmPowspecFilter:zi:
   *
   * The output initial time $z_i$ of the variance $\sigma^{2}(r,z)$.
   */
  g_object_class_install_property (object_class,
                                   PROP_ZI,
                                   g_param_spec_double ("zi",
                                                        NULL,
                                                        "Output initial time",
                                                        -G_MAXDOUBLE, G_MAXDOUBLE, 0.0,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmPowspecFilter:zf:
   *
   * The output final time $z_f$ of the variance $\sigma^{2}(r,z)$.
   */
  g_object_class_install_property (object_class,
                                   PROP_ZF,
                                   g_param_spec_double ("zf",
                                                        NULL,
                                                        "Output final time",
                                                        -G_MAXDOUBLE, G_MAXDOUBLE, 1.0,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmPowspecFilter:reltol:
   *
   * The relative tolerance for calibration in the distance direction.
   */
  g_object_class_install_property (object_class,
                                   PROP_RELTOL,
                                   g_param_spec_double ("reltol",
                                                        NULL,
                                                        "Relative tolerance for calibration",
                                                        GSL_DBL_EPSILON, 1.0, 1.0e-3,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmPowspecFilter:reltol-z:
   *
   * The relative tolerance for calibration in the redshift direction.
   */
  g_object_class_install_property (object_class,
                                   PROP_RELTOL_Z,
                                   g_param_spec_double ("reltol-z",
                                                        NULL,
                                                        "Relative tolerance for calibration in the redshift direction",
                                                        GSL_DBL_EPSILON, 1.0, 1.0e-6,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmPowspecFilter:max-k-knots:
   *
   * The maximum number of knots in the k direction.
   */
  g_object_class_install_property (object_class,
                                   PROP_MAX_K_KNOTS,
                                   g_param_spec_uint ("max-k-knots",
                                                      NULL,
                                                      "Maximum number of knots in the k direction",
                                                      0, G_MAXUINT, 10000,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmPowspecFilter:max-z-knots:
   *
   * The maximum number of knots in the redshift direction.
   */
  g_object_class_install_property (object_class,
                                   PROP_MAX_Z_KNOTS,
                                   g_param_spec_uint ("max-z-knots",
                                                      NULL,
                                                      "Maximum number of knots in the redshift direction",
                                                      0, G_MAXUINT, 1000,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmPowspecFilter:type:
   *
   * The type of fliter used $W(k,r)$.
   */
  g_object_class_install_property (object_class,
                                   PROP_TYPE,
                                   g_param_spec_enum ("type",
                                                      NULL,
                                                      "Filter type",
                                                      NCM_TYPE_POWSPEC_FILTER_TYPE, NCM_POWSPEC_FILTER_TYPE_TOPHAT,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmPowspecFilter:powerspectrum:
   *
   * The #NcmPowspec object to be used to compute the variance $\sigma^{2}(r,z)$.
   */
  g_object_class_install_property (object_class,
                                   PROP_POWERSPECTRUM,
                                   g_param_spec_object ("powerspectrum",
                                                        NULL,
                                                        "NcmPowspec object",
                                                        NCM_TYPE_POWSPEC,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
}

/**
 * ncm_powspec_filter_new:
 * @ps: a #NcmPowspec
 * @type: a type from #NcmPowspecFilterType
 *
 * Creates a new #NcmPowspecFilter from the power spectrum @ps.
 *
 * Returns: (transfer full): the newly created #NcmPowspecFilter.
 */
NcmPowspecFilter *
ncm_powspec_filter_new (NcmPowspec *ps, NcmPowspecFilterType type)
{
  NcmPowspecFilter *psf = g_object_new (NCM_TYPE_POWSPEC_FILTER,
                                        "type", type,
                                        "powerspectrum", ps,
                                        NULL);

  return psf;
}

/**
 * ncm_powspec_filter_ref:
 * @psf: a #NcmPowspecFilter
 *
 * Increases the reference count of @psf by one atomically.
 *
 * Returns: (transfer full): @psf
 */
NcmPowspecFilter *
ncm_powspec_filter_ref (NcmPowspecFilter *psf)
{
  return g_object_ref (psf);
}

/**
 * ncm_powspec_filter_free:
 * @psf: a #NcmPowspecFilter
 *
 * Atomically decrements the reference count of @psf by one.
 * If the reference count drops to 0, all memory allocated by @psf is released.
 *
 */
void
ncm_powspec_filter_free (NcmPowspecFilter *psf)
{
  g_object_unref (psf);
}

/**
 * ncm_powspec_filter_clear:
 * @psf: a #NcmPowspecFilter
 *
 * If @psf is different from NULL,
 * atomically decrements the reference count of @psf by one.
 * If the reference count drops to 0,
 * all memory allocated by @psf is released and @psf is set to NULL.
 *
 */
void
ncm_powspec_filter_clear (NcmPowspecFilter **psf)
{
  g_clear_object (psf);
}

/**
 * ncm_powspec_filter_set_type:
 * @psf: a #NcmPowspecFilter
 * @type: a type from #NcmPowspecFilterType
 *
 * Sets the @type of the #NcmPowspecFilter to be used.
 *
 */
void
ncm_powspec_filter_set_type (NcmPowspecFilter *psf, NcmPowspecFilterType type)
{
  if (!psf->constructed)
  {
    psf->type = type;
  }
  else if (type != psf->type)
  {
    const gdouble lnk_min = log (ncm_powspec_get_kmin (psf->ps));
    const gdouble lnk_max = log (ncm_powspec_get_kmax (psf->ps));

    psf->lnk0 = 0.5 * (lnk_max + lnk_min);
    psf->Lk   = (lnk_max - lnk_min);

    ncm_fftlog_clear (&psf->fftlog);
    psf->type = type;

    switch (psf->type)
    {
      case NCM_POWSPEC_FILTER_TYPE_TOPHAT:
        psf->fftlog = NCM_FFTLOG (ncm_fftlog_tophatwin2_new (psf->lnr0, psf->lnk0, psf->Lk, 100));
        break;
      case NCM_POWSPEC_FILTER_TYPE_GAUSS:
        psf->fftlog = NCM_FFTLOG (ncm_fftlog_gausswin2_new (psf->lnr0, psf->lnk0, psf->Lk, 100));
        break;
      default:
        g_assert_not_reached ();
        break;
    }

    ncm_fftlog_set_padding (psf->fftlog, 1.0);
    ncm_fftlog_set_nderivs (psf->fftlog, 1);

    ncm_powspec_filter_set_best_lnr0 (psf);

    ncm_model_ctrl_force_update (psf->ctrl);
    psf->calibrated = FALSE;
  }
}

typedef struct _NcmPowspecFilterArg
{
  NcmPowspecFilter *psf;
  NcmModel *model;
  gdouble z;
} NcmPowspecFilterArg;

static gdouble
_ncm_powspec_filter_k2Pk (gdouble k, gpointer userdata)
{
  NcmPowspecFilterArg *arg = (NcmPowspecFilterArg *) userdata;
  const gdouble k2         = k * k;
  const gdouble Pk         = ncm_powspec_eval (arg->psf->ps, arg->model, arg->z, k);
  const gdouble f          = Pk * k2 / ncm_c_two_pi_2 ();

  return f;
}

static gdouble
_ncm_powspec_filter_dummy_z (gdouble z, gpointer userdata)
{
  NcmPowspecFilterArg *arg = (NcmPowspecFilterArg *) userdata;
  gsl_function F;

  F.function = &_ncm_powspec_filter_k2Pk;
  F.params   = arg;

  arg->z = z;
  ncm_fftlog_eval_by_gsl_function (arg->psf->fftlog, &F);

  /*printf ("# z-knots % 20.15g % 20.15g\n", z, ncm_vector_get (ncm_fftlog_peek_output_vector (arg->psf->fftlog, 0), 0));*/
  return ncm_vector_get (ncm_fftlog_peek_output_vector (arg->psf->fftlog, 0), 0);
}

/**
 * ncm_powspec_filter_prepare:
 * @psf: a #NcmPowspecFilter
 * @model: a #NcmModel
 *
 * Prepares the object applying the filter to the power spectrum.
 *
 */
void
ncm_powspec_filter_prepare (NcmPowspecFilter *psf, NcmModel *model)
{
  NcmPowspecFilterArg arg;
  gsl_function F;

  F.function = &_ncm_powspec_filter_k2Pk;
  F.params   = &arg;

  arg.psf   = psf;
  arg.model = model;
  arg.z     = 0.0;

  ncm_powspec_prepare_if_needed (psf->ps, model);

  {
    const gdouble lnk_min = log (ncm_powspec_get_kmin (psf->ps));
    const gdouble lnk_max = log (ncm_powspec_get_kmax (psf->ps));

    psf->lnk0 = 0.5 * (lnk_max + lnk_min);
    psf->Lk   = (lnk_max - lnk_min);

    if ((psf->lnk0 != ncm_fftlog_get_lnk0 (psf->fftlog)) || (psf->Lk != ncm_fftlog_get_length (psf->fftlog)))
    {
      ncm_fftlog_set_lnk0 (psf->fftlog, psf->lnk0);
      ncm_fftlog_set_length (psf->fftlog, psf->Lk);
      psf->calibrated = FALSE;
    }
  }

  if (!psf->calibrated)
  {
    NcmMatrix *lnvar, *dlnvar;
    NcmVector *z_vec, *lnr_vec;
    guint N_k = 0, N_z = 0;
    guint i;

    ncm_powspec_get_nknots (psf->ps, &N_z, &N_k);

    ncm_fftlog_set_max_size (psf->fftlog, psf->max_k_knots);
    ncm_fftlog_calibrate_size_gsl (psf->fftlog, &F, psf->reltol);
    N_k = ncm_fftlog_get_size (psf->fftlog);

    {
      NcmSpline *dummy_z = NCM_SPLINE (ncm_spline_cubic_notaknot_new ());
      gsl_function Fdummy_z;

      Fdummy_z.function = &_ncm_powspec_filter_dummy_z;
      Fdummy_z.params   = &arg;

      ncm_spline_set_func (dummy_z, NCM_SPLINE_FUNCTION_SPLINE, &Fdummy_z, psf->zi, psf->zf, 0, psf->reltol_z);

      z_vec = ncm_spline_get_xv (dummy_z);
      N_z   = ncm_vector_len (z_vec);

      ncm_spline_clear (&dummy_z);
    }

    g_assert_cmpuint (N_z, >, 0);
    g_assert_cmpuint (N_k, >, 0);

/*
 *   printf ("# Calibrating in zmin % 20.15g zmax % 20.15g, rmin % 20.15g rmax % 20.15g, N_z = %u, N_k = %u\n",
 *           psf->zi, psf->zf,
 *           ncm_powspec_filter_get_r_min (psf),
 *           ncm_powspec_filter_get_r_max (psf),
 *           N_z, N_k);
 */
    lnvar   = ncm_matrix_new (N_z, N_k);
    dlnvar  = ncm_matrix_new (N_z, N_k);
    lnr_vec = ncm_fftlog_get_vector_lnr (psf->fftlog);

    for (i = 0; i < N_z; i++)
    {
      NcmVector *var_z  = ncm_matrix_get_row (lnvar, i);
      NcmVector *dvar_z = ncm_matrix_get_row (dlnvar, i);

      arg.z = ncm_vector_get (z_vec, i);
      ncm_fftlog_eval_by_gsl_function (psf->fftlog, &F);

      ncm_vector_memcpy (var_z, ncm_fftlog_peek_output_vector (psf->fftlog, 0));
      ncm_vector_memcpy (dvar_z, ncm_fftlog_peek_output_vector (psf->fftlog, 1));

      /*ncm_vector_log_vals (var_z, "NADA: ", "% 11.5e", TRUE);*/

      ncm_vector_free (var_z);
      ncm_vector_free (dvar_z);
    }

    ncm_spline2d_set (psf->var, lnr_vec, z_vec, lnvar, TRUE);
    ncm_spline2d_set (psf->dvar, lnr_vec, z_vec, dlnvar, TRUE);

    ncm_vector_free (z_vec);
    ncm_vector_free (lnr_vec);
    ncm_matrix_free (lnvar);
    ncm_matrix_free (dlnvar);

    psf->calibrated = TRUE;
  }
  else
  {
    NcmMatrix *lnvar  = ncm_spline2d_peek_zm (psf->var);
    NcmMatrix *dlnvar = ncm_spline2d_peek_zm (psf->dvar);
    NcmVector *var_yv = ncm_spline2d_peek_yv (psf->var);

    guint N_z = ncm_matrix_nrows (lnvar);
    guint i;

    for (i = 0; i < N_z; i++)
    {
      NcmVector *var_z  = ncm_matrix_get_row (lnvar, i);
      NcmVector *dvar_z = ncm_matrix_get_row (dlnvar, i);

      arg.z = ncm_vector_get (var_yv, i);
      ncm_fftlog_eval_by_gsl_function (psf->fftlog, &F);

      ncm_vector_memcpy (var_z, ncm_fftlog_peek_output_vector (psf->fftlog, 0));
      ncm_vector_memcpy (dvar_z, ncm_fftlog_peek_output_vector (psf->fftlog, 1));

      ncm_vector_free (var_z);
      ncm_vector_free (dvar_z);
    }

    ncm_spline2d_prepare (psf->var);
    ncm_spline2d_prepare (psf->dvar);
  }
}

/**
 * ncm_powspec_filter_prepare_if_needed:
 * @psf: a #NcmPowspecFilter
 * @model: a #NcmModel
 *
 * Prepares (if necessary) the object applying the filter to the power spectrum.
 *
 */
void
ncm_powspec_filter_prepare_if_needed (NcmPowspecFilter *psf, NcmModel *model)
{
  gboolean model_up = ncm_model_ctrl_update (psf->ctrl, model);

  if (model_up)
    ncm_powspec_filter_prepare (psf, model);
}

/**
 * ncm_powspec_filter_set_lnr0:
 * @psf: a #NcmPowspecFilter
 * @lnr0: the output center value $\ln(r_0)$
 *
 * Sets the center of the transform output $\ln(r_0)$ (see ncm_fftlog_set_lnr0()).
 *
 */
void
ncm_powspec_filter_set_lnr0 (NcmPowspecFilter *psf, gdouble lnr0)
{
  if (psf->lnr0 != lnr0)
  {
    const gdouble lnk_min = log (ncm_powspec_get_kmin (psf->ps));
    const gdouble lnk_max = log (ncm_powspec_get_kmax (psf->ps));

    psf->lnk0 = 0.5 * (lnk_max + lnk_min);
    psf->Lk   = (lnk_max - lnk_min);

    ncm_fftlog_set_lnk0 (psf->fftlog, psf->lnk0);
    ncm_fftlog_set_length (psf->fftlog, psf->Lk);

    psf->lnr0 = lnr0;
    ncm_model_ctrl_force_update (psf->ctrl);
    psf->calibrated = FALSE;
    ncm_fftlog_set_lnr0 (psf->fftlog, psf->lnr0);

    if (psf->lnr0 < -psf->lnk0)
      g_warning ("ncm_powspec_filter_set_lnr0: the requested center of the output does not satisfy r0k0 > 1.");
  }
}

/**
 * ncm_powspec_filter_set_best_lnr0:
 * @psf: a #NcmPowspecFilter
 *
 * Sets the value of $\ln(r_0)$ which gives the best results for
 * the transformation based on the current value of $\ln(k_0)$.
 *
 */
void
ncm_powspec_filter_set_best_lnr0 (NcmPowspecFilter *psf)
{
  const gdouble lnk_min = log (ncm_powspec_get_kmin (psf->ps));
  const gdouble lnk_max = log (ncm_powspec_get_kmax (psf->ps));

  psf->lnk0 = 0.5 * (lnk_max + lnk_min);
  psf->Lk   = (lnk_max - lnk_min);

  ncm_fftlog_set_lnk0 (psf->fftlog, psf->lnk0);
  ncm_fftlog_set_length (psf->fftlog, psf->Lk);

  ncm_powspec_filter_set_lnr0 (psf, -psf->lnk0);
}

/**
 * ncm_powspec_filter_set_reltol:
 * @psf: a #NcmPowspecFilter
 * @reltol: the relative tolerance for calibration
 *
 * Sets the relative tolerance for calibration in the distance direction.
 *
 */
void
ncm_powspec_filter_set_reltol (NcmPowspecFilter *psf, const gdouble reltol)
{
  psf->reltol = reltol;
}

/**
 * ncm_powspec_filter_set_reltol_z:
 * @psf: a #NcmPowspecFilter
 * @reltol_z: the relative tolerance for calibration in the redshift direction
 *
 * Sets the relative tolerance for calibration in the redshift direction.
 *
 */
void
ncm_powspec_filter_set_reltol_z (NcmPowspecFilter *psf, const gdouble reltol_z)
{
  psf->reltol_z = reltol_z;
}

/**
 * ncm_powspec_filter_set_zi:
 * @psf: a #NcmPowspecFilter
 * @zi: the output initial time $z_i$
 *
 * Sets the inital time $z_i$.
 *
 */
void
ncm_powspec_filter_set_zi (NcmPowspecFilter *psf, gdouble zi)
{
  if (psf->zi != zi)
  {
    psf->zi = zi;
    ncm_model_ctrl_force_update (psf->ctrl);
    psf->calibrated = FALSE;
    ncm_powspec_require_zi (psf->ps, zi);
  }
}

/**
 * ncm_powspec_filter_set_zf:
 * @psf: a #NcmPowspecFilter
 * @zf: the output final time $z_f$
 *
 * Sets the final time $z_f$.
 *
 */
void
ncm_powspec_filter_set_zf (NcmPowspecFilter *psf, gdouble zf)
{
  if (psf->zf != zf)
  {
    psf->zf = zf;
    ncm_model_ctrl_force_update (psf->ctrl);
    psf->calibrated = FALSE;
    ncm_powspec_require_zf (psf->ps, zf);
  }
}

/**
 * ncm_powspec_filter_require_zi:
 * @psf: a #NcmPowspecFilter
 * @zi: the output initial time $z_i$
 *
 * Require the initial time of at least $z_i$.
 *
 */
void
ncm_powspec_filter_require_zi (NcmPowspecFilter *psf, gdouble zi)
{
  if (psf->zi > zi)
    ncm_powspec_filter_set_zi (psf, zi);
}

/**
 * ncm_powspec_filter_require_zf:
 * @psf: a #NcmPowspecFilter
 * @zf: the output final time $z_f$
 *
 * Requires the final time of at least $z_f$.
 *
 */
void
ncm_powspec_filter_require_zf (NcmPowspecFilter *psf, gdouble zf)
{
  if (psf->zf < zf)
    ncm_powspec_filter_set_zf (psf, zf);
}

/**
 * ncm_powspec_filter_get_filter_type:
 * @psf: a #NcmPowspecFilter
 *
 * Gets the type of filter used.
 *
 * Returns: the type of filter used.
 */
NcmPowspecFilterType
ncm_powspec_filter_get_filter_type (NcmPowspecFilter *psf)
{
  return psf->type;
}

/**
 * ncm_powspec_filter_get_reltol:
 * @psf: a #NcmPowspecFilter
 *
 * Gets the relative tolerance for calibration in the distance direction.
 *
 * Returns: the relative tolerance for calibration in the distance direction.
 */
gdouble
ncm_powspec_filter_get_reltol (NcmPowspecFilter *psf)
{
  return psf->reltol;
}

/**
 * ncm_powspec_filter_get_reltol_z:
 * @psf: a #NcmPowspecFilter
 *
 * Gets the relative tolerance for calibration in the redshift direction.
 *
 * Returns: the relative tolerance for calibration in the redshift direction.
 */
gdouble
ncm_powspec_filter_get_reltol_z (NcmPowspecFilter *psf)
{
  return psf->reltol_z;
}

/**
 * ncm_powspec_filter_get_r_min:
 * @psf: a #NcmPowspecFilter
 *
 * This function returns $\sigma^2(r, z)$'s minimum evaluated distance.
 *
 * Returns: the minimum distance $r_{\mathrm{min}}$.
 */
gdouble
ncm_powspec_filter_get_r_min (NcmPowspecFilter *psf)
{
  return exp (psf->lnr0 - psf->Lk * 0.5);
}

/**
 * ncm_powspec_filter_get_r_max:
 * @psf: a #NcmPowspecFilter
 *
 * This function returns $\sigma^2(r, z)$'s maximum evaluated distance.
 *
 * Returns: the maximum distance $r_{\mathrm{max}}$.
 */
gdouble
ncm_powspec_filter_get_r_max (NcmPowspecFilter *psf)
{
  return exp (psf->lnr0 + psf->Lk * 0.5);
}

/**
 * ncm_powspec_filter_eval_lnvar_lnr:
 * @psf: a #NcmPowspecFilter
 * @z: redshift $z$
 * @lnr: logarithm base e of $r$
 *
 * Evaluates the logarithm base e of the filtered power spectrum at @lnr and @z.
 *
 * Returns: $\ln \left[ \sigma^2(\ln r, z)  \right]$.
 */
gdouble
ncm_powspec_filter_eval_lnvar_lnr (NcmPowspecFilter *psf, const gdouble z, const gdouble lnr)
{
  return log (ncm_spline2d_eval (psf->var, lnr, z));
}

/**
 * ncm_powspec_filter_eval_var_lnr:
 * @psf: a #NcmPowspecFilter
 * @z: redshift $z$
 * @lnr: logarithm base e of $r$
 *
 * Evaluates the filtered power spectrum at @lnr and @z.
 *
 * Returns: $\sigma^2(\ln r, z)$.
 */
gdouble
ncm_powspec_filter_eval_var_lnr (NcmPowspecFilter *psf, const gdouble z, const gdouble lnr)
{
  return ncm_spline2d_eval (psf->var, lnr, z);
}

/**
 * ncm_powspec_filter_eval_var:
 * @psf: a #NcmPowspecFilter
 * @z: redshift $z$
 * @r: distance $r$
 *
 * Evaluate the filtered variance at $r$.
 *
 * Returns: $\sigma^2(r, z)$.
 */
gdouble
ncm_powspec_filter_eval_var (NcmPowspecFilter *psf, const gdouble z, const gdouble r)
{
  return ncm_powspec_filter_eval_var_lnr (psf, z, log (r));
}

/**
 * ncm_powspec_filter_eval_sigma_lnr:
 * @psf: a #NcmPowspecFilter
 * @z: redshift $z$
 * @lnr: logarithm base e of $r$
 *
 * Evaluate the square root of the filtered power spectrum at @lnr and @z.
 *
 * Returns: $\sqrt{ \sigma^2(\ln r, z) }$.
 */
gdouble
ncm_powspec_filter_eval_sigma_lnr (NcmPowspecFilter *psf, const gdouble z, const gdouble lnr)
{
  return sqrt (ncm_powspec_filter_eval_var_lnr (psf, z, lnr));
}

/**
 * ncm_powspec_filter_eval_sigma:
 * @psf: a #NcmPowspecFilter
 * @z: redshift $z$
 * @r: distance $r$
 *
 * Evaluates the square root of the filtered power spectrum at @r and @z.
 *
 * Returns: $\sqrt{ \sigma^2(r, z) }$.
 */
gdouble
ncm_powspec_filter_eval_sigma (NcmPowspecFilter *psf, const gdouble z, const gdouble r)
{
  return ncm_powspec_filter_eval_sigma_lnr (psf, z, log (r));
}

/**
 * ncm_powspec_filter_eval_dvar_dlnr:
 * @psf: a #NcmPowspecFilter
 * @z: redshift $z$
 * @lnr: logarithm base e of $r$
 *
 * Evaluates the first derivative of the filtered
 * variance with respect to $\ln r$ at @lnr and @z.
 *
 * Returns: $\frac{\mathrm{d} \sigma^2(\ln r, z) }{\mathrm{d} \ln r }$.
 */
gdouble
ncm_powspec_filter_eval_dvar_dlnr (NcmPowspecFilter *psf, const gdouble z, const gdouble lnr)
{
  return ncm_spline2d_eval (psf->dvar, lnr, z);
}

/**
 * ncm_powspec_filter_eval_dlnvar_dlnr:
 * @psf: a #NcmPowspecFilter
 * @z: redshift $z$
 * @lnr: logarithm base e of $r$
 *
 * Evaluates the first derivative of the logarithm of the filtered
 * variance with respect to $\ln r$ at @lnr and @z.
 *
 * Returns:  $\frac{\mathrm{d} \left[ \ln \sigma^2(\ln r, z) \right] }{\mathrm{d} \ln r }$.
 */
gdouble
ncm_powspec_filter_eval_dlnvar_dlnr (NcmPowspecFilter *psf, const gdouble z, const gdouble lnr)
{
  return ncm_spline2d_eval (psf->dvar, lnr, z) / ncm_spline2d_eval (psf->var, lnr, z);
}

/**
 * ncm_powspec_filter_eval_dlnvar_dr:
 * @psf: a #NcmPowspecFilter
 * @z: redshift $z$
 * @lnr: logarithm base e of $r$
 *
 * Evaluates the first derivative of the logarithm of the filtered
 * variance with respect to $r$ at @lnr and @z.
 *
 * Returns:  $\frac{\mathrm{d} \left[ \ln \sigma^2(\ln r, z) \right] }{\mathrm{d} r }$.
 */
gdouble
ncm_powspec_filter_eval_dlnvar_dr (NcmPowspecFilter *psf, const gdouble z, const gdouble lnr)
{
  return ncm_powspec_filter_eval_dlnvar_dlnr (psf, z, lnr) * exp (-lnr);
}

/**
 * ncm_powspec_filter_eval_dnvar_dlnrn:
 * @psf: a #NcmPowspecFilter
 * @z: redshift $z$
 * @lnr: logarithm base e of $r$
 * @n: number of derivatives $n$
 *
 * Evaluates the derivatives of the filtered variance at @lnr and @z, namely:
 * - $n = 0 \rightarrow \sigma(r, z)^2$,
 * - $n = 1 \rightarrow \frac{\mathrm{d}\sigma^2}{\mathrm{d} \ln r}$,
 * - $n = 2 \rightarrow \frac{\mathrm{d}^2\sigma^2}{\mathrm{d}(\ln r)^2}$,
 * - $n = 3 \rightarrow \frac{\mathrm{d}^3\sigma^2}{\mathrm{d}(\ln r)^3}$.
 *
 * Returns: one of the four derivatives described above.
 */
gdouble
ncm_powspec_filter_eval_dnvar_dlnrn (NcmPowspecFilter *psf, const gdouble z, const gdouble lnr, guint n)
{
  switch (n)
  {
    case 0:

      return ncm_spline2d_eval (psf->var, lnr, z);

      break;
    case 1:

      return ncm_spline2d_eval (psf->dvar, lnr, z);

      break;
    case 2:

      return ncm_spline2d_deriv_dzdx (psf->dvar, lnr, z);

      break;
    case 3:

      return ncm_spline2d_deriv_d2zdx2 (psf->dvar, lnr, z);

      break;
    default:
      g_error ("ncm_powspec_filter_eval_dnvar_dlnrn: %u derivative not implemented.", n);

      return 0.0;

      break;
  }
}

/**
 * ncm_powspec_filter_eval_dnlnvar_dlnrn:
 * @psf: a #NcmPowspecFilter
 * @z: redshift $z$
 * @lnr: logarithm base e of $r$
 * @n: number of derivatives $n$
 *
 * Evaluates the derivatives of the logarithm of the filtered variance at @lnr and @z, namely:
 * - $n = 0 \rightarrow \ln \left[ \sigma(r, z)^2 \right]$,
 * - $n = 1 \rightarrow \frac{\mathrm{d}\ln \left( \sigma^2 \right)}{\mathrm{d} \ln r}$,
 * - $n = 2 \rightarrow \frac{\mathrm{d}^2 \ln \left( \sigma^2 \right)}{\mathrm{d}(\ln r)^2}$,
 * - $n = 3 \rightarrow \frac{\mathrm{d}^3 \ln \left( \sigma^2 \right)}{\mathrm{d}(\ln r)^3}$.
 *
 * Returns: one of the four derivatives described above.
 */
gdouble
ncm_powspec_filter_eval_dnlnvar_dlnrn (NcmPowspecFilter *psf, const gdouble z, const gdouble lnr, guint n)
{
  switch (n)
  {
    case 0:

      return ncm_powspec_filter_eval_lnvar_lnr (psf, z, lnr);

      break;
    case 1:

      return ncm_powspec_filter_eval_dlnvar_dlnr (psf, z, lnr);

      break;
    case 2:
    {
      const gdouble var   = ncm_spline2d_eval (psf->var, lnr, z);
      const gdouble dvar  = ncm_spline2d_eval (psf->dvar, lnr, z);
      const gdouble d2var = ncm_spline2d_deriv_dzdx (psf->dvar, lnr, z);

      const gdouble dlnvar = dvar / var;

      return d2var / var - dlnvar * dlnvar;

      break;
    }
    default:
      g_error ("ncm_powspec_filter_eval_dnlnvar_dlnrn: %u derivative not implemented.", n);

      return 0.0;

      break;
  }
}

/**
 * ncm_powspec_filter_volume_rm3:
 * @psf: a #NcmPowspecFilter
 *
 * Calculates the volume of the filter over $r^3$.
 *
 * Returns: Filter's volume over the radius squared $V r^{-3}$.
 */
gdouble
ncm_powspec_filter_volume_rm3 (NcmPowspecFilter *psf)
{
  const gdouble tophat_volumeRm3 = 4.0 * M_PI / 3.0;
  const gdouble gauss_volumeRm3  = sqrt (2.0 * M_PI) * sqrt (2.0 * M_PI) * sqrt (2.0 * M_PI);

  switch (psf->type)
  {
    case NCM_POWSPEC_FILTER_TYPE_TOPHAT:

      return tophat_volumeRm3;

      break;
    case NCM_POWSPEC_FILTER_TYPE_GAUSS:

      return gauss_volumeRm3;

      break;
    default:
      g_assert_not_reached ();

      return 0.0;

      break;
  }
}

/**
 * ncm_powspec_filter_peek_powspec:
 * @psf: a #NcmPowspecFilter
 *
 * Gets the #NcmPowspec object used to compute the filtered variance $\sigma^{2}(r,z)$.
 *
 * Returns: (transfer none): the #NcmPowspec object.
 */
NcmPowspec *
ncm_powspec_filter_peek_powspec (NcmPowspecFilter *psf)
{
  return psf->ps;
}

