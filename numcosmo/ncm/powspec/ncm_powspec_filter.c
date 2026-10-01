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
 * Variance of a power spectrum smoothed by a window, on a grid in $(r, z)$.
 *
 * Computes
 * $$\sigma^2(r, z) = \frac{1}{2\pi^2} \int k^2 P(k, z) \, W^2(kr) \, \mathrm{d}k,$$
 * with $W$ the top-hat or Gaussian window (#NcmPowspecFilterType), and its derivatives
 * with respect to $\ln r$ (ncm_powspec_filter_eval_dnvar_dlnrn()), by an #NcmFftlog
 * transform, and interpolates them with bicubic splines in $(\ln r, z)$.
 *
 * The transform runs over the $k$ range of the power spectrum, continued beyond it by
 * ncm_fftlog_use_smooth_padding(), with the bias that ncm_fftlog_get_best_bias() chooses
 * from the log-slopes of $k^2 P(k, z)$ at the ends of the table. Within a few e-foldings
 * of $R = 1/k_\mathrm{max}$ and $R = 1/k_\mathrm{min}$ the result depends on that
 * continuation, an extrapolation of the table; ncm_powspec_filter_get_r_min() and
 * ncm_powspec_filter_get_r_max() return the whole grid, those edges included.
 *
 * The evaluation functions do not check their arguments: they are valid for $r$ in
 * [ncm_powspec_filter_get_r_min(), ncm_powspec_filter_get_r_max()] and $z$ in
 * [NcmPowspecFilter:zi, NcmPowspecFilter:zf], and outside they return the extrapolation
 * of the splines.
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "ncm/powspec/ncm_powspec_filter.h"
#include "ncm/spline/ncm_spline_cubic_notaknot.h"
#include "ncm/spline/ncm_spline2d_bicubic.h"
#include "ncm/fftlog/ncm_fftlog_tophatwin2.h"
#include "ncm/fftlog/ncm_fftlog_gausswin2.h"
#include "ncm/core/ncm_c.h"
#include "ncm_enum_types.h"

/* Size the calibration of the k grid starts from. */
#define NCM_POWSPEC_FILTER_START_K_KNOTS 100

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
  PROP_NDERIVS,
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
  guint nderivs;
  GPtrArray *dnvar;
  NcmModelCtrl *ctrl;
  gboolean constructed;
};

/*
 * Ensures that dnvar holds one spline per derivative order in [0, nderivs].
 * Orders are only ever added, so splines already prepared are preserved.
 */
static void
_ncm_powspec_filter_alloc_dnvar (NcmPowspecFilter *psf)
{
  while (psf->dnvar->len < psf->nderivs + 1)
    g_ptr_array_add (psf->dnvar, ncm_spline2d_bicubic_notaknot_new ());

  if (psf->dnvar->len > psf->nderivs + 1)
    g_ptr_array_set_size (psf->dnvar, psf->nderivs + 1);
}

static NcmSpline2d *
_ncm_powspec_filter_peek_dnvar (NcmPowspecFilter *psf, const guint n)
{
  g_assert_cmpuint (n, <, psf->dnvar->len);

  return g_ptr_array_index (psf->dnvar, n);
}

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
  psf->nderivs     = 0;
  psf->dnvar       = g_ptr_array_new_with_free_func ((GDestroyNotify) & ncm_spline2d_free);
  psf->ctrl        = ncm_model_ctrl_new (NULL);
  psf->constructed = FALSE;
}

/* A calibration depends on the tolerances and knot limits: changing one after a prepare
 * must calibrate again at the next prepare, not keep the old grid. */
static void
_ncm_powspec_filter_invalidate (NcmPowspecFilter *psf)
{
  psf->calibrated = FALSE;

  if (psf->ctrl != NULL)
    ncm_model_ctrl_force_update (psf->ctrl);
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
      ncm_powspec_filter_set_reltol (psf, g_value_get_double (value));
      break;
    case PROP_RELTOL_Z:
      ncm_powspec_filter_set_reltol_z (psf, g_value_get_double (value));
      break;
    case PROP_MAX_K_KNOTS:
      ncm_powspec_filter_set_max_k_knots (psf, g_value_get_uint (value));
      break;
    case PROP_MAX_Z_KNOTS:
      ncm_powspec_filter_set_max_z_knots (psf, g_value_get_uint (value));
      break;
    case PROP_NDERIVS:
      ncm_powspec_filter_set_nderivs (psf, g_value_get_uint (value));
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
    case PROP_NDERIVS:
      g_value_set_uint (value, psf->nderivs);
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

  g_clear_pointer (&psf->dnvar, g_ptr_array_unref);

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
   * Center $\ln r_0$ of the output grid, $r$ in Mpc.
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
   * Lower end of the redshift range of $\sigma^2(r, z)$.
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
   * Upper end of the redshift range of $\sigma^2(r, z)$.
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
   * Tolerance of the calibration in the distance direction: the number of knots grows
   * until $\sigma^2$ and its derivatives change by less than this between sizes, relative
   * to each one's peak over the whole $R$ grid, not to its value at each $R$. Where
   * $\sigma^2$ is far below its peak, at large $R$, the relative accuracy is lower.
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
   * Relative tolerance of the calibration of the $z$ knots, on $\sigma^2$ at the
   * smallest $r$.
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
   * The maximum number of knots in the k direction; ncm_powspec_filter_prepare() aborts if
   * the calibration would need more (see ncm_fftlog_calibrate_size_gsl()).
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
   * The maximum number of knots in the redshift direction, zero for no limit;
   * ncm_powspec_filter_prepare() aborts if the grid would need more to reach
   * #NcmPowspecFilter:reltol-z.
   */
  g_object_class_install_property (object_class,
                                   PROP_MAX_Z_KNOTS,
                                   g_param_spec_uint ("max-z-knots",
                                                      NULL,
                                                      "Maximum number of knots in the redshift direction",
                                                      0, G_MAXUINT, 1000,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmPowspecFilter:nderivs:
   *
   * Highest order $n$ of $\mathrm{d}^n\sigma^2/\mathrm{d}(\ln r)^n$ computed by the
   * transform itself. Each additional order costs one extra Mellin kernel, output
   * vector and two-dimensional spline, so the default of one covers the common case.
   *
   * Do not set this directly to declare a requirement; use
   * ncm_powspec_filter_require_nderivs() instead, so that independent users of the
   * same filter cannot lower an order somebody else still needs.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_NDERIVS,
                                   g_param_spec_uint ("nderivs",
                                                      NULL,
                                                      "Number of derivatives computed by the transform",
                                                      1, G_MAXUINT, 1,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmPowspecFilter:type:
   *
   * The window $W$.
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
   * The #NcmPowspec whose variance is computed.
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
 * Returns: (transfer full): a new #NcmPowspecFilter for @ps with the window @type
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
 * If *@psf is not %NULL, decrements its reference count and sets *@psf to %NULL.
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
 * Sets the window to @type; a change rebuilds the transform and the next prepare
 * recalibrates.
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
        psf->fftlog = NCM_FFTLOG (ncm_fftlog_tophatwin2_new (psf->lnr0, psf->lnk0, psf->Lk, NCM_POWSPEC_FILTER_START_K_KNOTS));
        break;
      case NCM_POWSPEC_FILTER_TYPE_GAUSS:
        psf->fftlog = NCM_FFTLOG (ncm_fftlog_gausswin2_new (psf->lnr0, psf->lnk0, psf->Lk, NCM_POWSPEC_FILTER_START_K_KNOTS));
        break;
      default:
        g_assert_not_reached ();
        break;
    }

    ncm_fftlog_set_padding (psf->fftlog, 1.0);

    /* k^2 P(k) does not vanish at the ends of the table; zeros there ring at R near the
     * inverse ends and never converge. */
    ncm_fftlog_use_smooth_padding (psf->fftlog, TRUE);
    ncm_fftlog_set_nderivs (psf->fftlog, psf->nderivs);

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

  return ncm_vector_get (ncm_fftlog_peek_output_vector (arg->psf->fftlog, 0), 0);
}

/**
 * ncm_powspec_filter_prepare:
 * @psf: a #NcmPowspecFilter
 * @model: (allow-none): a #NcmModel
 *
 * Prepares the power spectrum if needed, calibrates the grid when it is not
 * calibrated, and fills $\sigma^2(r, z)$ and its derivatives on it.
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
  arg.z     = psf->zi;

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

  ncm_fftlog_set_nderivs (psf->fftlog, psf->nderivs);
  _ncm_powspec_filter_alloc_dnvar (psf);

  if (!psf->calibrated)
  {
    NcmMatrix **dnvar;
    NcmVector *z_vec, *lnr_vec;
    guint N_k, N_z;
    guint i, nd;

    /* Every calibration starts from the same size: starting from the last one, which the
     * calibration always passes, grew the grid at each recalibration. The bias comes from
     * the end slopes of the table at zi and holds for every
     * redshift, so that all share one transform. */
    ncm_fftlog_set_size (psf->fftlog, NCM_POWSPEC_FILTER_START_K_KNOTS);
    ncm_fftlog_eval_by_gsl_function (psf->fftlog, &F);
    ncm_fftlog_set_bias (psf->fftlog, ncm_fftlog_get_best_bias (psf->fftlog));

    ncm_fftlog_set_max_size (psf->fftlog, psf->max_k_knots);
    ncm_fftlog_calibrate_size_gsl (psf->fftlog, &F, psf->reltol);
    N_k = ncm_fftlog_get_size (psf->fftlog);

    {
      NcmSpline *dummy_z = NCM_SPLINE (ncm_spline_cubic_notaknot_new ());
      gsl_function Fdummy_z;

      Fdummy_z.function = &_ncm_powspec_filter_dummy_z;
      Fdummy_z.params   = &arg;

      ncm_spline_set_func (dummy_z, NCM_SPLINE_FUNCTION_SPLINE, &Fdummy_z, psf->zi, psf->zf, psf->max_z_knots, psf->reltol_z);

      z_vec = ncm_spline_get_xv (dummy_z);
      N_z   = ncm_vector_len (z_vec);

      /* The spline stops with a warning past its limit; the filter must not go on with a
       * grid that missed reltol-z. */
      if ((psf->max_z_knots > 0) && (N_z > psf->max_z_knots))
        g_error ("ncm_powspec_filter_prepare: the redshift grid needs more than %u knots (max-z-knots) "
                 "to reach the relative tolerance %e (reltol-z).", psf->max_z_knots, psf->reltol_z);

      ncm_spline_clear (&dummy_z);
    }

    g_assert_cmpuint (N_z, >, 0);
    g_assert_cmpuint (N_k, >, 0);

    dnvar = g_new0 (NcmMatrix *, psf->nderivs + 1);

    for (nd = 0; nd <= psf->nderivs; nd++)
      dnvar[nd] = ncm_matrix_new (N_z, N_k);

    lnr_vec = ncm_fftlog_get_vector_lnr (psf->fftlog);

    for (i = 0; i < N_z; i++)
    {
      arg.z = ncm_vector_get (z_vec, i);
      ncm_fftlog_eval_by_gsl_function (psf->fftlog, &F);

      for (nd = 0; nd <= psf->nderivs; nd++)
      {
        NcmVector *dnvar_z = ncm_matrix_get_row (dnvar[nd], i);

        ncm_vector_memcpy (dnvar_z, ncm_fftlog_peek_output_vector (psf->fftlog, nd));
        ncm_vector_free (dnvar_z);
      }
    }

    for (nd = 0; nd <= psf->nderivs; nd++)
    {
      ncm_spline2d_set (_ncm_powspec_filter_peek_dnvar (psf, nd), lnr_vec, z_vec, dnvar[nd], TRUE);
      ncm_matrix_free (dnvar[nd]);
    }

    g_free (dnvar);

    ncm_vector_free (z_vec);
    ncm_vector_free (lnr_vec);

    psf->calibrated = TRUE;
  }
  else
  {
    NcmSpline2d *var  = _ncm_powspec_filter_peek_dnvar (psf, 0);
    NcmVector *var_yv = ncm_spline2d_peek_yv (var);

    guint N_z = ncm_matrix_nrows (ncm_spline2d_peek_zm (var));
    guint i, nd;

    for (i = 0; i < N_z; i++)
    {
      arg.z = ncm_vector_get (var_yv, i);
      ncm_fftlog_eval_by_gsl_function (psf->fftlog, &F);

      for (nd = 0; nd <= psf->nderivs; nd++)
      {
        NcmVector *dnvar_z = ncm_matrix_get_row (ncm_spline2d_peek_zm (_ncm_powspec_filter_peek_dnvar (psf, nd)), i);

        ncm_vector_memcpy (dnvar_z, ncm_fftlog_peek_output_vector (psf->fftlog, nd));
        ncm_vector_free (dnvar_z);
      }
    }

    for (nd = 0; nd <= psf->nderivs; nd++)
      ncm_spline2d_prepare (_ncm_powspec_filter_peek_dnvar (psf, nd));
  }

  /* ncm_model_ctrl_update() dereferences its model. */
  if (model != NULL)
    ncm_model_ctrl_update (psf->ctrl, model);
}

/**
 * ncm_powspec_filter_prepare_if_needed:
 * @psf: a #NcmPowspecFilter
 * @model: (allow-none): a #NcmModel
 *
 * Calls ncm_powspec_filter_prepare() if @model or a setting of @psf changed since
 * the last preparation, and always when @model is %NULL.
 */
void
ncm_powspec_filter_prepare_if_needed (NcmPowspecFilter *psf, NcmModel *model)
{
  if ((model == NULL) || ncm_model_ctrl_update (psf->ctrl, model))
    ncm_powspec_filter_prepare (psf, model);
}

/**
 * ncm_powspec_filter_set_lnr0:
 * @psf: a #NcmPowspecFilter
 * @lnr0: output center $\ln r_0$
 *
 * Sets the center of the output grid (see ncm_fftlog_set_lnr0()); the next prepare
 * recalibrates. Warns when $k_0 r_0 < 1$.
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
 * Sets the transform over the $k$ range of the power spectrum with $k_0 r_0 = 1$, so the
 * output grid is $r \in [1/k_\mathrm{max}, 1/k_\mathrm{min}]$.
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
 * @reltol: relative tolerance
 *
 * Sets NcmPowspecFilter:reltol; a change makes the next prepare recalibrate.
 */
void
ncm_powspec_filter_set_reltol (NcmPowspecFilter *psf, const gdouble reltol)
{
  if (psf->reltol != reltol)
  {
    psf->reltol = reltol;
    _ncm_powspec_filter_invalidate (psf);
  }
}

/**
 * ncm_powspec_filter_set_reltol_z:
 * @psf: a #NcmPowspecFilter
 * @reltol_z: relative tolerance
 *
 * Sets NcmPowspecFilter:reltol-z; a change makes the next prepare recalibrate.
 */
void
ncm_powspec_filter_set_reltol_z (NcmPowspecFilter *psf, const gdouble reltol_z)
{
  if (psf->reltol_z != reltol_z)
  {
    psf->reltol_z = reltol_z;
    _ncm_powspec_filter_invalidate (psf);
  }
}

/**
 * ncm_powspec_filter_set_zi:
 * @psf: a #NcmPowspecFilter
 * @zi: lowest redshift $z_i$
 *
 * Sets NcmPowspecFilter:zi and requires it of the power spectrum (see
 * ncm_powspec_require_zi()); the next prepare recalibrates.
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
 * @zf: highest redshift $z_f$
 *
 * Sets NcmPowspecFilter:zf and requires it of the power spectrum (see
 * ncm_powspec_require_zf()); the next prepare recalibrates.
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
 * @zi: lowest redshift $z_i$
 *
 * Lowers NcmPowspecFilter:zi to @zi if @zi is below it.
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
 * @zf: highest redshift $z_f$
 *
 * Raises NcmPowspecFilter:zf to @zf if @zf is above it.
 */
void
ncm_powspec_filter_require_zf (NcmPowspecFilter *psf, gdouble zf)
{
  if (psf->zf < zf)
    ncm_powspec_filter_set_zf (psf, zf);
}

/**
 * ncm_powspec_filter_set_nderivs:
 * @psf: a #NcmPowspecFilter
 * @nderivs: the highest derivative order $n$
 *
 * Sets the highest order $n$ of $\mathrm{d}^n\sigma^2/\mathrm{d}(\ln r)^n$ obtained
 * from the transform itself. Lowering the order discards the extra tables, so
 * prefer ncm_powspec_filter_require_nderivs() whenever the filter may be shared.
 */
void
ncm_powspec_filter_set_nderivs (NcmPowspecFilter *psf, guint nderivs)
{
  g_assert_cmpuint (nderivs, >, 0);

  if (psf->nderivs != nderivs)
  {
    psf->nderivs = nderivs;

    _ncm_powspec_filter_alloc_dnvar (psf);

    if (psf->fftlog != NULL)
      ncm_fftlog_set_nderivs (psf->fftlog, nderivs);

    ncm_model_ctrl_force_update (psf->ctrl);
    psf->calibrated = FALSE;
  }
}

/**
 * ncm_powspec_filter_require_nderivs:
 * @psf: a #NcmPowspecFilter
 * @nderivs: the required derivative order $n$
 *
 * Requires derivatives up to at least order $n$. Requests at or below the order
 * already in use do nothing, so several users of the same filter may each state
 * their own minimum without any of them lowering an order another one needs.
 */
void
ncm_powspec_filter_require_nderivs (NcmPowspecFilter *psf, guint nderivs)
{
  if (psf->nderivs < nderivs)
    ncm_powspec_filter_set_nderivs (psf, nderivs);
}

/**
 * ncm_powspec_filter_get_nderivs:
 * @psf: a #NcmPowspecFilter
 *
 * Returns: NcmPowspecFilter:nderivs, the highest derivative order computed
 */
guint
ncm_powspec_filter_get_nderivs (NcmPowspecFilter *psf)
{
  return psf->nderivs;
}

/**
 * ncm_powspec_filter_get_filter_type:
 * @psf: a #NcmPowspecFilter
 *
 * Returns: the window, NcmPowspecFilter:type
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
 * Returns: NcmPowspecFilter:reltol
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
 * Returns: NcmPowspecFilter:reltol-z
 */
gdouble
ncm_powspec_filter_get_reltol_z (NcmPowspecFilter *psf)
{
  return psf->reltol_z;
}

/**
 * ncm_powspec_filter_set_max_k_knots:
 * @psf: a #NcmPowspecFilter
 * @max_k_knots: the maximum number of knots in $k$
 *
 * Sets #NcmPowspecFilter:max-k-knots, the most knots the calibration may use. A change
 * makes the next ncm_powspec_filter_prepare() calibrate again.
 */
void
ncm_powspec_filter_set_max_k_knots (NcmPowspecFilter *psf, guint max_k_knots)
{
  if (psf->max_k_knots != max_k_knots)
  {
    psf->max_k_knots = max_k_knots;
    _ncm_powspec_filter_invalidate (psf);
  }
}

/**
 * ncm_powspec_filter_get_max_k_knots:
 * @psf: a #NcmPowspecFilter
 *
 * Returns: the #NcmPowspecFilter:max-k-knots
 */
guint
ncm_powspec_filter_get_max_k_knots (NcmPowspecFilter *psf)
{
  return psf->max_k_knots;
}

/**
 * ncm_powspec_filter_set_max_z_knots:
 * @psf: a #NcmPowspecFilter
 * @max_z_knots: the maximum number of knots in $z$
 *
 * Sets #NcmPowspecFilter:max-z-knots. A change makes the next ncm_powspec_filter_prepare()
 * calibrate again.
 */
void
ncm_powspec_filter_set_max_z_knots (NcmPowspecFilter *psf, guint max_z_knots)
{
  if (psf->max_z_knots != max_z_knots)
  {
    psf->max_z_knots = max_z_knots;
    _ncm_powspec_filter_invalidate (psf);
  }
}

/**
 * ncm_powspec_filter_get_max_z_knots:
 * @psf: a #NcmPowspecFilter
 *
 * Returns: the #NcmPowspecFilter:max-z-knots
 */
guint
ncm_powspec_filter_get_max_z_knots (NcmPowspecFilter *psf)
{
  return psf->max_z_knots;
}

/**
 * ncm_powspec_filter_get_nknots:
 * @psf: a #NcmPowspecFilter
 * @N_k: (out): number of knots in $\ln r$
 * @N_z: (out): number of knots in $z$
 *
 * Gets the size of the grid the last calibration chose, both zero before the first
 * ncm_powspec_filter_prepare().
 */
void
ncm_powspec_filter_get_nknots (NcmPowspecFilter *psf, guint *N_k, guint *N_z)
{
  NcmSpline2d *var = (psf->dnvar->len > 0) ? g_ptr_array_index (psf->dnvar, 0) : NULL;

  if ((var == NULL) || !psf->calibrated)
  {
    *N_k = 0;
    *N_z = 0;

    return;
  }

  *N_k = ncm_vector_len (ncm_spline2d_peek_xv (var));
  *N_z = ncm_vector_len (ncm_spline2d_peek_yv (var));
}

/**
 * ncm_powspec_filter_get_r_min:
 * @psf: a #NcmPowspecFilter
 *
 * The minimum distance at which $\sigma^2(r, z)$ is tabulated: the first knot of the
 * calibrated grid, or its estimate from #NcmPowspecFilter:lnr0 before the first
 * ncm_powspec_filter_prepare().
 *
 * Returns: $r_\mathrm{min}$ in Mpc
 */
gdouble
ncm_powspec_filter_get_r_min (NcmPowspecFilter *psf)
{
  /* Once calibrated, the grid's own end: the no-ringing shift moves it by up to a knot,
   * and the knots stop one spacing short of ln r0 +- L / 2. */
  if (psf->calibrated && (psf->dnvar->len > 0))
  {
    NcmVector *lnr = ncm_spline2d_peek_xv (g_ptr_array_index (psf->dnvar, 0));

    return exp (ncm_vector_get (lnr, 0));
  }

  return exp (psf->lnr0 - psf->Lk * 0.5);
}

/**
 * ncm_powspec_filter_get_r_max:
 * @psf: a #NcmPowspecFilter
 *
 * The maximum distance at which $\sigma^2(r, z)$ is tabulated: the last knot of the
 * calibrated grid, or its estimate from #NcmPowspecFilter:lnr0 before the first
 * ncm_powspec_filter_prepare().
 *
 * Returns: $r_\mathrm{max}$ in Mpc
 */
gdouble
ncm_powspec_filter_get_r_max (NcmPowspecFilter *psf)
{
  /* Once calibrated, the grid's own end: the no-ringing shift moves it by up to a knot,
   * and the knots stop one spacing short of ln r0 +- L / 2. */
  if (psf->calibrated && (psf->dnvar->len > 0))
  {
    NcmVector *lnr = ncm_spline2d_peek_xv (g_ptr_array_index (psf->dnvar, 0));

    return exp (ncm_vector_get (lnr, ncm_vector_len (lnr) - 1));
  }

  return exp (psf->lnr0 + psf->Lk * 0.5);
}

/**
 * ncm_powspec_filter_eval_lnvar_lnr:
 * @psf: a #NcmPowspecFilter
 * @z: redshift
 * @lnr: $\ln r$, $r$ in Mpc
 *
 * Returns: $\ln \sigma^2(r, z)$
 */
gdouble
ncm_powspec_filter_eval_lnvar_lnr (NcmPowspecFilter *psf, const gdouble z, const gdouble lnr)
{
  return log (ncm_spline2d_eval (_ncm_powspec_filter_peek_dnvar (psf, 0), lnr, z));
}

/**
 * ncm_powspec_filter_eval_var_lnr:
 * @psf: a #NcmPowspecFilter
 * @z: redshift
 * @lnr: $\ln r$, $r$ in Mpc
 *
 * Returns: $\sigma^2(r, z)$
 */
gdouble
ncm_powspec_filter_eval_var_lnr (NcmPowspecFilter *psf, const gdouble z, const gdouble lnr)
{
  return ncm_spline2d_eval (_ncm_powspec_filter_peek_dnvar (psf, 0), lnr, z);
}

/**
 * ncm_powspec_filter_eval_var:
 * @psf: a #NcmPowspecFilter
 * @z: redshift
 * @r: radius in Mpc
 *
 * Returns: $\sigma^2(r, z)$
 */
gdouble
ncm_powspec_filter_eval_var (NcmPowspecFilter *psf, const gdouble z, const gdouble r)
{
  return ncm_powspec_filter_eval_var_lnr (psf, z, log (r));
}

/**
 * ncm_powspec_filter_eval_sigma_lnr:
 * @psf: a #NcmPowspecFilter
 * @z: redshift
 * @lnr: $\ln r$, $r$ in Mpc
 *
 * Returns: $\sigma(r, z)$
 */
gdouble
ncm_powspec_filter_eval_sigma_lnr (NcmPowspecFilter *psf, const gdouble z, const gdouble lnr)
{
  return sqrt (ncm_powspec_filter_eval_var_lnr (psf, z, lnr));
}

/**
 * ncm_powspec_filter_eval_sigma:
 * @psf: a #NcmPowspecFilter
 * @z: redshift
 * @r: radius in Mpc
 *
 * Returns: $\sigma(r, z)$
 */
gdouble
ncm_powspec_filter_eval_sigma (NcmPowspecFilter *psf, const gdouble z, const gdouble r)
{
  return ncm_powspec_filter_eval_sigma_lnr (psf, z, log (r));
}

/**
 * ncm_powspec_filter_eval_dvar_dlnr:
 * @psf: a #NcmPowspecFilter
 * @z: redshift
 * @lnr: $\ln r$, $r$ in Mpc
 *
 * Returns: $\mathrm{d}\sigma^2 / \mathrm{d}\ln r$ at $(r, z)$
 */
gdouble
ncm_powspec_filter_eval_dvar_dlnr (NcmPowspecFilter *psf, const gdouble z, const gdouble lnr)
{
  return ncm_spline2d_eval (_ncm_powspec_filter_peek_dnvar (psf, 1), lnr, z);
}

/**
 * ncm_powspec_filter_eval_dlnvar_dlnr:
 * @psf: a #NcmPowspecFilter
 * @z: redshift
 * @lnr: $\ln r$, $r$ in Mpc
 *
 * Returns: $\mathrm{d}\ln\sigma^2 / \mathrm{d}\ln r$ at $(r, z)$
 */
gdouble
ncm_powspec_filter_eval_dlnvar_dlnr (NcmPowspecFilter *psf, const gdouble z, const gdouble lnr)
{
  return ncm_spline2d_eval (_ncm_powspec_filter_peek_dnvar (psf, 1), lnr, z) /
         ncm_spline2d_eval (_ncm_powspec_filter_peek_dnvar (psf, 0), lnr, z);
}

/**
 * ncm_powspec_filter_eval_dlnvar_dr:
 * @psf: a #NcmPowspecFilter
 * @z: redshift
 * @lnr: $\ln r$, $r$ in Mpc
 *
 * Returns: $\mathrm{d}\ln\sigma^2 / \mathrm{d}r$ at $(r, z)$, in $\mathrm{Mpc}^{-1}$
 */
gdouble
ncm_powspec_filter_eval_dlnvar_dr (NcmPowspecFilter *psf, const gdouble z, const gdouble lnr)
{
  return ncm_powspec_filter_eval_dlnvar_dlnr (psf, z, lnr) * exp (-lnr);
}

/**
 * ncm_powspec_filter_eval_dnvar_dlnrn:
 * @psf: a #NcmPowspecFilter
 * @z: redshift
 * @lnr: $\ln r$, $r$ in Mpc
 * @n: derivative order
 *
 * Evaluates $\mathrm{d}^n\sigma^2 / \mathrm{d}(\ln r)^n$ at $(r, z)$, with
 * $n = 0$ giving $\sigma^2(r, z)$ itself.
 *
 * Every order comes from the transform itself rather than from differentiating an
 * interpolation, so @n must not exceed the order the filter was prepared for; see
 * ncm_powspec_filter_require_nderivs().
 *
 * Returns: $\mathrm{d}^n\sigma^2 / \mathrm{d}(\ln r)^n$
 */
gdouble
ncm_powspec_filter_eval_dnvar_dlnrn (NcmPowspecFilter *psf, const gdouble z, const gdouble lnr, guint n)
{
  if (n > psf->nderivs)
    g_error ("ncm_powspec_filter_eval_dnvar_dlnrn: derivative %u requested but the transform "
             "computes up to %u, call ncm_powspec_filter_require_nderivs () before preparing.", n, psf->nderivs);

  return ncm_spline2d_eval (_ncm_powspec_filter_peek_dnvar (psf, n), lnr, z);
}

/**
 * ncm_powspec_filter_eval_dnlnvar_dlnrn:
 * @psf: a #NcmPowspecFilter
 * @z: redshift
 * @lnr: $\ln r$, $r$ in Mpc
 * @n: derivative order, 0, 1 or 2
 *
 * Evaluates $\mathrm{d}^n\ln\sigma^2 / \mathrm{d}(\ln r)^n$ at $(r, z)$ from the
 * derivatives of $\sigma^2$; $n = 2$ needs NcmPowspecFilter:nderivs of at least 2,
 * and $n > 2$ aborts.
 *
 * Returns: $\mathrm{d}^n\ln\sigma^2 / \mathrm{d}(\ln r)^n$
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
      const gdouble var   = ncm_powspec_filter_eval_dnvar_dlnrn (psf, z, lnr, 0);
      const gdouble dvar  = ncm_powspec_filter_eval_dnvar_dlnrn (psf, z, lnr, 1);
      const gdouble d2var = ncm_powspec_filter_eval_dnvar_dlnrn (psf, z, lnr, 2);

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
 * Returns: the volume of the window divided by $r^3$: $4\pi/3$ for the top-hat and
 * $(2\pi)^{3/2}$ for the Gaussian
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
 * Returns: (transfer none): the #NcmPowspec, NcmPowspecFilter:powerspectrum
 */
NcmPowspec *
ncm_powspec_filter_peek_powspec (NcmPowspecFilter *psf)
{
  return psf->ps;
}

