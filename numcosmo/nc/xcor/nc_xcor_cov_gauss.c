/***************************************************************************
 *            nc_xcor_cov_gauss.c
 *
 *  Tue November 03 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) 2026 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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
 * NcXcorCovGauss:
 *
 * Gaussian covariance of band-power angular power spectra on a fraction of the sky.
 *
 * For tracers $a, b, c, d$ and band powers $L$, $L'$ with windows $w_L$, $w_{L'}$,
 * $$\mathrm{Cov}\left[C^{ab}_L, C^{cd}_{L'}\right] = \sum_{\ell} w_L(\ell)\, w_{L'}(\ell)\,
 * \frac{\tilde C^{ac}_\ell \tilde C^{bd}_\ell + \tilde C^{ad}_\ell \tilde C^{bc}_\ell}{(2\ell + 1) f_{\rm sky}},$$
 * where $\tilde C^{ab}_\ell = C^{ab}_\ell + \delta_{ab} N^a_\ell$ and the sum runs over the
 * common multipole range [#NcXcorCovGauss:ell-min, #NcXcorCovGauss:ell-max].
 *
 * The data vector is a list of spectra, each a tracer pair with its window, one row per band
 * power; nc_xcor_cov_gauss_add_spectrum() appends one. The window weights are used as given.
 * nc_xcor_cov_gauss_compute() takes the signal $C^{ab}_\ell$ of every tracer pair, in the
 * rows given by nc_xcor_cov_gauss_get_pair_index(), and the noise $N^a_\ell$ of every
 * tracer. The object holds no model; it writes its own workspaces, so one instance must not
 * be used by two threads at once.
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "nc/xcor/nc_xcor_cov_gauss.h"
#include "ncm/algebra/ncm_vector.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_blas.h>
#endif /* NUMCOSMO_GIR_SCAN */

enum
{
  PROP_0,
  PROP_N_TRACERS,
  PROP_ELL_MIN,
  PROP_ELL_MAX,
  PROP_F_SKY,
  PROP_SPECTRA,
  PROP_SIZE,
};

/* One spectrum of the data vector: its tracer pair and the rows it occupies. */
typedef struct _NcXcorCovGaussSpectrum
{
  guint a;
  guint b;
  guint offset;
  guint n_bp;
} NcXcorCovGaussSpectrum;

struct _NcXcorCovGauss
{
  /*< private > */
  GObject parent_instance;

  guint n_tracers;
  guint ell_min;
  guint ell_max;
  gdouble f_sky;

  GArray *spectra;    /* element-type NcXcorCovGaussSpectrum */
  GPtrArray *windows; /* element-type NcmMatrix*, n_bp x n_ell, owned */
  guint n_data;
  guint max_bp;

  /* Workspaces of nc_xcor_cov_gauss_compute(). */
  NcmMatrix *cl_tot;   /* signal plus noise, n_pairs x n_ell */
  NcmVector *g_ell;    /* 1 / ((2 ell + 1) f_sky) */
  NcmVector *ell_term; /* the per-ell covariance of one pair of spectra */
  NcmMatrix *scaled;   /* max_bp x n_ell */
};

G_DEFINE_TYPE (NcXcorCovGauss, nc_xcor_cov_gauss, G_TYPE_OBJECT)

static void
nc_xcor_cov_gauss_init (NcXcorCovGauss *cov)
{
  cov->n_tracers = 0;
  cov->ell_min   = 0;
  cov->ell_max   = 0;
  cov->f_sky     = 0.0;

  cov->spectra = NULL;
  cov->windows = NULL;
  cov->n_data  = 0;
  cov->max_bp  = 0;

  cov->cl_tot   = NULL;
  cov->g_ell    = NULL;
  cov->ell_term = NULL;
  cov->scaled   = NULL;
}

static void _nc_xcor_cov_gauss_set_spectra (NcXcorCovGauss *cov, GVariant *var);
static void _nc_xcor_cov_gauss_ell_term (NcXcorCovGauss *cov, const guint a, const guint b, const guint c, const guint d);
static guint _nc_xcor_cov_gauss_pair_row (const guint n, const guint a, const guint b);
static GVariant *_nc_xcor_cov_gauss_get_spectra (NcXcorCovGauss *cov);

static void
_nc_xcor_cov_gauss_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcXcorCovGauss *cov = NC_XCOR_COV_GAUSS (object);

  g_return_if_fail (NC_IS_XCOR_COV_GAUSS (object));

  switch (prop_id)
  {
    case PROP_N_TRACERS:
      cov->n_tracers = g_value_get_uint (value);
      break;
    case PROP_ELL_MIN:
      cov->ell_min = g_value_get_uint (value);
      break;
    case PROP_ELL_MAX:
      cov->ell_max = g_value_get_uint (value);
      break;
    case PROP_F_SKY:
      nc_xcor_cov_gauss_set_f_sky (cov, g_value_get_double (value));
      break;
    case PROP_SPECTRA:
      _nc_xcor_cov_gauss_set_spectra (cov, g_value_get_variant (value));
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_nc_xcor_cov_gauss_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcXcorCovGauss *cov = NC_XCOR_COV_GAUSS (object);

  g_return_if_fail (NC_IS_XCOR_COV_GAUSS (object));

  switch (prop_id)
  {
    case PROP_N_TRACERS:
      g_value_set_uint (value, cov->n_tracers);
      break;
    case PROP_ELL_MIN:
      g_value_set_uint (value, cov->ell_min);
      break;
    case PROP_ELL_MAX:
      g_value_set_uint (value, cov->ell_max);
      break;
    case PROP_F_SKY:
      g_value_set_double (value, cov->f_sky);
      break;
    case PROP_SPECTRA:
      g_value_take_variant (value, _nc_xcor_cov_gauss_get_spectra (cov));
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

/*
 * Allocates the spectrum list and the workspaces that depend only on the tracer count and
 * the multipole range. An instance built with no properties (the serialization registry
 * does this for every registered type) has no tracers and gets no workspaces;
 * nc_xcor_cov_gauss_add_spectrum() rejects it.
 */
static void
_nc_xcor_cov_gauss_constructed (GObject *object)
{
  NcXcorCovGauss *cov = NC_XCOR_COV_GAUSS (object);

  /* Chain up : start */
  G_OBJECT_CLASS (nc_xcor_cov_gauss_parent_class)->constructed (object);

  cov->spectra = g_array_new (FALSE, FALSE, sizeof (NcXcorCovGaussSpectrum));
  cov->windows = g_ptr_array_new_with_free_func ((GDestroyNotify) ncm_matrix_free);

  if (cov->n_tracers == 0)
    return;

  if (cov->ell_max < cov->ell_min)
    g_error ("_nc_xcor_cov_gauss_constructed: ell-max (%u) is smaller than ell-min (%u).",
             cov->ell_max, cov->ell_min);

  {
    const guint n_ell = nc_xcor_cov_gauss_get_n_ell (cov);

    cov->cl_tot   = ncm_matrix_new (nc_xcor_cov_gauss_get_n_pairs (cov), n_ell);
    cov->g_ell    = ncm_vector_new (n_ell);
    cov->ell_term = ncm_vector_new (n_ell);
  }
}

static void
_nc_xcor_cov_gauss_dispose (GObject *object)
{
  NcXcorCovGauss *cov = NC_XCOR_COV_GAUSS (object);

  g_clear_pointer (&cov->windows, g_ptr_array_unref);
  ncm_matrix_clear (&cov->cl_tot);
  ncm_vector_clear (&cov->g_ell);
  ncm_vector_clear (&cov->ell_term);
  ncm_matrix_clear (&cov->scaled);

  /* Chain up : end */
  G_OBJECT_CLASS (nc_xcor_cov_gauss_parent_class)->dispose (object);
}

static void
_nc_xcor_cov_gauss_finalize (GObject *object)
{
  NcXcorCovGauss *cov = NC_XCOR_COV_GAUSS (object);

  g_clear_pointer (&cov->spectra, g_array_unref);

  /* Chain up : end */
  G_OBJECT_CLASS (nc_xcor_cov_gauss_parent_class)->finalize (object);
}

static void
nc_xcor_cov_gauss_class_init (NcXcorCovGaussClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);

  object_class->set_property = &_nc_xcor_cov_gauss_set_property;
  object_class->get_property = &_nc_xcor_cov_gauss_get_property;
  object_class->constructed  = &_nc_xcor_cov_gauss_constructed;
  object_class->dispose      = &_nc_xcor_cov_gauss_dispose;
  object_class->finalize     = &_nc_xcor_cov_gauss_finalize;

  /**
   * NcXcorCovGauss:n-tracers:
   *
   * Number of tracers.
   */
  g_object_class_install_property (object_class,
                                   PROP_N_TRACERS,
                                   g_param_spec_uint ("n-tracers",
                                                      NULL,
                                                      "Number of tracers",
                                                      0, G_MAXUINT, 0,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcXcorCovGauss:ell-min:
   *
   * Smallest multipole of the common range.
   */
  g_object_class_install_property (object_class,
                                   PROP_ELL_MIN,
                                   g_param_spec_uint ("ell-min",
                                                      NULL,
                                                      "Smallest multipole",
                                                      0, G_MAXUINT, 0,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcXcorCovGauss:ell-max:
   *
   * Largest multipole of the common range.
   */
  g_object_class_install_property (object_class,
                                   PROP_ELL_MAX,
                                   g_param_spec_uint ("ell-max",
                                                      NULL,
                                                      "Largest multipole",
                                                      0, G_MAXUINT, 0,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcXcorCovGauss:f-sky:
   *
   * Observed fraction of the sky, in (0, 1].
   */
  g_object_class_install_property (object_class,
                                   PROP_F_SKY,
                                   g_param_spec_double ("f-sky",
                                                        NULL,
                                                        "Observed fraction of the sky",
                                                        G_MINDOUBLE, 1.0, 1.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcXcorCovGauss:spectra:
   *
   * The spectra of the data vector in order, each the tracer pair $(a, b)$ and the window
   * rows, one per band power over the common multipole range.
   */
  g_object_class_install_property (object_class,
                                   PROP_SPECTRA,
                                   g_param_spec_variant ("spectra",
                                                         NULL,
                                                         "Tracer pairs and windows of the spectra",
                                                         G_VARIANT_TYPE ("a(uuaad)"), NULL,
                                                         G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
}

/* Replaces the spectrum list by the one in @var, of type a(uuaad). */
static void
_nc_xcor_cov_gauss_set_spectra (NcXcorCovGauss *cov, GVariant *var)
{
  GVariantIter iter;
  GVariant *window_var;
  guint a, b;

  g_array_set_size (cov->spectra, 0);
  g_ptr_array_set_size (cov->windows, 0);
  cov->n_data = 0;

  if (var == NULL)
    return;

  g_variant_iter_init (&iter, var);

  while (g_variant_iter_next (&iter, "(uu@aad)", &a, &b, &window_var))
  {
    NcmMatrix *window = ncm_matrix_new_variant (window_var);

    nc_xcor_cov_gauss_add_spectrum (cov, a, b, window);

    ncm_matrix_free (window);
    g_variant_unref (window_var);
  }
}

/* Returns the spectrum list as a new #GVariant of type a(uuaad). */
static GVariant *
_nc_xcor_cov_gauss_get_spectra (NcXcorCovGauss *cov)
{
  GVariantBuilder builder;
  guint i;

  g_variant_builder_init (&builder, G_VARIANT_TYPE ("a(uuaad)"));

  for (i = 0; i < cov->spectra->len; i++)
  {
    const NcXcorCovGaussSpectrum *sp = &g_array_index (cov->spectra, NcXcorCovGaussSpectrum, i);
    GVariant *window_var             = ncm_matrix_get_variant (g_ptr_array_index (cov->windows, i));

    g_variant_builder_add (&builder, "(uu@aad)", sp->a, sp->b, window_var);
    g_variant_unref (window_var);
  }

  return g_variant_ref_sink (g_variant_builder_end (&builder));
}

/**
 * nc_xcor_cov_gauss_new:
 * @n_tracers: number of tracers
 * @ell_min: smallest multipole of the common range
 * @ell_max: largest multipole of the common range
 * @f_sky: observed fraction of the sky, in (0, 1]
 *
 * Creates a new #NcXcorCovGauss with no spectra.
 *
 * Returns: (transfer full): a new #NcXcorCovGauss.
 */
NcXcorCovGauss *
nc_xcor_cov_gauss_new (guint n_tracers, guint ell_min, guint ell_max, gdouble f_sky)
{
  return g_object_new (NC_TYPE_XCOR_COV_GAUSS,
                       "n-tracers", n_tracers,
                       "ell-min", ell_min,
                       "ell-max", ell_max,
                       "f-sky", f_sky,
                       NULL);
}

/**
 * nc_xcor_cov_gauss_ref:
 * @cov: a #NcXcorCovGauss
 *
 * Increases the reference count of @cov by one.
 *
 * Returns: (transfer full): @cov.
 */
NcXcorCovGauss *
nc_xcor_cov_gauss_ref (NcXcorCovGauss *cov)
{
  return g_object_ref (cov);
}

/**
 * nc_xcor_cov_gauss_free:
 * @cov: a #NcXcorCovGauss
 *
 * Decreases the reference count of @cov by one.
 */
void
nc_xcor_cov_gauss_free (NcXcorCovGauss *cov)
{
  g_object_unref (cov);
}

/**
 * nc_xcor_cov_gauss_clear:
 * @cov: a #NcXcorCovGauss
 *
 * Decreases the reference count of *@cov by one and sets *@cov to %NULL.
 */
void
nc_xcor_cov_gauss_clear (NcXcorCovGauss **cov)
{
  g_clear_object (cov);
}

/**
 * nc_xcor_cov_gauss_get_n_tracers:
 * @cov: a #NcXcorCovGauss
 *
 * Returns: the number of tracers.
 */
guint
nc_xcor_cov_gauss_get_n_tracers (NcXcorCovGauss *cov)
{
  return cov->n_tracers;
}

/**
 * nc_xcor_cov_gauss_get_ell_min:
 * @cov: a #NcXcorCovGauss
 *
 * Returns: the smallest multipole of the common range.
 */
guint
nc_xcor_cov_gauss_get_ell_min (NcXcorCovGauss *cov)
{
  return cov->ell_min;
}

/**
 * nc_xcor_cov_gauss_get_ell_max:
 * @cov: a #NcXcorCovGauss
 *
 * Returns: the largest multipole of the common range.
 */
guint
nc_xcor_cov_gauss_get_ell_max (NcXcorCovGauss *cov)
{
  return cov->ell_max;
}

/**
 * nc_xcor_cov_gauss_get_n_ell:
 * @cov: a #NcXcorCovGauss
 *
 * Returns: the number of multipoles of the common range, the number of columns of every
 * window, signal and noise row.
 */
guint
nc_xcor_cov_gauss_get_n_ell (NcXcorCovGauss *cov)
{
  return cov->ell_max - cov->ell_min + 1;
}

/**
 * nc_xcor_cov_gauss_set_f_sky:
 * @cov: a #NcXcorCovGauss
 * @f_sky: observed fraction of the sky
 *
 * Sets the observed fraction of the sky. Aborts unless 0 < @f_sky <= 1.
 */
void
nc_xcor_cov_gauss_set_f_sky (NcXcorCovGauss *cov, gdouble f_sky)
{
  if (!((f_sky > 0.0) && (f_sky <= 1.0)))
    g_error ("nc_xcor_cov_gauss_set_f_sky: f-sky must be in (0, 1], got %g.", f_sky);

  cov->f_sky = f_sky;
}

/**
 * nc_xcor_cov_gauss_get_f_sky:
 * @cov: a #NcXcorCovGauss
 *
 * Returns: the observed fraction of the sky.
 */
gdouble
nc_xcor_cov_gauss_get_f_sky (NcXcorCovGauss *cov)
{
  return cov->f_sky;
}

/**
 * nc_xcor_cov_gauss_get_n_pairs:
 * @cov: a #NcXcorCovGauss
 *
 * Returns: the number of tracer pairs, $n(n+1)/2$ for $n$ tracers, the number of rows of
 * the signal matrix of nc_xcor_cov_gauss_compute().
 */
guint
nc_xcor_cov_gauss_get_n_pairs (NcXcorCovGauss *cov)
{
  return cov->n_tracers * (cov->n_tracers + 1) / 2;
}

/**
 * nc_xcor_cov_gauss_get_pair_index:
 * @cov: a #NcXcorCovGauss
 * @a: first tracer
 * @b: second tracer
 *
 * Gives the row of the tracer pair (@a, @b) in the signal matrix of
 * nc_xcor_cov_gauss_compute(). The pairs (@a, @b) and (@b, @a) share a row; the rows run
 * over the upper triangle row by row: (0, 0), (0, 1), ..., (0, n - 1), (1, 1), ... Aborts
 * if a tracer index is out of range.
 *
 * Returns: the row of the pair.
 */
guint
nc_xcor_cov_gauss_get_pair_index (NcXcorCovGauss *cov, guint a, guint b)
{
  if ((a >= cov->n_tracers) || (b >= cov->n_tracers))
    g_error ("nc_xcor_cov_gauss_get_pair_index: tracer pair (%u, %u) out of range, %u tracers.",
             a, b, cov->n_tracers);

  return (a <= b) ? _nc_xcor_cov_gauss_pair_row (cov->n_tracers, a, b) : _nc_xcor_cov_gauss_pair_row (cov->n_tracers, b, a);
}

/**
 * nc_xcor_cov_gauss_add_spectrum:
 * @cov: a #NcXcorCovGauss
 * @a: first tracer
 * @b: second tracer
 * @window: window weights, one row per band power, one column per multipole
 *
 * Appends the spectrum of the tracer pair (@a, @b) to the data vector, with a copy of
 * @window. Its band powers take the next nc_xcor_cov_gauss_get_n_data() rows. Aborts if a
 * tracer index is out of range or @window does not have nc_xcor_cov_gauss_get_n_ell()
 * columns.
 */
void
nc_xcor_cov_gauss_add_spectrum (NcXcorCovGauss *cov, guint a, guint b, NcmMatrix *window)
{
  const guint n_ell = nc_xcor_cov_gauss_get_n_ell (cov);
  NcXcorCovGaussSpectrum sp;

  if ((a >= cov->n_tracers) || (b >= cov->n_tracers))
    g_error ("nc_xcor_cov_gauss_add_spectrum: tracer pair (%u, %u) out of range, %u tracers.",
             a, b, cov->n_tracers);

  if (ncm_matrix_ncols (window) != n_ell)
    g_error ("nc_xcor_cov_gauss_add_spectrum: the window has %u columns, the multipole range "
             "[%u, %u] has %u.", ncm_matrix_ncols (window), cov->ell_min, cov->ell_max, n_ell);

  sp.a      = a;
  sp.b      = b;
  sp.offset = cov->n_data;
  sp.n_bp   = ncm_matrix_nrows (window);

  g_array_append_val (cov->spectra, sp);
  g_ptr_array_add (cov->windows, ncm_matrix_dup (window));
  cov->n_data += sp.n_bp;

  if (sp.n_bp > cov->max_bp)
  {
    cov->max_bp = sp.n_bp;
    ncm_matrix_clear (&cov->scaled);
    cov->scaled = ncm_matrix_new (cov->max_bp, n_ell);
  }
}

/**
 * nc_xcor_cov_gauss_get_n_spectra:
 * @cov: a #NcXcorCovGauss
 *
 * Returns: the number of spectra of the data vector.
 */
guint
nc_xcor_cov_gauss_get_n_spectra (NcXcorCovGauss *cov)
{
  return cov->spectra->len;
}

/**
 * nc_xcor_cov_gauss_get_n_data:
 * @cov: a #NcXcorCovGauss
 *
 * Returns: the number of band powers of the data vector, the size of the covariance.
 */
guint
nc_xcor_cov_gauss_get_n_data (NcXcorCovGauss *cov)
{
  return cov->n_data;
}

/**
 * nc_xcor_cov_gauss_compute:
 * @cov: a #NcXcorCovGauss
 * @cl: signal $C^{ab}_\ell$, one row per tracer pair, one column per multipole
 * @noise: noise $N^a_\ell$, one row per tracer, one column per multipole
 * @cov_mat: output covariance, nc_xcor_cov_gauss_get_n_data() square
 *
 * Computes the covariance of the band powers of the data vector into @cov_mat. The rows of
 * @cl follow nc_xcor_cov_gauss_get_pair_index(); @noise is added to the auto spectra only.
 * Aborts if a matrix has the wrong shape.
 */
void
nc_xcor_cov_gauss_compute (NcXcorCovGauss *cov, NcmMatrix *cl, NcmMatrix *noise, NcmMatrix *cov_mat)
{
  const guint n_ell   = nc_xcor_cov_gauss_get_n_ell (cov);
  const guint n_pairs = nc_xcor_cov_gauss_get_n_pairs (cov);
  guint s, t, a, l;

  if ((ncm_matrix_nrows (cl) != n_pairs) || (ncm_matrix_ncols (cl) != n_ell))
    g_error ("nc_xcor_cov_gauss_compute: the signal matrix is %u x %u, expected %u x %u.",
             ncm_matrix_nrows (cl), ncm_matrix_ncols (cl), n_pairs, n_ell);

  if ((ncm_matrix_nrows (noise) != cov->n_tracers) || (ncm_matrix_ncols (noise) != n_ell))
    g_error ("nc_xcor_cov_gauss_compute: the noise matrix is %u x %u, expected %u x %u.",
             ncm_matrix_nrows (noise), ncm_matrix_ncols (noise), cov->n_tracers, n_ell);

  if ((ncm_matrix_nrows (cov_mat) != cov->n_data) || (ncm_matrix_ncols (cov_mat) != cov->n_data))
    g_error ("nc_xcor_cov_gauss_compute: the covariance matrix is %u x %u, expected %u x %u.",
             ncm_matrix_nrows (cov_mat), ncm_matrix_ncols (cov_mat), cov->n_data, cov->n_data);

  ncm_matrix_memcpy (cov->cl_tot, cl);

  for (a = 0; a < cov->n_tracers; a++)
  {
    const guint r_aa = _nc_xcor_cov_gauss_pair_row (cov->n_tracers, a, a);

    for (l = 0; l < n_ell; l++)
      ncm_matrix_addto (cov->cl_tot, r_aa, l, ncm_matrix_get (noise, a, l));
  }

  for (l = 0; l < n_ell; l++)
    ncm_vector_set (cov->g_ell, l, 1.0 / ((2.0 * (cov->ell_min + l) + 1.0) * cov->f_sky));

  for (s = 0; s < cov->spectra->len; s++)
  {
    const NcXcorCovGaussSpectrum *sp_s = &g_array_index (cov->spectra, NcXcorCovGaussSpectrum, s);
    NcmMatrix *window_s                = g_ptr_array_index (cov->windows, s);

    for (t = s; t < cov->spectra->len; t++)
    {
      const NcXcorCovGaussSpectrum *sp_t = &g_array_index (cov->spectra, NcXcorCovGaussSpectrum, t);
      NcmMatrix *window_t                = g_ptr_array_index (cov->windows, t);
      gsl_matrix_view scaled             = gsl_matrix_submatrix (ncm_matrix_gsl (cov->scaled), 0, 0, sp_t->n_bp, n_ell);
      gsl_matrix_view block              = gsl_matrix_submatrix (ncm_matrix_gsl (cov_mat), sp_s->offset, sp_t->offset, sp_s->n_bp, sp_t->n_bp);
      guint i, j;

      _nc_xcor_cov_gauss_ell_term (cov, sp_s->a, sp_s->b, sp_t->a, sp_t->b);

      for (i = 0; i < sp_t->n_bp; i++)
        for (l = 0; l < n_ell; l++)
          gsl_matrix_set (&scaled.matrix, i, l, ncm_matrix_get (window_t, i, l) * ncm_vector_get (cov->ell_term, l));

      gsl_blas_dgemm (CblasNoTrans, CblasTrans, 1.0, ncm_matrix_gsl (window_s), &scaled.matrix, 0.0, &block.matrix);

      if (t != s)
        for (i = 0; i < sp_s->n_bp; i++)
          for (j = 0; j < sp_t->n_bp; j++)
            ncm_matrix_set (cov_mat, sp_t->offset + j, sp_s->offset + i, gsl_matrix_get (&block.matrix, i, j));
    }
  }
}

/*
 * Fills cov->ell_term with the per-multipole covariance of the spectra (a, b) and (c, d),
 * g_ell (Ct^ac Ct^bd + Ct^ad Ct^bc), from cov->cl_tot and cov->g_ell. Summed with the two
 * windows it gives the band-power covariance.
 */
static void
_nc_xcor_cov_gauss_ell_term (NcXcorCovGauss *cov, const guint a, const guint b, const guint c, const guint d)
{
  const guint n_ell = nc_xcor_cov_gauss_get_n_ell (cov);
  const guint r_ac  = nc_xcor_cov_gauss_get_pair_index (cov, a, c);
  const guint r_bd  = nc_xcor_cov_gauss_get_pair_index (cov, b, d);
  const guint r_ad  = nc_xcor_cov_gauss_get_pair_index (cov, a, d);
  const guint r_bc  = nc_xcor_cov_gauss_get_pair_index (cov, b, c);
  guint l;

  for (l = 0; l < n_ell; l++)
  {
    const gdouble ac = ncm_matrix_get (cov->cl_tot, r_ac, l);
    const gdouble bd = ncm_matrix_get (cov->cl_tot, r_bd, l);
    const gdouble ad = ncm_matrix_get (cov->cl_tot, r_ad, l);
    const gdouble bc = ncm_matrix_get (cov->cl_tot, r_bc, l);

    ncm_vector_set (cov->ell_term, l, ncm_vector_get (cov->g_ell, l) * (ac * bd + ad * bc));
  }
}

/* Row of the pair (a, b), a <= b, in the upper triangle taken row by row. */
static guint
_nc_xcor_cov_gauss_pair_row (const guint n, const guint a, const guint b)
{
  return a * n - a * (a - 1) / 2 + (b - a);
}

