/***************************************************************************
 *            nc_hicosmo_de_wpoly_library.c
 *
 *  Mon Sep 28 16:00:00 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2026 <vitenti@uel.br>
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
 * NcHICosmoDEWPolyLibrary:
 *
 * Library of fiducial $w(z)$ models grouped by their number of turning points.
 *
 * Draws $\mathrm{d}w/\mathrm{d}z$ as a degree-$d$ Legendre series on $[0, z_f]$ with
 * independent standard Gaussian coefficients in the orthonormal basis, from an
 * #NcmRNG seeded with #NcHICosmoDEWPolyLibrary:seed, and integrates it to a zero-mean
 * $w + 1$. Each draw is counted by its number of turning points $m$
 * (nc_hicosmo_de_wpoly_turning_points()); draws with $m$ above
 * #NcHICosmoDEWPolyLibrary:max-turning-points are discarded. For $m \ge 1$ a draw is
 * accepted with the probability that makes the ordered turning points uniform on
 * $[0, z_f]$: the ratio of the target density to the draw density, the latter estimated
 * on a first pass by a histogram with #NcHICosmoDEWPolyLibrary:density-bins bins per
 * axis and #NcHICosmoDEWPolyLibrary:density-draws-per-cell draws per cell on average.
 * For every $m$ from $0$ to the maximum the library collects
 * #NcHICosmoDEWPolyLibrary:n-fiducials fiducials. An accepted $w + 1$ is rescaled so that
 * its largest excursion between consecutive stationary points ($0$, the turning points,
 * $z_f$) equals an amplitude drawn uniformly on $[0, $ #NcHICosmoDEWPolyLibrary:a-max $]$,
 * and the fiducials with the same $m$ are shuffled. Acceptance, amplitude and shuffle use
 * a second #NcmRNG seeded with seed + 1, so the coefficient stream is the same in both
 * passes and a fiducial is fixed by its draw index.
 *
 * Accessors prepare the library on first use.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "nc/background/nc_hicosmo_de_wpoly_library.h"
#include "ncm/core/ncm_rng.h"
#include "ncm/model/ncm_model.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_math.h>
#endif /* NUMCOSMO_GIR_SCAN */

typedef struct _NcHICosmoDEWPolyLibraryFiducial
{
  gulong draw;
  gdouble amplitude;
  GArray *turning_points;
  NcmVector *coeffs;
} NcHICosmoDEWPolyLibraryFiducial;

struct _NcHICosmoDEWPolyLibrary
{
  /*< private >*/
  GObject parent_instance;
  gulong seed;
  guint degree;
  gdouble z_f;
  guint n_fiducials;
  guint max_tp;
  gdouble a_max;
  guint density_bins;
  guint density_draws;
  gboolean prepared;
  gulong n_draws;
  NcHICosmoDEWPoly *work;
  GPtrArray *by_tp; /* GArray of NcHICosmoDEWPolyLibraryFiducial per number of turning points */
};

G_DEFINE_TYPE (NcHICosmoDEWPolyLibrary, nc_hicosmo_de_wpoly_library, G_TYPE_OBJECT)

enum
{
  PROP_0,
  PROP_SEED,
  PROP_DEGREE,
  PROP_Z_F,
  PROP_N_FIDUCIALS,
  PROP_MAX_TP,
  PROP_A_MAX,
  PROP_DENSITY_BINS,
  PROP_DENSITY_DRAWS,
};

static void
_nc_hicosmo_de_wpoly_library_fiducial_clear (gpointer data)
{
  NcHICosmoDEWPolyLibraryFiducial *fid = data;

  g_clear_pointer (&fid->turning_points, g_array_unref);
  ncm_vector_clear (&fid->coeffs);
}

static void
nc_hicosmo_de_wpoly_library_init (NcHICosmoDEWPolyLibrary *lib)
{
  lib->seed          = 0;
  lib->degree        = 0;
  lib->z_f           = 0.0;
  lib->n_fiducials   = 0;
  lib->max_tp        = 0;
  lib->a_max         = 0.0;
  lib->density_bins  = 0;
  lib->density_draws = 0;
  lib->prepared      = FALSE;
  lib->n_draws       = 0;
  lib->work          = NULL;
  lib->by_tp         = g_ptr_array_new_with_free_func ((GDestroyNotify) g_array_unref);
}

static void
_nc_hicosmo_de_wpoly_library_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcHICosmoDEWPolyLibrary *lib = NC_HICOSMO_DE_WPOLY_LIBRARY (object);

  switch (prop_id)
  {
    case PROP_SEED:
      lib->seed = g_value_get_ulong (value);
      break;
    case PROP_DEGREE:
      lib->degree = g_value_get_uint (value);
      break;
    case PROP_Z_F:
      lib->z_f = g_value_get_double (value);
      break;
    case PROP_N_FIDUCIALS:
      lib->n_fiducials = g_value_get_uint (value);
      break;
    case PROP_MAX_TP:
      lib->max_tp = g_value_get_uint (value);
      break;
    case PROP_A_MAX:
      lib->a_max = g_value_get_double (value);
      break;
    case PROP_DENSITY_BINS:
      lib->density_bins = g_value_get_uint (value);
      break;
    case PROP_DENSITY_DRAWS:
      lib->density_draws = g_value_get_uint (value);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_nc_hicosmo_de_wpoly_library_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcHICosmoDEWPolyLibrary *lib = NC_HICOSMO_DE_WPOLY_LIBRARY (object);

  switch (prop_id)
  {
    case PROP_SEED:
      g_value_set_ulong (value, lib->seed);
      break;
    case PROP_DEGREE:
      g_value_set_uint (value, lib->degree);
      break;
    case PROP_Z_F:
      g_value_set_double (value, lib->z_f);
      break;
    case PROP_N_FIDUCIALS:
      g_value_set_uint (value, lib->n_fiducials);
      break;
    case PROP_MAX_TP:
      g_value_set_uint (value, lib->max_tp);
      break;
    case PROP_A_MAX:
      g_value_set_double (value, lib->a_max);
      break;
    case PROP_DENSITY_BINS:
      g_value_set_uint (value, lib->density_bins);
      break;
    case PROP_DENSITY_DRAWS:
      g_value_set_uint (value, lib->density_draws);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_nc_hicosmo_de_wpoly_library_constructed (GObject *object)
{
  /* Chain up : start */
  G_OBJECT_CLASS (nc_hicosmo_de_wpoly_library_parent_class)->constructed (object);
  {
    NcHICosmoDEWPolyLibrary *lib = NC_HICOSMO_DE_WPOLY_LIBRARY (object);
    guint m;

    g_assert_cmpuint (lib->max_tp, <=, 4);

    lib->work = nc_hicosmo_de_wpoly_new (lib->degree + 2, lib->z_f);

    for (m = 0; m <= lib->max_tp; m++)
    {
      GArray *fids = g_array_new (FALSE, FALSE, sizeof (NcHICosmoDEWPolyLibraryFiducial));

      g_array_set_clear_func (fids, _nc_hicosmo_de_wpoly_library_fiducial_clear);
      g_ptr_array_add (lib->by_tp, fids);
    }
  }
}

static void
_nc_hicosmo_de_wpoly_library_dispose (GObject *object)
{
  NcHICosmoDEWPolyLibrary *lib = NC_HICOSMO_DE_WPOLY_LIBRARY (object);

  g_clear_object (&lib->work);
  g_clear_pointer (&lib->by_tp, g_ptr_array_unref);

  /* Chain up : end */
  G_OBJECT_CLASS (nc_hicosmo_de_wpoly_library_parent_class)->dispose (object);
}

static void
nc_hicosmo_de_wpoly_library_class_init (NcHICosmoDEWPolyLibraryClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);

  object_class->set_property = &_nc_hicosmo_de_wpoly_library_set_property;
  object_class->get_property = &_nc_hicosmo_de_wpoly_library_get_property;
  object_class->constructed  = &_nc_hicosmo_de_wpoly_library_constructed;
  object_class->dispose      = &_nc_hicosmo_de_wpoly_library_dispose;

  g_object_class_install_property (object_class,
                                   PROP_SEED,
                                   g_param_spec_ulong ("seed",
                                                       NULL,
                                                       "master seed of the coefficient stream",
                                                       0, G_MAXULONG, 0,
                                                       G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_DEGREE,
                                   g_param_spec_uint ("degree",
                                                      NULL,
                                                      "degree of dw/dz",
                                                      1, G_MAXUINT, 6,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_Z_F,
                                   g_param_spec_double ("zf",
                                                        NULL,
                                                        "final redshift",
                                                        1.0e-3, 1.0e10, 2.33,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_N_FIDUCIALS,
                                   g_param_spec_uint ("n-fiducials",
                                                      NULL,
                                                      "fiducials per number of turning points",
                                                      1, G_MAXUINT, 1000,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_MAX_TP,
                                   g_param_spec_uint ("max-turning-points",
                                                      NULL,
                                                      "largest number of turning points kept",
                                                      0, 4, 3,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_A_MAX,
                                   g_param_spec_double ("a-max",
                                                        NULL,
                                                        "upper end of the uniform amplitude",
                                                        0.0, 10.0, 0.5,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_DENSITY_BINS,
                                   g_param_spec_uint ("density-bins",
                                                      NULL,
                                                      "bins per axis of the turning-point density estimate",
                                                      1, 255, 10,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_DENSITY_DRAWS,
                                   g_param_spec_uint ("density-draws-per-cell",
                                                      NULL,
                                                      "first-pass draws per cell of the density estimate, on average",
                                                      1, G_MAXUINT, 50,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
}

/* Number of ordered cells of the histogram for @m turning points: C(bins + m - 1, m). */
static gulong
_nc_hicosmo_de_wpoly_library_n_cells (NcHICosmoDEWPolyLibrary *lib, guint m)
{
  gulong n = 1;
  guint i;

  for (i = 1; i <= m; i++)
    n = n * (lib->density_bins + i - 1) / i;

  return n;
}

/*
 * One draw of the coefficient stream: dw/dx as a Legendre series with
 * independent N(0, 1) coefficients in the orthonormal basis, integrated to
 * w + 1 with zero mean (a_0 = 0), written into @coeffs (length degree + 2).
 */
static void
_nc_hicosmo_de_wpoly_library_draw (NcHICosmoDEWPolyLibrary *lib, NcmRNG *rng, NcmVector *coeffs)
{
  const guint d = lib->degree;
  gdouble *c    = g_newa (gdouble, d + 1);
  guint l;

  for (l = 0; l <= d; l++)
    c[l] = ncm_rng_ugaussian_gen (rng) * sqrt (2.0 * l + 1.0);

  /* int P_l dx = (P_{l+1} - P_{l-1}) / (2 l + 1); the factor z_f / 2 from
   * dz = (z_f / 2) dx is absorbed by the amplitude normalization.
   */
  ncm_vector_set_zero (coeffs);

  for (l = 0; l <= d; l++)
  {
    const gdouble v = c[l] / (2.0 * l + 1.0);

    ncm_vector_addto (coeffs, l + 1, v);

    if (l >= 1)
      ncm_vector_addto (coeffs, l - 1, -v);
  }

  ncm_vector_set (coeffs, 0, 0.0);
}

static void
_nc_hicosmo_de_wpoly_library_set_work (NcHICosmoDEWPolyLibrary *lib, NcmVector *coeffs)
{
  ncm_model_orig_vparam_set_vector (NCM_MODEL (lib->work), NC_HICOSMO_DE_WPOLY_A, coeffs);
}

/* Largest |w(s_{i+1}) - w(s_i)| over the stationary set {0, z_1, ..., z_m, z_f}. */
static gdouble
_nc_hicosmo_de_wpoly_library_excursion (NcHICosmoDEWPolyLibrary *lib, GArray *tp)
{
  NcHICosmoDE *cosmo_de = NC_HICOSMO_DE (lib->work);
  gdouble w_prev        = nc_hicosmo_de_w_de (cosmo_de, 0.0);
  gdouble exc           = 0.0;
  guint i;

  for (i = 0; i <= tp->len; i++)
  {
    const gdouble z = (i < tp->len) ? g_array_index (tp, gdouble, i) : lib->z_f;
    const gdouble w = nc_hicosmo_de_w_de (cosmo_de, z);

    exc    = GSL_MAX (exc, fabs (w - w_prev));
    w_prev = w;
  }

  return exc;
}

/* Packed histogram key of the ordered turning-point bins, 8 bits per axis. */
static guint64
_nc_hicosmo_de_wpoly_library_cell (NcHICosmoDEWPolyLibrary *lib, GArray *tp)
{
  guint64 key = 0;
  guint i;

  for (i = 0; i < tp->len; i++)
  {
    const gdouble z = g_array_index (tp, gdouble, i);
    guint64 b       = (guint64) floor (lib->density_bins * z / lib->z_f);

    b    = MIN (b, (guint64) lib->density_bins - 1);
    key |= b << (8 * i);
  }

  return key;
}

/* Target density of a cell relative to a cell with distinct bins: 1 / prod k_j! over
 * the multiplicities k_j of repeated bins, the volume of the ordered cell.
 */
static gdouble
_nc_hicosmo_de_wpoly_library_cell_target (guint64 key, guint m)
{
  gdouble w = 1.0;
  guint i;

  for (i = 1; i < m; i++)
  {
    guint k = 1;
    guint j;

    for (j = 0; j < i; j++)
      if (((key >> (8 * j)) & 0xff) == ((key >> (8 * i)) & 0xff))
        k++;

    w /= k;
  }

  return w;
}

/**
 * nc_hicosmo_de_wpoly_library_prepare:
 * @lib: a #NcHICosmoDEWPolyLibrary
 *
 * Generates the library. A first pass over the coefficient stream estimates the
 * turning-point density for each $m \ge 1$ until it has
 * #NcHICosmoDEWPolyLibrary:density-draws-per-cell draws per histogram cell on
 * average; a second pass over the same stream accepts draws to the uniform target
 * until every $m$ has its fiducials, then
 * normalizes, draws amplitudes and shuffles. Does nothing when already prepared.
 */
void
nc_hicosmo_de_wpoly_library_prepare (NcHICosmoDEWPolyLibrary *lib)
{
  const guint nc     = lib->degree + 2;
  NcmVector *coeffs  = ncm_vector_new (nc);
  GHashTable **hist  = g_new0 (GHashTable *, lib->max_tp + 1);
  gdouble *max_ratio = g_new0 (gdouble, lib->max_tp + 1);
  gulong *n_draws_tp = g_new0 (gulong, lib->max_tp + 1);
  guint m;

  if (lib->prepared)
    return;

  for (m = 0; m <= lib->max_tp; m++)
    hist[m] = g_hash_table_new (g_int64_hash, g_int64_equal);

  /* Pass 1: histogram of the ordered turning-point bins for each m. */
  {
    NcmRNG *rng = ncm_rng_seeded_new (NULL, lib->seed);
    gboolean done;

    do {
      GArray *tp;

      _nc_hicosmo_de_wpoly_library_draw (lib, rng, coeffs);
      _nc_hicosmo_de_wpoly_library_set_work (lib, coeffs);
      tp = nc_hicosmo_de_wpoly_turning_points (lib->work);
      m  = tp->len;

      if ((m >= 1) && (m <= lib->max_tp))
      {
        guint64 *key = g_new (guint64, 1);
        gpointer val;

        *key = _nc_hicosmo_de_wpoly_library_cell (lib, tp);

        if (g_hash_table_lookup_extended (hist[m], key, NULL, &val))
        {
          /* the table keeps its own key */
          g_hash_table_insert (hist[m], key, GSIZE_TO_POINTER (GPOINTER_TO_SIZE (val) + 1));
          g_free (key);
        }
        else
        {
          g_hash_table_insert (hist[m], key, GSIZE_TO_POINTER (1));
        }

        n_draws_tp[m]++;
      }

      g_array_unref (tp);

      done = TRUE;

      for (m = 1; m <= lib->max_tp; m++)
        done = done && (n_draws_tp[m] >= lib->density_draws * _nc_hicosmo_de_wpoly_library_n_cells (lib, m));
    } while (!done);

    ncm_rng_free (rng);
  }

  /* Acceptance scale for each m: the largest target / density ratio over populated cells. */
  for (m = 1; m <= lib->max_tp; m++)
  {
    GHashTableIter iter;
    gpointer key, val;

    g_hash_table_iter_init (&iter, hist[m]);

    while (g_hash_table_iter_next (&iter, &key, &val))
    {
      const gdouble ratio = _nc_hicosmo_de_wpoly_library_cell_target (*(guint64 *) key, m) * n_draws_tp[m] / GPOINTER_TO_SIZE (val);

      max_ratio[m] = GSL_MAX (max_ratio[m], ratio);
    }
  }

  /* Pass 2: accept to the uniform target, normalize, draw the amplitude. */
  {
    NcmRNG *rng     = ncm_rng_seeded_new (NULL, lib->seed);
    NcmRNG *rng_acc = ncm_rng_seeded_new (NULL, lib->seed + 1);
    gulong draw     = 0;
    gboolean done;

    do {
      GArray *tp;
      gboolean accept;

      _nc_hicosmo_de_wpoly_library_draw (lib, rng, coeffs);
      _nc_hicosmo_de_wpoly_library_set_work (lib, coeffs);
      tp = nc_hicosmo_de_wpoly_turning_points (lib->work);
      m  = tp->len;
      draw++;

      if ((m > lib->max_tp) || (((GArray *) g_ptr_array_index (lib->by_tp, m))->len >= lib->n_fiducials))
      {
        accept = FALSE;
      }
      else if (m == 0)
      {
        accept = TRUE;
      }
      else
      {
        const guint64 key = _nc_hicosmo_de_wpoly_library_cell (lib, tp);
        gpointer val      = g_hash_table_lookup (hist[m], &key);
        const gdouble q   = (val == NULL) ? 1.0 : (gdouble) GPOINTER_TO_SIZE (val);
        const gdouble p   = _nc_hicosmo_de_wpoly_library_cell_target (key, m) * n_draws_tp[m] / q / max_ratio[m];

        accept = ncm_rng_uniform01_gen (rng_acc) < p;
      }

      if (accept)
      {
        NcHICosmoDEWPolyLibraryFiducial fid;
        const gdouble exc = _nc_hicosmo_de_wpoly_library_excursion (lib, tp);
        const gdouble A   = ncm_rng_uniform_gen (rng_acc, 0.0, lib->a_max);

        g_assert_cmpfloat (exc, >, 0.0);

        fid.draw           = draw;
        fid.amplitude      = A;
        fid.turning_points = tp;
        fid.coeffs         = ncm_vector_dup (coeffs);
        ncm_vector_scale (fid.coeffs, A / exc);

        g_array_append_val ((GArray *) g_ptr_array_index (lib->by_tp, m), fid);
      }
      else
      {
        g_array_unref (tp);
      }

      done = TRUE;

      for (m = 0; m <= lib->max_tp; m++)
        done = done && (((GArray *) g_ptr_array_index (lib->by_tp, m))->len >= lib->n_fiducials);
    } while (!done);

    lib->n_draws = draw;

    /* Shuffle the fiducials with the same m. */
    for (m = 0; m <= lib->max_tp; m++)
    {
      GArray *fids = g_ptr_array_index (lib->by_tp, m);
      guint i;

      for (i = fids->len - 1; i > 0; i--)
      {
        const guint j                       = ncm_rng_uniform_int_gen (rng_acc, i + 1);
        NcHICosmoDEWPolyLibraryFiducial tmp = g_array_index (fids, NcHICosmoDEWPolyLibraryFiducial, i);

        g_array_index (fids, NcHICosmoDEWPolyLibraryFiducial, i) = g_array_index (fids, NcHICosmoDEWPolyLibraryFiducial, j);
        g_array_index (fids, NcHICosmoDEWPolyLibraryFiducial, j) = tmp;
      }
    }

    ncm_rng_free (rng);
    ncm_rng_free (rng_acc);
  }

  for (m = 0; m <= lib->max_tp; m++)
  {
    GHashTableIter iter;
    gpointer key, val;

    g_hash_table_iter_init (&iter, hist[m]);

    while (g_hash_table_iter_next (&iter, &key, &val))
      g_free (key);

    g_hash_table_unref (hist[m]);
  }

  g_free (hist);
  g_free (max_ratio);
  g_free (n_draws_tp);
  ncm_vector_free (coeffs);

  lib->prepared = TRUE;
}

static NcHICosmoDEWPolyLibraryFiducial *
_nc_hicosmo_de_wpoly_library_peek (NcHICosmoDEWPolyLibrary *lib, guint n_tp, guint i)
{
  GArray *fids;

  nc_hicosmo_de_wpoly_library_prepare (lib);

  g_assert_cmpuint (n_tp, <=, lib->max_tp);
  fids = g_ptr_array_index (lib->by_tp, n_tp);
  g_assert_cmpuint (i, <, fids->len);

  return &g_array_index (fids, NcHICosmoDEWPolyLibraryFiducial, i);
}

/**
 * nc_hicosmo_de_wpoly_library_new:
 * @seed: master seed
 * @degree: degree $d$ of $\mathrm{d}w/\mathrm{d}z$
 * @z_f: final redshift $z_f$
 * @n_fiducials: fiducials per number of turning points
 *
 * Creates a library with the default #NcHICosmoDEWPolyLibrary:max-turning-points,
 * #NcHICosmoDEWPolyLibrary:a-max and #NcHICosmoDEWPolyLibrary:density-bins.
 *
 * Returns: a new #NcHICosmoDEWPolyLibrary
 */
NcHICosmoDEWPolyLibrary *
nc_hicosmo_de_wpoly_library_new (gulong seed, guint degree, const gdouble z_f, guint n_fiducials)
{
  return g_object_new (NC_TYPE_HICOSMO_DE_WPOLY_LIBRARY,
                       "seed", seed,
                       "degree", degree,
                       "zf", z_f,
                       "n-fiducials", n_fiducials,
                       NULL);
}

/**
 * nc_hicosmo_de_wpoly_library_ref:
 * @lib: a #NcHICosmoDEWPolyLibrary
 *
 * Returns: (transfer full): @lib
 */
NcHICosmoDEWPolyLibrary *
nc_hicosmo_de_wpoly_library_ref (NcHICosmoDEWPolyLibrary *lib)
{
  return g_object_ref (lib);
}

/**
 * nc_hicosmo_de_wpoly_library_free:
 * @lib: a #NcHICosmoDEWPolyLibrary
 *
 * Decreases the reference count of @lib by one.
 */
void
nc_hicosmo_de_wpoly_library_free (NcHICosmoDEWPolyLibrary *lib)
{
  g_object_unref (lib);
}

/**
 * nc_hicosmo_de_wpoly_library_clear:
 * @lib: a #NcHICosmoDEWPolyLibrary
 *
 * Decreases the reference count of *@lib by one and sets *@lib to %NULL.
 */
void
nc_hicosmo_de_wpoly_library_clear (NcHICosmoDEWPolyLibrary **lib)
{
  g_clear_object (lib);
}

/**
 * nc_hicosmo_de_wpoly_library_get_max_turning_points:
 * @lib: a #NcHICosmoDEWPolyLibrary
 *
 * Returns: the largest number of turning points kept
 */
guint
nc_hicosmo_de_wpoly_library_get_max_turning_points (NcHICosmoDEWPolyLibrary *lib)
{
  return lib->max_tp;
}

/**
 * nc_hicosmo_de_wpoly_library_len:
 * @lib: a #NcHICosmoDEWPolyLibrary
 * @n_tp: number of turning points
 *
 * Returns: the number of fiducials with @n_tp turning points
 */
guint
nc_hicosmo_de_wpoly_library_len (NcHICosmoDEWPolyLibrary *lib, guint n_tp)
{
  nc_hicosmo_de_wpoly_library_prepare (lib);
  g_assert_cmpuint (n_tp, <=, lib->max_tp);

  return ((GArray *) g_ptr_array_index (lib->by_tp, n_tp))->len;
}

/**
 * nc_hicosmo_de_wpoly_library_get_n_draws:
 * @lib: a #NcHICosmoDEWPolyLibrary
 *
 * Returns: the number of draws of the coefficient stream the second pass used
 */
gulong
nc_hicosmo_de_wpoly_library_get_n_draws (NcHICosmoDEWPolyLibrary *lib)
{
  nc_hicosmo_de_wpoly_library_prepare (lib);

  return lib->n_draws;
}

/**
 * nc_hicosmo_de_wpoly_library_get_draw_index:
 * @lib: a #NcHICosmoDEWPolyLibrary
 * @n_tp: number of turning points
 * @i: index among the fiducials with @n_tp turning points
 *
 * Returns: the 1-based index of the fiducial in the coefficient stream
 */
gulong
nc_hicosmo_de_wpoly_library_get_draw_index (NcHICosmoDEWPolyLibrary *lib, guint n_tp, guint i)
{
  return _nc_hicosmo_de_wpoly_library_peek (lib, n_tp, i)->draw;
}

/**
 * nc_hicosmo_de_wpoly_library_get_amplitude:
 * @lib: a #NcHICosmoDEWPolyLibrary
 * @n_tp: number of turning points
 * @i: index among the fiducials with @n_tp turning points
 *
 * Returns: the largest excursion of $w$ between consecutive stationary points
 */
gdouble
nc_hicosmo_de_wpoly_library_get_amplitude (NcHICosmoDEWPolyLibrary *lib, guint n_tp, guint i)
{
  return _nc_hicosmo_de_wpoly_library_peek (lib, n_tp, i)->amplitude;
}

/**
 * nc_hicosmo_de_wpoly_library_peek_turning_points:
 * @lib: a #NcHICosmoDEWPolyLibrary
 * @n_tp: number of turning points
 * @i: index among the fiducials with @n_tp turning points
 *
 * Returns: (element-type gdouble) (transfer none): the turning points in increasing $z$
 */
GArray *
nc_hicosmo_de_wpoly_library_peek_turning_points (NcHICosmoDEWPolyLibrary *lib, guint n_tp, guint i)
{
  return _nc_hicosmo_de_wpoly_library_peek (lib, n_tp, i)->turning_points;
}

/**
 * nc_hicosmo_de_wpoly_library_peek_coeffs:
 * @lib: a #NcHICosmoDEWPolyLibrary
 * @n_tp: number of turning points
 * @i: index among the fiducials with @n_tp turning points
 *
 * Returns: (transfer none): the Legendre coefficients of $w + 1$, $a_0 = 0$
 */
NcmVector *
nc_hicosmo_de_wpoly_library_peek_coeffs (NcHICosmoDEWPolyLibrary *lib, guint n_tp, guint i)
{
  return _nc_hicosmo_de_wpoly_library_peek (lib, n_tp, i)->coeffs;
}

/**
 * nc_hicosmo_de_wpoly_library_get_fiducial:
 * @lib: a #NcHICosmoDEWPolyLibrary
 * @n_tp: number of turning points
 * @i: index among the fiducials with @n_tp turning points
 *
 * Creates a #NcHICosmoDEWPoly carrying the fiducial, with the other cosmological
 * parameters at their defaults.
 *
 * Returns: (transfer full): a new #NcHICosmoDEWPoly
 */
NcHICosmoDEWPoly *
nc_hicosmo_de_wpoly_library_get_fiducial (NcHICosmoDEWPolyLibrary *lib, guint n_tp, guint i)
{
  NcHICosmoDEWPolyLibraryFiducial *fid = _nc_hicosmo_de_wpoly_library_peek (lib, n_tp, i);
  NcHICosmoDEWPoly *cosmo              = nc_hicosmo_de_wpoly_new (lib->degree + 2, lib->z_f);

  ncm_model_orig_vparam_set_vector (NCM_MODEL (cosmo), NC_HICOSMO_DE_WPOLY_A, fid->coeffs);

  return cosmo;
}

