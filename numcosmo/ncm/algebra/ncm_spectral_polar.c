/***************************************************************************
 *            ncm_spectral_polar.c
 *
 *  Wed Nov 11 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 *  Copyright  2026  Caio Lima de Oliveira
 *  <caiolimadeoliveira@pm.me>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2026 <vitenti@uel.br>
 * Copyright (C) 2026 Caio Lima de Oliveira <caiolimadeoliveira@pm.me>
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
 * NcmSpectralPolar:
 *
 * Adaptive integration on a disk in polar coordinates, with power laws at the radial ends.
 *
 * ncm_spectral_polar_integrate() computes
 * $$I = \int\!\!\int r^a (R - r)^b \, G (r, \theta) \, \mathrm{d}r \, \mathrm{d}\theta$$
 * for $a, b > -1$ on a domain given by its edges: the radial edges
 * $0 \le r_0 < r_1 < \dots < r_{n_r} \le R$ and the angular edges
 * $\theta_0 < \theta_1 < \dots < \theta_{n_\theta}$ split it into the $n_r n_\theta$ panels
 * $[r_i, r_{i+1}] \times [\theta_j, \theta_{j+1}]$, so an edge is a value of $r$ or of
 * $\theta$ where one panel ends and the next starts. The power laws are in the weight, so
 * $G$ is the factor that is smooth at $r = 0$ and $r = R$; the polar Jacobian is part of
 * $a$.
 *
 * On a panel, $F$ is $G$ times the part of the weight that is regular there: $r^a$ unless
 * the panel starts at $r = 0$, $(R - r)^b$ unless it ends at $r = R$. $F$ is sampled on a
 * tensor grid of Chebyshev-Lobatto points, $(2^{l_r} + 1) \times (2^{l_\theta} + 1)$, from
 * #NcmSpectralPolar:level-min in both directions. With $c^{kl}$ its Chebyshev coefficients,
 * $D_r = 2 \sum_{k > 2^{l_r - 1}} \sum_l |c^{kl}|$ and
 * $D_\theta = 2 \sum_{l > 2^{l_\theta - 1}} \sum_k |c^{kl}|$ are the $\ell_1$ changes of the
 * coefficients from the grids with every other point in $r$ and in $\theta$; each bounds
 * the largest change of the interpolant on the panel. The panel is accepted when
 * $$D_r + D_\theta \le \max (\varepsilon_\mathrm{rel} V, \varepsilon_F), \qquad
 * \varepsilon_F = \max (\varepsilon_\mathrm{rel} |I|, \varepsilon_\mathrm{abs}) / A,$$
 * with $V$ the largest $|F|$ on the grid of the panel, $\varepsilon_\mathrm{rel}$ the
 * relative and $\varepsilon_\mathrm{abs}$ the absolute tolerance, and $A$ the integral over
 * all panels of the weight left out of $F$. Otherwise the direction with the larger change
 * doubles, keeping its samples. At #NcmSpectralPolar:level-cap, or one level below when the
 * doubling is predicted to fail from the decay of the coefficients, the panel is bisected
 * in that direction. The value of a panel is the sum of its coefficients against the
 * moments of the weight left out of $F$.
 *
 * $|I|$ in $\varepsilon_F$ is first the sum of the panel values on the starting grids; when
 * the final sum is smaller, the accepted panels are tested again with it and those that
 * fail are expanded again.
 *
 * ncm_spectral_polar_peak_edges() gives the edges for an integrand with one peak. The
 * position and widths given to it change the number of evaluations, not the test.
 *
 * A call takes its scratch space from a pool, so one object may be used by calls in
 * several threads; its properties must not change during a call.
 *
 * See <a href="../../theory/ncm/algebra/spectral_polar.html">Spectral Integration on a
 * Disk</a> for the formulas.
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "ncm/algebra/ncm_spectral_polar.h"
#include "ncm/algebra/ncm_spectral_polar_private.h"
#include "ncm/algebra/ncm_spectral.h"
#include "ncm/algebra/ncm_spectral_private.h"
#include "ncm/core/ncm_cfg.h"
#include "ncm/core/ncm_memory_pool.h"
#include "ncm/core/ncm_util.h"
#include "ncm/core/ncm_vmath_private.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <math.h>
#include <string.h>
#include <fftw3.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_gamma.h>
#endif /* NUMCOSMO_GIR_SCAN */

/* Number of equal angular panels of ncm_spectral_polar_peak_edges() on the full circle */
#define NCM_SPECTRAL_POLAR_FULL_CIRCLE_PANELS (4)

/* Offset of the level-l table, of 2^l + 1 entries, in an array holding levels 0, 1, ... */
#define NCM_SPECTRAL_POLAR_LEVEL_OFFSET(l) ((1u << (l)) - 1 + (l))

enum
{
  PROP_0,
  PROP_LEVEL_MIN,
  PROP_LEVEL_CAP,
  PROP_MAX_DEPTH,
  PROP_EDGE_RATIO,
  PROP_MIN_RADIAL_PANELS,
};

struct _NcmSpectralPolar
{
  /*< private >*/
  GObject parent_instance;
  guint level_min;
  guint level_cap;
  guint max_depth;
  gdouble edge_ratio;
  guint min_radial_panels;
  NcmMemoryPool *work_pool; /* NcmSpectralPolarWork, one per concurrent call */
};

/* A panel and the number of bisections that led to it */
typedef struct _NcmSpectralPolarPanel
{
  gdouble ra;
  gdouble rb;
  gdouble ta;
  gdouble tb;
  guint depth;
} NcmSpectralPolarPanel;

/* A panel with its value, its change D_r + D_theta and its largest |F|, V */
typedef struct _NcmSpectralPolarEntry
{
  NcmSpectralPolarPanel panel;
  gdouble val;
  gdouble D;
  gdouble vs;
  guint lx; /* the levels it was accepted at */
  guint ly;
  gboolean replaced; /* expanded again by the second pass: its panels are later entries */
  gdouble Dr;        /* the radial and the angular changes of an initial panel */
  gdouble Dt;
  gsize series_off; /* its radial series in work->series_coeffs, when the call keeps G */
} NcmSpectralPolarEntry;

/*
 * The moments against T_k of (1 + s)^ca (1 - s)^cb on [-1, 1], for every level; have has
 * bit l set when level l is computed.
 */
typedef struct _NcmSpectralPolarMom
{
  gdouble ca;
  gdouble cb;
  GArray *m;
  guint have;
} NcmSpectralPolarMom;

/* The scratch space of one call, sized for level_cap */
typedef struct _NcmSpectralPolarWork
{
  guint level_cap;
  guint np;        /* 2^level_cap + 1 */
  gdouble *in;     /* the grid of the current level (lx, ly), (n + 1) x (m + 1), n = 2^lx
                    * and m = 2^ly, contiguous with row i = radius i: the input of the
                    * radial transform as it is; a doubling spreads it in place */
  gdouble *tmp;    /* its radial transform, same layout */
  gdouble *tmpT;   /* the transpose of tmp */
  gdouble *coeffs; /* the Chebyshev coefficients, (m + 1) x (n + 1): row l = angular index */
  gdouble *in_G;   /* the values of G without the weight, in the layout of in, when the call keeps them */
  gdouble *sv;     /* the angular integral of G at the radii of an accepted panel, np */
  gdouble *sc;     /* its Chebyshev coefficients, np */
  gdouble vs;      /* the largest |F| of the samples in in */
  GArray *lines;   /* the np-long arrays below, one after the other */
  gdouble *br;     /* l1 mass of the coefficients per radial index */
  gdouble *bt;     /* and per angular index */
  gdouble *rr;     /* radii of the panel at the current level, rr[i], i = 0..n */
  gdouble *wr;     /* regular weight at rr */
  gdouble *th;     /* angles of the panel at the current level, th[j], j = 0..m */
  gdouble *ct;     /* their cosines */
  gdouble *st;     /* their sines */
  gdouble *gth;    /* the angles of one callback, contiguous */
  gdouble *gct;
  gdouble *gst;
  gdouble *gf;                /* the values it returns */
  gdouble *xc;                /* cos (pi q / 2^level_cap), the Lobatto points of every level */
  gdouble *cr;                /* radii of one pass of _ncm_spectral_polar_nodes_r(), contiguous */
  gdouble *cw;                /* and their weights */
  gdouble *crho;              /* and the radial factor there (NcmSpectralPolarRadialF) */
  GArray *mt;                 /* angular moments, every level */
  GArray *cc;                 /* quadrature weights of the Lobatto points, every level */
  NcmSpectral *spectral;      /* the transforms from node values to Chebyshev coefficients */
  NcmSpectralPolarMom mom[4]; /* 0: weight 1, 1: r^a, 2: (R - r)^b, 3: both */
  GArray *initial;            /* NcmSpectralPolarEntry, the initial panels on their starting grids */
  GArray *initial_grid;       /* their samples */
  GArray *initial_grid_G;     /* their values of G, when the call keeps them */
  GArray *accepted;           /* NcmSpectralPolarEntry, the accepted panels */
  GArray *series_coeffs;      /* the radial series of the panels accepted so far, one after the other */
} NcmSpectralPolarWork;

/* The arguments and the state of one call of ncm_spectral_polar_integrate() */
typedef struct _NcmSpectralPolarCall
{
  NcmSpectralPolarF G;
  NcmSpectralPolarRadialF radial; /* the factor of G that depends only on r, or NULL */
  gpointer user_data;
  gdouble R;
  gdouble a;
  gdouble b;
  gdouble reltol;
  gdouble abstol_F; /* the absolute term of the test, MAX (reltol |I|, abstol) / A */
  gdouble A;        /* integral over the panels of the weight left out of F */
  guint level_min;
  guint level_cap;
  guint max_depth;
  guint64 n_eval;
  gboolean failed;
  gboolean keep_G; /* the values of G are kept and each accepted panel gets its radial series */
  NcmSpectralPolarPanel failed_panel;
  NcmSpectralPolarWork *work;
} NcmSpectralPolarCall;

/* *INDENT-OFF* */
G_DEFINE_QUARK (ncm-spectral-polar-error, ncm_spectral_polar_error)
/* *INDENT-ON* */
G_DEFINE_TYPE (NcmSpectralPolar, ncm_spectral_polar, G_TYPE_OBJECT)

static gpointer _ncm_spectral_polar_work_new (gpointer userdata);
static void _ncm_spectral_polar_work_free (gpointer p);

static void
ncm_spectral_polar_init (NcmSpectralPolar *sp)
{
  sp->level_min         = 0;
  sp->level_cap         = 0;
  sp->max_depth         = 0;
  sp->edge_ratio        = 0.0;
  sp->min_radial_panels = 0;
  sp->work_pool         = ncm_memory_pool_new (_ncm_spectral_polar_work_new, NULL,
                                               _ncm_spectral_polar_work_free);
}

static void
_ncm_spectral_polar_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmSpectralPolar *sp = NCM_SPECTRAL_POLAR (object);

  g_return_if_fail (NCM_IS_SPECTRAL_POLAR (object));

  switch (prop_id)
  {
    case PROP_LEVEL_MIN:
      ncm_spectral_polar_set_level_min (sp, g_value_get_uint (value));
      break;
    case PROP_LEVEL_CAP:
      ncm_spectral_polar_set_level_cap (sp, g_value_get_uint (value));
      break;
    case PROP_MAX_DEPTH:
      ncm_spectral_polar_set_max_depth (sp, g_value_get_uint (value));
      break;
    case PROP_EDGE_RATIO:
      ncm_spectral_polar_set_edge_ratio (sp, g_value_get_double (value));
      break;
    case PROP_MIN_RADIAL_PANELS:
      ncm_spectral_polar_set_min_radial_panels (sp, g_value_get_uint (value));
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_ncm_spectral_polar_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmSpectralPolar *sp = NCM_SPECTRAL_POLAR (object);

  g_return_if_fail (NCM_IS_SPECTRAL_POLAR (object));

  switch (prop_id)
  {
    case PROP_LEVEL_MIN:
      g_value_set_uint (value, sp->level_min);
      break;
    case PROP_LEVEL_CAP:
      g_value_set_uint (value, sp->level_cap);
      break;
    case PROP_MAX_DEPTH:
      g_value_set_uint (value, sp->max_depth);
      break;
    case PROP_EDGE_RATIO:
      g_value_set_double (value, sp->edge_ratio);
      break;
    case PROP_MIN_RADIAL_PANELS:
      g_value_set_uint (value, sp->min_radial_panels);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_ncm_spectral_polar_finalize (GObject *object)
{
  NcmSpectralPolar *sp = NCM_SPECTRAL_POLAR (object);

  ncm_memory_pool_free (sp->work_pool, TRUE);

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_spectral_polar_parent_class)->finalize (object);
}

static void
ncm_spectral_polar_class_init (NcmSpectralPolarClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);

  object_class->set_property = &_ncm_spectral_polar_set_property;
  object_class->get_property = &_ncm_spectral_polar_get_property;
  object_class->finalize     = &_ncm_spectral_polar_finalize;

  /**
   * NcmSpectralPolar:level-min:
   *
   * The level at which the expansion of a panel starts in both directions, below
   * #NcmSpectralPolar:level-cap; level $l$ has $2^l + 1$ points.
   */
  g_object_class_install_property (object_class,
                                   PROP_LEVEL_MIN,
                                   g_param_spec_uint ("level-min",
                                                      NULL,
                                                      "Starting level of a panel",
                                                      2, 9, 3,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmSpectralPolar:level-cap:
   *
   * The highest level of a direction of a panel; a panel that needs more is bisected in that
   * direction.
   */
  g_object_class_install_property (object_class,
                                   PROP_LEVEL_CAP,
                                   g_param_spec_uint ("level-cap",
                                                      NULL,
                                                      "Highest level of a panel",
                                                      3, 10, 7,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmSpectralPolar:max-depth:
   *
   * The largest number of bisections leading to a panel. A panel at this depth that needs
   * another bisection is the error %NCM_SPECTRAL_POLAR_ERROR_MAX_DEPTH.
   */
  g_object_class_install_property (object_class,
                                   PROP_MAX_DEPTH,
                                   g_param_spec_uint ("max-depth",
                                                      NULL,
                                                      "Largest number of bisections of a panel",
                                                      0, G_MAXUINT, 30,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmSpectralPolar:edge-ratio:
   *
   * The ratio $\omega > 1$ of ncm_spectral_polar_peak_edges(): its edges are at the
   * distances $L$, $\omega L$, $\omega^2 L$, ... from the peak, on each side of it, $L$ being
   * the half width of the peak in that direction.
   */
  g_object_class_install_property (object_class,
                                   PROP_EDGE_RATIO,
                                   g_param_spec_double ("edge-ratio",
                                                        NULL,
                                                        "Ratio between consecutive edge distances",
                                                        1.0, G_MAXDOUBLE, 4.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmSpectralPolar:min-radial-panels:
   *
   * The number $n \ge 1$ of equal radial intervals whose inner edges $k R / n$,
   * $k = 1, \dots, n - 1$, ncm_spectral_polar_peak_edges() adds to the edges of the peak.
   * They bound the width of the panels far from the peak, whose convergence test is
   * otherwise passed by wide panels with a larger error; a larger $n$ gives more panels of
   * lower degree.
   */
  g_object_class_install_property (object_class,
                                   PROP_MIN_RADIAL_PANELS,
                                   g_param_spec_uint ("min-radial-panels",
                                                      NULL,
                                                      "Number of equal radial intervals of the peak edges",
                                                      1, G_MAXUINT, 2,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
}

static gpointer
_ncm_spectral_polar_work_new (gpointer userdata)
{
  NcmSpectralPolarWork *work = g_new0 (NcmSpectralPolarWork, 1);
  guint i;

  work->level_cap    = 0;
  work->np           = 0;
  work->in           = NULL;
  work->tmp          = NULL;
  work->tmpT         = NULL;
  work->coeffs       = NULL;
  work->in_G         = NULL;
  work->sv           = NULL;
  work->sc           = NULL;
  work->lines        = g_array_new (FALSE, FALSE, sizeof (gdouble));
  work->mt           = g_array_new (FALSE, FALSE, sizeof (gdouble));
  work->cc           = g_array_new (FALSE, FALSE, sizeof (gdouble));
  work->spectral     = ncm_spectral_new ();
  work->initial      = g_array_new (FALSE, FALSE, sizeof (NcmSpectralPolarEntry));
  work->initial_grid = g_array_new (FALSE, FALSE, sizeof (gdouble));
  work->accepted     = g_array_new (FALSE, FALSE, sizeof (NcmSpectralPolarEntry));

  work->initial_grid_G = g_array_new (FALSE, FALSE, sizeof (gdouble));
  work->series_coeffs  = g_array_new (FALSE, FALSE, sizeof (gdouble));

  for (i = 0; i < 4; i++)
  {
    work->mom[i].ca   = 0.0;
    work->mom[i].cb   = 0.0;
    work->mom[i].m    = g_array_new (FALSE, FALSE, sizeof (gdouble));
    work->mom[i].have = 0;
  }

  return work;
}

static void
_ncm_spectral_polar_work_free (gpointer p)
{
  NcmSpectralPolarWork *work = p;
  guint i;

  g_clear_pointer (&work->in, fftw_free);
  g_clear_pointer (&work->tmp, fftw_free);
  g_clear_pointer (&work->tmpT, fftw_free);
  g_clear_pointer (&work->coeffs, fftw_free);
  g_clear_pointer (&work->in_G, fftw_free);
  g_clear_pointer (&work->sv, fftw_free);
  g_clear_pointer (&work->sc, fftw_free);
  g_array_unref (work->lines);
  g_array_unref (work->mt);
  g_array_unref (work->cc);
  ncm_spectral_free (work->spectral);
  g_array_unref (work->initial);
  g_array_unref (work->initial_grid);
  g_array_unref (work->initial_grid_G);
  g_array_unref (work->accepted);
  g_array_unref (work->series_coeffs);

  for (i = 0; i < 4; i++)
    g_array_unref (work->mom[i].m);

  g_free (work);
}

/**
 * ncm_spectral_polar_new:
 *
 * Creates an integrator with the default properties.
 *
 * Returns: (transfer full): a new #NcmSpectralPolar.
 */
NcmSpectralPolar *
ncm_spectral_polar_new (void)
{
  return g_object_new (NCM_TYPE_SPECTRAL_POLAR, NULL);
}

/**
 * ncm_spectral_polar_ref:
 * @sp: a #NcmSpectralPolar
 *
 * Increases the reference count of @sp by one.
 *
 * Returns: (transfer full): @sp.
 */
NcmSpectralPolar *
ncm_spectral_polar_ref (NcmSpectralPolar *sp)
{
  return g_object_ref (sp);
}

/**
 * ncm_spectral_polar_free:
 * @sp: a #NcmSpectralPolar
 *
 * Decreases the reference count of @sp by one.
 */
void
ncm_spectral_polar_free (NcmSpectralPolar *sp)
{
  g_object_unref (sp);
}

/**
 * ncm_spectral_polar_clear:
 * @sp: a #NcmSpectralPolar
 *
 * Decreases the reference count of *@sp by one and sets *@sp to %NULL.
 */
void
ncm_spectral_polar_clear (NcmSpectralPolar **sp)
{
  g_clear_object (sp);
}

/**
 * ncm_spectral_polar_set_level_min:
 * @sp: a #NcmSpectralPolar
 * @level_min: starting level of a panel
 *
 * Sets #NcmSpectralPolar:level-min.
 */
void
ncm_spectral_polar_set_level_min (NcmSpectralPolar *sp, guint level_min)
{
  sp->level_min = level_min;
}

/**
 * ncm_spectral_polar_set_level_cap:
 * @sp: a #NcmSpectralPolar
 * @level_cap: highest level of a panel
 *
 * Sets #NcmSpectralPolar:level-cap.
 */
void
ncm_spectral_polar_set_level_cap (NcmSpectralPolar *sp, guint level_cap)
{
  sp->level_cap = level_cap;
}

/**
 * ncm_spectral_polar_set_max_depth:
 * @sp: a #NcmSpectralPolar
 * @max_depth: largest number of bisections of a panel
 *
 * Sets #NcmSpectralPolar:max-depth.
 */
void
ncm_spectral_polar_set_max_depth (NcmSpectralPolar *sp, guint max_depth)
{
  sp->max_depth = max_depth;
}

/**
 * ncm_spectral_polar_set_edge_ratio:
 * @sp: a #NcmSpectralPolar
 * @edge_ratio: ratio between consecutive edge distances
 *
 * Sets #NcmSpectralPolar:edge-ratio, which must be above 1.
 */
void
ncm_spectral_polar_set_edge_ratio (NcmSpectralPolar *sp, gdouble edge_ratio)
{
  if (!(edge_ratio > 1.0))
    g_error ("ncm_spectral_polar_set_edge_ratio: edge-ratio %g is not above 1.", edge_ratio);

  sp->edge_ratio = edge_ratio;
}

/**
 * ncm_spectral_polar_set_min_radial_panels:
 * @sp: a #NcmSpectralPolar
 * @n: number of equal radial intervals
 *
 * Sets #NcmSpectralPolar:min-radial-panels, which must be at least 1.
 */
void
ncm_spectral_polar_set_min_radial_panels (NcmSpectralPolar *sp, guint n)
{
  if (n < 1)
    g_error ("ncm_spectral_polar_set_min_radial_panels: min-radial-panels must be at least 1.");

  sp->min_radial_panels = n;
}

/**
 * ncm_spectral_polar_get_level_min:
 * @sp: a #NcmSpectralPolar
 *
 * Returns: #NcmSpectralPolar:level-min.
 */
guint
ncm_spectral_polar_get_level_min (NcmSpectralPolar *sp)
{
  return sp->level_min;
}

/**
 * ncm_spectral_polar_get_level_cap:
 * @sp: a #NcmSpectralPolar
 *
 * Returns: #NcmSpectralPolar:level-cap.
 */
guint
ncm_spectral_polar_get_level_cap (NcmSpectralPolar *sp)
{
  return sp->level_cap;
}

/**
 * ncm_spectral_polar_get_max_depth:
 * @sp: a #NcmSpectralPolar
 *
 * Returns: #NcmSpectralPolar:max-depth.
 */
guint
ncm_spectral_polar_get_max_depth (NcmSpectralPolar *sp)
{
  return sp->max_depth;
}

/**
 * ncm_spectral_polar_get_edge_ratio:
 * @sp: a #NcmSpectralPolar
 *
 * Returns: #NcmSpectralPolar:edge-ratio.
 */
gdouble
ncm_spectral_polar_get_edge_ratio (NcmSpectralPolar *sp)
{
  return sp->edge_ratio;
}

/**
 * ncm_spectral_polar_get_min_radial_panels:
 * @sp: a #NcmSpectralPolar
 *
 * Returns: #NcmSpectralPolar:min-radial-panels.
 */
guint
ncm_spectral_polar_get_min_radial_panels (NcmSpectralPolar *sp)
{
  return sp->min_radial_panels;
}

static void _ncm_spectral_polar_check (NcmSpectralPolar *sp, const gchar *name, gdouble R, gdouble a, gdouble b, GArray *r_edges, GArray *theta_edges, gdouble reltol, gdouble abstol);
static gdouble _ncm_spectral_polar_integrate (NcmSpectralPolar *sp, const gchar *name, NcmSpectralPolarF G, NcmSpectralPolarRadialF radial, gpointer user_data, gdouble R, gdouble a, gdouble b, GArray *r_edges, GArray *theta_edges, gdouble reltol, gdouble abstol, guint64 *n_eval, NcmMatrix **panels, GArray *series, GArray *coeffs, GError **error);
static void _ncm_spectral_polar_check_edges (NcmSpectralPolar *sp, const gchar *func, GArray *r_edges, GArray *theta_edges);
static void _ncm_spectral_polar_work_prepare (NcmSpectralPolarWork *work, guint level_cap);
static void _ncm_spectral_polar_work_set_exponents (NcmSpectralPolarWork *work, gdouble a, gdouble b);
static gdouble _ncm_spectral_polar_first_pass (NcmSpectralPolarCall *call, GArray *r_edges, GArray *theta_edges);
static void _ncm_spectral_polar_expand_initial (NcmSpectralPolarCall *call);
static gdouble _ncm_spectral_polar_retest (NcmSpectralPolarCall *call, gdouble abstol);
static NcmMatrix *_ncm_spectral_polar_panels (NcmSpectralPolarWork *work);
static void _ncm_spectral_polar_series (NcmSpectralPolarWork *work, GArray *series, GArray *coeffs);

/**
 * ncm_spectral_polar_integrate:
 * @sp: a #NcmSpectralPolar
 * @G: (scope call) (closure user_data): the function $G$
 * @user_data: user data for @G
 * @R: radius of the disk
 * @a: exponent of the weight at $r = 0$
 * @b: exponent of the weight at $r = R$
 * @r_edges: (element-type gdouble): the radial edges $r_0, \dots, r_{n_r}$
 * @theta_edges: (element-type gdouble): the angular edges $\theta_0, \dots, \theta_{n_\theta}$
 * @reltol: relative tolerance
 * @abstol: absolute tolerance, in units of the integral
 * @n_eval: (out) (optional): number of values of $G$ computed
 * @panels: (out) (optional) (transfer full): the accepted panels
 * @error: a #GError
 *
 * Computes the integral of $r^a (R - r)^b G (r, \theta)$ over the panels of @r_edges times
 * @theta_edges, as described in #NcmSpectralPolar. @r_edges is ascending in $[0, @R]$,
 * @theta_edges is ascending and spans at most $2\pi$, each with two entries or more, and
 * @a and @b are above $-1$. The power law at $r = 0$ is integrated exactly only on panels
 * whose first radial edge is 0, and the one at $r = @R$ only on panels whose last radial edge
 * is @R.
 *
 * @G is called with one radius and the angles of a panel where it has no value yet, with
 * their cosines and sines; the three arrays are valid during the call only.
 *
 * With @panels not %NULL, the accepted panels are returned in a new #NcmMatrix, one row
 * per panel with the columns $r_a$, $r_b$, $\theta_a$, $\theta_b$, $l_r$, $l_\theta$, the
 * value of the panel and its change $D_r + D_\theta$; the values add up to the integral.
 * It is built after the integration and is meant for inspecting where the panels were
 * refined.
 *
 * A panel that needs more than #NcmSpectralPolar:max-depth bisections sets @error to
 * %NCM_SPECTRAL_POLAR_ERROR_MAX_DEPTH, the function returns NaN and *@panels is %NULL.
 *
 * Returns: the integral.
 */
gdouble
ncm_spectral_polar_integrate (NcmSpectralPolar *sp, NcmSpectralPolarF G, gpointer user_data,
                              gdouble R, gdouble a, gdouble b, GArray *r_edges, GArray *theta_edges,
                              gdouble reltol, gdouble abstol, guint64 *n_eval, NcmMatrix **panels, GError **error)
{
  return _ncm_spectral_polar_integrate (sp, "ncm_spectral_polar_integrate", G, NULL, user_data, R, a, b, r_edges, theta_edges,
                                        reltol, abstol, n_eval, panels, NULL, NULL, error);
}

/*
 * _ncm_spectral_polar_integrate_radial:
 *
 * As ncm_spectral_polar_integrate(), with the integrand r^a (R - r)^b rho (r) G (r, theta),
 * rho given by radial (NcmSpectralPolarRadialF). radial is called with the radii of a panel
 * where the integrand has no value yet and sets rho at them, which the weights are then
 * multiplied by; the integral and the refinement are those of
 * ncm_spectral_polar_integrate() with rho in G, up to rounding. A factor of G that depends
 * only on the radius is computed once per radius this way, for all of them in one call,
 * instead of in every call of G. With radial NULL it is ncm_spectral_polar_integrate().
 *
 * Declared in ncm_spectral_polar_private.h, for the other classes of the library.
 */
gdouble
_ncm_spectral_polar_integrate_radial (NcmSpectralPolar *sp, NcmSpectralPolarF G, NcmSpectralPolarRadialF radial, gpointer user_data,
                                      gdouble R, gdouble a, gdouble b, GArray *r_edges, GArray *theta_edges,
                                      gdouble reltol, gdouble abstol, guint64 *n_eval, NcmMatrix **panels, GError **error)
{
  return _ncm_spectral_polar_integrate (sp, "_ncm_spectral_polar_integrate_radial", G, radial, user_data, R, a, b, r_edges, theta_edges,
                                        reltol, abstol, n_eval, panels, NULL, NULL, error);
}

/*
 * _ncm_spectral_polar_integrate_radial_series:
 *
 * As _ncm_spectral_polar_integrate_radial(), and for each accepted panel the Chebyshev series
 * in r of the angular integral of G over the panel, G without the weight and without rho. The
 * integral, the panels and the number of values of G are those of
 * _ncm_spectral_polar_integrate_radial(), bit for bit: the values of G are stored a second
 * time as the callback returns them, and the series of a panel is computed from them when
 * the panel is accepted, by the quadrature of the angular Lobatto points and one transform
 * in r. series receives one NcmSpectralPolarRadialSeries per accepted panel and coeffs their
 * coefficients; both are arrays of the caller, emptied first, of elements
 * NcmSpectralPolarRadialSeries and gdouble. The integral of w (r) G over the domain, for a
 * radial function w, is the sum over the panels of the integral of w times the series. When
 * the integral fails both arrays are empty.
 *
 * Declared in ncm_spectral_polar_private.h, for the other classes of the library.
 */
gdouble
_ncm_spectral_polar_integrate_radial_series (NcmSpectralPolar *sp, NcmSpectralPolarF G, NcmSpectralPolarRadialF radial, gpointer user_data,
                                             gdouble R, gdouble a, gdouble b, GArray *r_edges, GArray *theta_edges,
                                             gdouble reltol, gdouble abstol, guint64 *n_eval, GArray *series, GArray *coeffs, GError **error)
{
  g_assert_nonnull (series);
  g_assert_nonnull (coeffs);
  g_assert_cmpuint (g_array_get_element_size (series), ==, sizeof (NcmSpectralPolarRadialSeries));
  g_assert_cmpuint (g_array_get_element_size (coeffs), ==, sizeof (gdouble));

  return _ncm_spectral_polar_integrate (sp, "_ncm_spectral_polar_integrate_radial_series", G, radial, user_data, R, a, b, r_edges, theta_edges,
                                        reltol, abstol, n_eval, NULL, series, coeffs, error);
}

/*
 * The integral of _ncm_spectral_polar_integrate_radial(), with name the function called, for
 * its messages. With series and coeffs not NULL, also the radial series of the accepted
 * panels (_ncm_spectral_polar_integrate_radial_series()).
 */
static gdouble
_ncm_spectral_polar_integrate (NcmSpectralPolar *sp, const gchar *name, NcmSpectralPolarF G, NcmSpectralPolarRadialF radial,
                               gpointer user_data, gdouble R, gdouble a, gdouble b, GArray *r_edges, GArray *theta_edges,
                               gdouble reltol, gdouble abstol, guint64 *n_eval, NcmMatrix **panels,
                               GArray *series, GArray *coeffs, GError **error)
{
  NcmSpectralPolarWork **ptr = ncm_memory_pool_get (sp->work_pool);
  NcmSpectralPolarWork *work = *ptr;
  NcmSpectralPolarCall call  = {
    .G            = G,
    .radial       = radial,
    .user_data    = user_data,
    .R            = R,
    .a            = a,
    .b            = b,
    .reltol       = reltol,
    .abstol_F     = 0.0,
    .A            = 0.0,
    .level_min    = sp->level_min,
    .level_cap    = sp->level_cap,
    .max_depth    = sp->max_depth,
    .n_eval       = 0,
    .failed       = FALSE,
    .keep_G       = (series != NULL),
    .failed_panel = {0.0, 0.0, 0.0, 0.0, 0},
    .work         = work
  };
  gdouble I0, res;

  _ncm_spectral_polar_check (sp, name, R, a, b, r_edges, theta_edges, reltol, abstol);
  _ncm_spectral_polar_work_prepare (work, sp->level_cap);
  _ncm_spectral_polar_work_set_exponents (work, a, b);
  g_array_set_size (work->series_coeffs, 0);

  I0            = _ncm_spectral_polar_first_pass (&call, r_edges, theta_edges);
  call.abstol_F = MAX (reltol * fabs (I0), abstol) / call.A;
  _ncm_spectral_polar_expand_initial (&call);
  res = _ncm_spectral_polar_retest (&call, abstol);

  if (panels != NULL)
    *panels = call.failed ? NULL : _ncm_spectral_polar_panels (work);

  if (series != NULL)
  {
    g_array_set_size (series, 0);
    g_array_set_size (coeffs, 0);

    if (!call.failed)
      _ncm_spectral_polar_series (work, series, coeffs);
  }

  ncm_memory_pool_return (ptr);

  if (n_eval != NULL)
    *n_eval = call.n_eval;

  if (call.failed)
  {
    ncm_util_set_or_call_error (error, NCM_SPECTRAL_POLAR_ERROR, NCM_SPECTRAL_POLAR_ERROR_MAX_DEPTH,
                                "%s: the panel r in [% 22.15g, % 22.15g], "
                                "theta in [% 22.15g, % 22.15g] does not pass the test after %u bisections.",
                                name, call.failed_panel.ra, call.failed_panel.rb, call.failed_panel.ta, call.failed_panel.tb,
                                call.failed_panel.depth);

    return GSL_NAN;
  }

  return res;
}

static void _ncm_spectral_polar_expand (NcmSpectralPolarCall *call, NcmSpectralPolarPanel panel, const NcmSpectralPolarEntry *start);

/* The arguments of an integral, checked before anything is evaluated */
static void
_ncm_spectral_polar_check (NcmSpectralPolar *sp, const gchar *name, gdouble R, gdouble a, gdouble b, GArray *r_edges, GArray *theta_edges,
                           gdouble reltol, gdouble abstol)
{
  if (!(R > 0.0))
    g_error ("%s: the radius %g is not positive.", name, R);

  if (!(a > -1.0) || !(b > -1.0))
    g_error ("%s: the exponents a = %g and b = %g must be above -1.", name, a, b);

  if (!(reltol > 0.0) || !(abstol >= 0.0))
    g_error ("%s: reltol %g must be positive and abstol %g non-negative.", name, reltol, abstol);

  _ncm_spectral_polar_check_edges (sp, name, r_edges, theta_edges);

  if (g_array_index (r_edges, gdouble, r_edges->len - 1) > R)
    g_error ("%s: the radial edges [%g, %g] are not inside [0, %g].", name,
             g_array_index (r_edges, gdouble, 0), g_array_index (r_edges, gdouble, r_edges->len - 1), R);
}

/*
 * The levels of sp and the edges of a call of func: two or more edges in each direction,
 * ascending, the radial ones from 0 or above and the angular ones spanning at most 2 pi.
 */
static void
_ncm_spectral_polar_check_edges (NcmSpectralPolar *sp, const gchar *func, GArray *r_edges, GArray *theta_edges)
{
  guint i;

  if (sp->level_min >= sp->level_cap)
    g_error ("%s: level-min %u is not below level-cap %u.", func, sp->level_min, sp->level_cap);

  if ((r_edges->len < 2) || (theta_edges->len < 2))
    g_error ("%s: %u radial and %u angular edges, two or more of each are needed.",
             func, r_edges->len, theta_edges->len);

  for (i = 0; i + 1 < r_edges->len; i++)
  {
    if (!(g_array_index (r_edges, gdouble, i) < g_array_index (r_edges, gdouble, i + 1)))
      g_error ("%s: the radial edges are not ascending at index %u.", func, i);
  }

  for (i = 0; i + 1 < theta_edges->len; i++)
  {
    if (!(g_array_index (theta_edges, gdouble, i) < g_array_index (theta_edges, gdouble, i + 1)))
      g_error ("%s: the angular edges are not ascending at index %u.", func, i);
  }

  if (g_array_index (r_edges, gdouble, 0) < 0.0)
    g_error ("%s: the first radial edge %g is negative.", func, g_array_index (r_edges, gdouble, 0));

  if (g_array_index (theta_edges, gdouble, theta_edges->len - 1) - g_array_index (theta_edges, gdouble, 0) > 2.0 * M_PI * (1.0 + GSL_DBL_EPSILON))
    g_error ("%s: the angular edges span %g, more than 2 pi.", func,
             g_array_index (theta_edges, gdouble, theta_edges->len - 1) - g_array_index (theta_edges, gdouble, 0));
}

/* Sizes the scratch space for level_cap and fills the tables that depend on it only */
static void
_ncm_spectral_polar_work_prepare (NcmSpectralPolarWork *work, guint level_cap)
{
  const guint np        = (1u << level_cap) + 1;
  const guint n_lines   = 15;
  const gsize table_len = NCM_SPECTRAL_POLAR_LEVEL_OFFSET (level_cap + 1);
  gdouble *lines, *mt, *cc;
  guint l, k;

  if (work->level_cap == level_cap)
    return;

  work->level_cap = level_cap;
  work->np        = np;

  g_array_set_size (work->lines, (gsize) n_lines * np);
  g_array_set_size (work->mt, table_len);
  g_array_set_size (work->cc, table_len);

  g_clear_pointer (&work->in, fftw_free);
  g_clear_pointer (&work->tmp, fftw_free);
  g_clear_pointer (&work->tmpT, fftw_free);
  g_clear_pointer (&work->coeffs, fftw_free);
  g_clear_pointer (&work->in_G, fftw_free);
  g_clear_pointer (&work->sv, fftw_free);
  g_clear_pointer (&work->sc, fftw_free);
  work->in     = fftw_malloc (sizeof (gdouble) * np * np);
  work->tmp    = fftw_malloc (sizeof (gdouble) * np * np);
  work->tmpT   = fftw_malloc (sizeof (gdouble) * np * np);
  work->coeffs = fftw_malloc (sizeof (gdouble) * np * np);
  work->in_G   = fftw_malloc (sizeof (gdouble) * np * np);
  work->sv     = fftw_malloc (sizeof (gdouble) * np);
  work->sc     = fftw_malloc (sizeof (gdouble) * np);

  lines      = &g_array_index (work->lines, gdouble, 0);
  work->br   = lines + 0 * np;
  work->bt   = lines + 1 * np;
  work->rr   = lines + 2 * np;
  work->wr   = lines + 3 * np;
  work->th   = lines + 4 * np;
  work->ct   = lines + 5 * np;
  work->st   = lines + 6 * np;
  work->gth  = lines + 7 * np;
  work->gct  = lines + 8 * np;
  work->gst  = lines + 9 * np;
  work->gf   = lines + 10 * np;
  work->xc   = lines + 11 * np;
  work->cr   = lines + 12 * np;
  work->cw   = lines + 13 * np;
  work->crho = lines + 14 * np;

  for (k = 0; k < np; k++)
    work->xc[k] = cos (M_PI * k / (np - 1.0));

  work->xc[(np - 1) / 2] = 0.0;

  /* Moments of weight 1 against T_l, 2 / (1 - l^2) for even l */
  mt = &g_array_index (work->mt, gdouble, 0);

  for (l = 0; l <= level_cap; l++)
  {
    const guint m = 1u << l;

    for (k = 0; k <= m; k++)
      mt[NCM_SPECTRAL_POLAR_LEVEL_OFFSET (l) + k] = (k % 2 == 0) ? 2.0 / (1.0 - (gdouble) k * k) : 0.0;
  }

  /*
   * The weights of the integral over [-1, 1] of the interpolant at the m + 1 Lobatto points
   * cos (pi j / m) (Clenshaw-Curtis): w_j = (c_j / m) (1 - sum_{k = 1}^{m / 2} b_k
   * cos (2 pi j k / m) / (4 k^2 - 1)), with c_j = 1 at the two ends and 2 otherwise, and
   * b_k = 1 for k = m / 2 and 2 otherwise
   */
  cc = &g_array_index (work->cc, gdouble, 0);

  for (l = 0; l <= level_cap; l++)
  {
    const guint m = 1u << l;
    guint j;

    for (j = 0; j <= m; j++)
    {
      gdouble s = 1.0;

      for (k = 1; k <= m / 2; k++)
        s -= ((2 * k == m) ? 1.0 : 2.0) * cos (2.0 * M_PI * j * k / m) / (4.0 * k * k - 1.0);

      cc[NCM_SPECTRAL_POLAR_LEVEL_OFFSET (l) + j] = (((j == 0) || (j == m)) ? 1.0 : 2.0) * s / m;
    }
  }

  for (k = 0; k < 4; k++)
  {
    g_array_set_size (work->mom[k].m, table_len);
    work->mom[k].have = 0;
  }
}

/*
 * Sets the exponent pairs of the four moment tables, (0, 0), (a, 0), (0, b) and (a, b),
 * discarding the levels computed for other exponents.
 */
static void
_ncm_spectral_polar_work_set_exponents (NcmSpectralPolarWork *work, gdouble a, gdouble b)
{
  guint k;

  for (k = 0; k < 4; k++)
  {
    const gdouble ca = (k & 1) ? a : 0.0;
    const gdouble cb = (k & 2) ? b : 0.0;

    if ((work->mom[k].ca != ca) || (work->mom[k].cb != cb))
    {
      work->mom[k].ca   = ca;
      work->mom[k].cb   = cb;
      work->mom[k].have = 0;
    }
  }
}

static void _ncm_spectral_polar_sample_start (NcmSpectralPolarCall *call, const NcmSpectralPolarPanel *panel);
static gsize _ncm_spectral_polar_series_append (NcmSpectralPolarWork *work, const NcmSpectralPolarPanel *panel, guint lx, guint ly);
static gdouble _ncm_spectral_polar_coeffs (NcmSpectralPolarWork *work, guint lx, guint ly);
static gdouble _ncm_spectral_polar_change (NcmSpectralPolarWork *work, guint lx, guint ly, gdouble *Dr, gdouble *Dt);
static gdouble _ncm_spectral_polar_panel_value (NcmSpectralPolarCall *call, const NcmSpectralPolarPanel *panel, guint lx, guint ly);
static guint _ncm_spectral_polar_panel_kind (const NcmSpectralPolarCall *call, const NcmSpectralPolarPanel *panel);

/*
 * Samples every initial panel on its starting grid and keeps its samples, value, change and
 * scale in work->initial and work->initial_grid. Sets call->A and returns the sum of the
 * values.
 */
static gdouble
_ncm_spectral_polar_first_pass (NcmSpectralPolarCall *call, GArray *r_edges, GArray *theta_edges)
{
  NcmSpectralPolarWork *work = call->work;
  const guint ns             = (1u << call->level_min) + 1;
  const gsize nb             = (gsize) (r_edges->len - 1) * (theta_edges->len - 1);
  gdouble I0                 = 0.0;
  gdouble A                  = 0.0;
  gsize ib                   = 0;
  guint i, j;

  g_array_set_size (work->initial, nb);
  g_array_set_size (work->initial_grid, nb * ns * ns);
  g_array_set_size (work->initial_grid_G, call->keep_G ? nb * ns * ns : 0);
  g_array_set_size (work->accepted, 0);

  for (i = 0; i + 1 < r_edges->len; i++)
  {
    for (j = 0; j + 1 < theta_edges->len; j++, ib++)
    {
      NcmSpectralPolarEntry *entry = &g_array_index (work->initial, NcmSpectralPolarEntry, ib);
      gdouble *start               = &g_array_index (work->initial_grid, gdouble, ib * ns * ns);
      gdouble Dr, Dt;
      guint kind;

      entry->panel.ra    = g_array_index (r_edges, gdouble, i);
      entry->panel.rb    = g_array_index (r_edges, gdouble, i + 1);
      entry->panel.ta    = g_array_index (theta_edges, gdouble, j);
      entry->panel.tb    = g_array_index (theta_edges, gdouble, j + 1);
      entry->panel.depth = 0;
      entry->lx          = call->level_min;
      entry->ly          = call->level_min;
      entry->replaced    = FALSE;

      _ncm_spectral_polar_sample_start (call, &entry->panel);
      memcpy (start, work->in, sizeof (gdouble) * ns * ns);

      entry->vs  = _ncm_spectral_polar_coeffs (work, call->level_min, call->level_min);
      entry->D   = _ncm_spectral_polar_change (work, call->level_min, call->level_min, &Dr, &Dt);
      entry->Dr  = Dr;
      entry->Dt  = Dt;
      entry->val = _ncm_spectral_polar_panel_value (call, &entry->panel, call->level_min, call->level_min);
      I0        += entry->val;

      entry->series_off = 0;

      if (call->keep_G)
      {
        memcpy (&g_array_index (work->initial_grid_G, gdouble, ib * ns * ns), work->in_G, sizeof (gdouble) * ns * ns);
        entry->series_off = _ncm_spectral_polar_series_append (work, &entry->panel, call->level_min, call->level_min);
      }

      /* The weight left out of F on this panel, integrated over it */
      kind = _ncm_spectral_polar_panel_kind (call, &entry->panel);
      {
        const gdouble ca = (kind & 1) ? call->a : 0.0;
        const gdouble cb = (kind & 2) ? call->b : 0.0;

        A += pow (entry->panel.rb - entry->panel.ra, 1.0 + ca + cb) * exp (gsl_sf_lnbeta (1.0 + ca, 1.0 + cb)) *
             (entry->panel.tb - entry->panel.ta);
      }
    }
  }

  call->A = A;

  return I0;
}

static void _ncm_spectral_polar_nodes_r (NcmSpectralPolarCall *call, const NcmSpectralPolarPanel *panel, guint lx, guint i0, guint step);
static void _ncm_spectral_polar_nodes_theta (NcmSpectralPolarCall *call, const NcmSpectralPolarPanel *panel, guint ly, guint j0, guint step);
static void _ncm_spectral_polar_rows (NcmSpectralPolarCall *call, guint lx, guint ly, guint i0, guint istep, guint j0, guint jstep);

/* Samples F on the starting grid of the panel; work->vs restarts from these samples */
static void
_ncm_spectral_polar_sample_start (NcmSpectralPolarCall *call, const NcmSpectralPolarPanel *panel)
{
  _ncm_spectral_polar_nodes_r (call, panel, call->level_min, 0, 1);
  _ncm_spectral_polar_nodes_theta (call, panel, call->level_min, 0, 1);
  call->work->vs = 0.0;
  _ncm_spectral_polar_rows (call, call->level_min, call->level_min, 0, 1, 0, 1);
}

/* Points of a pass of _ncm_spectral_polar_nodes_r() from which its weights are computed in
 * loops that vectorize (pow(), see ncm_vmath_private.h); below it they are computed point by point */
#define NCM_SPECTRAL_POLAR_SIMD_MIN (8)

/*
 * Radii of the panel at level lx and the regular weight there, r^a unless the panel starts at
 * 0 and (R - r)^b unless it ends at R: every point (i0 = 0, step 1) or the points a doubling
 * adds (i0 = 1, step 2), after the points of level lx - 1 are moved from i to 2 i.
 */
static void
_ncm_spectral_polar_nodes_r (NcmSpectralPolarCall *call, const NcmSpectralPolarPanel *panel, guint lx, guint i0, guint step)
{
  NcmSpectralPolarWork *work = call->work;
  const guint n              = 1u << lx;
  const guint sh             = call->level_cap - lx;
  const gdouble mid          = 0.5 * (panel->ra + panel->rb);
  const gdouble half         = 0.5 * (panel->rb - panel->ra);
  gdouble *rr                = work->rr;
  gdouble *wr                = work->wr;
  const gboolean wa          = (panel->ra != 0.0) && (call->a != 0.0);
  const gboolean wb          = (panel->rb != call->R) && (call->b != 0.0);
  const gdouble a            = call->a;
  const gdouble b            = call->b;
  const gdouble R            = call->R;
  gdouble *cr                = work->cr;
  gdouble *cw                = work->cw;
  guint i, c;

  if (step == 2)
  {
    for (i = n / 2; i > 0; i--)
    {
      rr[2 * i] = rr[i];
      wr[2 * i] = wr[i];
    }
  }

  /* The radii contiguous, the weights, then both at their positions */
  for (c = 0, i = i0; i <= n; i += step, c++)
    cr[c] = mid + half * work->xc[i << sh];

  if (c < NCM_SPECTRAL_POLAR_SIMD_MIN)
  {
    /*
     * A pass of few points: the weights point by point with pow(), as the vectorized
     * loops below cost more than they save under NCM_SPECTRAL_POLAR_SIMD_MIN points
     * (passes of 1 to 5 points occur from level_min 0)
     */
    for (i = 0; i < c; i++)
    {
      gdouble w = 1.0;

      if (wa)
        w *= (a == 1.0) ? cr[i] : pow (cr[i], a);

      if (wb)
        w *= pow (R - cr[i], b);

      cw[i] = w;
    }
  }
  else
  {
    for (i = 0; i < c; i++)
      cw[i] = 1.0;

    /* r >= ra > 0 here, and R - r >= R - rb > 0 below: the bases are positive */
    if (wa && (a == 1.0))
    {
      for (i = 0; i < c; i++)
        cw[i] = cr[i];
    }
    else if (wa)
    {
      #pragma omp simd

      for (i = 0; i < c; i++)
        cw[i] = pow (cr[i], a);
    }

    if (wb)
    {
      #pragma omp simd

      for (i = 0; i < c; i++)
        cw[i] *= pow (R - cr[i], b);
    }
  }

  if (call->radial != NULL)
  {
    gdouble *crho = work->crho;

    call->radial (call->user_data, c, cr, crho);

    #pragma omp simd

    for (i = 0; i < c; i++)
      cw[i] *= crho[i];
  }

  for (c = 0, i = i0; i <= n; i += step, c++)
  {
    rr[i] = cr[c];
    wr[i] = cw[c];
  }
}

/* Angles of the panel at level ly, with their cosines and sines, as the radii above */
static void
_ncm_spectral_polar_nodes_theta (NcmSpectralPolarCall *call, const NcmSpectralPolarPanel *panel, guint ly, guint j0, guint step)
{
  NcmSpectralPolarWork *work = call->work;
  const guint m              = 1u << ly;
  const guint sh             = call->level_cap - ly;
  const gdouble mid          = 0.5 * (panel->ta + panel->tb);
  const gdouble half         = 0.5 * (panel->tb - panel->ta);
  gdouble *th                = work->th;
  gdouble *ct                = work->ct;
  gdouble *st                = work->st;
  guint j;

  if (step == 2)
  {
    for (j = m / 2; j > 0; j--)
    {
      th[2 * j] = th[j];
      ct[2 * j] = ct[j];
      st[2 * j] = st[j];
    }
  }

  #pragma omp simd

  for (j = j0; j <= m; j += step)
  {
    const gdouble t = mid + half * work->xc[j << sh];

    th[j] = t;
    st[j] = sin (t);
    ct[j] = cos (t);
  }
}

/*
 * Calls G at the rows i0, i0 + istep, ... and the angles j0, j0 + jstep, ... of level
 * (lx, ly) and stores F, G times the regular weight, in the grid of that level. Updates
 * work->vs, the largest |F| of the samples in the grid, with the values it stores. Four
 * running maxima per row: with one, each comparison waits for the one before.
 */
static inline __attribute__ ((always_inline)) void

_ncm_spectral_polar_rows_body (NcmSpectralPolarCall *call, guint lx, guint ly, guint i0, guint istep, guint j0, guint jstep, const gboolean keep_G)
{
  NcmSpectralPolarWork *work = call->work;
  const guint n              = 1u << lx;
  const guint m              = 1u << ly;
  guint c                    = 0;
  guint i, j, q;

  for (j = j0; j <= m; j += jstep, c++)
  {
    work->gth[c] = work->th[j];
    work->gct[c] = work->ct[j];
    work->gst[c] = work->st[j];
  }

  for (i = i0; i <= n; i += istep)
  {
    const gdouble w = work->wr[i];
    gdouble *row    = work->in + (gsize) i * (m + 1);
    gdouble v0      = work->vs;
    gdouble v1      = 0.0;
    gdouble v2      = 0.0;
    gdouble v3      = 0.0;

    call->G (call->user_data, work->rr[i], c, work->gth, work->gct, work->gst, work->gf);

    if (keep_G)
    {
      gdouble *row_G = work->in_G + (gsize) i * (m + 1);

      for (q = 0, j = j0; q < c; q++, j += jstep)
        row_G[j] = work->gf[q];
    }

    for (q = 0, j = j0; q + 4 <= c; q += 4, j += 4 * jstep)
    {
      const gdouble f0 = w * work->gf[q];
      const gdouble f1 = w * work->gf[q + 1];
      const gdouble f2 = w * work->gf[q + 2];
      const gdouble f3 = w * work->gf[q + 3];

      row[j]             = f0;
      row[j + jstep]     = f1;
      row[j + 2 * jstep] = f2;
      row[j + 3 * jstep] = f3;

      v0 = MAX (v0, fabs (f0));
      v1 = MAX (v1, fabs (f1));
      v2 = MAX (v2, fabs (f2));
      v3 = MAX (v3, fabs (f3));
    }

    for ( ; q < c; q++, j += jstep)
    {
      const gdouble f = w * work->gf[q];

      row[j] = f;
      v0     = MAX (v0, fabs (f));
    }

    work->vs = MAX (MAX (v0, v1), MAX (v2, v3));
  }

  call->n_eval += (guint64) c * ((n - i0) / istep + 1);
}

/* _ncm_spectral_polar_rows_body() for a call that keeps the values of G and for one that does not, one copy for each */
static void
_ncm_spectral_polar_rows (NcmSpectralPolarCall *call, guint lx, guint ly, guint i0, guint istep, guint j0, guint jstep)
{
  if (call->keep_G)
    _ncm_spectral_polar_rows_body (call, lx, ly, i0, istep, j0, jstep, TRUE);
  else
    _ncm_spectral_polar_rows_body (call, lx, ly, i0, istep, j0, jstep, FALSE);
}

/*
 * Appends to work->series_coeffs the radial series of the panel at level (lx, ly), from the
 * values of G in work->in_G: at each radius the angular integral over [ta, tb] by the
 * quadrature of the Lobatto points, then its Chebyshev coefficients in r. Returns the index
 * of its first coefficient.
 */
static gsize
_ncm_spectral_polar_series_append (NcmSpectralPolarWork *work, const NcmSpectralPolarPanel *panel, guint lx, guint ly)
{
  const guint n        = 1u << lx;
  const guint m        = 1u << ly;
  const gdouble half_t = 0.5 * (panel->tb - panel->ta);
  const gdouble *w     = &g_array_index (work->cc, gdouble, NCM_SPECTRAL_POLAR_LEVEL_OFFSET (ly));
  const gsize off      = work->series_coeffs->len;
  guint i, j;

  for (i = 0; i <= n; i++)
  {
    const gdouble *row = work->in_G + (gsize) i * (m + 1);
    gdouble s          = 0.0;

    for (j = 0; j <= m; j++)
      s += w[j] * row[j];

    work->sv[i] = half_t * s;
  }

  _ncm_spectral_transform (work->spectral, n + 1, 1, work->sv, work->sc);
  g_array_append_vals (work->series_coeffs, work->sc, n + 1);

  return off;
}

/*
 * Transforms the grid of level (lx, ly) into its Chebyshev coefficients in work->coeffs,
 * with the angular index first, and fills work->br and work->bt with the l1 masses of the
 * coefficients per radial and per angular index. Returns work->vs, the largest |F| on the
 * grid, which the functions that write the samples keep.
 */
static gdouble
_ncm_spectral_polar_coeffs (NcmSpectralPolarWork *work, guint lx, guint ly)
{
  const guint n = 1u << lx;
  const guint m = 1u << ly;
  guint j, k, l;

  /*
   * Radial transform of the m + 1 columns of in, its transpose, then angular transform
   * of the n + 1 columns of the transpose. The transforms leave in unchanged (the kernels
   * and the FFTW plans, out of place, preserve their input): the next doubling spreads
   * these samples in place.
   */
  _ncm_spectral_transform (work->spectral, n + 1, m + 1, work->in, work->tmp);

  for (k = 0; k <= n; k++)
  {
    const gdouble *row = work->tmp + (gsize) k * (m + 1);

    for (j = 0; j <= m; j++)
      work->tmpT[(gsize) j * (n + 1) + k] = row[j];
  }

  _ncm_spectral_transform (work->spectral, m + 1, n + 1, work->tmpT, work->coeffs);

  memset (work->br, 0, sizeof (gdouble) * (n + 1));

  for (l = 0; l <= m; l++)
  {
    const gdouble *row = work->coeffs + (gsize) l * (n + 1);
    gdouble acc        = 0.0;

    for (k = 0; k <= n; k++)
    {
      const gdouble ac = fabs (row[k]);

      acc         += ac;
      work->br[k] += ac;
    }

    work->bt[l] = acc;
  }

  return work->vs;
}

/*
 * The changes D_r and D_theta of the coefficients of level (lx, ly), from work->br and
 * work->bt; returns their sum.
 */
static gdouble
_ncm_spectral_polar_change (NcmSpectralPolarWork *work, guint lx, guint ly, gdouble *Dr, gdouble *Dt)
{
  const guint n = 1u << lx;
  const guint m = 1u << ly;
  gdouble dr    = 0.0;
  gdouble dt    = 0.0;
  guint k;

  for (k = n / 2 + 1; k <= n; k++)
    dr += 2.0 * work->br[k];

  for (k = m / 2 + 1; k <= m; k++)
    dt += 2.0 * work->bt[k];

  *Dr = dr;
  *Dt = dt;

  return dr + dt;
}

static const gdouble *_ncm_spectral_polar_moments (NcmSpectralPolarMom *mom, guint l);

/*
 * The value of the panel from work->coeffs at level (lx, ly): the coefficients against the
 * radial moments of the weight left out of F and the angular moments of weight 1.
 */
static gdouble
_ncm_spectral_polar_panel_value (NcmSpectralPolarCall *call, const NcmSpectralPolarPanel *panel, guint lx, guint ly)
{
  NcmSpectralPolarWork *work = call->work;
  const guint n              = 1u << lx;
  const guint m              = 1u << ly;
  const guint kind           = _ncm_spectral_polar_panel_kind (call, panel);
  const gdouble *mr          = _ncm_spectral_polar_moments (&work->mom[kind], lx);
  const gdouble *mt          = &g_array_index (work->mt, gdouble, NCM_SPECTRAL_POLAR_LEVEL_OFFSET (ly));
  const gdouble ex           = 1.0 + ((kind & 1) ? call->a : 0.0) + ((kind & 2) ? call->b : 0.0);
  gdouble acc                = 0.0;
  guint k, l;

  for (l = 0; l <= m; l += 2)
  {
    const gdouble *row = work->coeffs + (gsize) l * (n + 1);
    gdouble s          = 0.0;

    for (k = 0; k <= n; k++)
      s += row[k] * mr[k];

    acc += mt[l] * s;
  }

  return pow (0.5 * (panel->rb - panel->ra), ex) * 0.5 * (panel->tb - panel->ta) * acc;
}

/* Which power laws are left out of F on the panel: bit 0 for r^a, bit 1 for (R - r)^b */
static guint
_ncm_spectral_polar_panel_kind (const NcmSpectralPolarCall *call, const NcmSpectralPolarPanel *panel)
{
  return (((panel->ra == 0.0) && (call->a != 0.0)) ? 1 : 0) | (((panel->rb == call->R) && (call->b != 0.0)) ? 2 : 0);
}

/*
 * The level-l table of mom, computed on first use: 2 / (1 - k^2) for even k when both
 * exponents are 0, otherwise a Gauss-Jacobi rule of 2^(l - 1) + 2 nodes, exact for the
 * degree 2^l of T_k.
 */
static const gdouble *
_ncm_spectral_polar_moments (NcmSpectralPolarMom *mom, guint l)
{
  const guint n = 1u << l;
  gdouble *m    = &g_array_index (mom->m, gdouble, NCM_SPECTRAL_POLAR_LEVEL_OFFSET (l));
  guint k, q;

  if (mom->have & (1u << l))
    return m;

  memset (m, 0, sizeof (gdouble) * (n + 1));

  if ((mom->ca == 0.0) && (mom->cb == 0.0))
  {
    for (k = 0; k <= n; k += 2)
      m[k] = 2.0 / (1.0 - (gdouble) k * k);
  }
  else
  {
    const guint nq = n / 2 + 2;
    /* GSL's weight is (b - x)^alpha (x - a)^beta on [a, b] */
    gsl_integration_fixed_workspace *ws = gsl_integration_fixed_alloc (gsl_integration_fixed_jacobi, nq, -1.0, 1.0, mom->cb, mom->ca);
    const gdouble *xn                   = gsl_integration_fixed_nodes (ws);
    const gdouble *wn                   = gsl_integration_fixed_weights (ws);

    for (q = 0; q < nq; q++)
    {
      const gdouble x = xn[q];
      gdouble t0      = 1.0;
      gdouble t1      = x;

      m[0] += wn[q];
      m[1] += wn[q] * x;

      for (k = 2; k <= n; k++)
      {
        const gdouble t2 = 2.0 * x * t1 - t0;

        m[k] += wn[q] * t2;
        t0    = t1;
        t1    = t2;
      }
    }

    gsl_integration_fixed_free (ws);
  }

  mom->have |= 1u << l;

  return m;
}

/*
 * Accepts the initial panels that pass the test on their starting grids and expands the
 * others from those grids.
 */
static void
_ncm_spectral_polar_expand_initial (NcmSpectralPolarCall *call)
{
  NcmSpectralPolarWork *work = call->work;
  const guint ns             = (1u << call->level_min) + 1;
  gsize ib;

  for (ib = 0; ib < work->initial->len; ib++)
  {
    const NcmSpectralPolarEntry *entry = &g_array_index (work->initial, NcmSpectralPolarEntry, ib);

    if (entry->D <= MAX (call->reltol * entry->vs, call->abstol_F))
    {
      g_array_append_val (work->accepted, *entry);
    }
    else
    {
      const gdouble *start = &g_array_index (work->initial_grid, gdouble, ib * ns * ns);

      _ncm_spectral_polar_nodes_r (call, &entry->panel, call->level_min, 0, 1);
      _ncm_spectral_polar_nodes_theta (call, &entry->panel, call->level_min, 0, 1);
      memcpy (work->in, start, sizeof (gdouble) * ns * ns);

      if (call->keep_G)
        memcpy (work->in_G, &g_array_index (work->initial_grid_G, gdouble, ib * ns * ns), sizeof (gdouble) * ns * ns);

      work->vs = entry->vs;

      _ncm_spectral_polar_expand (call, entry->panel, entry);
    }
  }
}

/*
 * The grid of level (lx, ly), n = 2^lx and m = 2^ly, moved in place to its positions in
 * the grid of level (lx + 1, ly): row i to row 2 i, from the last row, so that no row is
 * overwritten before it is moved.
 */
static void
_ncm_spectral_polar_spread_r (gdouble *grid, guint n, guint m)
{
  const gsize ld = m + 1;
  guint i;

  for (i = n; i > 0; i--)
    memcpy (grid + 2 * i * ld, grid + i * ld, sizeof (gdouble) * ld);
}

/*
 * The grid of level (lx, ly) moved in place to its positions in the grid of level
 * (lx, ly + 1): point (i, j) from i (m + 1) + j to i (2 m + 1) + 2 j, from the last point,
 * since every point moves to a position at or after its own.
 */
static void
_ncm_spectral_polar_spread_theta (gdouble *grid, guint n, guint m)
{
  gint i, j;

  for (i = (gint) n; i >= 0; i--)
  {
    const gdouble *src = grid + (gsize) i * (m + 1);
    gdouble *dst       = grid + (gsize) i * (2 * m + 1);

    for (j = (gint) m; j >= 0; j--)
      dst[2 * j] = src[j];
  }
}

/*
 * Doubles the level of the panel in r (in_r) or in theta, sampling F at the new nodes. The
 * grid of G of a call that keeps it is moved as the grid of F.
 */
static void
_ncm_spectral_polar_double (NcmSpectralPolarCall *call, const NcmSpectralPolarPanel *panel, guint *lx, guint *ly, gboolean in_r)
{
  NcmSpectralPolarWork *work = call->work;

  if (in_r)
  {
    _ncm_spectral_polar_spread_r (work->in, 1u << *lx, 1u << *ly);

    if (call->keep_G)
      _ncm_spectral_polar_spread_r (work->in_G, 1u << *lx, 1u << *ly);

    (*lx)++;
    _ncm_spectral_polar_nodes_r (call, panel, *lx, 1, 2);
    _ncm_spectral_polar_rows (call, *lx, *ly, 1, 2, 0, 1);
  }
  else
  {
    _ncm_spectral_polar_spread_theta (work->in, 1u << *lx, 1u << *ly);

    if (call->keep_G)
      _ncm_spectral_polar_spread_theta (work->in_G, 1u << *lx, 1u << *ly);

    (*ly)++;
    _ncm_spectral_polar_nodes_theta (call, panel, *ly, 1, 2);
    _ncm_spectral_polar_rows (call, *lx, *ly, 0, 1, 1, 2);
  }
}

static gdouble _ncm_spectral_polar_predict (const gdouble *b, guint M);

/*
 * Expands F on the panel from the starting level until the test passes, doubling the
 * direction with the larger change or, at the level cap, bisecting the panel in it. An
 * accepted panel is appended to work->accepted. With start, an initial panel that failed its test in
 * _ncm_spectral_polar_first_pass(), the starting grid and its nodes are in the workspace
 * already, and so are its changes: its starting level is not transformed again. A panel at max_depth that needs a bisection sets
 * call->failed, after which nothing more is evaluated.
 */
static void
_ncm_spectral_polar_expand (NcmSpectralPolarCall *call, NcmSpectralPolarPanel panel, const NcmSpectralPolarEntry *start)
{
  NcmSpectralPolarWork *work = call->work;
  guint lx                   = call->level_min;
  guint ly                   = call->level_min;

  if (call->failed)
    return;

  if (start == NULL)
    _ncm_spectral_polar_sample_start (call, &panel);
  else if (call->level_min + 1 < call->level_cap)
    /*
     * The first pass of the loop below would transform the starting grid again and get
     * the D, scale and tolerance the initial panel failed with, so it would not accept it,
     * and below level-cap - 1 it would not split it either: it would double it, in the
     * direction of the larger change. Doing that directly gives the same panels.
     */
    _ncm_spectral_polar_double (call, &panel, &lx, &ly, start->Dr >= start->Dt);

  while (TRUE)
  {
    const gdouble vs  = _ncm_spectral_polar_coeffs (work, lx, ly);
    const gdouble tol = MAX (call->reltol * vs, call->abstol_F);
    gboolean split    = FALSE;
    gdouble Dr, Dt, D;

    D = _ncm_spectral_polar_change (work, lx, ly, &Dr, &Dt);

    if (D <= tol)
    {
      NcmSpectralPolarEntry entry = { panel, _ncm_spectral_polar_panel_value (call, &panel, lx, ly), D, vs, lx, ly, FALSE, Dr, Dt, 0 };

      if (call->keep_G)
        entry.series_off = _ncm_spectral_polar_series_append (work, &panel, lx, ly);

      g_array_append_val (work->accepted, entry);

      return;
    }

    if (Dr >= Dt)
      split = (lx >= call->level_cap) ||
              ((lx + 1 == call->level_cap) && (_ncm_spectral_polar_predict (work->br, 1u << lx) > tol));
    else
      split = (ly >= call->level_cap) ||
              ((ly + 1 == call->level_cap) && (_ncm_spectral_polar_predict (work->bt, 1u << ly) > tol));

    if (split)
    {
      NcmSpectralPolarPanel h1 = panel;
      NcmSpectralPolarPanel h2 = panel;

      if (panel.depth >= call->max_depth)
      {
        call->failed       = TRUE;
        call->failed_panel = panel;

        return;
      }

      h1.depth = h2.depth = panel.depth + 1;

      if (Dr >= Dt)
        h1.rb = h2.ra = 0.5 * (panel.ra + panel.rb);
      else
        h1.tb = h2.ta = 0.5 * (panel.ta + panel.tb);

      _ncm_spectral_polar_expand (call, h1, NULL);
      _ncm_spectral_polar_expand (call, h2, NULL);

      return;
    }

    _ncm_spectral_polar_double (call, &panel, &lx, &ly, Dr >= Dt);
  }
}

/*
 * The value the test of the next doubling is predicted to have, from the l1 masses
 * b[0..M] of the coefficients along one direction: with B_1 and B_2 the masses of the
 * bands (M/4, M/2] and (M/2, M], geometric decay gives B_2 / B_1 = q (1 + q) and a next
 * band of B_2 q^2 (1 + q^2), as _ncm_spectral_batch_cannot_converge() in ncm_spectral.c.
 */
static gdouble
_ncm_spectral_polar_predict (const gdouble *b, guint M)
{
  gdouble B1 = 0.0;
  gdouble B2 = 0.0;
  guint i;

  for (i = M / 4 + 1; i <= M / 2; i++)
    B1 += b[i];

  for (i = M / 2 + 1; i <= M; i++)
    B2 += b[i];

  if (B2 == 0.0)
    return 0.0;

  if (B1 == 0.0)
    return GSL_POSINF;

  {
    const gdouble ratio = B2 / B1;
    const gdouble q     = 0.5 * (sqrt (1.0 + 4.0 * ratio) - 1.0);

    return 2.0 * B2 * q * q * (1.0 + q * q);
  }
}

/*
 * Sums the accepted panels. While that sum gives a smaller abstol_F than the one the panels
 * were accepted with, the panels that fail the test with the smaller one are expanded
 * again from their starting level. Returns the final sum.
 */
static gdouble
_ncm_spectral_polar_retest (NcmSpectralPolarCall *call, gdouble abstol)
{
  NcmSpectralPolarWork *work = call->work;

  while (TRUE)
  {
    const guint n_accepted = work->accepted->len;
    gboolean expanded      = FALSE;
    gdouble sum            = 0.0;
    gdouble abstol_F;
    guint i;

    for (i = 0; i < n_accepted; i++)
      sum += g_array_index (work->accepted, NcmSpectralPolarEntry, i).val;

    abstol_F = MAX (call->reltol * fabs (sum), abstol) / call->A;

    if (call->failed || (abstol_F >= call->abstol_F))
      return sum;

    call->abstol_F = abstol_F;

    for (i = 0; i < n_accepted; i++)
    {
      NcmSpectralPolarEntry *entry = &g_array_index (work->accepted, NcmSpectralPolarEntry, i);

      if (entry->D > MAX (call->reltol * entry->vs, abstol_F))
      {
        const NcmSpectralPolarPanel panel = entry->panel;

        /* The expansion appends to the array, which may move it */
        entry->val      = 0.0;
        entry->D        = 0.0;
        entry->replaced = TRUE;
        expanded        = TRUE;

        _ncm_spectral_polar_expand (call, panel, NULL);
      }
    }

    if (!expanded)
      return sum;
  }
}

/*
 * The radial series of the accepted panels that the second pass did not replace, in the two
 * arrays of _ncm_spectral_polar_integrate_radial_series().
 */
static void
_ncm_spectral_polar_series (NcmSpectralPolarWork *work, GArray *series, GArray *coeffs)
{
  guint i;

  for (i = 0; i < work->accepted->len; i++)
  {
    const NcmSpectralPolarEntry *entry = &g_array_index (work->accepted, NcmSpectralPolarEntry, i);

    if (!entry->replaced)
    {
      const guint n_coeffs                 = (1u << entry->lx) + 1;
      const NcmSpectralPolarRadialSeries s = {
        entry->panel.ra, entry->panel.rb, entry->panel.ta, entry->panel.tb, n_coeffs, coeffs->len
      };

      g_array_append_val (series, s);
      g_array_append_vals (coeffs, &g_array_index (work->series_coeffs, gdouble, entry->series_off), n_coeffs);
    }
  }
}

/* The accepted panels as a matrix, one row per panel, see ncm_spectral_polar_integrate() */
static NcmMatrix *
_ncm_spectral_polar_panels (NcmSpectralPolarWork *work)
{
  NcmMatrix *panels;
  guint n_rows = 0;
  guint i, row;

  for (i = 0; i < work->accepted->len; i++)
  {
    if (!g_array_index (work->accepted, NcmSpectralPolarEntry, i).replaced)
      n_rows++;
  }

  panels = ncm_matrix_new (n_rows, 8);

  for (i = 0, row = 0; i < work->accepted->len; i++)
  {
    const NcmSpectralPolarEntry *entry = &g_array_index (work->accepted, NcmSpectralPolarEntry, i);

    if (entry->replaced)
      continue;

    ncm_matrix_set (panels, row, 0, entry->panel.ra);
    ncm_matrix_set (panels, row, 1, entry->panel.rb);
    ncm_matrix_set (panels, row, 2, entry->panel.ta);
    ncm_matrix_set (panels, row, 3, entry->panel.tb);
    ncm_matrix_set (panels, row, 4, entry->lx);
    ncm_matrix_set (panels, row, 5, entry->ly);
    ncm_matrix_set (panels, row, 6, entry->val);
    ncm_matrix_set (panels, row, 7, entry->D);
    row++;
  }

  return panels;
}

static gint _ncm_spectral_polar_cmp_double (gconstpointer a, gconstpointer b);

/**
 * ncm_spectral_polar_peak_edges:
 * @sp: a #NcmSpectralPolar
 * @R: radius of the disk
 * @r_peak: radius of the peak, in $[0, @R]$
 * @theta_peak: angle of the peak
 * @L_r: radial half width of the peak
 * @L_theta: angular half width of the peak
 * @r_edges: (element-type gdouble): array to fill with the radial edges
 * @theta_edges: (element-type gdouble): array to fill with the angular edges
 *
 * Fills @r_edges and @theta_edges with the edges of panels that cover the disk for an
 * integrand with one peak, to be passed to ncm_spectral_polar_integrate(). The radial
 * edges are 0, @R, $k @R / n$ for $k = 1, \dots, n - 1$ and
 * $@r_peak \pm @L_r \omega^\nu$, $\nu \ge 0$, inside $(0, @R)$, with
 * $n$ = #NcmSpectralPolar:min-radial-panels and $\omega$ = #NcmSpectralPolar:edge-ratio; after
 * sorting, an edge closer than $@L_r / 4$ to the one before is left out, except @R. The angular edges are
 * $@theta_peak \pm @L_theta \omega^\nu$ up to $@theta_peak \pm \pi$; when
 * $@L_theta \ge \pi / 4$ they are those of four equal panels with @theta_peak at the
 * center of one. @L_r and @L_theta are positive.
 */
void
ncm_spectral_polar_peak_edges (NcmSpectralPolar *sp, gdouble R, gdouble r_peak, gdouble theta_peak,
                               gdouble L_r, gdouble L_theta, GArray *r_edges, GArray *theta_edges)
{
  const gdouble rho = sp->edge_ratio;
  const gdouble end = R;
  const gdouble org = 0.0;
  gdouble w;
  guint i, m;

  if (!(R > 0.0) || !(r_peak >= 0.0) || !(r_peak <= R))
    g_error ("ncm_spectral_polar_peak_edges: the peak radius %g is not inside [0, %g].", r_peak, R);

  if (!(L_r > 0.0) || !(L_theta > 0.0))
    g_error ("ncm_spectral_polar_peak_edges: the half widths %g and %g must be positive.", L_r, L_theta);

  g_array_set_size (r_edges, 0);
  g_array_append_val (r_edges, org);
  g_array_append_val (r_edges, end);

  for (i = 1; i < sp->min_radial_panels; i++)
  {
    const gdouble x = R * i / sp->min_radial_panels;

    g_array_append_val (r_edges, x);
  }

  for (w = L_r; w < 2.0 * R; w *= rho)
  {
    const gdouble x1 = r_peak - w;
    const gdouble x2 = r_peak + w;

    if ((x1 > 0.0) && (x1 < R))
      g_array_append_val (r_edges, x1);

    if ((x2 > 0.0) && (x2 < R))
      g_array_append_val (r_edges, x2);
  }

  g_array_sort (r_edges, &_ncm_spectral_polar_cmp_double);

  for (i = 1, m = 1; i < r_edges->len; i++)
  {
    const gdouble x = g_array_index (r_edges, gdouble, i);

    if ((x - g_array_index (r_edges, gdouble, m - 1) > 0.25 * L_r) || (i + 1 == r_edges->len))
      g_array_index (r_edges, gdouble, m++) = x;
  }

  g_array_set_size (r_edges, m);

  g_array_set_size (theta_edges, 0);

  if (L_theta >= M_PI / NCM_SPECTRAL_POLAR_FULL_CIRCLE_PANELS)
  {
    for (i = 0; i <= NCM_SPECTRAL_POLAR_FULL_CIRCLE_PANELS; i++)
    {
      const gdouble t = theta_peak + M_PI * (2.0 * i - 1.0) / NCM_SPECTRAL_POLAR_FULL_CIRCLE_PANELS;

      g_array_append_val (theta_edges, t);
    }
  }
  else
  {
    const gdouble t_min = theta_peak - M_PI;
    const gdouble t_max = theta_peak + M_PI;
    guint n_below;

    /* theta_peak - w from the largest w down, then theta_peak + w from the smallest up */
    g_array_append_val (theta_edges, t_min);

    for (w = L_theta; M_PI - w > 0.25 * L_theta; w *= rho)
    {
      const gdouble t = theta_peak - w;

      g_array_append_val (theta_edges, t);
    }

    n_below = theta_edges->len;

    for (i = 1; 2 * i < n_below; i++)
    {
      const gdouble t = g_array_index (theta_edges, gdouble, i);

      g_array_index (theta_edges, gdouble, i)           = g_array_index (theta_edges, gdouble, n_below - i);
      g_array_index (theta_edges, gdouble, n_below - i) = t;
    }

    for (w = L_theta; M_PI - w > 0.25 * L_theta; w *= rho)
    {
      const gdouble t = theta_peak + w;

      g_array_append_val (theta_edges, t);
    }

    g_array_append_val (theta_edges, t_max);
  }
}

static gint
_ncm_spectral_polar_cmp_double (gconstpointer a, gconstpointer b)
{
  const gdouble x = *(const gdouble *) a;
  const gdouble y = *(const gdouble *) b;

  return (x > y) - (x < y);
}

