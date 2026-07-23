/***************************************************************************
 *            nc_data_planck_commander.c
 *
 *  Wed July 23 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_data_planck_commander.c
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
 * NcDataPlanckCommander:
 *
 * Planck low-ℓ commander (gibbs_gauss) temperature likelihood, native port.
 *
 * A Gaussianized Blackwell-Rao likelihood over the low-ℓ TT band. Each theory
 * $D_\ell = C_\ell \ell(\ell+1)/2\pi$ (with $C_\ell$ divided by the calibration
 * $A_\mathrm{planck}^2$) is mapped through a per-$\ell$ cubic spline to a
 * Gaussianized variable $x_\ell$, and
 * $$\ln L = -\tfrac12 (x-\mu)^\mathsf{T} P (x-\mu) + \sum_\ell \ln\frac{dx_\ell}{dD_\ell} - \mathrm{offset},$$
 * where $P$ is the inverse of the (band-limited) covariance and the offset fixes
 * $\ln L(\mu_\sigma)=0$. This #NcmData reports $-2\ln L$; unlike the clik wrapper
 * it can also resample.
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "nc/data/nc_data_planck_commander.h"
#include "nc/cmb/nc_planck_fi.h"
#include "ncm/core/ncm_serialize.h"

#include <math.h>

struct _NcDataPlanckCommander
{
  /*< private >*/
  NcmData parent_instance;
  NcHIPertBoltzmann *pb;
  guint lmin;
  guint lmax;
  guint nl; /* lmax - lmin + 1 */
  guint nbin;
  guint delta_l;
  gchar *calib_name;
  NcmVector *cl2x_xa;  /* nl*nbin: per-ell spline Cl (Dl) abscissa */
  NcmVector *cl2x_ya;  /* nl*nbin: per-ell spline x ordinate */
  NcmVector *cl2x_y2a; /* nl*nbin: per-ell spline 2nd derivatives */
  NcmVector *mu;       /* nl: Gaussianized mean */
  NcmMatrix *cov;      /* nl x nl: covariance (inverted at construction) */
  NcmVector *mu_sigma; /* nl: mean sigma_l (Dl), sets the offset */
  /* runtime (constructed) */
  NcmMatrix *P;      /* nl x nl: inverse of band-limited cov */
  NcmMatrix *cov_ll; /* nl x nl: Cholesky (lower) of band-limited cov, for resample */
  gdouble offset;
  gdouble *prior_lo; /* nl */
  gdouble *prior_hi; /* nl */
  gdouble *xa;       /* cached data pointers (stride-1) */
  gdouble *ya;
  gdouble *y2a;
  gdouble *x;    /* nl scratch */
  gdouble *d;    /* nl scratch */
  gint calib_pi; /* resolved calib param index (-1 => unresolved) */
  NcmVector *cl_TT;
};

enum
{
  PROP_0,
  PROP_PB,
  PROP_LMIN,
  PROP_LMAX,
  PROP_NBIN,
  PROP_DELTA_L,
  PROP_CALIB_NAME,
  PROP_CL2X_XA,
  PROP_CL2X_YA,
  PROP_CL2X_Y2A,
  PROP_MU,
  PROP_COV,
  PROP_MU_SIGMA,
  PROP_SIZE,
};

G_DEFINE_TYPE (NcDataPlanckCommander, nc_data_planck_commander, NCM_TYPE_DATA)

static void
nc_data_planck_commander_init (NcDataPlanckCommander *cmd)
{
  cmd->pb         = NULL;
  cmd->lmin       = 0;
  cmd->lmax       = 0;
  cmd->nl         = 0;
  cmd->nbin       = 0;
  cmd->delta_l    = 0;
  cmd->calib_name = NULL;
  cmd->cl2x_xa    = NULL;
  cmd->cl2x_ya    = NULL;
  cmd->cl2x_y2a   = NULL;
  cmd->mu         = NULL;
  cmd->cov        = NULL;
  cmd->mu_sigma   = NULL;
  cmd->P          = NULL;
  cmd->cov_ll     = NULL;
  cmd->offset     = 0.0;
  cmd->prior_lo   = NULL;
  cmd->prior_hi   = NULL;
  cmd->xa         = NULL;
  cmd->ya         = NULL;
  cmd->y2a        = NULL;
  cmd->x          = NULL;
  cmd->d          = NULL;
  cmd->calib_pi   = -1;
  cmd->cl_TT      = NULL;
}

/* Locate klo in [0, n-2] with xa[klo] <= x < xa[klo+1] (clamped). */
static guint
commander_locate (const gdouble *xa, guint n, gdouble x)
{
  guint lo = 0, hi = n - 1;

  while (hi - lo > 1)
  {
    const guint mid = (lo + hi) / 2;

    if (xa[mid] <= x)
      lo = mid;
    else
      hi = mid;
  }

  return lo;
}

/* Natural cubic-spline evaluation using the stored second derivatives; also
 * returns the derivative dx/dDl when @deriv is non-NULL. */
static gdouble
commander_splint (const gdouble *xa, const gdouble *ya, const gdouble *y2a, guint n, gdouble x, gdouble *deriv)
{
  const guint klo = commander_locate (xa, n, x);
  const guint khi = klo + 1;
  const gdouble h = xa[khi] - xa[klo];
  const gdouble a = (xa[khi] - x) / h;
  const gdouble b = (x - xa[klo]) / h;

  if (deriv != NULL)
    *deriv = (ya[khi] - ya[klo]) / h -
             (3.0 * a * a - 1.0) / 6.0 * h * y2a[klo] +
             (3.0 * b * b - 1.0) / 6.0 * h * y2a[khi];

  return a * ya[klo] + b * ya[khi] +
         ((a * a * a - a) * y2a[klo] + (b * b * b - b) * y2a[khi]) * h * h / 6.0;
}

/* lnL without the offset subtraction for a band-power vector @dl (length nl,
 * in Dl units). Sets *bad when a prior/monotonicity bound is violated. */
static gdouble
commander_lnL_raw (NcDataPlanckCommander *cmd, const gdouble *dl, gboolean *bad)
{
  const guint nl = cmd->nl;
  gdouble jac    = 0.0;
  gdouble chi2   = 0.0;
  guint i, j;

  *bad = FALSE;

  for (i = 0; i < nl; i++)
  {
    const gdouble *xa  = cmd->xa + i * cmd->nbin;
    const gdouble *ya  = cmd->ya + i * cmd->nbin;
    const gdouble *y2a = cmd->y2a + i * cmd->nbin;
    gdouble deriv;

    if ((dl[i] < cmd->prior_lo[i]) || (dl[i] > cmd->prior_hi[i]))
    {
      *bad = TRUE;

      return 0.0;
    }

    cmd->x[i] = commander_splint (xa, ya, y2a, cmd->nbin, dl[i], &deriv);

    if (deriv <= 0.0)
    {
      *bad = TRUE;

      return 0.0;
    }

    jac += log (deriv);
  }

  for (i = 0; i < nl; i++)
    cmd->d[i] = cmd->x[i] - ncm_vector_get (cmd->mu, i);

  for (i = 0; i < nl; i++)
  {
    gdouble u = 0.0;

    for (j = 0; j < nl; j++)
      u += ncm_matrix_get (cmd->P, i, j) * cmd->d[j];

    chi2 += cmd->d[i] * u;
  }

  return -0.5 * chi2 + jac;
}

static void
nc_data_planck_commander_constructed (GObject *object)
{
  NcDataPlanckCommander *cmd = NC_DATA_PLANCK_COMMANDER (object);
  guint i, j;
  gint ret;
  gboolean bad;

  g_assert_cmpuint (cmd->lmax, >=, cmd->lmin);
  g_assert_cmpuint (cmd->nbin, >, 1);

  cmd->nl = cmd->lmax - cmd->lmin + 1;

  g_assert_nonnull (cmd->cl2x_xa);
  g_assert_nonnull (cmd->cl2x_ya);
  g_assert_nonnull (cmd->cl2x_y2a);
  g_assert_cmpuint (ncm_vector_len (cmd->cl2x_xa), ==, cmd->nl * cmd->nbin);
  g_assert_cmpuint (ncm_vector_len (cmd->cl2x_ya), ==, cmd->nl * cmd->nbin);
  g_assert_cmpuint (ncm_vector_len (cmd->cl2x_y2a), ==, cmd->nl * cmd->nbin);
  g_assert_nonnull (cmd->mu);
  g_assert_nonnull (cmd->mu_sigma);
  g_assert_cmpuint (ncm_vector_len (cmd->mu), ==, cmd->nl);
  g_assert_cmpuint (ncm_vector_len (cmd->mu_sigma), ==, cmd->nl);
  g_assert_nonnull (cmd->cov);
  g_assert_cmpuint (ncm_matrix_nrows (cmd->cov), ==, cmd->nl);
  g_assert_cmpuint (ncm_matrix_ncols (cmd->cov), ==, cmd->nl);

  cmd->xa  = ncm_vector_data (cmd->cl2x_xa);
  cmd->ya  = ncm_vector_data (cmd->cl2x_ya);
  cmd->y2a = ncm_vector_data (cmd->cl2x_y2a);
  g_assert_cmpuint (ncm_vector_stride (cmd->cl2x_xa), ==, 1);
  g_assert_cmpuint (ncm_vector_stride (cmd->cl2x_ya), ==, 1);
  g_assert_cmpuint (ncm_vector_stride (cmd->cl2x_y2a), ==, 1);

  cmd->x        = g_new0 (gdouble, cmd->nl);
  cmd->d        = g_new0 (gdouble, cmd->nl);
  cmd->prior_lo = g_new0 (gdouble, cmd->nl);
  cmd->prior_hi = g_new0 (gdouble, cmd->nl);
  cmd->cl_TT    = ncm_vector_new (cmd->lmax + 1);

  /* Band-limit the covariance (|i-j| > delta_l => 0), keep a Cholesky for
   * resampling, and invert it (P) for the likelihood quadratic form. */
  {
    NcmMatrix *cov_bl = ncm_matrix_new (cmd->nl, cmd->nl);

    for (i = 0; i < cmd->nl; i++)
      for (j = 0; j < cmd->nl; j++)
      {
        const gdouble v = ((guint) ABS ((gint) i - (gint) j) > cmd->delta_l) ? 0.0 : ncm_matrix_get (cmd->cov, i, j);

        ncm_matrix_set (cov_bl, i, j, v);
      }

    cmd->cov_ll = ncm_matrix_dup (cov_bl);
    ret         = ncm_matrix_cholesky_decomp (cmd->cov_ll, 'L');

    if (ret != 0)
      g_error ("nc_data_planck_commander_constructed: covariance not positive definite (%d).", ret);

    cmd->P = ncm_matrix_dup (cov_bl);
    ret    = ncm_matrix_cholesky_decomp (cmd->P, 'L');

    if (ret != 0)
      g_error ("nc_data_planck_commander_constructed: covariance Cholesky failed (%d).", ret);

    ret = ncm_matrix_cholesky_inverse (cmd->P, 'L');

    if (ret != 0)
      g_error ("nc_data_planck_commander_constructed: covariance inversion failed (%d).", ret);

    /* Symmetrize the lower-only inverse into the full matrix. */
    for (i = 0; i < cmd->nl; i++)
      for (j = i + 1; j < cmd->nl; j++)
        ncm_matrix_set (cmd->P, i, j, ncm_matrix_get (cmd->P, j, i));

    ncm_matrix_free (cov_bl);
  }

  /* Priors: two nodes inside where the tabulated x crosses -/+5. */
  for (i = 0; i < cmd->nl; i++)
  {
    const gdouble *xa = cmd->xa + i * cmd->nbin;
    const gdouble *ya = cmd->ya + i * cmd->nbin;

    j = 0;

    while (fabs (ya[j] + 5.0) < 1.0e-4)
      j++;

    cmd->prior_lo[i] = xa[j + 2];

    j = cmd->nbin - 1;

    while (fabs (ya[j] - 5.0) < 1.0e-4)
      j--;

    cmd->prior_hi[i] = xa[j - 2];
  }

  /* Offset so that lnL(mu_sigma) = 0. */
  {
    gdouble *ms = ncm_vector_data (cmd->mu_sigma);

    cmd->offset = 0.0;
    cmd->offset = commander_lnL_raw (cmd, ms, &bad);

    if (bad)
      g_error ("nc_data_planck_commander_constructed: mu_sigma outside the tabulated range.");
  }

  /* Chain up : end */
  G_OBJECT_CLASS (nc_data_planck_commander_parent_class)->constructed (object);
}

static void
nc_data_planck_commander_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcDataPlanckCommander *cmd = NC_DATA_PLANCK_COMMANDER (object);

  g_return_if_fail (NC_IS_DATA_PLANCK_COMMANDER (object));

  switch (prop_id)
  {
    case PROP_PB:
      nc_data_planck_commander_set_hipert_boltzmann (cmd, g_value_get_object (value));
      break;
    case PROP_LMIN:
      cmd->lmin = g_value_get_uint (value);
      break;
    case PROP_LMAX:
      cmd->lmax = g_value_get_uint (value);
      break;
    case PROP_NBIN:
      cmd->nbin = g_value_get_uint (value);
      break;
    case PROP_DELTA_L:
      cmd->delta_l = g_value_get_uint (value);
      break;
    case PROP_CALIB_NAME:
      g_free (cmd->calib_name);
      cmd->calib_name = g_value_dup_string (value);
      break;
    case PROP_CL2X_XA:
      ncm_vector_substitute (&cmd->cl2x_xa, g_value_get_object (value), TRUE);
      break;
    case PROP_CL2X_YA:
      ncm_vector_substitute (&cmd->cl2x_ya, g_value_get_object (value), TRUE);
      break;
    case PROP_CL2X_Y2A:
      ncm_vector_substitute (&cmd->cl2x_y2a, g_value_get_object (value), TRUE);
      break;
    case PROP_MU:
      ncm_vector_substitute (&cmd->mu, g_value_get_object (value), TRUE);
      break;
    case PROP_COV:
      ncm_matrix_substitute (&cmd->cov, g_value_get_object (value), TRUE);
      break;
    case PROP_MU_SIGMA:
      ncm_vector_substitute (&cmd->mu_sigma, g_value_get_object (value), TRUE);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
nc_data_planck_commander_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcDataPlanckCommander *cmd = NC_DATA_PLANCK_COMMANDER (object);

  g_return_if_fail (NC_IS_DATA_PLANCK_COMMANDER (object));

  switch (prop_id)
  {
    case PROP_PB:
      g_value_set_object (value, cmd->pb);
      break;
    case PROP_LMIN:
      g_value_set_uint (value, cmd->lmin);
      break;
    case PROP_LMAX:
      g_value_set_uint (value, cmd->lmax);
      break;
    case PROP_NBIN:
      g_value_set_uint (value, cmd->nbin);
      break;
    case PROP_DELTA_L:
      g_value_set_uint (value, cmd->delta_l);
      break;
    case PROP_CALIB_NAME:
      g_value_set_string (value, cmd->calib_name);
      break;
    case PROP_CL2X_XA:
      g_value_set_object (value, cmd->cl2x_xa);
      break;
    case PROP_CL2X_YA:
      g_value_set_object (value, cmd->cl2x_ya);
      break;
    case PROP_CL2X_Y2A:
      g_value_set_object (value, cmd->cl2x_y2a);
      break;
    case PROP_MU:
      g_value_set_object (value, cmd->mu);
      break;
    case PROP_COV:
      g_value_set_object (value, cmd->cov);
      break;
    case PROP_MU_SIGMA:
      g_value_set_object (value, cmd->mu_sigma);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
nc_data_planck_commander_dispose (GObject *object)
{
  NcDataPlanckCommander *cmd = NC_DATA_PLANCK_COMMANDER (object);

  nc_hipert_boltzmann_clear (&cmd->pb);
  ncm_vector_clear (&cmd->cl2x_xa);
  ncm_vector_clear (&cmd->cl2x_ya);
  ncm_vector_clear (&cmd->cl2x_y2a);
  ncm_vector_clear (&cmd->mu);
  ncm_matrix_clear (&cmd->cov);
  ncm_vector_clear (&cmd->mu_sigma);
  ncm_matrix_clear (&cmd->P);
  ncm_matrix_clear (&cmd->cov_ll);
  ncm_vector_clear (&cmd->cl_TT);

  /* Chain up : end */
  G_OBJECT_CLASS (nc_data_planck_commander_parent_class)->dispose (object);
}

static void
nc_data_planck_commander_finalize (GObject *object)
{
  NcDataPlanckCommander *cmd = NC_DATA_PLANCK_COMMANDER (object);

  g_clear_pointer (&cmd->calib_name, g_free);
  g_clear_pointer (&cmd->prior_lo, g_free);
  g_clear_pointer (&cmd->prior_hi, g_free);
  g_clear_pointer (&cmd->x, g_free);
  g_clear_pointer (&cmd->d, g_free);

  /* Chain up : end */
  G_OBJECT_CLASS (nc_data_planck_commander_parent_class)->finalize (object);
}

static guint _nc_data_planck_commander_get_length (NcmData *data);
static void _nc_data_planck_commander_prepare (NcmData *data, NcmMSet *mset);
static void _nc_data_planck_commander_m2lnL_val (NcmData *data, NcmMSet *mset, gdouble *m2lnL);
static void _nc_data_planck_commander_resample (NcmData *data, NcmMSet *mset, NcmRNG *rng);

static void
nc_data_planck_commander_class_init (NcDataPlanckCommanderClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);
  NcmDataClass *data_class   = NCM_DATA_CLASS (klass);

  object_class->set_property = &nc_data_planck_commander_set_property;
  object_class->get_property = &nc_data_planck_commander_get_property;
  object_class->constructed  = &nc_data_planck_commander_constructed;
  object_class->dispose      = &nc_data_planck_commander_dispose;
  object_class->finalize     = &nc_data_planck_commander_finalize;

  g_object_class_install_property (object_class,
                                   PROP_PB,
                                   g_param_spec_object ("hipert-boltzmann",
                                                        NULL,
                                                        "Perturbations (Cls source)",
                                                        NC_TYPE_HIPERT_BOLTZMANN,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_LMIN,
                                   g_param_spec_uint ("lmin",
                                                      NULL,
                                                      "Minimum multipole",
                                                      0, G_MAXUINT, 2,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_LMAX,
                                   g_param_spec_uint ("lmax",
                                                      NULL,
                                                      "Maximum multipole",
                                                      0, G_MAXUINT, 29,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_NBIN,
                                   g_param_spec_uint ("nbin",
                                                      NULL,
                                                      "Number of spline nodes per multipole",
                                                      0, G_MAXUINT, 1000,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_DELTA_L,
                                   g_param_spec_uint ("delta-l",
                                                      NULL,
                                                      "Covariance band-limit width",
                                                      0, G_MAXUINT, 1000,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_CALIB_NAME,
                                   g_param_spec_string ("calib-name",
                                                        NULL,
                                                        "Free-calibration parameter name (NULL => none)",
                                                        "A_planck",
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_CL2X_XA,
                                   g_param_spec_object ("cl2x-xa",
                                                        NULL,
                                                        "Per-ell spline Dl abscissa (flattened nl*nbin)",
                                                        NCM_TYPE_VECTOR,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_CL2X_YA,
                                   g_param_spec_object ("cl2x-ya",
                                                        NULL,
                                                        "Per-ell spline x ordinate (flattened nl*nbin)",
                                                        NCM_TYPE_VECTOR,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_CL2X_Y2A,
                                   g_param_spec_object ("cl2x-y2a",
                                                        NULL,
                                                        "Per-ell spline second derivatives (flattened nl*nbin)",
                                                        NCM_TYPE_VECTOR,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_MU,
                                   g_param_spec_object ("mu",
                                                        NULL,
                                                        "Gaussianized mean vector",
                                                        NCM_TYPE_VECTOR,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_COV,
                                   g_param_spec_object ("cov",
                                                        NULL,
                                                        "Covariance matrix (band-limited and inverted internally)",
                                                        NCM_TYPE_MATRIX,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_MU_SIGMA,
                                   g_param_spec_object ("mu-sigma",
                                                        NULL,
                                                        "Mean sigma_l (Dl) defining the offset normalization",
                                                        NCM_TYPE_VECTOR,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  data_class->get_length = &_nc_data_planck_commander_get_length;
  data_class->prepare    = &_nc_data_planck_commander_prepare;
  data_class->m2lnL_val  = &_nc_data_planck_commander_m2lnL_val;
  data_class->resample   = &_nc_data_planck_commander_resample;
}

static guint
_nc_data_planck_commander_get_length (NcmData *data)
{
  NcDataPlanckCommander *cmd = NC_DATA_PLANCK_COMMANDER (data);

  return cmd->nl;
}

static void
_nc_data_planck_commander_prepare (NcmData *data, NcmMSet *mset)
{
  NcDataPlanckCommander *cmd = NC_DATA_PLANCK_COMMANDER (data);
  NcHICosmo *cosmo           = NC_HICOSMO (ncm_mset_peek (mset, nc_hicosmo_id ()));

  if (cmd->pb == NULL)
    g_error ("_nc_data_planck_commander_prepare: no NcHIPertBoltzmann set.");

  if (cosmo == NULL)
    g_error ("_nc_data_planck_commander_prepare: no NcHICosmo in mset.");

  nc_hipert_boltzmann_prepare_if_needed (cmd->pb, cosmo);
}

/* Fill @dl (length nl, Dl units) from the theory Cls, applying the calibration. */
static void
commander_model_dl (NcDataPlanckCommander *cmd, NcmMSet *mset, gdouble *dl)
{
  const gdouble twopi = 2.0 * M_PI;
  gdouble calib2      = 1.0;
  guint l;

  nc_hipert_boltzmann_get_TT_Cls (cmd->pb, cmd->cl_TT);

  if (cmd->calib_name != NULL)
  {
    NcPlanckFI *pfi = NC_PLANCK_FI (ncm_mset_peek (mset, nc_planck_fi_id ()));
    gdouble cal;

    if (pfi == NULL)
      g_error ("commander_model_dl: calibration '%s' requested but no NcPlanckFI in mset.", cmd->calib_name);

    if (cmd->calib_pi < 0)
    {
      guint idx;

      if (!ncm_model_param_index_from_name (NCM_MODEL (pfi), cmd->calib_name, &idx, NULL))
        g_error ("commander_model_dl: calibration parameter '%s' not in model.", cmd->calib_name);

      cmd->calib_pi = idx;
    }

    cal    = ncm_model_param_get (NCM_MODEL (pfi), cmd->calib_pi);
    calib2 = cal * cal;
  }

  for (l = cmd->lmin; l <= cmd->lmax; l++)
    dl[l - cmd->lmin] = ncm_vector_get (cmd->cl_TT, l) / calib2 * (l * (l + 1.0)) / twopi;
}

static void
_nc_data_planck_commander_m2lnL_val (NcmData *data, NcmMSet *mset, gdouble *m2lnL)
{
  NcDataPlanckCommander *cmd = NC_DATA_PLANCK_COMMANDER (data);
  gdouble *dl                = g_newa (gdouble, cmd->nl);
  gboolean bad;
  gdouble lnL;

  commander_model_dl (cmd, mset, dl);
  lnL = commander_lnL_raw (cmd, dl, &bad);

  if (bad)
    *m2lnL = 1.0e30;
  else
    *m2lnL = -2.0 * (lnL - cmd->offset);
}

static void
_nc_data_planck_commander_resample (NcmData *data, NcmMSet *mset, NcmRNG *rng)
{
  NcDataPlanckCommander *cmd = NC_DATA_PLANCK_COMMANDER (data);
  gdouble *dl                = g_newa (gdouble, cmd->nl);
  gboolean bad;
  guint i, j;

  /* Draw x_obs ~ N(x(theory), cov); the Gaussianized mean mu becomes the mock
   * observation, and the offset is recomputed so lnL(mu_sigma) stays 0. */
  commander_model_dl (cmd, mset, dl);

  for (i = 0; i < cmd->nl; i++)
  {
    gdouble deriv;
    const gdouble *xa  = cmd->xa + i * cmd->nbin;
    const gdouble *ya  = cmd->ya + i * cmd->nbin;
    const gdouble *y2a = cmd->y2a + i * cmd->nbin;

    cmd->x[i] = commander_splint (xa, ya, y2a, cmd->nbin, dl[i], &deriv);
    cmd->d[i] = ncm_rng_gaussian_gen (rng, 0.0, 1.0);
  }

  /* mu = x(theory) + L z, with L the lower Cholesky of the band-limited cov. */
  for (i = 0; i < cmd->nl; i++)
  {
    gdouble s = cmd->x[i];

    for (j = 0; j <= i; j++)
      s += ncm_matrix_get (cmd->cov_ll, i, j) * cmd->d[j];

    ncm_vector_set (cmd->mu, i, s);
  }

  {
    gdouble *ms = ncm_vector_data (cmd->mu_sigma);

    cmd->offset = 0.0;
    cmd->offset = commander_lnL_raw (cmd, ms, &bad);
  }
}

/**
 * nc_data_planck_commander_new:
 * @lmin: lowest multipole
 * @lmax: highest multipole
 * @nbin: number of spline nodes per multipole
 * @delta_l: covariance band-limit width
 * @cl2x_xa: (transfer none): per-ell spline Dl abscissa (flattened nl*nbin)
 * @cl2x_ya: (transfer none): per-ell spline x ordinate (flattened nl*nbin)
 * @cl2x_y2a: (transfer none): per-ell spline second derivatives (flattened nl*nbin)
 * @mu: (transfer none): Gaussianized mean vector (length nl)
 * @cov: (transfer none): covariance matrix (nl x nl), inverted internally
 * @mu_sigma: (transfer none): mean sigma_l defining the offset (length nl)
 *
 * Builds a fully-specified #NcDataPlanckCommander (calibration defaults to
 * A_planck). The theory $C_\ell$ source is set separately.
 *
 * Returns: (transfer full): the new #NcDataPlanckCommander.
 */
NcDataPlanckCommander *
nc_data_planck_commander_new (guint     lmin,
                              guint     lmax,
                              guint     nbin,
                              guint     delta_l,
                              NcmVector *cl2x_xa,
                              NcmVector *cl2x_ya,
                              NcmVector *cl2x_y2a,
                              NcmVector *mu,
                              NcmMatrix *cov,
                              NcmVector *mu_sigma)
{
  return g_object_new (NC_TYPE_DATA_PLANCK_COMMANDER,
                       "lmin", lmin,
                       "lmax", lmax,
                       "nbin", nbin,
                       "delta-l", delta_l,
                       "cl2x-xa", cl2x_xa,
                       "cl2x-ya", cl2x_ya,
                       "cl2x-y2a", cl2x_y2a,
                       "mu", mu,
                       "cov", cov,
                       "mu-sigma", mu_sigma,
                       NULL);
}

/**
 * nc_data_planck_commander_new_from_file:
 * @filename: a serialized #NcDataPlanckCommander file.
 *
 * Returns: (transfer full): the new #NcDataPlanckCommander.
 */
NcDataPlanckCommander *
nc_data_planck_commander_new_from_file (const gchar *filename)
{
  NcDataPlanckCommander *cmd = NC_DATA_PLANCK_COMMANDER (ncm_serialize_global_from_file (filename));

  g_assert (NC_IS_DATA_PLANCK_COMMANDER (cmd));

  return cmd;
}

/**
 * nc_data_planck_commander_set_hipert_boltzmann:
 * @cmd: a #NcDataPlanckCommander
 * @pb: a #NcHIPertBoltzmann
 *
 * Sets the theory $C_\ell$ source.
 */
void
nc_data_planck_commander_set_hipert_boltzmann (NcDataPlanckCommander *cmd, NcHIPertBoltzmann *pb)
{
  nc_hipert_boltzmann_clear (&cmd->pb);

  if (pb != NULL)
    cmd->pb = nc_hipert_boltzmann_ref (pb);
}

/**
 * nc_data_planck_commander_peek_hipert_boltzmann:
 * @cmd: a #NcDataPlanckCommander
 *
 * Returns: (transfer none): the perturbations object, or %NULL.
 */
NcHIPertBoltzmann *
nc_data_planck_commander_peek_hipert_boltzmann (NcDataPlanckCommander *cmd)
{
  return cmd->pb;
}

