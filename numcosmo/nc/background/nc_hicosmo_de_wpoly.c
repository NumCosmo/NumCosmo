/***************************************************************************
 *            nc_hicosmo_de_wpoly.c
 *
 *  Mon Sep 28 10:00:00 2026
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
 * NcHICosmoDEWPoly:
 *
 * Dark energy equation of state as a Legendre polynomial in redshift.
 *
 * On $0 \le z \le z_f$, with $x = 2z/z_f - 1$,
 * $$w(z) = -1 + \sum_{\ell=0}^{L} a_\ell P_\ell(x),$$
 * so $a_\ell = 0$ is $\Lambda$CDM and $a_0$ is the mean of $w + 1$ over the range.
 * For $z > z_f$, $w$ is held at $w(z_f)$. The dark energy density
 * $\Omega_{de}(z) = \Omega_x \exp\left[3\int_0^z \frac{1+w(z')}{1+z'}\mathrm{d}z'\right]$
 * is evaluated by fixed-order Gauss-Legendre quadrature whose order is set at
 * construction from $L$ and $z_f$ so that the rule is exact to double precision.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "nc/background/nc_hicosmo_de_wpoly.h"
#include "nc_enum_types.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_sort.h>
#endif /* NUMCOSMO_GIR_SCAN */

struct _NcHICosmoDEWPolyPrivate
{
  guint ncoeffs;
  guint lmax;
  gdouble z_f;
  gdouble alpha_f;
  gdouble w_f;
  gdouble int_f;
  gsize gl_n;
  gdouble *gl_x;
  gdouble *gl_w;
};

G_DEFINE_TYPE_WITH_PRIVATE (NcHICosmoDEWPoly, nc_hicosmo_de_wpoly, NC_TYPE_HICOSMO_DE)

enum
{
  PROP_0,
  PROP_Z_F,
  PROP_SIZE,
};

static void
nc_hicosmo_de_wpoly_init (NcHICosmoDEWPoly *wpoly)
{
  NcHICosmoDEWPolyPrivate * const self = wpoly->priv = nc_hicosmo_de_wpoly_get_instance_private (wpoly);

  self->ncoeffs = 0;
  self->lmax    = 0;
  self->z_f     = 0.0;
  self->alpha_f = 0.0;
  self->w_f     = 0.0;
  self->int_f   = 0.0;
  self->gl_n    = 0;
  self->gl_x    = NULL;
  self->gl_w    = NULL;
}

static void
_nc_hicosmo_de_wpoly_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcHICosmoDEWPoly *wpoly              = NC_HICOSMO_DE_WPOLY (object);
  NcHICosmoDEWPolyPrivate * const self = wpoly->priv;

  switch (prop_id)
  {
    case PROP_Z_F:
      g_value_set_double (value, self->z_f);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_nc_hicosmo_de_wpoly_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcHICosmoDEWPoly *wpoly              = NC_HICOSMO_DE_WPOLY (object);
  NcHICosmoDEWPolyPrivate * const self = wpoly->priv;

  switch (prop_id)
  {
    case PROP_Z_F:
      self->z_f = g_value_get_double (value);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_nc_hicosmo_de_wpoly_constructed (GObject *object)
{
  /* Chain up : start */
  G_OBJECT_CLASS (nc_hicosmo_de_wpoly_parent_class)->constructed (object);
  {
    NcHICosmoDEWPoly *wpoly              = NC_HICOSMO_DE_WPOLY (object);
    NcHICosmoDEWPolyPrivate * const self = wpoly->priv;
    NcmModel *model                      = NCM_MODEL (wpoly);

    self->ncoeffs = ncm_model_vparam_len (model, NC_HICOSMO_DE_WPOLY_A);
    self->lmax    = self->ncoeffs - 1;
    self->alpha_f = log1p (self->z_f);

    g_assert_cmpuint (self->ncoeffs, >, 0);

    {
      /* The integrand w(z)/(1+z) on [0, z] is analytic with its nearest
       * singularity at z = -1. On the Bernstein ellipse of [0, z_f] that
       * singularity sits at parameter rho, and an n-point Gauss-Legendre rule
       * has error of order rho^(L - 2n); n is chosen so that this is below
       * double precision, plus a margin for the prefactor.
       */
      const gdouble x_s = -1.0 - 2.0 / self->z_f;
      const gdouble rho = -x_s + sqrt (x_s * x_s - 1.0);
      const gdouble n   = 0.5 * (self->lmax - log (GSL_DBL_EPSILON) / log (rho));

      self->gl_n = (gsize) ceil (n) + 4;
    }

    {
      /* Nodes from GSL, weights from 2 / ((1 - x^2) P_n'(x)^2): the weights
       * GSL tabulates for some n are accurate only to 1e-12.
       */
      gsl_integration_glfixed_table *table = gsl_integration_glfixed_table_alloc (self->gl_n);
      gdouble *Pl                          = g_new (gdouble, self->gl_n + 1);
      gdouble *dPl                         = g_new (gdouble, self->gl_n + 1);
      gsize i;

      self->gl_x = g_new (gdouble, self->gl_n);
      self->gl_w = g_new (gdouble, self->gl_n);

      for (i = 0; i < self->gl_n; i++)
      {
        gdouble xi, wi;

        gsl_integration_glfixed_point (-1.0, 1.0, i, &xi, &wi, table);
        gsl_sf_legendre_Pl_deriv_array (self->gl_n, xi, Pl, dPl);

        self->gl_x[i] = xi;
        self->gl_w[i] = 2.0 / ((1.0 - xi * xi) * dPl[self->gl_n] * dPl[self->gl_n]);
      }

      g_free (Pl);
      g_free (dPl);
      gsl_integration_glfixed_table_free (table);
    }
  }
}

static void
_nc_hicosmo_de_wpoly_finalize (GObject *object)
{
  NcHICosmoDEWPoly *wpoly              = NC_HICOSMO_DE_WPOLY (object);
  NcHICosmoDEWPolyPrivate * const self = wpoly->priv;

  g_clear_pointer (&self->gl_x, g_free);
  g_clear_pointer (&self->gl_w, g_free);

  /* Chain up : end */
  G_OBJECT_CLASS (nc_hicosmo_de_wpoly_parent_class)->finalize (object);
}

static gdouble _nc_hicosmo_de_wpoly_E2Omega_de (NcHICosmoDE *cosmo_de, gdouble z);
static gdouble _nc_hicosmo_de_wpoly_dE2Omega_de_dz (NcHICosmoDE *cosmo_de, gdouble z);
static gdouble _nc_hicosmo_de_wpoly_d2E2Omega_de_dz2 (NcHICosmoDE *cosmo_de, gdouble z);
static gdouble _nc_hicosmo_de_wpoly_w_de (NcHICosmoDE *cosmo_de, gdouble z);

static void
nc_hicosmo_de_wpoly_class_init (NcHICosmoDEWPolyClass *klass)
{
  GObjectClass *object_class     = G_OBJECT_CLASS (klass);
  NcHICosmoDEClass *parent_class = NC_HICOSMO_DE_CLASS (klass);
  NcmModelClass *model_class     = NCM_MODEL_CLASS (klass);

  model_class->set_property = &_nc_hicosmo_de_wpoly_set_property;
  model_class->get_property = &_nc_hicosmo_de_wpoly_get_property;
  object_class->constructed = &_nc_hicosmo_de_wpoly_constructed;
  object_class->finalize    = &_nc_hicosmo_de_wpoly_finalize;

  ncm_model_class_set_name_nick (model_class, "WPoly - DE EOS Legendre polynomial w(z)", "WPoly");
  ncm_model_class_add_params (model_class, 0,
                              NC_HICOSMO_DE_WPOLY_VPARAM_LEN - NC_HICOSMO_DE_VPARAM_LEN,
                              PROP_SIZE);

  g_object_class_install_property (object_class,
                                   PROP_Z_F,
                                   g_param_spec_double ("zf",
                                                        NULL,
                                                        "final redshift",
                                                        1.0e-3, 1.0e10, 2.33,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  ncm_model_class_set_vparam (model_class, NC_HICOSMO_DE_WPOLY_A, NC_HICOSMO_DE_WPOLY_DEFAULT_A_LEN, "a", "a",
                              -5.0, 5.0, 1.0e-2,
                              NC_HICOSMO_DEFAULT_PARAMS_ABSTOL, NC_HICOSMO_DE_WPOLY_DEFAULT_A,
                              NCM_PARAM_TYPE_FREE);

  ncm_model_class_check_params_info (model_class);

  nc_hicosmo_de_set_E2Omega_de_impl (parent_class, &_nc_hicosmo_de_wpoly_E2Omega_de);
  nc_hicosmo_de_set_dE2Omega_de_dz_impl (parent_class, &_nc_hicosmo_de_wpoly_dE2Omega_de_dz);
  nc_hicosmo_de_set_d2E2Omega_de_dz2_impl (parent_class, &_nc_hicosmo_de_wpoly_d2E2Omega_de_dz2);
  nc_hicosmo_de_set_w_de_impl (parent_class, &_nc_hicosmo_de_wpoly_w_de);
}

#define VECTOR (NCM_MODEL (cosmo_de))
#define OMEGA_X (ncm_model_orig_param_get (VECTOR, NC_HICOSMO_DE_OMEGA_X))

/* w(z) + 1 on [0, z_f]. */
static gdouble
_nc_hicosmo_de_wpoly_w1 (NcHICosmoDEWPoly *wpoly, const gdouble z)
{
  NcHICosmoDEWPolyPrivate * const self = wpoly->priv;
  NcmModel *model                      = NCM_MODEL (wpoly);
  gdouble *Pl                          = g_newa (gdouble, self->ncoeffs);
  const gdouble x                      = 2.0 * z / self->z_f - 1.0;
  gdouble w1                           = 0.0;
  guint l;

  gsl_sf_legendre_Pl_array (self->lmax, x, Pl);

  for (l = 0; l < self->ncoeffs; l++)
    w1 += ncm_model_orig_vparam_get (model, NC_HICOSMO_DE_WPOLY_A, l) * Pl[l];

  return w1;
}

/* dw/dz on [0, z_f]. */
static gdouble
_nc_hicosmo_de_wpoly_dw1 (NcHICosmoDEWPoly *wpoly, const gdouble z)
{
  NcHICosmoDEWPolyPrivate * const self = wpoly->priv;
  NcmModel *model                      = NCM_MODEL (wpoly);
  gdouble *Pl                          = g_newa (gdouble, self->ncoeffs);
  gdouble *dPl                         = g_newa (gdouble, self->ncoeffs);
  const gdouble x                      = 2.0 * z / self->z_f - 1.0;
  gdouble dw                           = 0.0;
  guint l;

  gsl_sf_legendre_Pl_deriv_array (self->lmax, x, Pl, dPl);

  for (l = 0; l < self->ncoeffs; l++)
    dw += ncm_model_orig_vparam_get (model, NC_HICOSMO_DE_WPOLY_A, l) * dPl[l];

  return dw * 2.0 / self->z_f;
}

/* int_0^z w(z')/(1+z') dz' on [0, z_f]. */
static gdouble
_nc_hicosmo_de_wpoly_int_w (NcHICosmoDEWPoly *wpoly, const gdouble z)
{
  NcHICosmoDEWPolyPrivate * const self = wpoly->priv;
  gdouble res                          = 0.0;
  gsize i;

  for (i = 0; i < self->gl_n; i++)
  {
    const gdouble zi = 0.5 * z * (self->gl_x[i] + 1.0);

    res += self->gl_w[i] * (_nc_hicosmo_de_wpoly_w1 (wpoly, zi) - 1.0) / (1.0 + zi);
  }

  return 0.5 * z * res;
}

static void
_nc_hicosmo_de_wpoly_prepare (NcHICosmoDEWPoly *wpoly)
{
  NcHICosmoDEWPolyPrivate * const self = wpoly->priv;

  if (!ncm_model_lstate_is_update (NCM_MODEL (wpoly), 0))
  {
    self->w_f   = _nc_hicosmo_de_wpoly_w1 (wpoly, self->z_f) - 1.0;
    self->int_f = _nc_hicosmo_de_wpoly_int_w (wpoly, self->z_f);

    ncm_model_lstate_set_update (NCM_MODEL (wpoly), 0);
  }
}

/* ln (E2Omega_de / Omega_x) = 3 int_0^z (1+w)/(1+z') dz', both branches. */
static gdouble
_nc_hicosmo_de_wpoly_ln_E2Omega_de (NcHICosmoDEWPoly *wpoly, const gdouble z, gdouble *w, gdouble *dw)
{
  NcHICosmoDEWPolyPrivate * const self = wpoly->priv;
  const gdouble alpha                  = log1p (z);

  _nc_hicosmo_de_wpoly_prepare (wpoly);

  if (z < self->z_f)
  {
    *w  = _nc_hicosmo_de_wpoly_w1 (wpoly, z) - 1.0;
    *dw = _nc_hicosmo_de_wpoly_dw1 (wpoly, z);

    return 3.0 * (alpha + _nc_hicosmo_de_wpoly_int_w (wpoly, z));
  }
  else
  {
    *w  = self->w_f;
    *dw = 0.0;

    return 3.0 * (alpha + self->int_f + self->w_f * (alpha - self->alpha_f));
  }
}

static gdouble
_nc_hicosmo_de_wpoly_E2Omega_de (NcHICosmoDE *cosmo_de, gdouble z)
{
  NcHICosmoDEWPoly *wpoly = NC_HICOSMO_DE_WPOLY (cosmo_de);
  gdouble w, dw;

  return OMEGA_X * exp (_nc_hicosmo_de_wpoly_ln_E2Omega_de (wpoly, z, &w, &dw));
}

static gdouble
_nc_hicosmo_de_wpoly_dE2Omega_de_dz (NcHICosmoDE *cosmo_de, gdouble z)
{
  NcHICosmoDEWPoly *wpoly = NC_HICOSMO_DE_WPOLY (cosmo_de);
  gdouble w, dw;
  const gdouble E2Omega_de = OMEGA_X * exp (_nc_hicosmo_de_wpoly_ln_E2Omega_de (wpoly, z, &w, &dw));

  return 3.0 * (1.0 + w) * E2Omega_de / (1.0 + z);
}

static gdouble
_nc_hicosmo_de_wpoly_d2E2Omega_de_dz2 (NcHICosmoDE *cosmo_de, gdouble z)
{
  NcHICosmoDEWPoly *wpoly = NC_HICOSMO_DE_WPOLY (cosmo_de);
  gdouble w, dw;
  const gdouble E2Omega_de = OMEGA_X * exp (_nc_hicosmo_de_wpoly_ln_E2Omega_de (wpoly, z, &w, &dw));
  const gdouble onepz      = 1.0 + z;

  return 3.0 * (2.0 + w * (5.0 + 3.0 * w) + dw * onepz) * E2Omega_de / (onepz * onepz);
}

static gdouble
_nc_hicosmo_de_wpoly_w_de (NcHICosmoDE *cosmo_de, gdouble z)
{
  NcHICosmoDEWPoly *wpoly              = NC_HICOSMO_DE_WPOLY (cosmo_de);
  NcHICosmoDEWPolyPrivate * const self = wpoly->priv;

  _nc_hicosmo_de_wpoly_prepare (wpoly);

  if (z < self->z_f)
    return _nc_hicosmo_de_wpoly_w1 (wpoly, z) - 1.0;
  else
    return self->w_f;
}

/**
 * nc_hicosmo_de_wpoly_new:
 * @ncoeffs: number of Legendre coefficients, $L + 1$
 * @z_f: final redshift $z_f$
 *
 * Creates a #NcHICosmoDEWPoly with all coefficients zero, that is $w = -1$.
 *
 * Returns: a new #NcHICosmoDEWPoly
 */
NcHICosmoDEWPoly *
nc_hicosmo_de_wpoly_new (gsize ncoeffs, const gdouble z_f)
{
  NcHICosmoDEWPoly *wpoly = g_object_new (NC_TYPE_HICOSMO_DE_WPOLY,
                                          "zf", z_f,
                                          "a-length", ncoeffs,
                                          NULL);

  return wpoly;
}

/**
 * nc_hicosmo_de_wpoly_get_z_f:
 * @wpoly: a #NcHICosmoDEWPoly
 *
 * Returns: the final redshift $z_f$
 */
gdouble
nc_hicosmo_de_wpoly_get_z_f (NcHICosmoDEWPoly *wpoly)
{
  NcHICosmoDEWPolyPrivate * const self = wpoly->priv;

  return self->z_f;
}

/**
 * nc_hicosmo_de_wpoly_get_ncoeffs:
 * @wpoly: a #NcHICosmoDEWPoly
 *
 * Returns: the number of Legendre coefficients, $L + 1$
 */
guint
nc_hicosmo_de_wpoly_get_ncoeffs (NcHICosmoDEWPoly *wpoly)
{
  NcHICosmoDEWPolyPrivate * const self = wpoly->priv;

  return self->ncoeffs;
}

/**
 * nc_hicosmo_de_wpoly_dw_dz:
 * @wpoly: a #NcHICosmoDEWPoly
 * @z: redshift, $0 \le z$
 *
 * Computes $\mathrm{d}w/\mathrm{d}z$, zero for $z \ge z_f$.
 *
 * Returns: $\mathrm{d}w/\mathrm{d}z$ at @z
 */
gdouble
nc_hicosmo_de_wpoly_dw_dz (NcHICosmoDEWPoly *wpoly, const gdouble z)
{
  NcHICosmoDEWPolyPrivate * const self = wpoly->priv;

  if (z < self->z_f)
    return _nc_hicosmo_de_wpoly_dw1 (wpoly, z);
  else
    return 0.0;
}

/**
 * nc_hicosmo_de_wpoly_turning_points:
 * @wpoly: a #NcHICosmoDEWPoly
 *
 * Finds the turning points of $w(z)$ inside $(0, z_f)$: the real roots of
 * $\mathrm{d}w/\mathrm{d}z$, from the eigenvalues of the Legendre comrade matrix
 * of the derivative series. An eigenvalue counts as real when its imaginary part
 * is below $10^{-8}$; a repeated root is then reported twice or not at all.
 *
 * Returns: (element-type gdouble) (transfer container): the turning points in increasing $z$, possibly empty
 */
GArray *
nc_hicosmo_de_wpoly_turning_points (NcHICosmoDEWPoly *wpoly)
{
  NcHICosmoDEWPolyPrivate * const self = wpoly->priv;
  NcmModel *model                      = NCM_MODEL (wpoly);
  GArray *tp                           = g_array_new (FALSE, FALSE, sizeof (gdouble));
  const guint nc                       = self->ncoeffs;
  gdouble *c                           = g_new (gdouble, nc);
  gdouble *d                           = g_new0 (gdouble, nc);
  guint nd, j;

  for (j = 0; j < nc; j++)
    c[j] = ncm_model_orig_vparam_get (model, NC_HICOSMO_DE_WPOLY_A, j);

  /* Legendre series of dw/dx, x = 2z/z_f - 1, degree nc - 2. */
  for (j = nc - 1; j >= 1; j--)
  {
    d[j - 1] = (2.0 * j - 1.0) * c[j];

    if (j >= 2)
      c[j - 2] += c[j];
  }

  nd = nc - 1;

  while ((nd > 0) && (d[nd - 1] == 0.0))
    nd--;

  if (nd == 2)
  {
    const gdouble x = -d[0] / d[1];

    if ((x > -1.0) && (x < 1.0))
    {
      const gdouble z = 0.5 * self->z_f * (x + 1.0);

      g_array_append_val (tp, z);
    }
  }
  else if (nd > 2)
  {
    /* Comrade matrix of the derivative series, symmetrized as in numpy's
     * legcompanion; the last column carries the coefficients.
     */
    const guint n                   = nd - 1;
    gsl_matrix *mat                 = gsl_matrix_calloc (n, n);
    gsl_vector_complex *eval        = gsl_vector_complex_alloc (n);
    gsl_eigen_nonsymm_workspace *ws = gsl_eigen_nonsymm_alloc (n);
    gdouble *scl                    = g_new (gdouble, n);
    guint i;

    for (i = 0; i < n; i++)
      scl[i] = 1.0 / sqrt (2.0 * i + 1.0);

    for (i = 0; i + 1 < n; i++)
    {
      const gdouble v = (i + 1.0) * scl[i] * scl[i + 1];

      gsl_matrix_set (mat, i, i + 1, v);
      gsl_matrix_set (mat, i + 1, i, v);
    }

    for (i = 0; i < n; i++)
    {
      const gdouble v = gsl_matrix_get (mat, i, n - 1) - (d[i] / d[n]) * (scl[i] / scl[n - 1]) * (n / (2.0 * n - 1.0));

      gsl_matrix_set (mat, i, n - 1, v);
    }

    gsl_eigen_nonsymm (mat, eval, ws);

    for (i = 0; i < n; i++)
    {
      const gsl_complex r = gsl_vector_complex_get (eval, i);
      const gdouble x     = GSL_REAL (r);

      if ((fabs (GSL_IMAG (r)) < 1.0e-8) && (x > -1.0) && (x < 1.0))
      {
        const gdouble z = 0.5 * self->z_f * (x + 1.0);

        g_array_append_val (tp, z);
      }
    }

    if (tp->len > 1)
      gsl_sort ((gdouble *) tp->data, 1, tp->len);

    g_free (scl);
    gsl_eigen_nonsymm_free (ws);
    gsl_vector_complex_free (eval);
    gsl_matrix_free (mat);
  }

  g_free (c);
  g_free (d);

  return tp;
}

