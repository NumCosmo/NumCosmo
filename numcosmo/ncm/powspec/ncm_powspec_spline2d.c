/***************************************************************************
 *            ncm_powspec_spline2d.c
 *
 *  Tue February 16 17:00:52 2016
 *  Copyright  2016  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * ncm_powspec_spline2d.c
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
 * NcmPowspecSpline2d:
 *
 * Power spectrum read from a table: an #NcmSpline2d of $\ln P$ on $(z, \ln k)$,
 * with $k$ in $\mathrm{Mpc}^{-1}$ and $P$ in $\mathrm{Mpc}^3$.
 *
 * The knots of the table set NcmPowspec:zi, NcmPowspec:zf, NcmPowspec:kmin and
 * NcmPowspec:kmax. Outside $[k_\mathrm{min}, k_\mathrm{max}]$, $\ln P$ continues
 * linearly in $\ln k$ with the value and slope of the spline at the nearest end,
 * so $P$ and its first derivatives are continuous there. Evaluating at a $z$
 * outside the table aborts, and so does preparing after NcmPowspec:zi or
 * NcmPowspec:zf was moved outside it.
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "ncm/powspec/ncm_powspec_spline2d.h"

typedef struct _NcmPowspecSpline2dPrivate
{
  NcmSpline2d *spline2d;
  gdouble zmin;
  gdouble zmax;
  gdouble lnkmin;
  gdouble lnkmax;
} NcmPowspecSpline2dPrivate;

struct _NcmPowspecSpline2d
{
  NcmPowspec parent_instance;
};

enum
{
  PROP_0,
  PROP_SPLINE2D,
};

G_DEFINE_TYPE_WITH_PRIVATE (NcmPowspecSpline2d, ncm_powspec_spline2d, NCM_TYPE_POWSPEC)

static void
ncm_powspec_spline2d_init (NcmPowspecSpline2d *ps_s2d)
{
  NcmPowspecSpline2dPrivate * const self = ncm_powspec_spline2d_get_instance_private (ps_s2d);

  self->spline2d = NULL;
  self->zmin     = 0.0;
  self->zmax     = 0.0;
  self->lnkmin   = 0.0;
  self->lnkmax   = 0.0;
}

static void
_ncm_powspec_spline2d_dispose (GObject *object)
{
  NcmPowspecSpline2d *ps_s2d             = NCM_POWSPEC_SPLINE2D (object);
  NcmPowspecSpline2dPrivate * const self = ncm_powspec_spline2d_get_instance_private (ps_s2d);

  ncm_spline2d_clear (&self->spline2d);

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_powspec_spline2d_parent_class)->dispose (object);
}

static void
_ncm_powspec_spline2d_finalize (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (ncm_powspec_spline2d_parent_class)->finalize (object);
}

static void
_ncm_powspec_spline2d_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmPowspecSpline2d *ps_s2d = NCM_POWSPEC_SPLINE2D (object);

  g_return_if_fail (NCM_IS_POWSPEC_SPLINE2D (object));

  switch (prop_id)
  {
    case PROP_SPLINE2D:
      ncm_powspec_spline2d_set_spline2d (ps_s2d, g_value_get_object (value));
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_ncm_powspec_spline2d_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmPowspecSpline2d *ps_s2d = NCM_POWSPEC_SPLINE2D (object);

  g_return_if_fail (NCM_IS_POWSPEC_SPLINE2D (object));

  switch (prop_id)
  {
    case PROP_SPLINE2D:
      g_value_set_object (value, ncm_powspec_spline2d_peek_spline2d (ps_s2d));
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void _ncm_powspec_spline2d_prepare (NcmPowspec *powspec, NcmModel *model);
static gdouble _ncm_powspec_spline2d_eval (NcmPowspec *powspec, NcmModel *model, const gdouble z, const gdouble k);
static void _ncm_powspec_spline2d_eval_vec (NcmPowspec *powspec, NcmModel *model, const gdouble z, NcmVector *k, NcmVector *Pk);
static gdouble _ncm_powspec_spline2d_deriv_z (NcmPowspec *powspec, NcmModel *model, const gdouble z, const gdouble k);
static gdouble _ncm_powspec_spline2d_deriv_k (NcmPowspec *powspec, NcmModel *model, const gdouble z, const gdouble k);
static void _ncm_powspec_spline2d_get_nknots (NcmPowspec *powspec, guint *Nz, guint *Nk);

static void
ncm_powspec_spline2d_class_init (NcmPowspecSpline2dClass *klass)
{
  GObjectClass *object_class     = G_OBJECT_CLASS (klass);
  NcmPowspecClass *powspec_class = NCM_POWSPEC_CLASS (klass);

  object_class->set_property = &_ncm_powspec_spline2d_set_property;
  object_class->get_property = &_ncm_powspec_spline2d_get_property;
  object_class->dispose      = &_ncm_powspec_spline2d_dispose;
  object_class->finalize     = &_ncm_powspec_spline2d_finalize;

  /**
   * NcmPowspecSpline2d:spline2d:
   *
   * The #NcmSpline2d of $\ln P$ on $(z, \ln k)$; required at construction.
   */
  g_object_class_install_property (object_class,
                                   PROP_SPLINE2D,
                                   g_param_spec_object ("spline2d",
                                                        NULL,
                                                        "Spline2d of ln P on (z, ln k)",
                                                        NCM_TYPE_SPLINE2D,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  powspec_class->prepare    = &_ncm_powspec_spline2d_prepare;
  powspec_class->eval       = &_ncm_powspec_spline2d_eval;
  powspec_class->eval_vec   = &_ncm_powspec_spline2d_eval_vec;
  powspec_class->deriv_z    = &_ncm_powspec_spline2d_deriv_z;
  powspec_class->deriv_k    = &_ncm_powspec_spline2d_deriv_k;
  powspec_class->get_nknots = &_ncm_powspec_spline2d_get_nknots;
}

static void
_ncm_powspec_spline2d_prepare (NcmPowspec *powspec, NcmModel *model)
{
  NcmPowspecSpline2d *ps_s2d             = NCM_POWSPEC_SPLINE2D (powspec);
  NcmPowspecSpline2dPrivate * const self = ncm_powspec_spline2d_get_instance_private (ps_s2d);
  const gdouble zi                       = ncm_powspec_get_zi (powspec);
  const gdouble zf                       = ncm_powspec_get_zf (powspec);

  if ((zi < self->zmin) || (zf > self->zmax))
    g_error ("_ncm_powspec_spline2d_prepare: the requested z range [%g, %g] is not inside the table [%g, %g].",
             zi, zf, self->zmin, self->zmax);

  if (!ncm_spline2d_is_init (self->spline2d))
    ncm_spline2d_prepare (self->spline2d);
}

static void _ncm_powspec_spline2d_check_z (NcmPowspecSpline2dPrivate * const self, const gdouble z);
static gdouble _ncm_powspec_spline2d_edge (NcmPowspecSpline2dPrivate * const self, const gdouble lnk);

static gdouble
_ncm_powspec_spline2d_eval (NcmPowspec *powspec, NcmModel *model, const gdouble z, const gdouble k)
{
  NcmPowspecSpline2d *ps_s2d             = NCM_POWSPEC_SPLINE2D (powspec);
  NcmPowspecSpline2dPrivate * const self = ncm_powspec_spline2d_get_instance_private (ps_s2d);
  const gdouble lnk                      = log (k);
  const gdouble lnk_e                    = _ncm_powspec_spline2d_edge (self, lnk);

  _ncm_powspec_spline2d_check_z (self, z);

  if (lnk_e == lnk)
    return exp (ncm_spline2d_eval (self->spline2d, z, lnk));
  else
    return exp (ncm_spline2d_eval (self->spline2d, z, lnk_e) + ncm_spline2d_deriv_dzdy (self->spline2d, z, lnk_e) * (lnk - lnk_e));
}

static void
_ncm_powspec_spline2d_eval_vec (NcmPowspec *powspec, NcmModel *model, const gdouble z, NcmVector *k, NcmVector *Pk)
{
  const guint n = ncm_vector_len (k);
  guint i;

  g_assert_cmpuint (n, ==, ncm_vector_len (Pk));

  for (i = 0; i < n; i++)
    ncm_vector_set (Pk, i, _ncm_powspec_spline2d_eval (powspec, model, z, ncm_vector_get (k, i)));
}

/*
 * Derivatives of the spline of ln P: dP/dz = P d(ln P)/dz and
 * dP/dk = P d(ln P)/d(ln k) / k. On the continuation ln P = s(z, e) + s_y(z, e) (ln k - e),
 * with e the nearest end, so d(ln P)/dz = s_x(z, e) + s_xy(z, e) (ln k - e) and
 * d(ln P)/d(ln k) = s_y(z, e).
 */
static gdouble
_ncm_powspec_spline2d_deriv_z (NcmPowspec *powspec, NcmModel *model, const gdouble z, const gdouble k)
{
  NcmPowspecSpline2d *ps_s2d             = NCM_POWSPEC_SPLINE2D (powspec);
  NcmPowspecSpline2dPrivate * const self = ncm_powspec_spline2d_get_instance_private (ps_s2d);
  const gdouble lnk                      = log (k);
  const gdouble lnk_e                    = _ncm_powspec_spline2d_edge (self, lnk);
  const gdouble P                        = _ncm_powspec_spline2d_eval (powspec, model, z, k);

  return P * (ncm_spline2d_deriv_dzdx (self->spline2d, z, lnk_e) + ncm_spline2d_deriv_d2zdxy (self->spline2d, z, lnk_e) * (lnk - lnk_e));
}

static gdouble
_ncm_powspec_spline2d_deriv_k (NcmPowspec *powspec, NcmModel *model, const gdouble z, const gdouble k)
{
  NcmPowspecSpline2d *ps_s2d             = NCM_POWSPEC_SPLINE2D (powspec);
  NcmPowspecSpline2dPrivate * const self = ncm_powspec_spline2d_get_instance_private (ps_s2d);
  const gdouble lnk_e                    = _ncm_powspec_spline2d_edge (self, log (k));
  const gdouble P                        = _ncm_powspec_spline2d_eval (powspec, model, z, k);

  return P * ncm_spline2d_deriv_dzdy (self->spline2d, z, lnk_e) / k;
}

static void
_ncm_powspec_spline2d_get_nknots (NcmPowspec *powspec, guint *Nz, guint *Nk)
{
  NcmPowspecSpline2d *ps_s2d             = NCM_POWSPEC_SPLINE2D (powspec);
  NcmPowspecSpline2dPrivate * const self = ncm_powspec_spline2d_get_instance_private (ps_s2d);

  Nz[0] = ncm_vector_len (ncm_spline2d_peek_xv (self->spline2d));
  Nk[0] = ncm_vector_len (ncm_spline2d_peek_yv (self->spline2d));
}

static void
_ncm_powspec_spline2d_check_z (NcmPowspecSpline2dPrivate * const self, const gdouble z)
{
  if ((z < self->zmin) || (z > self->zmax))
    g_error ("_ncm_powspec_spline2d_check_z: z = %g is outside the table [%g, %g].", z, self->zmin, self->zmax);
}

/* @lnk inside the table, or the nearest end of it. */
static gdouble
_ncm_powspec_spline2d_edge (NcmPowspecSpline2dPrivate * const self, const gdouble lnk)
{
  return GSL_MIN (GSL_MAX (lnk, self->lnkmin), self->lnkmax);
}

/**
 * ncm_powspec_spline2d_new:
 * @spline2d: a #NcmSpline2d of $\ln P$ on $(z, \ln k)$
 *
 * Returns: (transfer full): a new #NcmPowspecSpline2d
 */
NcmPowspecSpline2d *
ncm_powspec_spline2d_new (NcmSpline2d *spline2d)
{
  NcmPowspecSpline2d *ps_s2d = g_object_new (NCM_TYPE_POWSPEC_SPLINE2D,
                                             "spline2d", spline2d,
                                             NULL);

  return ps_s2d;
}

/**
 * ncm_powspec_spline2d_ref:
 * @ps_s2d: a #NcmPowspecSpline2d
 *
 * Increases the reference count of @ps_s2d by one atomically.
 *
 * Returns: (transfer full): @ps_s2d
 */
NcmPowspecSpline2d *
ncm_powspec_spline2d_ref (NcmPowspecSpline2d *ps_s2d)
{
  return g_object_ref (ps_s2d);
}

/**
 * ncm_powspec_spline2d_free:
 * @ps_s2d: a #NcmPowspecSpline2d
 *
 * Atomically decrements the reference count of @ps_s2d by one.
 * If the reference count drops to 0,
 * all memory allocated by @ps_s2d is released.
 */
void
ncm_powspec_spline2d_free (NcmPowspecSpline2d *ps_s2d)
{
  g_object_unref (ps_s2d);
}

/**
 * ncm_powspec_spline2d_clear:
 * @ps_s2d: a #NcmPowspecSpline2d
 *
 * If *@ps_s2d is not %NULL, decrements its reference count and sets
 * *@ps_s2d to %NULL.
 */
void
ncm_powspec_spline2d_clear (NcmPowspecSpline2d **ps_s2d)
{
  g_clear_object (ps_s2d);
}

/**
 * ncm_powspec_spline2d_set_spline2d:
 * @ps_s2d: a #NcmPowspecSpline2d
 * @spline2d: a #NcmSpline2d of $\ln P$ on $(z, \ln k)$
 *
 * Sets the table to @spline2d and the ranges of @ps_s2d to its knots, and forces
 * the next ncm_powspec_prepare_if_needed() to prepare.
 */
void
ncm_powspec_spline2d_set_spline2d (NcmPowspecSpline2d *ps_s2d, NcmSpline2d *spline2d)
{
  NcmPowspecSpline2dPrivate * const self = ncm_powspec_spline2d_get_instance_private (ps_s2d);
  NcmPowspec *powspec                    = NCM_POWSPEC (ps_s2d);
  NcmSpline2d *old                       = self->spline2d;

  g_assert_nonnull (spline2d);

  self->spline2d = ncm_spline2d_ref (spline2d);
  ncm_spline2d_clear (&old);

  {
    NcmVector *z_vec   = ncm_spline2d_peek_xv (self->spline2d);
    NcmVector *lnk_vec = ncm_spline2d_peek_yv (self->spline2d);

    self->zmin   = ncm_vector_get (z_vec, 0);
    self->zmax   = ncm_vector_get (z_vec, ncm_vector_len (z_vec) - 1);
    self->lnkmin = ncm_vector_get (lnk_vec, 0);
    self->lnkmax = ncm_vector_get (lnk_vec, ncm_vector_len (lnk_vec) - 1);
  }

  ncm_powspec_set_zi (powspec, self->zmin);
  ncm_powspec_set_zf (powspec, self->zmax);
  ncm_powspec_set_kmin (powspec, exp (self->lnkmin));
  ncm_powspec_set_kmax (powspec, exp (self->lnkmax));
  ncm_model_ctrl_force_update (ncm_powspec_peek_model_ctrl (powspec));
}

/**
 * ncm_powspec_spline2d_peek_spline2d:
 * @ps_s2d: a #NcmPowspecSpline2d
 *
 * Returns: (transfer none): the #NcmSpline2d of $\ln P$ on $(z, \ln k)$
 */
NcmSpline2d *
ncm_powspec_spline2d_peek_spline2d (NcmPowspecSpline2d *ps_s2d)
{
  NcmPowspecSpline2dPrivate * const self = ncm_powspec_spline2d_get_instance_private (ps_s2d);

  return self->spline2d;
}

