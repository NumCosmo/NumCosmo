/***************************************************************************
 *            ncm_sf_spherical_harmonics.c
 *
 *  Thu December 14 11:18:00 2017
 *  Copyright  2017  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * ncm_sf_spherical_harmonics.c
 * Copyright (C) 2017 Sandro Dias Pinto Vitenti <vitenti@uel.br>
 *
 * NumCosmo is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * NumCosmo is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * NcmSFSphericalHarmonics:
 *
 * Recursions for the spherical harmonics.
 *
 * Computes
 * $$
 * \bar{Y}_l^m(x) = \sqrt{\frac{2l+1}{4\pi}\frac{(l-m)!}{(l+m)!}}\,P_l^m(x), \qquad
 * x = \cos\theta,
 * $$
 * with the Condon-Shortley phase in $P_l^m$, so that $Y_l^m(\theta, \phi) =
 * \bar{Y}_l^m(\cos\theta)\,e^{im\phi}$. These are the basis of #NcmSphereMap.
 *
 * The object holds the recurrence coefficients up to #NcmSFSphericalHarmonics:lmax. A
 * #NcmSFSphericalHarmonicsY follows the recursion at one angle, and a
 * #NcmSFSphericalHarmonicsYArray at up to #NCM_SF_SPHERICAL_HARMONICS_MAX_LEN angles in
 * lockstep. Both start at $l = m = 0$ and step in $l$ with
 * $$
 * \bar{Y}_{l+2}^m = K^{l+1}_{lm}\,x\,\bar{Y}_{l+1}^m - K^l_{lm}\,\bar{Y}_l^m
 * $$
 * (see #NcmSFSphericalHarmonicsK). The step to $m + 1$ skips the leading orders whose
 * $|\bar{Y}_l^m|$ is below the absolute tolerance, so the first order must be read with
 * ncm_sf_spherical_harmonics_Y_get_l() after it. The seeds of the recursion in $l$ are
 * stored scaled by $10^{280}$, which keeps $\bar{Y}_m^m \propto \sin^m\theta$
 * representable at large $m$.
 *
 * The errors are small relative to the largest $|\bar{Y}_l^m|$ at the angle, as in any
 * recursion in $\cos\theta$. Near the poles, the step to $m + 1$ after skipped orders
 * subtracts nearly equal seeds, so the rows it reaches, whose values are tiny, are
 * accurate only at that global scale and not relative to their own peak.
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "ncm/specfunc/ncm_sf_spherical_harmonics.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <math.h>
#endif /* NUMCOSMO_GIR_SCAN */

enum
{
  PROP_0,
  PROP_LMAX,
};

G_DEFINE_BOXED_TYPE (NcmSFSphericalHarmonicsY, ncm_sf_spherical_harmonics_Y, ncm_sf_spherical_harmonics_Y_dup, ncm_sf_spherical_harmonics_Y_free)
G_DEFINE_BOXED_TYPE (NcmSFSphericalHarmonicsYArray, ncm_sf_spherical_harmonics_Y_array, ncm_sf_spherical_harmonics_Y_array_dup, ncm_sf_spherical_harmonics_Y_array_free)
G_DEFINE_TYPE (NcmSFSphericalHarmonics, ncm_sf_spherical_harmonics, G_TYPE_OBJECT)

static void
ncm_sf_spherical_harmonics_init (NcmSFSphericalHarmonics *spha)
{
  spha->lmax     = -1; /* so that set_lmax (0) builds the tables */
  spha->sqrt_n   = g_array_new (FALSE, FALSE, sizeof (gdouble));
  spha->sqrtm1_n = g_array_new (FALSE, FALSE, sizeof (gdouble));
  spha->K_array  = g_ptr_array_new ();

  g_ptr_array_set_free_func (spha->K_array, (GDestroyNotify) g_array_unref);
}

static void
_ncm_sf_spherical_harmonics_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmSFSphericalHarmonics *spha = NCM_SF_SPHERICAL_HARMONICS (object);

  g_return_if_fail (NCM_IS_SF_SPHERICAL_HARMONICS (object));

  switch (prop_id)
  {
    case PROP_LMAX:
      ncm_sf_spherical_harmonics_set_lmax (spha, g_value_get_int (value));
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_ncm_sf_spherical_harmonics_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmSFSphericalHarmonics *spha = NCM_SF_SPHERICAL_HARMONICS (object);

  g_return_if_fail (NCM_IS_SF_SPHERICAL_HARMONICS (object));

  switch (prop_id)
  {
    case PROP_LMAX:
      g_value_set_int (value, ncm_sf_spherical_harmonics_get_lmax (spha));
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_ncm_sf_spherical_harmonics_dispose (GObject *object)
{
  NcmSFSphericalHarmonics *spha = NCM_SF_SPHERICAL_HARMONICS (object);

  g_clear_pointer (&spha->sqrt_n, g_array_unref);
  g_clear_pointer (&spha->sqrtm1_n, g_array_unref);
  g_clear_pointer (&spha->K_array, g_ptr_array_unref);

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_sf_spherical_harmonics_parent_class)->dispose (object);
}

static void
_ncm_sf_spherical_harmonics_finalize (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (ncm_sf_spherical_harmonics_parent_class)->finalize (object);
}

static void
ncm_sf_spherical_harmonics_class_init (NcmSFSphericalHarmonicsClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);

  object_class->set_property = &_ncm_sf_spherical_harmonics_set_property;
  object_class->get_property = &_ncm_sf_spherical_harmonics_get_property;
  object_class->dispose      = &_ncm_sf_spherical_harmonics_dispose;
  object_class->finalize     = &_ncm_sf_spherical_harmonics_finalize;

  /**
   * NcmSFSphericalHarmonics:lmax:
   *
   * Largest multipole. Setting it precomputes the recurrence coefficients for every
   * $(l, m)$ with $m \le l < \ell_\mathrm{max}$, so memory and setup cost grow as
   * $\ell_\mathrm{max}^2$.
   */
  g_object_class_install_property (object_class,
                                   PROP_LMAX,
                                   g_param_spec_int ("lmax",
                                                     NULL,
                                                     "max l",
                                                     0, G_MAXINT, 1024,
                                                     G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
}

/**
 * ncm_sf_spherical_harmonics_Y_new:
 * @spha: a #NcmSFSphericalHarmonics
 * @abstol: absolute tolerance
 *
 * Creates a new #NcmSFSphericalHarmonicsY. ncm_sf_spherical_harmonics_Y_next_m() skips
 * the orders whose $|\bar{Y}_l^m|$ is below @abstol. Start it with
 * ncm_sf_spherical_harmonics_start_rec().
 *
 * Returns: (transfer full): a new #NcmSFSphericalHarmonicsY
 */
NcmSFSphericalHarmonicsY *
ncm_sf_spherical_harmonics_Y_new (NcmSFSphericalHarmonics *spha, const gdouble abstol)
{
  NcmSFSphericalHarmonicsY *sphaY = g_slice_new0 (NcmSFSphericalHarmonicsY);

  sphaY->x        = 0.0;
  sphaY->sqrt1mx2 = 0.0;

  sphaY->l      = 0;
  sphaY->l0     = 0;
  sphaY->m      = 0;
  sphaY->Klm    = NULL;
  sphaY->Pl0m   = 0.0;
  sphaY->Pl0p1m = 0.0;
  sphaY->Plm    = 0.0;
  sphaY->Plp1m  = 0.0;

  sphaY->spha   = ncm_sf_spherical_harmonics_ref (spha);
  sphaY->abstol = abstol;

  return sphaY;
}

/**
 * ncm_sf_spherical_harmonics_Y_dup:
 * @sphaY: a #NcmSFSphericalHarmonicsY
 *
 * Duplicates @sphaY, including the state of its recursion.
 *
 * Returns: (transfer full): a copy of @sphaY
 */
NcmSFSphericalHarmonicsY *
ncm_sf_spherical_harmonics_Y_dup (NcmSFSphericalHarmonicsY *sphaY)
{
  NcmSFSphericalHarmonicsY *sphaY_dup = ncm_sf_spherical_harmonics_Y_new (sphaY->spha, sphaY->abstol);

  /* The copy keeps the spha pointer whose reference the call above took. */
  sphaY_dup[0] = sphaY[0];

  return sphaY_dup;
}

/**
 * ncm_sf_spherical_harmonics_Y_free:
 * @sphaY: a #NcmSFSphericalHarmonicsY
 *
 * Frees @sphaY.
 */
void
ncm_sf_spherical_harmonics_Y_free (NcmSFSphericalHarmonicsY *sphaY)
{
  ncm_sf_spherical_harmonics_clear (&sphaY->spha);
  g_slice_free (NcmSFSphericalHarmonicsY, sphaY);
}

/**
 * ncm_sf_spherical_harmonics_Y_array_new:
 * @spha: a #NcmSFSphericalHarmonics
 * @len: number of angles, at most #NCM_SF_SPHERICAL_HARMONICS_MAX_LEN
 * @abstol: absolute tolerance
 *
 * Creates a new #NcmSFSphericalHarmonicsYArray for @len angles.
 * ncm_sf_spherical_harmonics_Y_array_next_m() skips the orders whose $|\bar{Y}_l^m|$ is
 * below @abstol at every angle. Start it with
 * ncm_sf_spherical_harmonics_start_rec_array().
 *
 * Returns: (transfer full): a new #NcmSFSphericalHarmonicsYArray
 */
NcmSFSphericalHarmonicsYArray *
ncm_sf_spherical_harmonics_Y_array_new (NcmSFSphericalHarmonics *spha, const gint len, const gdouble abstol)
{
  NcmSFSphericalHarmonicsYArray *sphaYa = g_slice_new0 (NcmSFSphericalHarmonicsYArray);

  sphaYa->l   = 0;
  sphaYa->l0  = 0;
  sphaYa->m   = 0;
  sphaYa->Klm = NULL;
  sphaYa->len = len;

  g_assert_cmpuint (len, <=, NCM_SF_SPHERICAL_HARMONICS_MAX_LEN);

  sphaYa->spha   = ncm_sf_spherical_harmonics_ref (spha);
  sphaYa->abstol = abstol;

  return sphaYa;
}

/**
 * ncm_sf_spherical_harmonics_Y_array_dup:
 * @sphaYa: a #NcmSFSphericalHarmonicsYArray
 *
 * Duplicates @sphaYa, including the state of its recursion.
 *
 * Returns: (transfer full): a copy of @sphaYa
 */
NcmSFSphericalHarmonicsYArray *
ncm_sf_spherical_harmonics_Y_array_dup (NcmSFSphericalHarmonicsYArray *sphaYa)
{
  NcmSFSphericalHarmonicsYArray *sphaYa_dup = ncm_sf_spherical_harmonics_Y_array_new (sphaYa->spha, sphaYa->len, sphaYa->abstol);

  sphaYa_dup[0] = sphaYa[0];

  return sphaYa_dup;
}

/**
 * ncm_sf_spherical_harmonics_Y_array_free:
 * @sphaYa: a #NcmSFSphericalHarmonicsYArray
 *
 * Frees @sphaYa.
 */
void
ncm_sf_spherical_harmonics_Y_array_free (NcmSFSphericalHarmonicsYArray *sphaYa)
{
  ncm_sf_spherical_harmonics_clear (&sphaYa->spha);

  g_slice_free (NcmSFSphericalHarmonicsYArray, sphaYa);
}

/**
 * ncm_sf_spherical_harmonics_Y_next_l:
 * @sphaY: a #NcmSFSphericalHarmonicsY
 *
 * Advances the recursion to $l + 1$.
 */
/**
 * ncm_sf_spherical_harmonics_Y_next_l2:
 * @sphaY: a #NcmSFSphericalHarmonicsY
 * @Yblm: (array fixed-size=2) (element-type gdouble) (out caller-allocates): $\bar{Y}_l^m$ and $\bar{Y}_{l+1}^m$
 *
 * Advances the recursion to $l + 2$, returning the two values it passes.
 */
/**
 * ncm_sf_spherical_harmonics_Y_next_l4:
 * @sphaY: a #NcmSFSphericalHarmonicsY
 * @Yblm: (array fixed-size=4) (element-type gdouble) (out caller-allocates): $\bar{Y}_l^m$ to $\bar{Y}_{l+3}^m$
 *
 * Advances the recursion to $l + 4$, returning the four values it passes.
 */
/**
 * ncm_sf_spherical_harmonics_Y_next_l2pn: (skip)
 * @sphaY: a #NcmSFSphericalHarmonicsY
 * @Yblm: output, $\bar{Y}_l^m$ to $\bar{Y}_{l+n+1}^m$
 * @n: number of steps beyond two
 *
 * Advances the recursion to $l + n + 2$, returning the $n + 2$ values it passes.
 */
/**
 * ncm_sf_spherical_harmonics_Y_next_m:
 * @sphaY: a #NcmSFSphericalHarmonicsY
 *
 * Moves the recursion to $m + 1$, restarting it at the lowest order not skipped. The
 * orders $l \geq m + 1$ whose $|\bar{Y}_l^{m+1}|$ is below the absolute tolerance are
 * skipped, up to $\ell_\mathrm{max} + 1$; read the order reached with
 * ncm_sf_spherical_harmonics_Y_get_l(). Once it exceeds $\ell_\mathrm{max}$ the
 * recursion is over and must not be advanced further.
 */

/**
 * ncm_sf_spherical_harmonics_Y_get_lm:
 * @sphaY: a #NcmSFSphericalHarmonicsY
 *
 * Returns: the current value of $\bar{Y}_l^m(x)$
 */
/**
 * ncm_sf_spherical_harmonics_Y_get_lp1m:
 * @sphaY: a #NcmSFSphericalHarmonicsY
 *
 * Returns: the current value of $\bar{Y}_{l+1}^m(x)$
 */
/**
 * ncm_sf_spherical_harmonics_Y_get_x:
 * @sphaY: a #NcmSFSphericalHarmonicsY
 *
 * Returns: the current value of $x$
 */
/**
 * ncm_sf_spherical_harmonics_Y_get_l:
 * @sphaY: a #NcmSFSphericalHarmonicsY
 *
 * Returns: the current value of $l$
 */
/**
 * ncm_sf_spherical_harmonics_Y_get_m:
 * @sphaY: a #NcmSFSphericalHarmonicsY
 *
 * Returns: the current value of $m$
 */

/**
 * ncm_sf_spherical_harmonics_Y_reset:
 * @sphaY: a #NcmSFSphericalHarmonicsY
 *
 * Restarts the recursion at $l = m = 0$, keeping the angle set by
 * ncm_sf_spherical_harmonics_start_rec().
 */

/**
 * ncm_sf_spherical_harmonics_Y_array_next_l:
 * @sphaYa: a #NcmSFSphericalHarmonicsYArray
 * @len: number of angles
 *
 * Advances the recursion at every angle to $l + 1$.
 */
/**
 * ncm_sf_spherical_harmonics_Y_array_next_l2: (skip)
 * @sphaYa: a #NcmSFSphericalHarmonicsYArray
 * @len: number of angles
 * @Yblm: output, $2 \times$ @len values, $\bar{Y}_l^m(x_i)$ and $\bar{Y}_{l+1}^m(x_i)$
 *
 * Advances the recursion at every angle to $l + 2$, returning the values it passes.
 * Index @Yblm with NCM_SF_SPHERICAL_HARMONICS_ARRAY_INDEX().
 */
/**
 * ncm_sf_spherical_harmonics_Y_array_next_l4: (skip)
 * @sphaYa: a #NcmSFSphericalHarmonicsYArray
 * @len: number of angles
 * @Yblm: output, $4 \times$ @len values, $\bar{Y}_l^m(x_i)$ to $\bar{Y}_{l+3}^m(x_i)$
 *
 * Advances the recursion at every angle to $l + 4$, returning the values it passes.
 * Index @Yblm with NCM_SF_SPHERICAL_HARMONICS_ARRAY_INDEX().
 */
/**
 * ncm_sf_spherical_harmonics_Y_array_next_l2pn: (skip)
 * @sphaYa: a #NcmSFSphericalHarmonicsYArray
 * @len: number of angles
 * @Yblm: output, $(n + 2) \times$ @len values, $\bar{Y}_l^m(x_i)$ to $\bar{Y}_{l+n+1}^m(x_i)$
 * @n: number of steps beyond two
 *
 * Advances the recursion at every angle to $l + n + 2$, returning the values it passes.
 * Index @Yblm with NCM_SF_SPHERICAL_HARMONICS_ARRAY_INDEX().
 */
/**
 * ncm_sf_spherical_harmonics_Y_array_next_m:
 * @sphaYa: a #NcmSFSphericalHarmonicsYArray
 * @len: number of angles
 *
 * Same as ncm_sf_spherical_harmonics_Y_next_m() at every angle. An order is skipped
 * only when $|\bar{Y}_l^{m+1}|$ is below the absolute tolerance at some angle, so the
 * angles advance together until the smallest value reaches it.
 */
/**
 * ncm_sf_spherical_harmonics_Y_array_reset:
 * @sphaYa: a #NcmSFSphericalHarmonicsYArray
 * @len: number of angles
 *
 * Restarts the recursion at $l = m = 0$, keeping the angles set by
 * ncm_sf_spherical_harmonics_start_rec_array().
 */

/**
 * ncm_sf_spherical_harmonics_Y_array_get_lm:
 * @sphaYa: a #NcmSFSphericalHarmonicsYArray
 * @len: number of angles
 * @i: angle index
 *
 * Returns: the current value of $\bar{Y}_l^m(x_i)$
 */
/**
 * ncm_sf_spherical_harmonics_Y_array_get_lp1m:
 * @sphaYa: a #NcmSFSphericalHarmonicsYArray
 * @len: number of angles
 * @i: angle index
 *
 * Returns: the current value of $\bar{Y}_{l+1}^m(x_i)$
 */
/**
 * ncm_sf_spherical_harmonics_Y_array_get_x:
 * @sphaYa: a #NcmSFSphericalHarmonicsYArray
 * @i: angle index
 *
 * Returns: the current value of $x_i$
 */
/**
 * ncm_sf_spherical_harmonics_Y_array_get_l:
 * @sphaYa: a #NcmSFSphericalHarmonicsYArray
 *
 * Returns: the current value of $l$
 */
/**
 * ncm_sf_spherical_harmonics_Y_array_get_m:
 * @sphaYa: a #NcmSFSphericalHarmonicsYArray
 *
 * Returns: the current value of $m$
 */

/**
 * ncm_sf_spherical_harmonics_new:
 * @lmax: largest multipole
 *
 * Creates a new #NcmSFSphericalHarmonics.
 *
 * Returns: (transfer full): a new #NcmSFSphericalHarmonics
 */
NcmSFSphericalHarmonics *
ncm_sf_spherical_harmonics_new (const gint lmax)
{
  NcmSFSphericalHarmonics *spha = g_object_new (NCM_TYPE_SF_SPHERICAL_HARMONICS,
                                                "lmax", lmax,
                                                NULL);

  return spha;
}

/**
 * ncm_sf_spherical_harmonics_ref:
 * @spha: a #NcmSFSphericalHarmonics
 *
 * Increases the reference count of @spha by one.
 *
 * Returns: (transfer full): @spha
 */
NcmSFSphericalHarmonics *
ncm_sf_spherical_harmonics_ref (NcmSFSphericalHarmonics *spha)
{
  return g_object_ref (spha);
}

/**
 * ncm_sf_spherical_harmonics_free:
 * @spha: a #NcmSFSphericalHarmonics
 *
 * Decreases the reference count of @spha by one.
 */
void
ncm_sf_spherical_harmonics_free (NcmSFSphericalHarmonics *spha)
{
  g_object_unref (spha);
}

/**
 * ncm_sf_spherical_harmonics_clear:
 * @spha: a #NcmSFSphericalHarmonics
 *
 * If *@spha is not NULL, decreases its reference count by one and sets *@spha to NULL.
 */
void
ncm_sf_spherical_harmonics_clear (NcmSFSphericalHarmonics **spha)
{
  g_clear_object (spha);
}

#define SN(n)   g_array_index (spha->sqrt_n, gdouble, (n))
#define SNM1(n) g_array_index (spha->sqrtm1_n, gdouble, (n))

/**
 * ncm_sf_spherical_harmonics_set_lmax:
 * @spha: a #NcmSFSphericalHarmonics
 * @lmax: largest multipole
 *
 * Sets #NcmSFSphericalHarmonics:lmax.
 */
void
ncm_sf_spherical_harmonics_set_lmax (NcmSFSphericalHarmonics *spha, const gint lmax)
{
  if (lmax != spha->lmax)
  {
    const gint nmax_old = spha->sqrt_n->len;
    const gint nmax     = 2 * lmax + 3;
    gint n;

    g_array_set_size (spha->sqrt_n,   nmax + 1);
    g_array_set_size (spha->sqrtm1_n, nmax + 1);

    for (n = nmax_old; n <= nmax; n++)
    {
      const gdouble sqrt_n = sqrt (n);

      g_array_index (spha->sqrt_n, gdouble, n)   = sqrt_n;
      g_array_index (spha->sqrtm1_n, gdouble, n) = 1.0 / sqrt_n;
    }

    for (n = 0; n <= lmax; n++)
    {
      GArray *Km_array = NULL;
      gint l;

      if (n < (gint) spha->K_array->len)
      {
        Km_array = g_ptr_array_index (spha->K_array, n);
      }
      else if (n == (gint) spha->K_array->len)
      {
        Km_array = g_array_new (TRUE, TRUE, sizeof (NcmSFSphericalHarmonicsK));
        g_ptr_array_add (spha->K_array, Km_array);
      }
      else
      {
        g_assert_not_reached ();
      }

      for (l = Km_array->len + n; l < lmax; l++)
      {
        NcmSFSphericalHarmonicsK Km;
        const gint m       = n;
        const gint twol    = 2 * l;
        const gint lmm     = l - m;
        const gint lpm     = l + m;
        const gdouble pref = SN (twol + 5) * SNM1 (lpm + 2) * SNM1 (lmm + 2);

        Km.lp1 = pref * SN (twol + 3);
        Km.l   = pref * SN (lmm + 1) * SN (lpm + 1) * SNM1 (twol + 1);

        g_array_append_val (Km_array, Km);
      }

      g_array_set_size (Km_array, lmax - n);
    }

    spha->lmax = lmax;
  }
}

/**
 * ncm_sf_spherical_harmonics_get_lmax:
 * @spha: a #NcmSFSphericalHarmonics
 *
 * Returns: the #NcmSFSphericalHarmonics:lmax
 */
guint
ncm_sf_spherical_harmonics_get_lmax (NcmSFSphericalHarmonics *spha)
{
  return spha->lmax;
}

/**
 * ncm_sf_spherical_harmonics_start_rec:
 * @spha: a #NcmSFSphericalHarmonics
 * @sphaY: a #NcmSFSphericalHarmonicsY
 * @theta: polar angle $\theta \in [0, \pi]$
 *
 * Starts the recursion of @sphaY at @theta and $l = m = 0$.
 */
/**
 * ncm_sf_spherical_harmonics_start_rec_array:
 * @spha: a #NcmSFSphericalHarmonics
 * @sphaYa: a #NcmSFSphericalHarmonicsYArray
 * @len: number of angles
 * @theta: (array length=len) (element-type gdouble): polar angles $\theta_i \in [0, \pi]$
 *
 * Starts the recursion of @sphaYa at the angles @theta and $l = m = 0$.
 */
/**
 * ncm_sf_spherical_harmonics_get_Klm: (skip)
 * @spha: a #NcmSFSphericalHarmonics
 * @l0: first order $l_0 \geq m$
 * @m: $m$
 *
 * Returns: (array) (element-type NcmSFSphericalHarmonicsK): the coefficients of the
 * steps from $l_0$ to $\ell_\mathrm{max} - 1$ at @m, in order
 */

