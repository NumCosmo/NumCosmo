/***************************************************************************
 *            ncm_sbessel_integrator_gl.c
 *
 *  Thu January 09 12:00:00 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * ncm_sbessel_integrator_gl.c
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
 * NcmSBesselIntegratorGL:
 *
 * Gauss-Legendre based spherical Bessel function integrator.
 *
 * This class implements integration of functions multiplied by spherical
 * Bessel functions using Gauss-Legendre quadrature.
 *
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_sbessel_integrator_gl.h"
#include "math/ncm_sf_sbessel.h"
#include "math/ncm_c.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_bessel.h>
#endif /* NUMCOSMO_GIR_SCAN */

typedef struct _NcmSBesselIntegratorGLParams
{
  NcmSBesselIntegratorGL *sbigl;
  NcmSBesselIntegratorF F;
  gpointer user_data;
  gint ell;
} NcmSBesselIntegratorGLParams;

struct _NcmSBesselIntegratorGL
{
  /*< private >*/
  NcmSBesselIntegrator parent_instance;
  gsl_integration_workspace *ws;
  gsl_integration_glfixed_table *gl_table;
  guint npts;
  gdouble margin;
  gdouble nosc;
  NcmSBesselIntegratorGLParams params;
};

enum
{
  PROP_0,
  PROP_NPTS,
  PROP_MARGIN,
  PROP_NOSC,
};

G_DEFINE_TYPE (NcmSBesselIntegratorGL, ncm_sbessel_integrator_gl, NCM_TYPE_SBESSEL_INTEGRATOR)

static void
ncm_sbessel_integrator_gl_init (NcmSBesselIntegratorGL *sbigl)
{
  sbigl->ws       = gsl_integration_workspace_alloc (1000);
  sbigl->npts     = 10;
  sbigl->gl_table = gsl_integration_glfixed_table_alloc (sbigl->npts);
  sbigl->margin   = 5.0;
  sbigl->nosc     = 6.0;

  sbigl->params.sbigl     = sbigl;
  sbigl->params.F         = NULL;
  sbigl->params.user_data = NULL;
  sbigl->params.ell       = 0;
}

static void
_ncm_sbessel_integrator_gl_dispose (GObject *object)
{
  NcmSBesselIntegratorGL *sbigl = NCM_SBESSEL_INTEGRATOR_GL (object);

  if (sbigl->ws != NULL)
  {
    gsl_integration_workspace_free (sbigl->ws);
    sbigl->ws = NULL;
  }

  if (sbigl->gl_table != NULL)
  {
    gsl_integration_glfixed_table_free (sbigl->gl_table);
    sbigl->gl_table = NULL;
  }

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_sbessel_integrator_gl_parent_class)->dispose (object);
}

static void
_ncm_sbessel_integrator_gl_finalize (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (ncm_sbessel_integrator_gl_parent_class)->finalize (object);
}

static void _ncm_sbessel_integrator_gl_prepare (NcmSBesselIntegrator *sbi);
static gdouble _ncm_sbessel_integrator_gl_integrate_ell (NcmSBesselIntegrator *sbi, NcmSBesselIntegratorF F, gpointer user_data, gdouble a, gdouble b, gint ell);

static void
_ncm_sbessel_integrator_gl_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmSBesselIntegratorGL *sbigl = NCM_SBESSEL_INTEGRATOR_GL (object);

  g_return_if_fail (NCM_IS_SBESSEL_INTEGRATOR_GL (object));

  switch (prop_id)
  {
    case PROP_NPTS:
      ncm_sbessel_integrator_gl_set_npts (sbigl, g_value_get_uint (value));
      break;
    case PROP_MARGIN:
      ncm_sbessel_integrator_gl_set_margin (sbigl, g_value_get_double (value));
      break;
    case PROP_NOSC:
      ncm_sbessel_integrator_gl_set_nosc (sbigl, g_value_get_double (value));
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_ncm_sbessel_integrator_gl_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmSBesselIntegratorGL *sbigl = NCM_SBESSEL_INTEGRATOR_GL (object);

  g_return_if_fail (NCM_IS_SBESSEL_INTEGRATOR_GL (object));

  switch (prop_id)
  {
    case PROP_NPTS:
      g_value_set_uint (value, ncm_sbessel_integrator_gl_get_npts (sbigl));
      break;
    case PROP_MARGIN:
      g_value_set_double (value, ncm_sbessel_integrator_gl_get_margin (sbigl));
      break;
    case PROP_NOSC:
      g_value_set_double (value, ncm_sbessel_integrator_gl_get_nosc (sbigl));
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
ncm_sbessel_integrator_gl_class_init (NcmSBesselIntegratorGLClass *klass)
{
  GObjectClass *object_class              = G_OBJECT_CLASS (klass);
  NcmSBesselIntegratorClass *parent_class = NCM_SBESSEL_INTEGRATOR_CLASS (klass);

  object_class->set_property = &_ncm_sbessel_integrator_gl_set_property;
  object_class->get_property = &_ncm_sbessel_integrator_gl_get_property;
  object_class->dispose      = &_ncm_sbessel_integrator_gl_dispose;
  object_class->finalize     = &_ncm_sbessel_integrator_gl_finalize;

  /**
   * NcmSBesselIntegratorGL:npts:
   *
   * Number of points for Gauss-Legendre quadrature in oscillatory region.
   */
  g_object_class_install_property (object_class,
                                   PROP_NPTS,
                                   g_param_spec_uint ("npts",
                                                      NULL,
                                                      "Number of GL quadrature points",
                                                      1, G_MAXUINT, 10,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmSBesselIntegratorGL:margin:
   *
   * Safety margin beyond turning point for non-oscillatory region.
   */
  g_object_class_install_property (object_class,
                                   PROP_MARGIN,
                                   g_param_spec_double ("margin",
                                                        NULL,
                                                        "Safety margin beyond turning point",
                                                        0.0, G_MAXDOUBLE, 5.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmSBesselIntegratorGL:nosc:
   *
   * Number of oscillations per panel width in oscillatory region.
   */
  g_object_class_install_property (object_class,
                                   PROP_NOSC,
                                   g_param_spec_double ("nosc",
                                                        NULL,
                                                        "Number of oscillations per panel",
                                                        0.1, G_MAXDOUBLE, 6.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  parent_class->prepare       = &_ncm_sbessel_integrator_gl_prepare;
  parent_class->integrate_ell = &_ncm_sbessel_integrator_gl_integrate_ell;
}

/* Integrand: f(x) * j_ell(x) */
static gdouble
_ncm_sbessel_integrator_gl_integrand (gdouble x, gpointer user_data)
{
  NcmSBesselIntegratorGLParams *params = (NcmSBesselIntegratorGLParams *) user_data;
  gdouble fx                           = params->F (params->user_data, x);
  gdouble jl                           = gsl_sf_bessel_jl (params->ell, x);

  return fx * jl;
}

static void
_ncm_sbessel_integrator_gl_prepare (NcmSBesselIntegrator *sbi)
{
  /* Preparation could reinitialize workspace if needed */
}

static gdouble
_ncm_sbessel_integrator_gl_integrate_ell (NcmSBesselIntegrator *sbi,
                                          NcmSBesselIntegratorF F,
                                          gpointer user_data,
                                          gdouble a, gdouble b,
                                          gint ell)
{
  NcmSBesselIntegratorGL *sbigl = NCM_SBESSEL_INTEGRATOR_GL (sbi);
  gdouble result                = 0.0;
  gdouble abserr;

  /* Update integrand params */
  sbigl->params.F         = F;
  sbigl->params.user_data = user_data;
  sbigl->params.ell       = ell;

  /* Turning point and split */
  const gdouble x_tp  = sqrt (ell * (ell + 1.0));
  const gdouble x_mid = x_tp + sbigl->margin;

  /* ----------------------------
   * Region 1: non-oscillatory
   * ---------------------------- */
  if (a < x_mid)
  {
    gsl_function G;

    G.function = &_ncm_sbessel_integrator_gl_integrand;
    G.params   = &sbigl->params;

    const gdouble a1 = a;
    const gdouble b1 = GSL_MIN (x_mid, b);

    gdouble I1 = 0.0;

    gsl_integration_qag (&G, a1, b1,
                         0.0, 1.0e-7,
                         1000,
                         GSL_INTEG_GAUSS61,
                         sbigl->ws, &I1, &abserr);

    result += I1;

    if (b1 == b)
      return result;
  }

  /* ----------------------------
   * Region 2: oscillatory
   * ---------------------------- */

  const gdouble a2 = GSL_MAX (x_mid, a);
  const gdouble b2 = b;

  gsl_function G;

  G.function = &_ncm_sbessel_integrator_gl_integrand;
  G.params   = &sbigl->params;

  /* Base oscillation-determined panel width */
  const gdouble dx0 = M_PI / sbigl->nosc;

  /* Controls how dx grows with x */
  const gdouble x_scale = 1000.0;

  /* Early-termination threshold */
  const gdouble eps_tail = 1.0e-8;

  gdouble x = a2;

  while (x < b2)
  {
    /* ----------------------------
     * Amplitude-based early exit
     * ---------------------------- */
    /* Amplitude-based early termination */
    gdouble fx = F (user_data, x);

    if (fabs (fx) / x < eps_tail)
    {
      gdouble jl1_x = ncm_sf_sbessel (ell + 1, x);
      gdouble fb    = F (user_data, b);
      gdouble jl1_b = ncm_sf_sbessel (ell + 1, b);
      gdouble tail  = fb * jl1_b - fx * jl1_x;

      result += tail;
      break;
    }

    /* ----------------------------
     * dx that grows smoothly w.r.t x
     * ---------------------------- */
    gdouble dx   = GSL_MIN (dx0 * (1.0 + x / x_scale), 4.0 * M_PI);
    gdouble x_hi = x + dx;

    if (x_hi > b2)
      x_hi = b2;

    /* One GL fixed panel */
    gdouble panel_result =
      gsl_integration_glfixed (&G, x, x_hi, sbigl->gl_table);

    result += panel_result;
    x       = x_hi;
  }

  return result;
}

/**
 * ncm_sbessel_integrator_gl_new:
 * @lmin: minimum multipole
 * @lmax: maximum multipole
 *
 * Creates a new #NcmSBesselIntegratorGL.
 *
 * Returns: (transfer full): a new #NcmSBesselIntegratorGL
 */
NcmSBesselIntegratorGL *
ncm_sbessel_integrator_gl_new (guint lmin, guint lmax)
{
  NcmSBesselIntegratorGL *sbigl = g_object_new (NCM_TYPE_SBESSEL_INTEGRATOR_GL,
                                                "lmin", lmin,
                                                "lmax", lmax,
                                                NULL);

  return sbigl;
}

/**
 * ncm_sbessel_integrator_gl_ref:
 * @sbigl: a #NcmSBesselIntegratorGL
 *
 * Increases the reference count of @sbigl by one.
 *
 * Returns: (transfer full): @sbigl
 */
NcmSBesselIntegratorGL *
ncm_sbessel_integrator_gl_ref (NcmSBesselIntegratorGL *sbigl)
{
  return g_object_ref (sbigl);
}

/**
 * ncm_sbessel_integrator_gl_free:
 * @sbigl: a #NcmSBesselIntegratorGL
 *
 * Decreases the reference count of @sbigl by one.
 *
 */
void
ncm_sbessel_integrator_gl_free (NcmSBesselIntegratorGL *sbigl)
{
  g_object_unref (sbigl);
}

/**
 * ncm_sbessel_integrator_gl_clear:
 * @sbigl: a #NcmSBesselIntegratorGL
 *
 * If @sbigl is different from NULL, decreases the reference count of
 * @sbigl by one and sets @sbigl to NULL.
 *
 */
void
ncm_sbessel_integrator_gl_clear (NcmSBesselIntegratorGL **sbigl)
{
  g_clear_object (sbigl);
}

/**
 * ncm_sbessel_integrator_gl_set_npts:
 * @sbigl: a #NcmSBesselIntegratorGL
 * @npts: number of Gauss-Legendre quadrature points
 *
 * Sets the number of Gauss-Legendre quadrature points used in the oscillatory region.
 * This will reallocate the GL table if the value changes.
 *
 */
void
ncm_sbessel_integrator_gl_set_npts (NcmSBesselIntegratorGL *sbigl, guint npts)
{
  g_return_if_fail (NCM_IS_SBESSEL_INTEGRATOR_GL (sbigl));
  g_return_if_fail (npts > 0);

  if (npts != sbigl->npts)
  {
    sbigl->npts = npts;

    if (sbigl->gl_table != NULL)
    {
      gsl_integration_glfixed_table_free (sbigl->gl_table);
      sbigl->gl_table = gsl_integration_glfixed_table_alloc (sbigl->npts);
    }
  }
}

/**
 * ncm_sbessel_integrator_gl_get_npts:
 * @sbigl: a #NcmSBesselIntegratorGL
 *
 * Gets the number of Gauss-Legendre quadrature points.
 *
 * Returns: the number of quadrature points
 */
guint
ncm_sbessel_integrator_gl_get_npts (NcmSBesselIntegratorGL *sbigl)
{
  g_return_val_if_fail (NCM_IS_SBESSEL_INTEGRATOR_GL (sbigl), 0);

  return sbigl->npts;
}

/**
 * ncm_sbessel_integrator_gl_set_margin:
 * @sbigl: a #NcmSBesselIntegratorGL
 * @margin: safety margin beyond turning point
 *
 * Sets the safety margin beyond the turning point for the non-oscillatory region.
 *
 */
void
ncm_sbessel_integrator_gl_set_margin (NcmSBesselIntegratorGL *sbigl, gdouble margin)
{
  g_return_if_fail (NCM_IS_SBESSEL_INTEGRATOR_GL (sbigl));
  g_return_if_fail (margin >= 0.0);

  sbigl->margin = margin;
}

/**
 * ncm_sbessel_integrator_gl_get_margin:
 * @sbigl: a #NcmSBesselIntegratorGL
 *
 * Gets the safety margin beyond the turning point.
 *
 * Returns: the margin value
 */
gdouble
ncm_sbessel_integrator_gl_get_margin (NcmSBesselIntegratorGL *sbigl)
{
  g_return_val_if_fail (NCM_IS_SBESSEL_INTEGRATOR_GL (sbigl), 0.0);

  return sbigl->margin;
}

/**
 * ncm_sbessel_integrator_gl_set_nosc:
 * @sbigl: a #NcmSBesselIntegratorGL
 * @nosc: number of oscillations per panel
 *
 * Sets the number of oscillations per panel width in the oscillatory region.
 *
 */
void
ncm_sbessel_integrator_gl_set_nosc (NcmSBesselIntegratorGL *sbigl, gdouble nosc)
{
  g_return_if_fail (NCM_IS_SBESSEL_INTEGRATOR_GL (sbigl));
  g_return_if_fail (nosc > 0.0);

  sbigl->nosc = nosc;
}

/**
 * ncm_sbessel_integrator_gl_get_nosc:
 * @sbigl: a #NcmSBesselIntegratorGL
 *
 * Gets the number of oscillations per panel.
 *
 * Returns: the nosc value
 */
gdouble
ncm_sbessel_integrator_gl_get_nosc (NcmSBesselIntegratorGL *sbigl)
{
  g_return_val_if_fail (NCM_IS_SBESSEL_INTEGRATOR_GL (sbigl), 0.0);

  return sbigl->nosc;
}

