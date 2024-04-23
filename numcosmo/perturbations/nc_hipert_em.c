/***************************************************************************
 *            nc_hipert_em.h
 *
 *  Sat March 16 10:53:30 2024
 *  Copyright  2024  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_hipert_em.c
 * Copyright (C) 2024 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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
 * SECTION:nc_hipert_em
 * @title: NcHIPertEM
 * @short_description: Perturbation object for electromagnetic mode.
 *
 * This object provides the computation of the electromagnetic wave mode for the cosmological
 * perturbations. It solves the equation of motion for the (cosmological) gauge invariant
 * variable $A$.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_spline_cubic_notaknot.h"
#include "perturbations/nc_hipert_em.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <cvode/cvode.h>
#include <nvector/nvector_serial.h>
#include <gsl/gsl_roots.h>
#endif /* NUMCOSMO_GIR_SCAN */

struct _NcHIPertEM
{
  NcmCSQ1D parent_instance;
  gdouble k;
};

G_DEFINE_INTERFACE (NcHIPertIEM, nc_hipert_iem, G_TYPE_OBJECT)
G_DEFINE_TYPE (NcHIPertEM, nc_hipert_em, NCM_TYPE_CSQ1D)

static gdouble _nc_hipert_iem_eval_unit (NcHIPertIEM *iem);

static void
nc_hipert_iem_default_init (NcHIPertIEMInterface *iface)
{
  iface->eval_xi   = NULL;
  iface->eval_F1   = NULL;
  iface->eval_nu   = NULL;
  iface->eval_m    = NULL;
  iface->eval_unit = NULL;
}

enum
{
  PROP_0,
  PROP_K,
  PROP_SIZE,
};

typedef struct _NcHIPertEMArg
{
  NcHICosmo *cosmo;
  NcHIPertEM *pem;
} NcHIPertEMArg;

static void
nc_hipert_em_init (NcHIPertEM *pem)
{
  pem->k = 0.0;
}

static void
_nc_hipert_em_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcHIPertEM *pem = NC_HIPERT_EM (object);

  g_return_if_fail (NC_IS_HIPERT_EM (object));

  switch (prop_id)
  {
    case PROP_K:
      nc_hipert_em_set_k (pem, g_value_get_double (value));
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_nc_hipert_em_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  /* NcHIPertEM *pem = NC_HIPERT_EM (object); */
  g_return_if_fail (NC_IS_HIPERT_EM (object));

  switch (prop_id)
  {
    case PROP_K:
      g_value_set_double (value, nc_hipert_em_get_k (NC_HIPERT_EM (object)));
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_nc_hipert_em_dispose (GObject *object)
{
  /*NcHIPertEM *pem = NC_HIPERT_EM (object);*/

  /* Chain up : end */
  G_OBJECT_CLASS (nc_hipert_em_parent_class)->dispose (object);
}

static void
_nc_hipert_em_finalize (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (nc_hipert_em_parent_class)->finalize (object);
}

static gdouble _nc_hipert_em_eval_xi (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t);
static gdouble _nc_hipert_em_eval_F1 (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t);
static gdouble _nc_hipert_em_eval_nu (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t);
static gdouble _nc_hipert_em_eval_m (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t);
static void _nc_hipert_em_prepare (NcmCSQ1D *csq1d, NcmModel *model);

static void
nc_hipert_em_class_init (NcHIPertEMClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);
  NcmCSQ1DClass *csq1d_class = NCM_CSQ1D_CLASS (klass);

  object_class->set_property = &_nc_hipert_em_set_property;
  object_class->get_property = &_nc_hipert_em_get_property;
  object_class->dispose      = &_nc_hipert_em_dispose;
  object_class->finalize     = &_nc_hipert_em_finalize;

  g_object_class_install_property (object_class,
                                   PROP_K,
                                   g_param_spec_double ("k",
                                                        NULL,
                                                        "Wave number",
                                                        0.0, G_MAXDOUBLE, 1.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  csq1d_class->eval_xi = &_nc_hipert_em_eval_xi;
  csq1d_class->eval_F1 = &_nc_hipert_em_eval_F1;
  csq1d_class->eval_nu = &_nc_hipert_em_eval_nu;
  csq1d_class->eval_m  = &_nc_hipert_em_eval_m;
  csq1d_class->prepare = &_nc_hipert_em_prepare;
}

static gdouble
_nc_hipert_em_eval_xi (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t)
{
  NcHIPertEM *pem = NC_HIPERT_EM (csq1d);
  const gdouble k = pem->k;

  return nc_hipert_iem_eval_xi (NC_HIPERT_IEM (model), t, k);
}

static gdouble
_nc_hipert_em_eval_F1 (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t)
{
  NcHIPertEM *pem = NC_HIPERT_EM (csq1d);
  const gdouble k = pem->k;

  return nc_hipert_iem_eval_F1 (NC_HIPERT_IEM (model), t, k);
}

static gdouble
_nc_hipert_em_eval_nu (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t)
{
  NcHIPertEM *pem = NC_HIPERT_EM (csq1d);
  const gdouble k = pem->k;

  return nc_hipert_iem_eval_nu (NC_HIPERT_IEM (model), t, k);
}

static gdouble
_nc_hipert_em_eval_m (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t)
{
  NcHIPertEM *pem = NC_HIPERT_EM (csq1d);
  const gdouble k = pem->k;

  return nc_hipert_iem_eval_m (NC_HIPERT_IEM (model), t, k);
}

static void
_nc_hipert_em_prepare (NcmCSQ1D *csq1d, NcmModel *model)
{
  g_assert (NC_IS_HIPERT_IEM (model));
}

/**
 * nc_hipert_iem_eval_xi:
 * @iem: a #NcHIPertIEM
 * @tau: $\tau$
 * @k: $k$
 *
 * Computes the value of $\xi = \ln(m\nu)$
 *
 * Returns: $\xi$.
 */

gdouble
nc_hipert_iem_eval_xi (NcHIPertIEM *iem, const gdouble tau, const gdouble k)
{
  return NC_HIPERT_IEM_GET_IFACE (iem)->eval_xi (iem, tau, k);
}

/**
 * nc_hipert_iem_eval_F1:
 * @iem: a #NcHIPertIEM
 * @tau: $\tau$
 * @k: $k$
 *
 * Computes the value of $F_1 = \dot{\xi}/(2\nu)$.
 *
 * Returns: $F_1$.
 */

gdouble
nc_hipert_iem_eval_F1 (NcHIPertIEM *iem, const gdouble tau, const gdouble k)
{
  return NC_HIPERT_IEM_GET_IFACE (iem)->eval_F1 (iem, tau, k);
}

/**
 * nc_hipert_iem_eval_nu:
 * @iem: a #NcHIPertIEM
 * @tau: $\tau$
 * @k: $k$
 *
 * Computes the value of $\nu$.
 *
 * Returns: $\nu$.
 */
gdouble
nc_hipert_iem_eval_nu (NcHIPertIEM *iem, const gdouble tau, const gdouble k)
{
  return NC_HIPERT_IEM_GET_IFACE (iem)->eval_nu (iem, tau, k);
}

/**
 * nc_hipert_iem_eval_m:
 * @iem: a #NcHIPertIEM
 * @tau: $\tau$
 * @k: $k$
 *
 * Computes the value of $m$.
 *
 * Returns: $m$.
 */
gdouble
nc_hipert_iem_eval_m (NcHIPertIEM *iem, const gdouble tau, const gdouble k)
{
  return NC_HIPERT_IEM_GET_IFACE (iem)->eval_m (iem, tau, k);
}

/**
 * nc_hipert_iem_eval_unit:
 * @iem: a #NcHIPertIEM
 *
 * FIXME
 *
 * Returns: FIXME.
 */
gdouble
nc_hipert_iem_eval_unit (NcHIPertIEM *iem)
{
  return NC_HIPERT_IEM_GET_IFACE (iem)->eval_unit (iem);
}

/**
 * nc_hipert_iem_eval_x:
 * @iem: a #NcHIPertIEM
 * @tau: $\tau$
 *
 * Evaluates the value of $x = a_0 / a$ at a given time $\tau$.
 *
 * Returns: $x$.
 */
gdouble
nc_hipert_iem_eval_x (NcHIPertIEM *iem, const gdouble tau)
{
  return NC_HIPERT_IEM_GET_IFACE (iem)->eval_x (iem, tau);
}

/**
 * nc_hipert_iem_eval_lapse:
 * @iem: a #NcHIPertIEM
 * @tau: $\tau$
 *
 * Evaluates the value of the lapse function
 * $N = \mathrm{d}t/\mathrm{d}\tau$ at a given time $\tau$.
 *
 * Returns: the lapse function $N(\tau)$.
 */
gdouble
nc_hipert_iem_eval_lapse (NcHIPertIEM *iem, const gdouble tau)
{
  return NC_HIPERT_IEM_GET_IFACE (iem)->eval_lapse (iem, tau);
}

/**
 * nc_hipert_em_new:
 *
 * Creates a new #NcHIPertEM object.
 *
 * Returns: (transfer full): a new #NcHIPertEM.
 */
NcHIPertEM *
nc_hipert_em_new (void)
{
  NcHIPertEM *pem = g_object_new (NC_TYPE_HIPERT_EM,
                                  NULL);

  return pem;
}

/**
 * nc_hipert_em_ref:
 * @pem: a #NcHIPertEM
 *
 * Increases the reference count of @pem.
 *
 * Returns: (transfer full): @pem.
 */
NcHIPertEM *
nc_hipert_em_ref (NcHIPertEM *pem)
{
  return g_object_ref (pem);
}

/**
 * nc_hipert_em_free:
 * @pem: a #NcHIPertEM
 *
 * Decreases the reference count of @pem.
 *
 */
void
nc_hipert_em_free (NcHIPertEM *pem)
{
  g_object_unref (pem);
}

/**
 * nc_hipert_em_clear:
 * @pem: a #NcHIPertEM
 *
 * Decreases the reference count of *@pem and sets *@pem to NULL.
 *
 */
void
nc_hipert_em_clear (NcHIPertEM **pem)
{
  g_clear_object (pem);
}

/**
 * nc_hipert_em_set_k:
 * @pem: a #NcHIPertEM
 * @k: the mode $k$
 *
 * Sets the mode $k$ for the gravitational wave perturbation mode.
 *
 */
void
nc_hipert_em_set_k (NcHIPertEM *pem, const gdouble k)
{
  pem->k = k;
}

/**
 * nc_hipert_em_get_k:
 * @pem: a #NcHIPertEM
 *
 * Returns the mode $k$ for the gravitational wave perturbation mode.
 *
 * Returns: the mode $k$.
 */
gdouble
nc_hipert_em_get_k (NcHIPertEM *pem)
{
  return pem->k;
}

/**
 * nc_hipert_em_eval_PE_PB:
 * @pem: a #NcHIPertEM
 * @tau: $\tau$
 * @PE: (out): the electric field power spectrum
 * @PB: (out): the magnetic field power spectrum
 *
 * Evaluates the electric and magnetic field power spectra in
 * units of Gauss squared $G_\mathrm{s}^2$.
 */
void
nc_hipert_em_eval_PE_PB (NcHIPertEM *pem, NcmModel *model, const gdouble tau, gdouble *PE, gdouble *PB)
{
  NcmCSQ1D *csq1d    = NCM_CSQ1D (pem);
  NcHIPertIEM *iem   = NC_HIPERT_IEM (model);
  const gdouble unit = nc_hipert_iem_eval_unit (iem);
  const gdouble m    = ncm_csq1d_eval_m (csq1d, model, tau);
  const gdouble k    = pem->k;
  const gdouble x    = nc_hipert_iem_eval_x (iem, tau);
  const gdouble N    = nc_hipert_iem_eval_lapse (iem, tau);
  const gdouble fact = ncm_c_two_pi_2 ();

  NcmCSQ1DState state;
  gdouble J11, J12, J22;

  ncm_csq1d_eval_at (csq1d, model, tau, &state);

  ncm_csq1d_state_get_J (&state, &J11, &J12, &J22);

  *PE = unit * unit * gsl_pow_4 (x) * gsl_pow_3 (k) * J22 / gsl_pow_2 (N * m) / fact;
  *PB = unit * unit * gsl_pow_4 (x) * gsl_pow_5 (k) * J11 / fact;
}

