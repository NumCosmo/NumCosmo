/***************************************************************************
 *            nc_hipert_gw.c
 *
 *  Fri December 09 11:25:16 2016
 *  Copyright  2016  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_hipert_gw.c
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
 * SECTION:nc_hipert_gw
 * @title: NcHIPertGW
 * @short_description: Perturbation object for gravitational wave mode.
 *
 * This object provides the computation of the gravitational wave mode for the
 * cosmological perturbations. It solves the equation of motion for the gauge invariant
 * variable $h$.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_spline_cubic_notaknot.h"
#include "perturbations/nc_hipert_gw.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <cvode/cvode.h>
#include <nvector/nvector_serial.h>
#include <gsl/gsl_roots.h>
#endif /* NUMCOSMO_GIR_SCAN */

struct _NcHIPertGW
{
  NcmCSQ1D parent_instance;
  gdouble k;
};

G_DEFINE_INTERFACE (NcHIPertIGW, nc_hipert_igw, G_TYPE_OBJECT)
G_DEFINE_TYPE (NcHIPertGW, nc_hipert_gw, NCM_TYPE_CSQ1D)

static gdouble _nc_hipert_igw_eval_unit (NcHIPertIGW *igw);

static void
nc_hipert_igw_default_init (NcHIPertIGWInterface *iface)
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

typedef struct _NcHIPertGWArg
{
  NcHICosmo *cosmo;
  NcHIPertGW *pgw;
} NcHIPertGWArg;

static void
nc_hipert_gw_init (NcHIPertGW *pgw)
{
  pgw->k = 0.0;
}

static void
_nc_hipert_gw_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcHIPertGW *pgw = NC_HIPERT_GW (object);

  g_return_if_fail (NC_IS_HIPERT_GW (object));

  switch (prop_id)
  {
    case PROP_K:
      nc_hipert_gw_set_k (pgw, g_value_get_double (value));
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_nc_hipert_gw_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  /* NcHIPertGW *pgw = NC_HIPERT_GW (object); */
  g_return_if_fail (NC_IS_HIPERT_GW (object));

  switch (prop_id)
  {
    case PROP_K:
      g_value_set_double (value, nc_hipert_gw_get_k (NC_HIPERT_GW (object)));
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_nc_hipert_gw_dispose (GObject *object)
{
  /*NcHIPertGW *pgw = NC_HIPERT_GW (object);*/

  /* Chain up : end */
  G_OBJECT_CLASS (nc_hipert_gw_parent_class)->dispose (object);
}

static void
_nc_hipert_gw_finalize (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (nc_hipert_gw_parent_class)->finalize (object);
}

static gdouble _nc_hipert_gw_eval_xi (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t);
static gdouble _nc_hipert_gw_eval_F1 (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t);
static gdouble _nc_hipert_gw_eval_nu (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t);
static gdouble _nc_hipert_gw_eval_m (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t);
static void _nc_hipert_gw_prepare (NcmCSQ1D *csq1d, NcmModel *model);

static void
nc_hipert_gw_class_init (NcHIPertGWClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);
  NcmCSQ1DClass *csq1d_class = NCM_CSQ1D_CLASS (klass);

  object_class->set_property = &_nc_hipert_gw_set_property;
  object_class->get_property = &_nc_hipert_gw_get_property;
  object_class->dispose      = &_nc_hipert_gw_dispose;
  object_class->finalize     = &_nc_hipert_gw_finalize;

  g_object_class_install_property (object_class,
                                   PROP_K,
                                   g_param_spec_double ("k",
                                                        NULL,
                                                        "Wave number",
                                                        0.0, G_MAXDOUBLE, 1.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  csq1d_class->eval_xi = &_nc_hipert_gw_eval_xi;
  csq1d_class->eval_F1 = &_nc_hipert_gw_eval_F1;
  csq1d_class->eval_nu = &_nc_hipert_gw_eval_nu;
  csq1d_class->eval_m  = &_nc_hipert_gw_eval_m;
  csq1d_class->prepare = &_nc_hipert_gw_prepare;
}

static gdouble
_nc_hipert_gw_eval_xi (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t)
{
  NcHIPertGW *pgw = NC_HIPERT_GW (csq1d);
  const gdouble k = pgw->k;

  return nc_hipert_igw_eval_xi (NC_HIPERT_IGW (model), t, k);
}

static gdouble
_nc_hipert_gw_eval_F1 (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t)
{
  NcHIPertGW *pgw = NC_HIPERT_GW (csq1d);
  const gdouble k = pgw->k;

  return nc_hipert_igw_eval_F1 (NC_HIPERT_IGW (model), t, k);
}

static gdouble
_nc_hipert_gw_eval_nu (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t)
{
  NcHIPertGW *pgw = NC_HIPERT_GW (csq1d);
  const gdouble k = pgw->k;

  return nc_hipert_igw_eval_nu (NC_HIPERT_IGW (model), t, k);
}

static gdouble
_nc_hipert_gw_eval_m (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t)
{
  NcHIPertGW *pgw = NC_HIPERT_GW (csq1d);
  const gdouble k = pgw->k;

  return nc_hipert_igw_eval_m (NC_HIPERT_IGW (model), t, k);
}

static void
_nc_hipert_gw_prepare (NcmCSQ1D *csq1d, NcmModel *model)
{
  g_assert (NC_IS_HIPERT_IGW (model));
}

/**
 * nc_hipert_igw_eval_xi:
 * @igw: a #NcHIPertIGW
 * @tau: $\tau$
 * @k: $k$
 *
 * Computes the value of $\xi = \ln(m\nu)$
 *
 * Returns: $\xi$.
 */

gdouble
nc_hipert_igw_eval_xi (NcHIPertIGW *igw, const gdouble tau, const gdouble k)
{
  return NC_HIPERT_IGW_GET_IFACE (igw)->eval_xi (igw, tau, k);
}

/**
 * nc_hipert_igw_eval_F1:
 * @igw: a #NcHIPertIGW
 * @tau: $\tau$
 * @k: $k$
 *
 * Computes the value of $F_1 = \dot{\xi}/(2\nu)$.
 *
 * Returns: $F_1$.
 */

gdouble
nc_hipert_igw_eval_F1 (NcHIPertIGW *igw, const gdouble tau, const gdouble k)
{
  return NC_HIPERT_IGW_GET_IFACE (igw)->eval_F1 (igw, tau, k);
}

/**
 * nc_hipert_igw_eval_nu:
 * @igw: a #NcHIPertIGW
 * @tau: $\tau$
 * @k: $k$
 *
 * Computes the value of $\nu$.
 *
 * Returns: $\nu$.
 */
gdouble
nc_hipert_igw_eval_nu (NcHIPertIGW *igw, const gdouble tau, const gdouble k)
{
  return NC_HIPERT_IGW_GET_IFACE (igw)->eval_nu (igw, tau, k);
}

/**
 * nc_hipert_igw_eval_m:
 * @igw: a #NcHIPertIGW
 * @tau: $\tau$
 * @k: $k$
 *
 * Computes the value of $m$.
 *
 * Returns: $m$.
 */
gdouble
nc_hipert_igw_eval_m (NcHIPertIGW *igw, const gdouble tau, const gdouble k)
{
  return NC_HIPERT_IGW_GET_IFACE (igw)->eval_m (igw, tau, k);
}

/**
 * nc_hipert_igw_eval_unit:
 * @igw: a #NcHIPertIGW
 *
 * FIXME
 *
 * Returns: FIXME.
 */
gdouble
nc_hipert_igw_eval_unit (NcHIPertIGW *igw)
{
  return NC_HIPERT_IGW_GET_IFACE (igw)->eval_unit (igw);
}

/**
 * nc_hipert_igw_eval_x:
 * @igw: a #NcHIPertIGW
 * @tau: $\tau$
 *
 * Evaluates the value of $x = a_0 / a$ at a given time $\tau$.
 */
gdouble
nc_hipert_igw_eval_x (NcHIPertIGW *igw, const gdouble tau)
{
  return NC_HIPERT_IGW_GET_IFACE (igw)->eval_x (igw, tau);
}

/**
 * nc_hipert_gw_new:
 *
 * Creates a new #NcHIPertGW object.
 *
 * Returns: (transfer full): a new #NcHIPertGW.
 */
NcHIPertGW *
nc_hipert_gw_new (void)
{
  NcHIPertGW *pgw = g_object_new (NC_TYPE_HIPERT_GW,
                                  NULL);

  return pgw;
}

/**
 * nc_hipert_gw_ref:
 * @pgw: a #NcHIPertGW
 *
 * Increases the reference count of @pgw.
 *
 * Returns: (transfer full): @pgw.
 */
NcHIPertGW *
nc_hipert_gw_ref (NcHIPertGW *pgw)
{
  return g_object_ref (pgw);
}

/**
 * nc_hipert_gw_free:
 * @pgw: a #NcHIPertGW
 *
 * Decreases the reference count of @pgw.
 *
 */
void
nc_hipert_gw_free (NcHIPertGW *pgw)
{
  g_object_unref (pgw);
}

/**
 * nc_hipert_gw_clear:
 * @pgw: a #NcHIPertGW
 *
 * Decreases the reference count of *@pgw and sets *@pgw to NULL.
 *
 */
void
nc_hipert_gw_clear (NcHIPertGW **pgw)
{
  g_clear_object (pgw);
}

/**
 * nc_hipert_gw_set_k:
 * @pgw: a #NcHIPertGW
 * @k: the mode $k$
 *
 * Sets the mode $k$ for the gravitational wave perturbation mode.
 *
 */
void
nc_hipert_gw_set_k (NcHIPertGW *pgw, const gdouble k)
{
  pgw->k = k;
}

/**
 * nc_hipert_gw_get_k:
 * @pgw: a #NcHIPertGW
 *
 * Returns the mode $k$ for the gravitational wave perturbation mode.
 *
 * Returns: the mode $k$.
 */
gdouble
nc_hipert_gw_get_k (NcHIPertGW *pgw)
{
  return pgw->k;
}

