/***************************************************************************
 *            nc_hipert_gw.c
 *
 *  Fri December 09 11:25:16 2016
 *  Copyright  2016  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_hipert_gw.c
 * Copyright (C) 2014 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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
 * @short_description: Perturbation object for gwatic mode only.
 *
 * This object provides the computation of the gwatic mode for the cosmological
 * perturbations. It solves the equation of motion for the gauge invariant variable
 * (see [Vitenti (2013)][XVitenti2013] for notation and details)
 * $$\zeta \equiv \Psi - \frac{2\bar{K}}{\kappa(\bar{\rho} + \bar{p})} + H\mathcal{V}.$$
 * Its conjugated momentum is give by
 * \begin{split}
 * P_\zeta &= \frac{2\bar{D}^2_\bar{K}\Psi}{x^3H},
 * \end{split}
 *
 * The equations of motion in their first order form are
 * \begin{align}
 * \zeta^\prime &= \frac{P_\zeta}{m_\zeta}, \\\\
 * P_\zeta^\prime &= -m_\zeta\mu_\zeta^2\zeta.
 * \end{align}
 * The mass $m_\zeta$ and the frequency $\mu_\zeta$ are defined by
 * \begin{align}
 * m_\zeta     &= \frac{3\Delta_\bar{K}(\bar{\rho} + \bar{p})}{\rho_\text{crit0} N x^3 c_s^2 E^2}, \\\\
 * \mu_\zeta^2 &= x^2N^2c_s^2k^2,
 * \end{align}
 * where $\bar{\rho} + \bar{p}$ is the background total energy density plus pressure,
 * $E^2 = H^2/H_0^2$ is the dimensionless Hubble function squared (nc_hicosmo_E2()), $c_s^2$ the speed of sound,
 * $N$ is the lapse function that in this case (using $\alpha$ as time variable) is $N \equiv \vert{}E\vert^{-1}$, $\rho_\text{crit0}$
 * is the critical density today defined by $\rho_\text{crit0} \equiv 3H_0^2/\kappa$ and $$\Delta_\bar{K} \equiv \frac{k^2}{k^2 + \Omega_{k0}}.$$
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
};

G_DEFINE_INTERFACE (NcHIPertIGW, nc_hipert_igw, G_TYPE_OBJECT)
G_DEFINE_TYPE (NcHIPertGW, nc_hipert_gw, NCM_TYPE_CSQ1D)

static gdouble _nc_hipert_igw_eval_powspec_factor (NcHIPertIGW *igw, const gdouble k);

static void
nc_hipert_igw_default_init (NcHIPertIGWInterface *iface)
{
  iface->eval_xi     = NULL;
  iface->eval_F1     = NULL;
  iface->eval_nu     = NULL;
  iface->eval_m      = NULL;
  iface->eval_system = NULL;

  iface->eval_powspec_factor = &_nc_hipert_igw_eval_powspec_factor;
}

enum
{
  PROP_0,
  PROP_SIZE,
};

typedef struct _NcHIPertGWArg
{
  NcHICosmo *cosmo;
  NcHIPertGW *pa;
} NcHIPertGWArg;

static void
nc_hipert_gw_init (NcHIPertGW *pa)
{
}

static void
_nc_hipert_gw_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  /* NcHIPertGW *pa = NC_HIPERT_GW (object); */
  g_return_if_fail (NC_IS_HIPERT_GW (object));

  switch (prop_id)
  {
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_hipert_gw_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  /* NcHIPertGW *pa = NC_HIPERT_GW (object); */
  g_return_if_fail (NC_IS_HIPERT_GW (object));

  switch (prop_id)
  {
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_hipert_gw_dispose (GObject *object)
{
  /*NcHIPertGW *pa = NC_HIPERT_GW (object);*/

  /* Chain up : end */
  G_OBJECT_CLASS (nc_hipert_gw_parent_class)->dispose (object);
}

static void
_nc_hipert_gw_finalize (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (nc_hipert_gw_parent_class)->finalize (object);
}

static gdouble _nc_hipert_gw_eval_xi (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t, const gdouble k);
static gdouble _nc_hipert_gw_eval_F1 (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t, const gdouble k);
static gdouble _nc_hipert_gw_eval_nu (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t, const gdouble k);
static gdouble _nc_hipert_gw_eval_m (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t, const gdouble k);
static void _nc_hipert_gw_eval_system (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t, const gdouble k, gdouble *nu, gdouble *dlnmnu, gdouble *Vnu);

static gdouble _nc_hipert_adiab_eval_powspec_factor (NcmCSQ1D *csq1d, NcmModel *model, const gdouble k);

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

  csq1d_class->eval_xi = &_nc_hipert_gw_eval_xi;
  csq1d_class->eval_F1 = &_nc_hipert_gw_eval_F1;
  csq1d_class->eval_nu = &_nc_hipert_gw_eval_nu;
  csq1d_class->eval_m  = &_nc_hipert_gw_eval_m;
  csq1d_class->prepare = &_nc_hipert_gw_prepare;

  csq1d_class->eval_powspec_factor = &_nc_hipert_adiab_eval_powspec_factor;
}

static gdouble
_nc_hipert_gw_eval_xi (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t, const gdouble k)
{
  return nc_hipert_igw_eval_xi (NC_HIPERT_IGW (model), t, k);
}

static gdouble
_nc_hipert_gw_eval_F1 (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t, const gdouble k)
{
  return nc_hipert_igw_eval_F1 (NC_HIPERT_IGW (model), t, k);
}

static gdouble
_nc_hipert_gw_eval_nu (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t, const gdouble k)
{
  return nc_hipert_igw_eval_nu (NC_HIPERT_IGW (model), t, k);
}

static gdouble
_nc_hipert_gw_eval_m (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t, const gdouble k)
{
  return nc_hipert_igw_eval_m (NC_HIPERT_IGW (model), t, k);
}

static void
_nc_hipert_gw_eval_system (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t, const gdouble k, gdouble *nu, gdouble *dlnmnu, gdouble *Vnu)
{
  Vnu[0] = 0.0;
  nc_hipert_igw_eval_system (NC_HIPERT_IGW (model), t, k, nu, dlnmnu);
}

static void
_nc_hipert_gw_prepare (NcmCSQ1D *csq1d, NcmModel *model)
{
  g_assert (NC_IS_HIPERT_IGW (model));
}

static gdouble
_nc_hipert_adiab_eval_powspec_factor (NcmCSQ1D *csq1d, NcmModel *model, const gdouble k)
{
  return nc_hipert_igw_eval_powspec_factor (NC_HIPERT_IGW (model), k);
}

static gdouble
_nc_hipert_igw_eval_powspec_factor (NcHIPertIGW *igw, const gdouble k)
{
  g_assert (NC_IS_HICOSMO (igw));
  {
    NcHICosmo *cosmo    = NC_HICOSMO (igw);
    const gdouble RH_lp = nc_hicosmo_RH_planck (cosmo);

    return 16.0 / (M_PI * gsl_pow_2 (RH_lp));
  }
}

/**
 * nc_hipert_igw_eval_xi:
 * @igw: a #NcHIPertIGW
 * @tau: $\tau$
 * @k: $k$
 *
 * FIXME
 *
 * Returns: FIXME.
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
 * FIXME
 *
 * Returns: FIXME.
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
 * FIXME
 *
 * Returns: FIXME.
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
 * FIXME
 *
 * Returns: FIXME.
 */
gdouble
nc_hipert_igw_eval_m (NcHIPertIGW *igw, const gdouble tau, const gdouble k)
{
  return NC_HIPERT_IGW_GET_IFACE (igw)->eval_m (igw, tau, k);
}

/**
 * nc_hipert_igw_eval_system:
 * @igw: a #NcHIPertIGW
 * @tau: $\tau$
 * @k: $k$
 * @nu: (out): $\nu$
 * @dlnmnu: (out): FIXME
 *
 * FIXME
 *
 */

void
nc_hipert_igw_eval_system (NcHIPertIGW *igw, const gdouble tau, const gdouble k, gdouble *nu, gdouble *dlnmnu)
{
  NC_HIPERT_IGW_GET_IFACE (igw)->eval_system (igw, tau, k, nu, dlnmnu);
}

/**
 * nc_hipert_igw_eval_powspec_factor:
 * @igw: a #NcHIPertIGW
 * @k: $k$
 *
 * FIXME
 *
 * Returns: FIXME.
 */
gdouble
nc_hipert_igw_eval_powspec_factor (NcHIPertIGW *igw, const gdouble k)
{
  return NC_HIPERT_IGW_GET_IFACE (igw)->eval_powspec_factor (igw, k);
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
  NcHIPertGW *pa = g_object_new (NC_TYPE_HIPERT_GW,
                                 NULL);

  return pa;
}

/**
 * nc_hipert_gw_ref:
 * @pa: a #NcHIPertGW.
 *
 * Increases the reference count of @pa.
 *
 * Returns: (transfer full): @pa.
 */
NcHIPertGW *
nc_hipert_gw_ref (NcHIPertGW *pa)
{
  return g_object_ref (pa);
}

/**
 * nc_hipert_gw_free:
 * @pa: a #NcHIPertGW.
 *
 * Decreases the reference count of @pa.
 *
 */
void
nc_hipert_gw_free (NcHIPertGW *pa)
{
  g_object_unref (pa);
}

/**
 * nc_hipert_gw_clear:
 * @pa: a #NcHIPertGW.
 *
 * Decreases the reference count of *@pa and sets *@pa to NULL.
 *
 */
void
nc_hipert_gw_clear (NcHIPertGW **pa)
{
  g_clear_object (pa);
}

