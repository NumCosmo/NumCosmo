/***************************************************************************
 *            nc_hipert_adiab.c
 *
 *  Tue June 03 17:20:42 2014
 *  Copyright  2014  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * nc_hipert_adiab.c
 * Copyright (C) 2014 Sandro Dias Pinto Vitenti <sandro@isoftware.com.br>
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
 * SECTION:nc_hipert_adiab
 * @title: NcHIPertAdiab
 * @short_description: Perturbation object for adiabatic mode only.
 *
 * This object provides the computation of the adiabatic mode for the cosmological
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
#include "perturbations/nc_hipert_adiab.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <cvode/cvode.h>
#include <nvector/nvector_serial.h>
#include <gsl/gsl_roots.h>
#endif /* NUMCOSMO_GIR_SCAN */

G_DEFINE_INTERFACE (NcHIPertIAdiab, nc_hipert_iadiab, G_TYPE_OBJECT);

static gdouble _nc_hipert_iadiab_eval_powspec_factor (NcHIPertIAdiab *iad);

static void
nc_hipert_iadiab_default_init (NcHIPertIAdiabInterface *iface)
{
  iface->eval_mnu            = NULL;
  iface->eval_nu             = NULL;
  iface->eval_dlnmnu         = NULL;
  iface->eval_system         = NULL;

	iface->eval_powspec_factor = &_nc_hipert_iadiab_eval_powspec_factor;
}

G_DEFINE_TYPE (NcHIPertAdiab, nc_hipert_adiab, NCM_TYPE_HOAA);

enum {
  PROP_0,
  PROP_SIZE,
};

typedef struct _NcHIPertAdiabArg
{
  NcHICosmo *cosmo;
  NcHIPertAdiab *pa;
} NcHIPertAdiabArg;

static void
nc_hipert_adiab_init (NcHIPertAdiab *pa)
{
}

static void
_nc_hipert_adiab_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  /* NcHIPertAdiab *pa = NC_HIPERT_ADIAB (object); */
  g_return_if_fail (NC_IS_HIPERT_ADIAB (object));

  switch (prop_id)
  {
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_hipert_adiab_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  /* NcHIPertAdiab *pa = NC_HIPERT_ADIAB (object); */
  g_return_if_fail (NC_IS_HIPERT_ADIAB (object));

  switch (prop_id)
  {
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_hipert_adiab_dispose (GObject *object)
{
  /*NcHIPertAdiab *pa = NC_HIPERT_ADIAB (object);*/

  /* Chain up : end */
  G_OBJECT_CLASS (nc_hipert_adiab_parent_class)->dispose (object);
}

static void
_nc_hipert_adiab_finalize (GObject *object)
{

  /* Chain up : end */
  G_OBJECT_CLASS (nc_hipert_adiab_parent_class)->finalize (object);
}

static gdouble _nc_hipert_adiab_eval_mnu (NcmHOAA *hoaa, NcmModel *model, const gdouble t, const gdouble k);
static gdouble _nc_hipert_adiab_eval_nu (NcmHOAA *hoaa, NcmModel *model, const gdouble t, const gdouble k);
static gdouble _nc_hipert_adiab_eval_dlnmnu (NcmHOAA *hoaa, NcmModel *model, const gdouble t, const gdouble k);
static void _nc_hipert_adiab_eval_system (NcmHOAA *hoaa, NcmModel *model, const gdouble t, const gdouble k, gdouble *nu, gdouble *dlnmnu, gdouble *Vnu);


static gdouble _nc_hipert_adiab_eval_powspec_factor (NcmHOAA *hoaa, NcmModel *model);

static void _nc_hipert_adiab_prepare (NcmHOAA *hoaa, NcmModel *model);

static void
nc_hipert_adiab_class_init (NcHIPertAdiabClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);
  NcmHOAAClass *hoaa_class   = NCM_HOAA_CLASS (klass);

  object_class->set_property = &_nc_hipert_adiab_set_property;
  object_class->get_property = &_nc_hipert_adiab_get_property;
  object_class->dispose      = &_nc_hipert_adiab_dispose;
  object_class->finalize     = &_nc_hipert_adiab_finalize;

  hoaa_class->eval_mnu    = &_nc_hipert_adiab_eval_mnu;
  hoaa_class->eval_nu     = &_nc_hipert_adiab_eval_nu;
  hoaa_class->eval_dlnmnu = &_nc_hipert_adiab_eval_dlnmnu;
  hoaa_class->eval_V      = NULL;
  hoaa_class->eval_system = &_nc_hipert_adiab_eval_system;
  hoaa_class->prepare     = &_nc_hipert_adiab_prepare;

	hoaa_class->eval_powspec_factor = &_nc_hipert_adiab_eval_powspec_factor;
}

static gdouble 
_nc_hipert_adiab_eval_mnu (NcmHOAA *hoaa, NcmModel *model, const gdouble t, const gdouble k)
{
  return nc_hipert_iadiab_eval_mnu (NC_HIPERT_IADIAB (model), t, k);  
}

static gdouble 
_nc_hipert_adiab_eval_nu (NcmHOAA *hoaa, NcmModel *model, const gdouble t, const gdouble k)
{
  return nc_hipert_iadiab_eval_nu (NC_HIPERT_IADIAB (model), t, k);
}

static gdouble 
_nc_hipert_adiab_eval_dlnmnu (NcmHOAA *hoaa, NcmModel *model, const gdouble t, const gdouble k)
{
  return nc_hipert_iadiab_eval_dlnmnu (NC_HIPERT_IADIAB (model), t, k);
}

static void 
_nc_hipert_adiab_eval_system (NcmHOAA *hoaa, NcmModel *model, const gdouble t, const gdouble k, gdouble *nu, gdouble *dlnmnu, gdouble *Vnu)
{
  Vnu[0] = 0.0;  
  nc_hipert_iadiab_eval_system (NC_HIPERT_IADIAB (model), t, k, nu, dlnmnu);
}

static void 
_nc_hipert_adiab_prepare (NcmHOAA *hoaa, NcmModel *model)
{
  g_assert (NC_IS_HIPERT_IADIAB (model));
}

static gdouble
_nc_hipert_adiab_eval_powspec_factor (NcmHOAA *hoaa, NcmModel *model)
{
	return nc_hipert_iadiab_eval_powspec_factor (NC_HIPERT_IADIAB (model));
}

static gdouble 
_nc_hipert_iadiab_eval_powspec_factor (NcHIPertIAdiab *iad) 
{ 
	g_assert (NC_IS_HICOSMO (iad));
	{
		NcHICosmo *cosmo    = NC_HICOSMO (iad);
		const gdouble RH_lp = nc_hicosmo_RH_planck (cosmo);
		return 4.0 / (3.0 * M_PI * gsl_pow_2 (RH_lp));
	}
}

/**
 * nc_hipert_iadiab_eval_mnu:
 * @iad: a #NcHIPertIAdiab
 * @tau: $\tau$
 * @k: $k$
 *
 * FIXME
 *
 * Returns: FIXME.
 */
/**
 * nc_hipert_iadiab_eval_nu:
 * @iad: a #NcHIPertIAdiab
 * @tau: $\tau$
 * @k: $k$
 *
 * FIXME
 *
 * Returns: FIXME.
 */
/**
 * nc_hipert_iadiab_eval_dlnmnu:
 * @iad: a #NcHIPertIAdiab
 * @tau: $\tau$
 * @k: $k$
 *
 * FIXME
 *
 * Returns: FIXME.
 */
/**
 * nc_hipert_iadiab_eval_system:
 * @iad: a #NcHIPertIAdiab
 * @tau: $\tau$
 * @k: $k$
 * @nu: (out): $\nu$
 * @dlnmnu: (out): FIXME
 *
 * FIXME
 *
 */

/**
 * nc_hipert_adiab_new:
 *
 * Creates a new #NcHIPertAdiab object.
 *
 * Returns: (transfer full): a new #NcHIPertAdiab.
 */
NcHIPertAdiab *
nc_hipert_adiab_new (void)
{
  NcHIPertAdiab *pa = g_object_new (NC_TYPE_HIPERT_ADIAB,
                                    "opt", NCM_HOAA_OPT_DLNMNU_ONLY,
                                    NULL);
  return pa;
}

/**
 * nc_hipert_adiab_ref:
 * @pa: a #NcHIPertAdiab.
 *
 * Increases the reference count of @pa.
 *
 * Returns: (transfer full): @pa.
 */
NcHIPertAdiab *
nc_hipert_adiab_ref (NcHIPertAdiab *pa)
{
  return g_object_ref (pa);
}

/**
 * nc_hipert_adiab_free:
 * @pa: a #NcHIPertAdiab.
 *
 * Decreases the reference count of @pa.
 *
 */
void
nc_hipert_adiab_free (NcHIPertAdiab *pa)
{
  g_object_unref (pa);
}

/**
 * nc_hipert_adiab_clear:
 * @pa: a #NcHIPertAdiab.
 *
 * Decreases the reference count of *@pa and sets *@pa to NULL.
 *
 */
void
nc_hipert_adiab_clear (NcHIPertAdiab **pa)
{
  g_clear_object (pa);
}
