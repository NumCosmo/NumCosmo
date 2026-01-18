/***************************************************************************
 *            ncm_sbessel_integrator.h
 *
 *  Thu January 09 12:00:00 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * ncm_sbessel_integrator.h
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

#ifndef _NCM_SBESSEL_INTEGRATOR_H_
#define _NCM_SBESSEL_INTEGRATOR_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_vector.h>

G_BEGIN_DECLS

/**
 * NcmSBesselIntegratorF:
 * @user_data: (closure): user data
 * @x: the value at which to evaluate the function
 *
 * Function to be integrated with spherical Bessel functions.
 *
 * Returns: the function value at @x
 */
typedef gdouble (*NcmSBesselIntegratorF) (gpointer user_data, gdouble x);

#define NCM_TYPE_SBESSEL_INTEGRATOR (ncm_sbessel_integrator_get_type ())

G_DECLARE_DERIVABLE_TYPE (NcmSBesselIntegrator, ncm_sbessel_integrator, NCM, SBESSEL_INTEGRATOR, GObject)

struct _NcmSBesselIntegratorClass
{
  /*< private >*/
  GObjectClass parent_class;

  void (*prepare) (NcmSBesselIntegrator *sbi);
  gdouble (*integrate_ell) (NcmSBesselIntegrator *sbi, NcmSBesselIntegratorF F, gdouble a, gdouble b, gint ell, gpointer user_data);
  void (*integrate) (NcmSBesselIntegrator *sbi, NcmSBesselIntegratorF F, gdouble a, gdouble b, NcmVector *result, gpointer user_data);

  /* Padding to allow 18 virtual functions without breaking ABI. */
  gpointer padding[15];
};

NcmSBesselIntegrator *ncm_sbessel_integrator_ref (NcmSBesselIntegrator *sbi);
void ncm_sbessel_integrator_free (NcmSBesselIntegrator *sbi);
void ncm_sbessel_integrator_clear (NcmSBesselIntegrator **sbi);

guint ncm_sbessel_integrator_get_lmin (NcmSBesselIntegrator *sbi);
guint ncm_sbessel_integrator_get_lmax (NcmSBesselIntegrator *sbi);
void ncm_sbessel_integrator_set_lmin (NcmSBesselIntegrator *sbi, guint lmin);
void ncm_sbessel_integrator_set_lmax (NcmSBesselIntegrator *sbi, guint lmax);

void ncm_sbessel_integrator_prepare (NcmSBesselIntegrator *sbi);
gdouble ncm_sbessel_integrator_integrate_ell (NcmSBesselIntegrator *sbi, NcmSBesselIntegratorF F, gdouble a, gdouble b, gint ell, gpointer user_data);
void ncm_sbessel_integrator_integrate (NcmSBesselIntegrator *sbi, NcmSBesselIntegratorF F, gdouble a, gdouble b, NcmVector *result, gpointer user_data);

gdouble ncm_sbessel_integrator_integrate_gaussian_ell (NcmSBesselIntegrator *sbi, gdouble center, gdouble std, gdouble k, gdouble a, gdouble b, gint ell);
void ncm_sbessel_integrator_integrate_gaussian (NcmSBesselIntegrator *sbi, gdouble center, gdouble std, gdouble k, gdouble a, gdouble b, NcmVector *result);

gdouble ncm_sbessel_integrator_integrate_rational_ell (NcmSBesselIntegrator *sbi, gdouble center, gdouble std, gdouble k, gdouble a, gdouble b, gint ell);
void ncm_sbessel_integrator_integrate_rational (NcmSBesselIntegrator *sbi, gdouble center, gdouble std, gdouble k, gdouble a, gdouble b, NcmVector *result);

G_END_DECLS

#endif /* _NCM_SBESSEL_INTEGRATOR_H_ */

