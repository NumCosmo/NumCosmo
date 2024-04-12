/***************************************************************************
 *            nc_hipert_adiab.h
 *
 *  Tue June 03 17:20:57 2014
 *  Copyright  2014  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_hipert_adiab.h
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

#ifndef _NC_HIPERT_ADIAB_H_
#define _NC_HIPERT_ADIAB_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_ode_spline.h>
#include <numcosmo/math/ncm_powspec_spline2d.h>
#include <numcosmo/nc_hicosmo.h>
#include <numcosmo/math/ncm_csq1d.h>

G_BEGIN_DECLS

#define NC_TYPE_HIPERT_ADIAB (nc_hipert_adiab_get_type ())
#define NC_TYPE_HIPERT_IADIAB (nc_hipert_iadiab_get_type ())

G_DECLARE_INTERFACE (NcHIPertIAdiab, nc_hipert_iadiab, NC, HIPERT_IADIAB, GObject)
G_DECLARE_FINAL_TYPE (NcHIPertAdiab, nc_hipert_adiab, NC, HIPERT_ADIAB, NcmCSQ1D)

struct _NcHIPertIAdiabInterface
{
  /*< private >*/
  GTypeInterface parent;
  gdouble (*eval_xi) (NcHIPertIAdiab *iad, const gdouble tau, const gdouble k);
  gdouble (*eval_F1) (NcHIPertIAdiab *iad, const gdouble tau, const gdouble k);
  gdouble (*eval_nu) (NcHIPertIAdiab *iad, const gdouble tau, const gdouble k);
  gdouble (*eval_m) (NcHIPertIAdiab *iad, const gdouble tau, const gdouble k);
  gdouble (*eval_unit) (NcHIPertIAdiab *iad);
  gdouble (*eval_x) (NcHIPertIAdiab *iad, const gdouble tau);
  gdouble (*eval_p2Psi) (NcHIPertIAdiab *iad, const gdouble tau, const gdouble k);
  gdouble (*eval_p2drho) (NcHIPertIAdiab *iad, const gdouble tau, const gdouble k);
  gdouble (*eval_lapse) (NcHIPertIAdiab *iad, const gdouble tau);
  gdouble (*eval_tau_hubble) (NcHIPertIAdiab *iad, const gdouble k);
  gdouble (*eval_tau_jeans) (NcHIPertIAdiab *iad, const gdouble k);
  gdouble (*eval_hubble) (NcHIPertIAdiab *iad, const gdouble tau);

  /* Padding to allow 18 virtual functions without breaking ABI. */
  gpointer padding[12];
};


/**
 * NcHIPertAdiabVars:
 * @NC_HIPERT_ADIAB_RE_ZETA: $\text{Re}(\zeta)$
 * @NC_HIPERT_ADIAB_IM_ZETA: $\text{Im}(\zeta)$
 * @NC_HIPERT_ADIAB_RE_PZETA: $\text{Re}(P_\zeta)$
 * @NC_HIPERT_ADIAB_IM_PZETA: $\text{Im}(P_\zeta)$
 *
 * Perturbation variables enumerator.
 *
 */
typedef enum /*< enum,underscore_name=NC_HIPERT_ADIAB_VARS  >*/
{
  NC_HIPERT_ADIAB_RE_ZETA = 0,
  NC_HIPERT_ADIAB_IM_ZETA,
  NC_HIPERT_ADIAB_RE_PZETA,
  NC_HIPERT_ADIAB_IM_PZETA,
} NcHIPertAdiabVars;

gdouble nc_hipert_iadiab_eval_xi (NcHIPertIAdiab *iad, const gdouble tau, const gdouble k);
gdouble nc_hipert_iadiab_eval_F1 (NcHIPertIAdiab *iad, const gdouble tau, const gdouble k);
gdouble nc_hipert_iadiab_eval_nu (NcHIPertIAdiab *iad, const gdouble tau, const gdouble k);
gdouble nc_hipert_iadiab_eval_m (NcHIPertIAdiab *iad, const gdouble tau, const gdouble k);
gdouble nc_hipert_iadiab_eval_unit (NcHIPertIAdiab *iad);
gdouble nc_hipert_iadiab_eval_x (NcHIPertIAdiab *iad, const gdouble tau);
gdouble nc_hipert_iadiab_eval_p2Psi (NcHIPertIAdiab *iad, const gdouble tau, const gdouble k);
gdouble nc_hipert_iadiab_eval_p2drho (NcHIPertIAdiab *iad, const gdouble tau, const gdouble k);
gdouble nc_hipert_iadiab_eval_lapse (NcHIPertIAdiab *iad, const gdouble tau);
gdouble nc_hipert_iadiab_eval_tau_hubble (NcHIPertIAdiab *iad, const gdouble k);
gdouble nc_hipert_iadiab_eval_tau_jeans (NcHIPertIAdiab *iad, const gdouble k);
gdouble nc_hipert_iadiab_eval_hubble (NcHIPertIAdiab *iad, const gdouble tau);

NcHIPertAdiab *nc_hipert_adiab_new (void);
NcHIPertAdiab *nc_hipert_adiab_ref (NcHIPertAdiab *pa);
void nc_hipert_adiab_free (NcHIPertAdiab *pa);
void nc_hipert_adiab_clear (NcHIPertAdiab **pa);

void nc_hipert_adiab_set_k (NcHIPertAdiab *adiab, const gdouble k);
gdouble nc_hipert_adiab_get_k (NcHIPertAdiab *adiab);

gdouble nc_hipert_adiab_eval_cosmic_time (NcHIPertAdiab *adiab, NcmModel *model, const gdouble tau);
gdouble nc_hipert_adiab_eval_delta_critial (NcHIPertAdiab *adiab, NcmModel *model, const gdouble tau);

gdouble nc_hipert_adiab_eval_powspec_zeta_at (NcHIPertAdiab *adiab, NcmModel *model, const gdouble tau);
gdouble nc_hipert_adiab_eval_powspec_Psi_at (NcHIPertAdiab *adiab, NcmModel *model, const gdouble tau);
gdouble nc_hipert_adiab_eval_powspec_drho_at (NcHIPertAdiab *adiab, NcmModel *model, const gdouble tau);

void nc_hipert_adiab_prepare_spectrum (NcHIPertAdiab *adiab, NcmModel *model, GArray *k_array, GArray *tau_array);

NcmPowspecSpline2d *nc_hipert_adiab_eval_powspec_zeta (NcHIPertAdiab *adiab, NcmModel *model);
NcmPowspecSpline2d *nc_hipert_adiab_eval_powspec_Psi (NcHIPertAdiab *adiab, NcmModel *model);
NcmPowspecSpline2d *nc_hipert_adiab_eval_powspec_drho (NcHIPertAdiab *adiab, NcmModel *model);

G_END_DECLS

#endif /* _NC_HIPERT_ADIAB_H_ */

