/***************************************************************************
 *            nc_hipert_em.h
 *
 *  Sat March 16 10:53:30 2024
 *  Copyright  2024  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_hipert_em.h
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

#ifndef _NC_HIPERT_EM_H_
#define _NC_HIPERT_EM_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_ode_spline.h>
#include <numcosmo/nc_hicosmo.h>
#include <numcosmo/math/ncm_csq1d.h>

G_BEGIN_DECLS

#define NC_TYPE_HIPERT_IEM (nc_hipert_iem_get_type ())
#define NC_TYPE_HIPERT_EM (nc_hipert_em_get_type ())

G_DECLARE_INTERFACE (NcHIPertIEM, nc_hipert_iem, NC, HIPERT_IEM, GObject)
G_DECLARE_FINAL_TYPE (NcHIPertEM, nc_hipert_em, NC, HIPERT_EM, NcmCSQ1D)

struct _NcHIPertIEMInterface
{
  /*< private >*/
  GTypeInterface parent;
  gdouble (*eval_xi) (NcHIPertIEM *iem, const gdouble tau, const gdouble k);
  gdouble (*eval_F1) (NcHIPertIEM *iem, const gdouble tau, const gdouble k);
  gdouble (*eval_nu) (NcHIPertIEM *iem, const gdouble tau, const gdouble k);
  gdouble (*eval_m) (NcHIPertIEM *iem, const gdouble tau, const gdouble k);
  gdouble (*eval_unit) (NcHIPertIEM *iem);
  gdouble (*eval_x) (NcHIPertIEM *iem, const gdouble tau);
  gdouble (*eval_lapse) (NcHIPertIEM *iem, const gdouble tau);

  /* Padding to allow 18 virtual functions without breaking ABI. */
  gpointer padding[11];
};

/**
 * NcHIPertEMVars:
 * @NC_HIPERT_EM_RE_H: $\text{Re}(\zeta)$
 * @NC_HIPERT_EM_IM_H: $\text{Im}(\zeta)$
 * @NC_HIPERT_EM_RE_PH: $\text{Re}(P_\zeta)$
 * @NC_HIPERT_EM_IM_PH: $\text{Im}(P_\zeta)$
 *
 * Perturbation variables enumerator.
 *
 */
typedef enum /*< enum,underscore_name=NC_HIPERT_EM_VARS  >*/
{
  NC_HIPERT_EM_RE_H = 0,
  NC_HIPERT_EM_IM_H,
  NC_HIPERT_EM_RE_PH,
  NC_HIPERT_EM_IM_PH,
} NcHIPertEMVars;


gdouble nc_hipert_iem_eval_xi (NcHIPertIEM *iem, const gdouble tau, const gdouble k);
gdouble nc_hipert_iem_eval_F1 (NcHIPertIEM *iem, const gdouble tau, const gdouble k);
gdouble nc_hipert_iem_eval_nu (NcHIPertIEM *iem, const gdouble tau, const gdouble k);
gdouble nc_hipert_iem_eval_m (NcHIPertIEM *iem, const gdouble tau, const gdouble k);
gdouble nc_hipert_iem_eval_unit (NcHIPertIEM *iem);
gdouble nc_hipert_iem_eval_x (NcHIPertIEM *iem, const gdouble tau);

NcHIPertEM *nc_hipert_em_new (void);
NcHIPertEM *nc_hipert_em_ref (NcHIPertEM *pem);
void nc_hipert_em_free (NcHIPertEM *pem);
void nc_hipert_em_clear (NcHIPertEM **pem);

void nc_hipert_em_set_k (NcHIPertEM *pem, const gdouble k);
gdouble nc_hipert_em_get_k (NcHIPertEM *pem);

void nc_hipert_em_eval_PE_PB (NcHIPertEM *pem, NcmModel *model, const gdouble tau, gdouble *PE, gdouble *PB);

G_END_DECLS

#endif /* _NC_HIPERT_EM_H_ */

