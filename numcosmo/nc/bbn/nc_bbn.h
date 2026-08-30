/***************************************************************************
 *            nc_bbn.h
 *
 *  Fri August 29 20:30:00 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_bbn.h
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

#ifndef _NC_BBN_H_
#define _NC_BBN_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/nc/background/nc_hicosmo.h>
#include <numcosmo/ncm/model/ncm_model.h>
#include <numcosmo/ncm/model/ncm_mset.h>

G_BEGIN_DECLS

#define NC_TYPE_BBN (nc_bbn_get_type ())

G_DECLARE_DERIVABLE_TYPE (NcBBN, nc_bbn, NC, BBN, NcmModel);

struct _NcBBNClass
{
  /*< private >*/
  NcmModelClass parent_class;
  gdouble (*Yp_4He) (NcBBN *bbn, NcHICosmo *cosmo);
  void (*get_domain) (NcBBN *bbn, gdouble *wb_lb, gdouble *wb_ub, gdouble *DNeff_lb, gdouble *DNeff_ub);
  gdouble (*DH) (NcBBN *bbn, NcHICosmo *cosmo);
  gdouble (*He3H) (NcBBN *bbn, NcHICosmo *cosmo);
  gdouble (*Li7H) (NcBBN *bbn, NcHICosmo *cosmo);

  /* Padding to allow adding up to 13 more virtual functions without breaking ABI. */
  gpointer padding[13];
};

/**
 * NcBBNImpl:
 * @NC_BBN_IMPL_Yp_4He: primordial Helium-4 mass fraction
 * @NC_BBN_IMPL_get_domain: range of $(\omega_b, \Delta N_\mathrm{eff})$ the model answers for
 * @NC_BBN_IMPL_DH: Deuterium abundance $\mathrm{D}/\mathrm{H}$
 * @NC_BBN_IMPL_He3H: Helium-3 abundance ${}^3\mathrm{He}/\mathrm{H}$
 * @NC_BBN_IMPL_Li7H: Lithium-7 abundance ${}^7\mathrm{Li}/\mathrm{H}$
 *
 * Implementation flags for #NcBBN. Only @NC_BBN_IMPL_Yp_4He and
 * @NC_BBN_IMPL_get_domain are required of every implementation: the tables
 * shipped with NumCosmo tabulate $Y_p$ alone, so the remaining abundances are
 * optional and a caller must check ncm_model_check_impl_opt() before asking for them.
 *
 */
typedef enum _NcBBNImpl
{
  NC_BBN_IMPL_Yp_4He = 0,
  NC_BBN_IMPL_get_domain,
  NC_BBN_IMPL_DH,
  NC_BBN_IMPL_He3H,
  NC_BBN_IMPL_Li7H,
  /* < private > */
} NcBBNImpl;

#define NC_BBN_IMPL_ALL NCM_MODEL_CLASS_IMPL_ALL

NCM_MSET_MODEL_DECLARE_ID (nc_bbn);

NcBBN *nc_bbn_ref (NcBBN *bbn);
void nc_bbn_free (NcBBN *bbn);
void nc_bbn_clear (NcBBN **bbn);

gdouble nc_bbn_Yp_4He (NcBBN *bbn, NcHICosmo *cosmo);
void nc_bbn_get_domain (NcBBN *bbn, gdouble *wb_lb, gdouble *wb_ub, gdouble *DNeff_lb, gdouble *DNeff_ub);
gdouble nc_bbn_DH (NcBBN *bbn, NcHICosmo *cosmo);
gdouble nc_bbn_He3H (NcBBN *bbn, NcHICosmo *cosmo);
gdouble nc_bbn_Li7H (NcBBN *bbn, NcHICosmo *cosmo);

void nc_bbn_check_domain (NcBBN *bbn, const gdouble wb, const gdouble DNeff);

#define NC_BBN_DEFAULT_PARAMS_ABSTOL (0.0)

G_END_DECLS

#endif /* _NC_BBN_H_ */

