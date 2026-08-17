/***************************************************************************
 *            nc_halo_bias_castro.h
 *
 *  Thu Aug 14 12:00:00 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * numcosmo
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

#ifndef _NC_HALO_BIAS_CASTRO_H_
#define _NC_HALO_BIAS_CASTRO_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/nc/lss/halo/nc_halo_bias.h>

G_BEGIN_DECLS

#define NC_TYPE_HALO_BIAS_CASTRO             (nc_halo_bias_castro_get_type ())

G_DECLARE_FINAL_TYPE (NcHaloBiasCastro, nc_halo_bias_castro, NC, HALO_BIAS_CASTRO, NcHaloBias)

NcHaloBiasCastro *nc_halo_bias_castro_new (NcHaloMassFunction * mfp);

NcHaloBiasCastro *nc_halo_bias_castro_ref (NcHaloBiasCastro *biasf);
void nc_halo_bias_castro_free (NcHaloBiasCastro *biasf);
void nc_halo_bias_castro_clear (NcHaloBiasCastro **biasf);

void nc_halo_bias_castro_set_A0 (NcHaloBiasCastro *biasf, gdouble A0);
gdouble nc_halo_bias_castro_get_A0 (NcHaloBiasCastro *biasf);

void nc_halo_bias_castro_set_a1 (NcHaloBiasCastro *biasf, gdouble a1);
gdouble nc_halo_bias_castro_get_a1 (NcHaloBiasCastro *biasf);

void nc_halo_bias_castro_set_b1 (NcHaloBiasCastro *biasf, gdouble b1);
gdouble nc_halo_bias_castro_get_b1 (NcHaloBiasCastro *biasf);

void nc_halo_bias_castro_set_b2 (NcHaloBiasCastro *biasf, gdouble b2);
gdouble nc_halo_bias_castro_get_b2 (NcHaloBiasCastro *biasf);

void nc_halo_bias_castro_set_c1 (NcHaloBiasCastro *biasf, gdouble c1);
gdouble nc_halo_bias_castro_get_c1 (NcHaloBiasCastro *biasf);

gdouble nc_halo_bias_castro_pbs (NcHaloBiasCastro *biasf, NcHICosmo *cosmo, gdouble sigma, gdouble dlnsigma_dlnR, gdouble z);
gdouble nc_halo_bias_castro_correction (NcHaloBiasCastro *biasf, NcHICosmo *cosmo, gdouble dlnsigma_dlnR, gdouble z);
gdouble nc_halo_bias_castro_S8 (NcHaloBiasCastro *biasf, NcHICosmo *cosmo);

G_END_DECLS

#endif /* _NC_HALO_BIAS_CASTRO_H_ */

