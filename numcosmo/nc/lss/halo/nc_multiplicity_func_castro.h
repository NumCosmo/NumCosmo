/***************************************************************************
 *            nc_multiplicity_func_castro.h
 *
 *  Thu Aug 14 10:00:00 2026
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

#ifndef _NC_MULTIPLICITY_FUNC_CASTRO_H_
#define _NC_MULTIPLICITY_FUNC_CASTRO_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/nc/lss/halo/nc_multiplicity_func.h>

G_BEGIN_DECLS

#define NC_TYPE_MULTIPLICITY_FUNC_CASTRO             (nc_multiplicity_func_castro_get_type ())

G_DECLARE_FINAL_TYPE (NcMultiplicityFuncCastro, nc_multiplicity_func_castro, NC, MULTIPLICITY_FUNC_CASTRO, NcMultiplicityFunc)

/**
 * NcMultiplicityFuncCastroModel:
 * @NC_MULTIPLICITY_FUNC_CASTRO_MODEL_C23: $\Lambda(\nu)$CDM calibration, Castro et al. (2023)
 * @NC_MULTIPLICITY_FUNC_CASTRO_MODEL_C25: dynamical dark energy calibration, Castro et al. (2025)
 *
 * Selects which of the two calibrations to evaluate. They share the whole
 * coefficient structure but are distinct fits: @NC_MULTIPLICITY_FUNC_CASTRO_MODEL_C25
 * does not reduce to @NC_MULTIPLICITY_FUNC_CASTRO_MODEL_C23 when the dark energy
 * term vanishes, since it carries an extra shape parameter.
 *
 * See [Castro et al. (2023)](https://arxiv.org/abs/2208.02174) and
 * [Castro et al. (2025)](https://arxiv.org/abs/2504.07608).
 *
 */
typedef enum _NcMultiplicityFuncCastroModel
{
  NC_MULTIPLICITY_FUNC_CASTRO_MODEL_C23 = 0,
  NC_MULTIPLICITY_FUNC_CASTRO_MODEL_C25,
  /* < private > */
  NC_MULTIPLICITY_FUNC_CASTRO_MODEL_LEN, /*< skip >*/
} NcMultiplicityFuncCastroModel;

/**
 * NcMultiplicityFuncCastroHaloFinder:
 * @NC_MULTIPLICITY_FUNC_CASTRO_HALO_FINDER_ROCKSTAR: ROCKSTAR
 * @NC_MULTIPLICITY_FUNC_CASTRO_HALO_FINDER_AHF: AHF
 * @NC_MULTIPLICITY_FUNC_CASTRO_HALO_FINDER_SUBFIND: SUBFIND
 * @NC_MULTIPLICITY_FUNC_CASTRO_HALO_FINDER_VELOCIRAPTOR: VELOCIraptor
 *
 * Selects the halo finder whose best-fit parameters are used. This applies to
 * @NC_MULTIPLICITY_FUNC_CASTRO_MODEL_C23 only; the Castro et al. (2025)
 * calibration provides a single parameter set and ignores this property.
 *
 */
typedef enum _NcMultiplicityFuncCastroHaloFinder
{
  NC_MULTIPLICITY_FUNC_CASTRO_HALO_FINDER_ROCKSTAR = 0,
  NC_MULTIPLICITY_FUNC_CASTRO_HALO_FINDER_AHF,
  NC_MULTIPLICITY_FUNC_CASTRO_HALO_FINDER_SUBFIND,
  NC_MULTIPLICITY_FUNC_CASTRO_HALO_FINDER_VELOCIRAPTOR,
  /* < private > */
  NC_MULTIPLICITY_FUNC_CASTRO_HALO_FINDER_LEN, /*< skip >*/
} NcMultiplicityFuncCastroHaloFinder;

NcMultiplicityFuncCastro *nc_multiplicity_func_castro_new (void);
NcMultiplicityFuncCastro *nc_multiplicity_func_castro_new_full (NcMultiplicityFuncCastroModel model, NcMultiplicityFuncCastroHaloFinder halo_finder);
NcMultiplicityFuncCastro *nc_multiplicity_func_castro_ref (NcMultiplicityFuncCastro *mc);

void nc_multiplicity_func_castro_free (NcMultiplicityFuncCastro *mc);
void nc_multiplicity_func_castro_clear (NcMultiplicityFuncCastro **mc);

void nc_multiplicity_func_castro_set_model (NcMultiplicityFuncCastro *mc, NcMultiplicityFuncCastroModel model);
NcMultiplicityFuncCastroModel nc_multiplicity_func_castro_get_model (NcMultiplicityFuncCastro *mc);

void nc_multiplicity_func_castro_set_halo_finder (NcMultiplicityFuncCastro *mc, NcMultiplicityFuncCastroHaloFinder halo_finder);
NcMultiplicityFuncCastroHaloFinder nc_multiplicity_func_castro_get_halo_finder (NcMultiplicityFuncCastro *mc);

gdouble nc_multiplicity_func_castro_eval_full (NcMultiplicityFuncCastro *mc, NcHICosmo *cosmo, gdouble sigma, gdouble dlnsigma_dlnR, gdouble z);
gdouble nc_multiplicity_func_castro_eval_lnf (NcMultiplicityFuncCastro *mc, NcHICosmo *cosmo, gdouble sigma, gdouble dlnsigma_dlnR, gdouble z);

gdouble nc_multiplicity_func_castro_delta_c (NcMultiplicityFuncCastro *mc, NcHICosmo *cosmo, gdouble z);
gdouble nc_multiplicity_func_castro_z_ta (NcMultiplicityFuncCastro *mc, gdouble z);

G_END_DECLS

#endif /* _NC_MULTIPLICITY_FUNC_CASTRO_H_ */

