/***************************************************************************
 *            nc_multiplicity_func_bhattacharya.h
 *
 *  Mon May 18 14:26:00 2026
 *  Copyright  2026  Cinthia Nunes de Lima / Henrique C. N. Lettieri
 *  <cinthia.nlima@gmail.com> <henrique.cnl@hotmail.com>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) 2026  Cinthia Nunes de Lima / Henrique C. N. Lettieri
 * <cinthia.nlima@gmail.com> <henrique.cnl@hotmail.com>
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

#ifndef _NC_MULTIPLICITY_FUNC_BHATTACHARYA_H_
#define _NC_MULTIPLICITY_FUNC_BHATTACHARYA_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/lss/nc_multiplicity_func.h>

G_BEGIN_DECLS

#define NC_TYPE_MULTIPLICITY_FUNC_BHATTACHARYA             (nc_multiplicity_func_bhattacharya_get_type ())

G_DECLARE_FINAL_TYPE (NcMultiplicityFuncBhattacharya, nc_multiplicity_func_bhattacharya, NC, MULTIPLICITY_FUNC_BHATTACHARYA, NcMultiplicityFunc)


/**
 * NcMultiplicityFuncBhattacharyaConvention:
 * @NC_MULTIPLICITY_FUNC_BHATTACHARYA_CONVENTION_BHATTACHARYA2011: original fit, $a(z) = a_0 (1+z)^{-0.01}$
 * @NC_MULTIPLICITY_FUNC_BHATTACHARYA_CONVENTION_HEITMANN2019: Outer Rim form, $a$ constant in redshift
 *
 * Selects the redshift evolution of the shape parameter $a$. The two papers state
 * the same functional form with the same $z = 0$ parameters; they differ only in
 * whether $a$ evolves with redshift. See arXiv:1005.2239 (Bhattacharya et al. 2011)
 * and arXiv:1904.11970 (Heitmann et al. 2019).
 *
 */
typedef enum _NcMultiplicityFuncBhattacharyaConvention
{
  NC_MULTIPLICITY_FUNC_BHATTACHARYA_CONVENTION_BHATTACHARYA2011 = 0,
  NC_MULTIPLICITY_FUNC_BHATTACHARYA_CONVENTION_HEITMANN2019,
  /* < private > */
  NC_MULTIPLICITY_FUNC_BHATTACHARYA_CONVENTION_LEN, /*< skip >*/
} NcMultiplicityFuncBhattacharyaConvention;

NcMultiplicityFuncBhattacharya *nc_multiplicity_func_bhattacharya_new (void);
NcMultiplicityFuncBhattacharya *nc_multiplicity_func_bhattacharya_new_full (NcMultiplicityFuncBhattacharyaConvention convention);
NcMultiplicityFuncBhattacharya *nc_multiplicity_func_bhattacharya_ref (NcMultiplicityFuncBhattacharya *mbt);

void nc_multiplicity_func_bhattacharya_free (NcMultiplicityFuncBhattacharya *mbt);
void nc_multiplicity_func_bhattacharya_clear (NcMultiplicityFuncBhattacharya **mbt);

void nc_multiplicity_func_bhattacharya_set_convention (NcMultiplicityFuncBhattacharya *mbt, NcMultiplicityFuncBhattacharyaConvention convention);
NcMultiplicityFuncBhattacharyaConvention nc_multiplicity_func_bhattacharya_get_convention (const NcMultiplicityFuncBhattacharya *mbt);

void nc_multiplicity_func_bhattacharya_set_A (NcMultiplicityFuncBhattacharya *mbt, gdouble A);
gdouble nc_multiplicity_func_bhattacharya_get_A (const NcMultiplicityFuncBhattacharya *mbt);

void nc_multiplicity_func_bhattacharya_set_a (NcMultiplicityFuncBhattacharya *mbt, gdouble a);
gdouble nc_multiplicity_func_bhattacharya_get_a (const NcMultiplicityFuncBhattacharya *mbt);

void nc_multiplicity_func_bhattacharya_set_p (NcMultiplicityFuncBhattacharya *mbt, gdouble p);
gdouble nc_multiplicity_func_bhattacharya_get_p (const NcMultiplicityFuncBhattacharya *mbt);

void nc_multiplicity_func_bhattacharya_set_q (NcMultiplicityFuncBhattacharya *mbt, gdouble q);
gdouble nc_multiplicity_func_bhattacharya_get_q (const NcMultiplicityFuncBhattacharya *mbt);

void nc_multiplicity_func_bhattacharya_set_delta_c (NcMultiplicityFuncBhattacharya *mbt, gdouble delta_c);
gdouble nc_multiplicity_func_bhattacharya_get_delta_c (const NcMultiplicityFuncBhattacharya *mbt);

G_END_DECLS

#endif /* _NC_MULTIPLICITY_FUNC_BHATTACHARYA_H_ */

