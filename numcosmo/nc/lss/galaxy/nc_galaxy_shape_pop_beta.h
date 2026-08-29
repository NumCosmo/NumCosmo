/***************************************************************************
 *            nc_galaxy_shape_pop_beta.h
 *
 *  Thu Jun 19 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_galaxy_shape_pop_beta.h
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
 * with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef _NC_GALAXY_SHAPE_POP_BETA_H_
#define _NC_GALAXY_SHAPE_POP_BETA_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/nc/lss/galaxy/nc_galaxy_shape_pop.h>

G_BEGIN_DECLS

#define NC_TYPE_GALAXY_SHAPE_POP_BETA (nc_galaxy_shape_pop_beta_get_type ())

G_DECLARE_FINAL_TYPE (NcGalaxyShapePopBeta, nc_galaxy_shape_pop_beta, NC, GALAXY_SHAPE_POP_BETA, NcGalaxyShapePop)

/**
 * NcGalaxyShapePopBetaParams:
 * @NC_GALAXY_SHAPE_POP_BETA_ALPHA: shape parameter $\alpha$ of the Beta distribution of $r = |\chi_I|$.
 * @NC_GALAXY_SHAPE_POP_BETA_BETA: shape parameter $\beta$ of the Beta distribution of $r = |\chi_I|$.
 *
 * Beta intrinsic ellipticity model parameters; $\alpha\ge1$, $\beta\ge1$
 * (see the class documentation).
 *
 */
typedef enum /*< enum,underscore_name=NC_GALAXY_SHAPE_POP_BETA_PARAMS,prefix=NC_GALAXY_SHAPE_POP_BETA >*/
{
  NC_GALAXY_SHAPE_POP_BETA_ALPHA = 0,
  NC_GALAXY_SHAPE_POP_BETA_BETA,
  /* < private > */
  NC_GALAXY_SHAPE_POP_BETA_SPARAM_LEN, /*< skip >*/
} NcGalaxyShapePopBetaParams;

#define NC_GALAXY_SHAPE_POP_BETA_DEFAULT_ALPHA (1.55)
#define NC_GALAXY_SHAPE_POP_BETA_DEFAULT_BETA (1.55)
#define NC_GALAXY_SHAPE_POP_BETA_DEFAULT_PARAMS_ABSTOL (0.0)

NcGalaxyShapePopBeta *nc_galaxy_shape_pop_beta_new (void);
NcGalaxyShapePopBeta *nc_galaxy_shape_pop_beta_ref (NcGalaxyShapePopBeta *gspb);

void nc_galaxy_shape_pop_beta_free (NcGalaxyShapePopBeta *gspb);
void nc_galaxy_shape_pop_beta_clear (NcGalaxyShapePopBeta **gspb);

gdouble nc_galaxy_shape_pop_beta_get_e_rms (NcGalaxyShapePopBeta *gspb);
gdouble nc_galaxy_shape_pop_beta_get_mode (NcGalaxyShapePopBeta *gspb);

G_END_DECLS

#endif /* _NC_GALAXY_SHAPE_POP_BETA_H_ */

