/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 2; tab-width: 2 -*-  */

/***************************************************************************
 *            nc_galaxy_wl_proj.h
 *
 *  Sun Aug 02 11:13:56 2020
 *  Copyright  2020  Sandro Dias Pinto Vitenti & Mariana Penna Lima
 *  <sandro@isoftware.com.br>, <pennalima@gmail.com>
 ****************************************************************************/
/*
 * nc_galaxy_wl_proj.h
 * Copyright (C) 2020 Sandro Dias Pinto Vitenti <sandro@isoftware.com.br>
 * Copyright (C) 2020 Mariana Penna Lima <pennalima@gmail.com>
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

#ifndef _NC_GALAXY_WL_PROJ_H_
#define _NC_GALAXY_WL_PROJ_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_rng.h>
#include <numcosmo/math/ncm_matrix.h>
#include <numcosmo/lss/nc_galaxy_wl_dist.h>

G_BEGIN_DECLS

#define NC_TYPE_GALAXY_WL_PROJ             (nc_galaxy_wl_proj_get_type ())
#define NC_GALAXY_WL_PROJ(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_GALAXY_WL_PROJ, NcGalaxyWLProj))
#define NC_GALAXY_WL_PROJ_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_GALAXY_WL_PROJ, NcGalaxyWLProjClass))
#define NC_IS_GALAXY_WL_PROJ(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_GALAXY_WL_PROJ))
#define NC_IS_GALAXY_WL_PROJ_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_GALAXY_WL_PROJ))
#define NC_GALAXY_WL_PROJ_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_GALAXY_WL_PROJ, NcGalaxyWLProjClass))

typedef struct _NcGalaxyWLProjClass NcGalaxyWLProjClass;
typedef struct _NcGalaxyWLProj NcGalaxyWLProj;
typedef struct _NcGalaxyWLProjPrivate NcGalaxyWLProjPrivate;

struct _NcGalaxyWLProjClass
{
  /*< private >*/
  NcGalaxyWLDistClass parent_class;
};

struct _NcGalaxyWLProj
{
  /*< private >*/
  NcGalaxyWLDist parent_instance;
  NcGalaxyWLProjPrivate *priv;
};

/**
 * NcGalaxyWLProjPos:
 * @NC_GALAXY_WL_PROJ_POS_ANG: FIXME
 * @NC_GALAXY_WL_PROJ_POS_R: FIXME
 *
 * FIXME
 *
 */
typedef enum _NcGalaxyWLProjPos
{
  NC_GALAXY_WL_PROJ_POS_ANG,
  NC_GALAXY_WL_PROJ_POS_R,
  /* < private > */
  NC_GALAXY_WL_PROJ_POS_LEN, /*< skip >*/
} NcGalaxyWLProjPos;

GType nc_galaxy_wl_proj_get_type (void) G_GNUC_CONST;

NcGalaxyWLProj *nc_galaxy_wl_proj_new (NcGalaxyWLProjPos pos);
NcGalaxyWLProj *nc_galaxy_wl_proj_ref (NcGalaxyWLProj *gwlp);

void nc_galaxy_wl_proj_free (NcGalaxyWLProj *gwlp);
void nc_galaxy_wl_proj_clear (NcGalaxyWLProj **gwlp);

void nc_galaxy_wl_proj_set_pos (NcGalaxyWLProj *gwlp, NcGalaxyWLProjPos pos);
NcGalaxyWLProjPos nc_galaxy_wl_proj_get_pos (NcGalaxyWLProj *gwlp);

void nc_galaxy_wl_proj_set_obs (NcGalaxyWLProj *gwlp, NcmMatrix *obs);
NcmMatrix *nc_galaxy_wl_proj_peek_obs (NcGalaxyWLProj *gwlp);

G_END_DECLS

#endif /* _NC_GALAXY_WL_PROJ_H_ */

