/* -*- Mode: C; indent-tabs-mode: nil; c-basic-offset: 2; tab-width: 2 -*-  */
/***************************************************************************
 *            nc_hipert_comp_pb.h
 *
 *  Fri October 13 11:10:39 2017
 *  Copyright  2017  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * nc_hipert_comp_pb.h
 * Copyright (C) 2017 Sandro Dias Pinto Vitenti <sandro@isoftware.com.br>
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

#ifndef _NC_HIPERT_COMP_PB_H_
#define _NC_HIPERT_COMP_PB_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/perturbations/nc_hipert_comp.h>

G_BEGIN_DECLS

#define NC_TYPE_HIPERT_COMP_PB             (nc_hipert_comp_pb_get_type ())
#define NC_HIPERT_COMP_PB(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_HIPERT_COMP_PB, NcHIPertCompPB))
#define NC_HIPERT_COMP_PB_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_HIPERT_COMP_PB, NcHIPertCompPBClass))
#define NC_IS_HIPERT_COMP_PB(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_HIPERT_COMP_PB))
#define NC_IS_HIPERT_COMP_PB_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_HIPERT_COMP_PB))
#define NC_HIPERT_COMP_PB_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_HIPERT_COMP_PB, NcHIPertCompPBClass))

typedef struct _NcHIPertCompPBClass NcHIPertCompPBClass;
typedef struct _NcHIPertCompPB NcHIPertCompPB;
typedef struct _NcHIPertCompPBPrivate NcHIPertCompPBPrivate;

struct _NcHIPertCompPBClass
{
  /*< private >*/
  NcHIPertCompClass parent_class;
};

struct _NcHIPertCompPB
{
  /*< private >*/
  NcHIPertComp parent_instance;
  NcHIPertCompPBPrivate *priv;
};

/**
 * NcHIPertCompPBVar:
 * @NC_HIPERT_COMP_PB_VAR_DELTA_B: $\delta_b+3\psi$
 * @NC_HIPERT_COMP_PB_VAR_V_B: $\mathcal{V}_b$
 * @NC_HIPERT_COMP_PB_VAR_DELTA_G: $\delta_\gamma+3\psi$
 * @NC_HIPERT_COMP_PB_VAR_V_G: $\mathcal{V}_\gamma$
 * @NC_HIPERT_COMP_PB_VAR_THETA_G: $\theta_\gamma$
 * @NC_HIPERT_COMP_PB_VAR_F_G3: $\mathcal{F}_{\gamma3}$
 * 
 * Photon-baryon component variables
 * 
 */
typedef enum /*< enum,underscore_name=NC_HIPERT_COMP_PB_VAR >*/
{
  NC_HIPERT_COMP_PB_VAR_DELTA_B,
  NC_HIPERT_COMP_PB_VAR_V_B,
  NC_HIPERT_COMP_PB_VAR_DELTA_G,
  NC_HIPERT_COMP_PB_VAR_V_G,
  NC_HIPERT_COMP_PB_VAR_THETA_G,
  NC_HIPERT_COMP_PB_VAR_F_G3,
  /* < private > */
  NC_HIPERT_COMP_PB_VAR_LEN, /*< skip >*/
} NcHIPertCompPBVar;

#define NC_HIPERT_COMP_PB_VAR_F_G(ell) (NC_HIPERT_COMP_PB_VAR_F_G3 + (ell) - 3)

GType nc_hipert_comp_pb_get_type (void) G_GNUC_CONST;

NC_HIPERT_BG_VAR_ID_FUNC_DECL (nc_hipert_comp_pb);

NcHIPertCompPB *nc_hipert_comp_pb_new (void);
NcHIPertCompPB *nc_hipert_comp_pb_ref (NcHIPertCompPB *pb);

void nc_hipert_comp_pb_free (NcHIPertCompPB *pb);
void nc_hipert_comp_pb_clear (NcHIPertCompPB **pb);

void nc_hipert_comp_pb_set_lmax (NcHIPertCompPB *pb, guint lmax);
guint nc_hipert_comp_pb_get_lmax (NcHIPertCompPB *pb);

G_END_DECLS

#endif /* _NC_HIPERT_COMP_PB_H_ */

