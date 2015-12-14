/***************************************************************************
 *            nc_hireion.c
 *
 *  Tue December 08 14:51:37 2015
 *  Copyright  2015  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * nc_hireion.c
 * Copyright (C) 2015 Sandro Dias Pinto Vitenti <sandro@isoftware.com.br>
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

/**
 * SECTION:nc_hireion
 * @title: NcHIReion
 * @short_description: Abstract class for implementing homogeneous and isotropic reionization models.
 *
 * FIXME
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "nc_hireion.h"
#include "math/ncm_serialize.h"
#include "math/ncm_cfg.h"

G_DEFINE_ABSTRACT_TYPE (NcHIReion, nc_hireion, NCM_TYPE_MODEL);

static void
nc_hireion_init (NcHIReion *nc_hireion)
{
}

static void
nc_hireion_finalize (GObject *object)
{

  /* Chain up : end */
  G_OBJECT_CLASS (nc_hireion_parent_class)->finalize (object);
}

NCM_MSET_MODEL_REGISTER_ID (nc_hireion, NC_TYPE_HIREION);

static void _nc_hireion_update_Xe (NcHIReion *reion, NcRecomb *recomb);

static void
nc_hireion_class_init (NcHIReionClass *klass)
{
  GObjectClass *object_class  = G_OBJECT_CLASS (klass);
  NcmModelClass *model_class  = NCM_MODEL_CLASS (klass);
  NcHIReionClass *reion_class = NC_HIREION_CLASS (klass);

  object_class->finalize = nc_hireion_finalize;
  
  ncm_model_class_set_name_nick (model_class, "Abstract class for H and He reionization models.", "NcHIReion");
  ncm_model_class_add_params (model_class, 0, 0, 1);

  ncm_mset_model_register_id (model_class,
                              "NcHIReion",
                              "Homogeneous and isotropic reionization models.",
                              NULL);

  ncm_model_class_check_params_info (model_class);

  reion_class->update_Xe = &_nc_hireion_update_Xe;
}

static void 
_nc_hireion_update_Xe (NcHIReion *reion, NcRecomb *recomb)
{
  NCM_UNUSED (reion);
  NCM_UNUSED (recomb);
  g_error ("_nc_hireion_update_Xe: error object `%s' do not implement this virtual function.", 
           g_type_name (G_OBJECT_TYPE (reion)));
}
