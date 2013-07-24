/***************************************************************************
 *            nc_data_hubble_bao.h
 *
 *  Thu Apr 22 14:35:37 2010
 *  Copyright  2010  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) 2012 Sandro Dias Pinto Vitenti <sandro@isoftware.com.br>
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

#ifndef _NC_DATA_HUBBLE_BAO_H_
#define _NC_DATA_HUBBLE_BAO_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_data_gauss_diag.h>
#include <numcosmo/nc_distance.h>

G_BEGIN_DECLS

#define NC_TYPE_DATA_HUBBLE_BAO             (nc_data_hubble_bao_get_type ())
#define NC_DATA_HUBBLE_BAO(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_DATA_HUBBLE_BAO, NcDataHubbleBao))
#define NC_DATA_HUBBLE_BAO_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_DATA_HUBBLE_BAO, NcDataHubbleBaoClass))
#define NC_IS_DATA_HUBBLE_BAO(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_DATA_HUBBLE_BAO))
#define NC_IS_DATA_HUBBLE_BAO_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_DATA_HUBBLE_BAO))
#define NC_DATA_HUBBLE_BAO_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_DATA_HUBBLE_BAO, NcDataHubbleBaoClass))

typedef struct _NcDataHubbleBaoClass NcDataHubbleBaoClass;
typedef struct _NcDataHubbleBao NcDataHubbleBao;

/**
 * NcDataHubbleBaoId:
 * @NC_DATA_HUBBLE_BAO_BUSCA2013: FIXME
 *
 * FIXME
 */
typedef enum _NcDataHubbleBaoId
{
  NC_DATA_HUBBLE_BAO_BUSCA2013, /*< private >*/
  NC_DATA_HUBBLE_BAO_NSAMPLES,  /*< skip >*/
} NcDataHubbleBaoId;

struct _NcDataHubbleBao
{
  /*< private >*/
  NcmDataGaussDiag parent_instance;
  NcmVector *x;
  NcDataHubbleBaoId id;
  NcDistance *dist;
};

struct _NcDataHubbleBaoClass
{
  /*< private >*/
  NcmDataGaussDiagClass parent_class;
};

GType nc_data_hubble_bao_get_type (void) G_GNUC_CONST;

NcmData *nc_data_hubble_bao_new (NcDistance *dist, NcDataHubbleBaoId id);

void nc_data_hubble_bao_set_size (NcDataHubbleBao *hubble_bao, guint np);
guint nc_data_hubble_bao_get_size (NcDataHubbleBao *hubble_bao);

void nc_data_hubble_bao_set_sample (NcDataHubbleBao *hubble_bao, NcDataHubbleBaoId id);
NcDataHubbleBaoId nc_data_hubble_bao_get_sample (NcDataHubbleBao *hubble_bao);

G_END_DECLS

#endif /* _NC_DATA_HUBBLE_BAO_H_ */

