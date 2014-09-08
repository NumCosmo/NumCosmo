/***************************************************************************
 *            ncm_mc_sampler_gauss.h
 *
 *  Wed September 03 14:55:28 2014
 *  Copyright  2014  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * ncm_mc_sampler_gauss.h
 * Copyright (C) 2014 Sandro Dias Pinto Vitenti <sandro@isoftware.com.br>
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

#ifndef _NCM_MC_SAMPLER_GAUSS_H_
#define _NCM_MC_SAMPLER_GAUSS_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_mc_sampler.h>

G_BEGIN_DECLS

#define NCM_TYPE_MC_SAMPLER_GAUSS             (ncm_mc_sampler_gauss_get_type ())
#define NCM_MC_SAMPLER_GAUSS(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NCM_TYPE_MC_SAMPLER_GAUSS, NcmMCSamplerGauss))
#define NCM_MC_SAMPLER_GAUSS_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NCM_TYPE_MC_SAMPLER_GAUSS, NcmMCSamplerGaussClass))
#define NCM_IS_MC_SAMPLER_GAUSS(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NCM_TYPE_MC_SAMPLER_GAUSS))
#define NCM_IS_MC_SAMPLER_GAUSS_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NCM_TYPE_MC_SAMPLER_GAUSS))
#define NCM_MC_SAMPLER_GAUSS_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NCM_TYPE_MC_SAMPLER_GAUSS, NcmMCSamplerGaussClass))

typedef struct _NcmMCSamplerGaussClass NcmMCSamplerGaussClass;
typedef struct _NcmMCSamplerGauss NcmMCSamplerGauss;

struct _NcmMCSamplerGaussClass
{
  /*< private >*/
  NcmMCSamplerClass parent_class;
};

struct _NcmMCSamplerGauss
{
  /*< private >*/
  NcmMCSampler parent_instance;
  guint len;
  NcmMatrix *cov;
  NcmMatrix *LLT;
  gboolean init;
};

GType ncm_mc_sampler_gauss_get_type (void) G_GNUC_CONST;

NcmMCSamplerGauss *ncm_mc_sampler_gauss_new (guint len);

void ncm_mc_sampler_gauss_set_size (NcmMCSamplerGauss *mcsg, guint len);
guint ncm_mc_sampler_gauss_get_size (NcmMCSamplerGauss *mcsg);

void ncm_mc_sampler_gauss_set_cov (NcmMCSamplerGauss *mcsg, const NcmMatrix *cov);
void ncm_mc_sampler_gauss_set_cov_variant (NcmMCSamplerGauss *mcsg, GVariant *cov);
NcmMatrix *ncm_mc_sampler_gauss_get_cov (NcmMCSamplerGauss *mcsg);

void ncm_mc_sampler_gauss_set_cov_from_scale (NcmMCSamplerGauss *mcsg);

G_END_DECLS

#endif /* _NCM_MC_SAMPLER_GAUSS_H_ */
