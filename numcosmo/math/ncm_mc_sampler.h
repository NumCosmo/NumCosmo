/***************************************************************************
 *            ncm_mc_sampler.h
 *
 *  Fri August 29 18:56:39 2014
 *  Copyright  2014  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * ncm_mc_sampler.h
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

#ifndef _NCM_MC_SAMPLER_H_
#define _NCM_MC_SAMPLER_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_mset.h>
#include <numcosmo/math/ncm_rng.h>

G_BEGIN_DECLS

#define NCM_TYPE_MC_SAMPLER             (ncm_mc_sampler_get_type ())
#define NCM_MC_SAMPLER(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NCM_TYPE_MC_SAMPLER, NcmMCSampler))
#define NCM_MC_SAMPLER_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NCM_TYPE_MC_SAMPLER, NcmMCSamplerClass))
#define NCM_IS_MC_SAMPLER(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NCM_TYPE_MC_SAMPLER))
#define NCM_IS_MC_SAMPLER_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NCM_TYPE_MC_SAMPLER))
#define NCM_MC_SAMPLER_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NCM_TYPE_MC_SAMPLER, NcmMCSamplerClass))

typedef struct _NcmMCSamplerClass NcmMCSamplerClass;
typedef struct _NcmMCSampler NcmMCSampler;

struct _NcmMCSampler
{
  /*< private >*/
  GObject parent_instance;
  NcmMSet *mset;
};

struct _NcmMCSamplerClass
{
  /*< private >*/
  GObjectClass parent_class;
  gboolean bernoulli_scheme;
  void (*set_mset) (NcmMCSampler *mcs, NcmMSet *mset);
  void (*generate) (NcmMCSampler *mcs, NcmVector *theta, NcmVector *thetastar, NcmRNG *rng);
  const gchar *(*get_name) (NcmMCSampler *mcs);
};

GType ncm_mc_sampler_get_type (void) G_GNUC_CONST;

NcmMCSampler *ncm_mc_sampler_ref (NcmMCSampler *mcs);
void ncm_mc_sampler_free (NcmMCSampler *mcs);
void ncm_mc_sampler_clear (NcmMCSampler **mcs);

void ncm_mc_sampler_set_mset (NcmMCSampler *mcs, NcmMSet *mset);
void ncm_mc_sampler_generate (NcmMCSampler *mcs, NcmVector *theta, NcmVector *thetastar, NcmRNG *rng);

const gchar *ncm_mc_sampler_get_name (NcmMCSampler *mcs);

G_END_DECLS

#endif /* _NCM_MC_SAMPLER_H_ */
