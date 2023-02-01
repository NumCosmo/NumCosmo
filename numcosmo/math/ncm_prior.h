/***************************************************************************
 *            ncm_prior.h
 *
 *  Wed August 03 10:08:38 2016
 *  Copyright  2016  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * ncm_prior.h
 * Copyright (C) 2016 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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

#ifndef _NCM_PRIOR_H_
#define _NCM_PRIOR_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_mset_func.h>

G_BEGIN_DECLS

#define NCM_TYPE_PRIOR             (ncm_prior_get_type ())
#define NCM_PRIOR(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NCM_TYPE_PRIOR, NcmPrior))
#define NCM_PRIOR_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NCM_TYPE_PRIOR, NcmPriorClass))
#define NCM_IS_PRIOR(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NCM_TYPE_PRIOR))
#define NCM_IS_PRIOR_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NCM_TYPE_PRIOR))
#define NCM_PRIOR_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NCM_TYPE_PRIOR, NcmPriorClass))

typedef struct _NcmPriorClass NcmPriorClass;
typedef struct _NcmPrior NcmPrior;

struct _NcmPriorClass
{
  /*< private >*/
  NcmMSetFuncClass parent_class;
  gboolean is_m2lnL;
};

struct _NcmPrior
{
  /*< private >*/  
  NcmMSetFunc parent_instance;
};

GType ncm_prior_get_type (void) G_GNUC_CONST;

NcmPrior *ncm_prior_ref (NcmPrior *prior);
void ncm_prior_free (NcmPrior *prior);
void ncm_prior_clear (NcmPrior **prior);

gboolean ncm_prior_is_m2lnL (NcmPrior *prior);

G_END_DECLS

#endif /* _NCM_PRIOR_H_ */
