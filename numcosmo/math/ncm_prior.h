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

#define NCM_TYPE_PRIOR (ncm_prior_get_type ())

G_DECLARE_DERIVABLE_TYPE (NcmPrior, ncm_prior, NCM, PRIOR, NcmMSetFunc)
struct _NcmPriorClass
{
  /*< private >*/
  NcmMSetFuncClass parent_class;
  gboolean is_m2lnL;

  /* Padding to allow 18 virtual functions without breaking ABI. */
  gpointer padding[18];
};

NcmPrior *ncm_prior_ref (NcmPrior *prior);
void ncm_prior_free (NcmPrior *prior);
void ncm_prior_clear (NcmPrior **prior);

gboolean ncm_prior_is_m2lnL (NcmPrior *prior);

G_END_DECLS

#endif /* _NCM_PRIOR_H_ */

