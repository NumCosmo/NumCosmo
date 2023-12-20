/***************************************************************************
 *            ncm_mset_func1.h
 *
 *  Sun May 20 21:32:22 2018
 *  Copyright  2018  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * ncm_mset_func1.h
 * Copyright (C) 2018 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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

#ifndef _NCM_MSET_FUNC1_H_
#define _NCM_MSET_FUNC1_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_mset_func.h>

G_BEGIN_DECLS

#define NCM_TYPE_MSET_FUNC1 (ncm_mset_func1_get_type ())

G_DECLARE_DERIVABLE_TYPE (NcmMSetFunc1, ncm_mset_func1, NCM, MSET_FUNC1, NcmMSetFunc)

typedef GArray *(*NcmMSetFunc1N) (NcmMSetFunc1 *f1, NcmMSet *mset, GArray *x);

struct _NcmMSetFunc1Class
{
  /*< private >*/
  NcmMSetFuncClass parent_class;
  NcmMSetFunc1N eval1;

  /* Padding to allow 18 virtual functions without breaking ABI. */
  gpointer padding[17];
};

NcmMSetFunc1 *ncm_mset_func1_ref (NcmMSetFunc1 *f1);

void ncm_mset_func1_free (NcmMSetFunc1 *f1);
void ncm_mset_func1_clear (NcmMSetFunc1 **f1);
GArray *ncm_mset_func1_eval1 (NcmMSetFunc1 *f1, NcmMSet *mset, GArray *x);

G_END_DECLS

#endif /* _NCM_MSET_FUNC1_H_ */

