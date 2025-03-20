/***************************************************************************
 *            nc_hipert_private.h
 *
 *  Wed January 23 16:23:21 2019
 *  Copyright  2019  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_hipert_private.h
 *
 * Copyright (C) 2019 - Sandro Dias Pinto Vitenti
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef _NC_HIPERT_PRIVATE_H_
#define _NC_HIPERT_PRIVATE_H_

#include <glib.h>

G_BEGIN_DECLS

typedef struct _NcHIPertPrivate
{
  /*< private >*/
  gdouble alpha0;
  N_Vector y;
  SUNMatrix A;
  SUNLinearSolver LS;
  guint sys_size;
  gpointer cvode;
  gboolean cvode_init;
  gboolean cvode_stiff;
  gdouble reltol;
  gdouble abstol;
  N_Vector vec_abstol;
  gdouble k;
  gboolean prepared;
} NcHIPertPrivate;

NcHIPertPrivate *nc_hipert_get_private (NcHIPert *pert);

G_END_DECLS

#endif /* _NC_HIPERT_PRIVATE_H_ */

