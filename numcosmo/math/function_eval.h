/***************************************************************************
 *            function_eval.h
 *
 *  Sun Jul 25 19:50:16 2010
 *  Copyright  2010  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2012 <sandro@lapsandro>
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
 * @file
 * @brief FIXME
 *
 * FIXME
 */

#ifndef _NC_FUNCTION_EVAL_H
#define _NC_FUNCTION_EVAL_H

#include <glib.h>

G_BEGIN_DECLS

typedef void (*NcFunctionLoop1) (gsize i, gpointer data, gpointer glob_data);

void ncm_function_eval_threaded (gsl_function *F, gdouble *x, gdouble *val, gulong n, guint x_stride, guint val_stride);
void ncm_function_eval_threaded_vec (gsl_function *F, gsl_vector *x, gsl_vector *val);
/* void ncm_function_eval_threaded_loop (NcFunctionLoop1 fl, gsize lp_size, gpointer glob_data); */

G_END_DECLS

#endif /* _NC_FUNCTION_EVAL_H */
