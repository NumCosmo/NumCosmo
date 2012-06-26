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

typedef void (*NcmLoopFunc) (glong i, glong f, gpointer data);

void ncm_function_eval_set_max_threads (gint mt);
void ncm_function_eval_threaded_loop (NcmLoopFunc lfunc, glong i, glong f, gpointer data);

G_END_DECLS

#endif /* _NC_FUNCTION_EVAL_H */
