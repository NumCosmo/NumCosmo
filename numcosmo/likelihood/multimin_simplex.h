/***************************************************************************
 *            multimin_simplex.h
 *
 *  Sat Aug 16 19:57:28 2008
 *  Copyright  2008  Sandro Dias Pinto Vitenti
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

#ifndef _NC_MULTIMIN_SIMPLEX_H
#define _NC_MULTIMIN_SIMPLEX_H

#include <glib.h>

G_BEGIN_DECLS

gboolean ncm_fit_run_sp (NcmFit *fit, gint niters, NcmFitRunMsgs mtype);

G_END_DECLS

#endif /* _NC_MULTIMIN_SIMPLEX_H */
