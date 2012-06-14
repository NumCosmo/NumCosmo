/***************************************************************************
 *            nc_nlopt.h
 *
 *  Sat Apr  3 16:07:17 2010
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

#ifndef _NC_NLOPT_H
#define _NC_NLOPT_H

#include <glib.h>
#include <nlopt.h>

G_BEGIN_DECLS

#ifdef NUMCOSMO_HAVE_NLOPT
gboolean ncm_fit_run_nlopt (NcmFit *fit, nlopt_algorithm algo, gint niters, NcmFitRunMsgs mtype);
#endif /* NUMCOSMO_HAVE_NLOPT */

G_END_DECLS

#endif /* _NC_NLOPT_H */
