/***************************************************************************
 *            print_data.h
 *
 *  Sat Sep 10 18:56:31 2011
 *  Copyright  2011  Mariana Penna Lima
 *  <pennalima@gmail.com>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Mariana Penna Lima 2012 <pennalima@gmail.com>
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

#ifndef _NC_PRINT_DATA_H
#define _NC_PRINT_DATA_H

G_BEGIN_DECLS

void nc_mass_function_print (NcData *ca_unbinned, NcHICosmo *model, FILE *out, gchar *header);

G_END_DECLS

#endif /* _NC_PRINT_DATA_H */
