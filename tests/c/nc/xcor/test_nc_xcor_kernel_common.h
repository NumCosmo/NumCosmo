/***************************************************************************
 *            test_nc_xcor_kernel_common.h
 *
 *  Tue Oct 13 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) 2026 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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

#ifndef _TEST_NC_XCOR_KERNEL_COMMON_H_
#define _TEST_NC_XCOR_KERNEL_COMMON_H_

#include <glib.h>

G_BEGIN_DECLS

gint test_nc_xcor_kernel_main (gint argc, gchar *argv[], const gchar *kernel);

G_END_DECLS

#endif /* _TEST_NC_XCOR_KERNEL_COMMON_H_ */

