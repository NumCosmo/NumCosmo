/***************************************************************************
 *            ncm_build_inline.h
 *
 *  Mon Oct 29 23:35:29 2012
 *  Copyright  2012  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2012 <sandro@isoftware.com.br>
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
/* Inspired on build.h from gsl */


/* Compile subsequent inline functions as static functions */

#ifdef _NCM_BUILD_INLINE_H_ 
#error ncm_build_inline.h must not be included multiple times
#endif
#define _NCM_BUILD_INLINE_H_

#ifdef G_INLINE_FUNC
#undef G_INLINE_FUNC
#endif /* G_INLINE_FUNC */
#define G_INLINE_FUNC          /* disable inline in declarations */

#ifndef NUMCOSMO_HAVE_INLINE /* enable compilation of definitions in .h files */
#define NUMCOSMO_HAVE_INLINE
#endif     
