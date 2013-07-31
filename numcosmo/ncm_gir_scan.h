/***************************************************************************
 *            ncm_gir_scan.h
 *
 *  Wed Oct 31 11:42:12 2012
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

/**
 * SECTION:ncm_gir_scan
 * @title: Gir Scanning Compatibility.
 * @short_description: Gir scanning types stubs
 *
 * Stubs to avoid warnings from gir scanning 
 * all functions/structs using these types
 * must be skipped using (skip).
 * 
 * These types do not represent anything, do not
 * use this documentation.
 * 
 */

/**
 * mpq_ptr: (skip)
 */
#define mpq_ptr gint

/**
 * mpq_t: (skip)
 */
#define mpq_t gint

/**
 * mp_rnd_t: (skip)
 */
#define mp_rnd_t gint

/**
 * mpz_t: (skip)
 */
#define mpz_t gint

/**
 * mpfr_ptr: (skip)
 */
#define mpfr_ptr gint

/**
 * mpfr_t: (skip)
 */
#define mpfr_t gint

/**
 * CVRhsFn: (skip)
 */
#define CVRhsFn gint

/**
 * CVDlsDenseJacFn: (skip)
 */
#define CVDlsDenseJacFn gint

/**
 * fftw_plan: (skip)
 */
#define fftw_plan gint

/**
 * fftw_complex: (skip)
 */
#define fftw_complex gint

/**
 * FILE: (skip)
 *
 */
#define FILE gint

/**
 * NCM_MUTEX_TYPE: (skip)
 *
 */
#define NCM_MUTEX_TYPE _NCM_MUTEX_TYPE
