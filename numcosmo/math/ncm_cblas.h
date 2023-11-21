/***************************************************************************
 *            ncm_cblas.h
 *
 *  Tue Nov 21 09:31:15 2023
 *  Copyright  2023  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2023 <vitenti@uel.br>
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

#ifndef _NCM_CBLAS_H_
#define _NCM_CBLAS_H_

#include <glib.h>

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */

#ifdef HAVE_MKL_CBLAS_H
#  include <mkl_cblas.h>
#  define CBLAS_H "mkl_cblas.h"
#else
#  ifdef HAVE_CBLAS_H
#    include <cblas.h>
#    define CBLAS_H "cblas.h"
#  else
#    ifdef HAVE_GSL_CBLAS_H
#      include <gsl/gsl_cblas.h>
#      define CBLAS_H "gsl/gsl_cblas.h"
#    endif
#  endif
#endif

G_BEGIN_DECLS

#if !defined (CBLAS_H)
#error "No CBLAS header found."
#endif

#ifdef BLAS_NOT_TYPEDEFED
typedef  enum CBLAS_ORDER     CBLAS_ORDER;
typedef  enum CBLAS_TRANSPOSE CBLAS_TRANSPOSE;
typedef  enum CBLAS_UPLO      CBLAS_UPLO;
typedef  enum CBLAS_DIAG      CBLAS_DIAG;
typedef  enum CBLAS_SIDE      CBLAS_SIDE;
#endif /* BLAS_NOT_TYPEDEFED */

typedef  CBLAS_INDEX       CBLAS_INDEX_t;
typedef  CBLAS_ORDER       CBLAS_ORDER_t;
typedef  CBLAS_TRANSPOSE   CBLAS_TRANSPOSE_t;
typedef  CBLAS_UPLO        CBLAS_UPLO_t;
typedef  CBLAS_DIAG        CBLAS_DIAG_t;
typedef  CBLAS_SIDE        CBLAS_SIDE_t;

G_END_DECLS

#endif /* _NCM_CBLAS_H_ */
