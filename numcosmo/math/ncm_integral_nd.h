/***************************************************************************
 *            ncm_integral_nd.h
 *
 *  Thu July 20 08:39:30 2023
 *  Copyright  2023 Eduardo José Barroso
 *  <eduardo.jsbarroso@uel.br>
 ****************************************************************************/
/*
 * ncm_integral_nd.h
 * Copyright (C) 2023 Eduardo José Barroso <eduardo.jsbarroso@uel.br>*
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

#ifndef _NCM_INTEGRAL_ND_H_
#define _NCM_INTEGRAL_ND_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>

#include <numcosmo/math/ncm_vector.h>

G_BEGIN_DECLS

#define NCM_TYPE_INTEGRAL_ND             (ncm_integral_nd_get_type ())
#define NCM_INTEGRAL_ND(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NCM_TYPE_INTEGRAL_ND, NcmIntegralND))
#define NCM_INTEGRAL_ND_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NCM_TYPE_INTEGRAL_ND, NcmIntegralNDClass))
#define NCM_IS_INTEGRAL_ND(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NCM_TYPE_INTEGRAL_ND))
#define NCM_IS_INTEGRAL_ND_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NCM_TYPE_INTEGRAL_ND))
#define NCM_INTEGRAL_ND_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NCM_TYPE_INTEGRAL_ND, NcmIntegralNDClass))

typedef struct _NcmIntegralNDClass NcmIntegralNDClass;
typedef struct _NcmIntegralND NcmIntegralND;
typedef struct _NcmIntegralNDPrivate NcmIntegralNDPrivate;

/**
 * NcmIntegralNDF:
 * @intnd: a #NcmIntegralND
 * @x: a #NcmVector containing the value of the variable of integration
 * @dim: the dimension of the integral argument
 * @npoints: the number of points in the array @x
 * @fdim: the dimension of the function to be integrated
 * @fval: a #NcmVector containing the @fdim values of the integrand at all points in @x
 *
 * The type of the function that must be implemented by a subclass of #NcmIntegralND.
 *
 * This function receives @npoints points in the array @x (size @dim * @npoints), and
 * returns an array (size @fdim * @npoints) of @npoints values of the integrand at all
 * points in @x. The @x is an array of doubles in row-major order
 * (i.e. the first @dim elements of @x are the coordinates of the first point, the next
 * @dim elements are the coordinates of the second point, and so on). The return value
 * is an array of @fdim values of the integrand at all points in @x (e.g. the first
 * @fdim elements are the values of the integrand at the first point, the next @fdim
 * elements are the values of the integrand at the second point, and so on).
 *
 */
typedef void (*NcmIntegralNDF) (NcmIntegralND *intnd, NcmVector *x, guint dim, guint npoints, guint fdim, NcmVector *fval);

/**
 * NcmIntegralNDGetDimensions:
 * @intnd: a #NcmIntegralND
 * @dim: (out): the dimension of the integral argument
 * @fdim: (out): the dimension of the function to be integrated
 *
 * The type of the function that must be implemented by a subclass of #NcmIntegralND.
 *
 * This function returns the dimension of the integral argument and the dimension of
 * the function to be integrated.
 *
 */
typedef void (*NcmIntegralNDGetDimensions) (NcmIntegralND *intnd, guint *dim, guint *fdim);

struct _NcmIntegralNDClass
{
  /*< private >*/
  GObjectClass parent_class;
  NcmIntegralNDF integrand;
  NcmIntegralNDGetDimensions get_dimensions;
  gpointer padding[16];
};

typedef enum _NcmIntegralNDMethod
{
  NCM_INTEGRAL_ND_METHOD_CUBATURE_H,
  NCM_INTEGRAL_ND_METHOD_CUBATURE_P,
  NCM_INTEGRAL_ND_METHOD_CUBATURE_H_V,
  NCM_INTEGRAL_ND_METHOD_CUBATURE_P_V,
  /* < private > */
  NCM_INTEGRAL_ND_METHOD_LEN, /*< skip >*/
} NcmIntegralNDMethod;

typedef enum _NcmIntegralNDError
{
  NCM_INTEGRAL_ND_ERROR_INDIVIDUAL,
  NCM_INTEGRAL_ND_ERROR_PAIRWISE,
  NCM_INTEGRAL_ND_ERROR_L2,
  NCM_INTEGRAL_ND_ERROR_L1,
  NCM_INTEGRAL_ND_ERROR_LINF,
  /* < private > */
  NCM_INTEGRAL_ND_ERROR_LEN, /*< skip >*/
} NcmIntegralNDError;

struct _NcmIntegralND
{
  /*< private >*/
  GObject parent_instance;
  NcmIntegralNDPrivate *priv;
};

GType ncm_integral_nd_get_type (void) G_GNUC_CONST;

NcmIntegralND *ncm_integral_nd_ref (NcmIntegralND *intnd);
void ncm_integral_nd_free (NcmIntegralND *intnd);
void ncm_integral_nd_clear (NcmIntegralND **intnd);

void ncm_integral_nd_set_method (NcmIntegralND *intnd, NcmIntegralNDMethod method);
void ncm_integral_nd_set_error (NcmIntegralND *intnd, NcmIntegralNDError error);
void ncm_integral_nd_set_maxeval (NcmIntegralND *intnd, guint maxeval);
void ncm_integral_nd_set_reltol (NcmIntegralND *intnd, gdouble reltol);
void ncm_integral_nd_set_abstol (NcmIntegralND *intnd, gdouble abstol);

NcmIntegralNDMethod ncm_integral_nd_get_method (NcmIntegralND *intnd);
NcmIntegralNDError ncm_integral_nd_get_error (NcmIntegralND *intnd);
guint ncm_integral_nd_get_maxeval (NcmIntegralND *intnd);
gdouble ncm_integral_nd_get_reltol (NcmIntegralND *intnd);
gdouble ncm_integral_nd_get_abstol (NcmIntegralND *intnd);

void ncm_integral_nd_eval (NcmIntegralND *intnd, const NcmVector *xi, const NcmVector *xf, NcmVector *res, NcmVector *err);

#define NCM_INTEGRAL_ND_DEFAULT_RELTOL 1e-7
#define NCM_INTEGRAL_ND_DEFAULT_ABSTOL 0.0

G_END_DECLS

#endif /* _NCM_INTEGRAL_ND_H_ */

