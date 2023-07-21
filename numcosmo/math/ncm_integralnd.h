/***************************************************************************
 *            ncm_integralnd.h
 *
 *  Thu July 20 08:39:30 2023
 *  Copyright  2023 Eduardo José Barroso
 *  <eduardo.jsbarroso@uel.br>
 ****************************************************************************/
/*
 * ncm_integralnd.h
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

#ifndef _NCM_INTEGRALND_H_
#define _NCM_INTEGRALND_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>

#include <numcosmo/math/ncm_vector.h>

G_BEGIN_DECLS

#define NCM_TYPE_INTEGRALND             (ncm_integralnd_get_type ())
#define NCM_INTEGRALND(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NCM_TYPE_INTEGRALND, NcmIntegralnd))
#define NCM_INTEGRALND_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NCM_TYPE_INTEGRALND, NcmIntegralndClass))
#define NCM_IS_INTEGRALND(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NCM_TYPE_INTEGRALND))
#define NCM_IS_INTEGRALND_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NCM_TYPE_INTEGRALND))
#define NCM_INTEGRALND_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NCM_TYPE_INTEGRALND, NcmIntegralndClass))

typedef struct _NcmIntegralndClass NcmIntegralndClass;
typedef struct _NcmIntegralnd NcmIntegralnd;
typedef struct _NcmIntegralndPrivate NcmIntegralndPrivate;

/**
 * NcmIntegralndF:
 * @intnd: a #NcmIntegralnd
 * @x: a #NcmVector containing the value of the variable of integration
 * @dim: the dimension of the integral argument
 * @npoints: the number of points in the array @x
 * @fdim: the dimension of the function to be integrated
 * @fval: a #NcmVector containing the @fdim values of the integrand at all points in @x
 *
 * The type of the function that must be implemented by a subclass of #NcmIntegralnd.
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
typedef void (*NcmIntegralndF) (NcmIntegralnd *intnd, NcmVector *x, guint dim, guint npoints, guint fdim, NcmVector *fval);

/**
 * NcmIntegralndGetDimensions:
 * @intnd: a #NcmIntegralnd
 * @dim: (out): the dimension of the integral argument
 * @fdim: (out): the dimension of the function to be integrated
 *
 * The type of the function that must be implemented by a subclass of #NcmIntegralnd.
 *
 * This function returns the dimension of the integral argument and the dimension of
 * the function to be integrated.
 *
 */
typedef void (*NcmIntegralndGetDimensions) (NcmIntegralnd *intnd, guint *dim, guint *fdim);

struct _NcmIntegralndClass
{
  /*< private >*/
  GObjectClass parent_class;
  NcmIntegralndF integrand;
  NcmIntegralndGetDimensions get_dimensions;
  gpointer padding[16];
};

typedef enum _NcmIntegralndMethod
{
  NCM_INTEGRALND_METHOD_CUBATURE_H,
  NCM_INTEGRALND_METHOD_CUBATURE_P,
  NCM_INTEGRALND_METHOD_CUBATURE_H_V,
  NCM_INTEGRALND_METHOD_CUBATURE_P_V,
  /* < private > */
  NCM_INTEGRALND_METHOD_LEN, /*< skip >*/
} NcmIntegralndMethod;

typedef enum _NcmIntegralndError
{
  NCM_INTEGRALND_ERROR_INDIVIDUAL,
  NCM_INTEGRALND_ERROR_PAIRWISE,
  NCM_INTEGRALND_ERROR_L2,
  NCM_INTEGRALND_ERROR_L1,
  NCM_INTEGRALND_ERROR_LINF,
  /* < private > */
  NCM_INTEGRALND_ERROR_LEN, /*< skip >*/
} NcmIntegralndError;

struct _NcmIntegralnd
{
  /*< private >*/
  GObject parent_instance;
  NcmIntegralndPrivate *priv;
};

GType ncm_integralnd_get_type (void) G_GNUC_CONST;

NcmIntegralnd *ncm_integralnd_ref (NcmIntegralnd *intnd);
void ncm_integralnd_free (NcmIntegralnd *intnd);
void ncm_integralnd_clear (NcmIntegralnd **intnd);

void ncm_integralnd_set_method (NcmIntegralnd *intnd, NcmIntegralndMethod method);
void ncm_integralnd_set_error (NcmIntegralnd *intnd, NcmIntegralndError error);
void ncm_integralnd_set_maxeval (NcmIntegralnd *intnd, guint maxeval);
void ncm_integralnd_set_reltol (NcmIntegralnd *intnd, gdouble reltol);
void ncm_integralnd_set_abstol (NcmIntegralnd *intnd, gdouble abstol);

NcmIntegralndMethod ncm_integralnd_get_method (NcmIntegralnd *intnd);
NcmIntegralndError ncm_integralnd_get_error (NcmIntegralnd *intnd);
guint ncm_integralnd_get_maxeval (NcmIntegralnd *intnd);
gdouble ncm_integralnd_get_reltol (NcmIntegralnd *intnd);
gdouble ncm_integralnd_get_abstol (NcmIntegralnd *intnd);

void ncm_integralnd_eval (NcmIntegralnd *intnd, const NcmVector *xi, const NcmVector *xf, NcmVector *res, NcmVector *err);

#define NCM_INTEGRALND_DEFAULT_RELTOL 1e-13
#define NCM_INTEGRALND_DEFAULT_ABSTOL 0.0

G_END_DECLS

#endif /* _NCM_INTEGRALND_H_ */

