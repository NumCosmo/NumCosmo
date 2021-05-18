/***************************************************************************
 *           ncm_spline_func_test.h
 *
 *  Wed May 20 16:30:36 2020
 *  Copyright  2020 Fernando de Simoni
 *  <fsimoni@id.uff.br>
 ****************************************************************************/
/*
 * ncm_spline_func_test.h
 * Copyright (C) 2020 Fernando de Simoni <fsimoni@id.uff.brr>
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

#ifndef _NCM_SPLINE_FUNC_TEST_H_
#define _NCM_SPLINE_FUNC_TEST_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_matrix.h>
#include <numcosmo/math/ncm_spline_func.h>

G_BEGIN_DECLS

#define NCM_TYPE_SPLINE_FUNC_TEST             (ncm_spline_func_test_get_type ())
#define NCM_SPLINE_FUNC_TEST(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NCM_TYPE_SPLINE_FUNC_TEST, NcmSplineFuncTest))
#define NCM_SPLINE_FUNC_TEST_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NCM_TYPE_SPLINE_FUNC_TEST, NcmSplineFuncTestClass))
#define NCM_IS_SPLINE_FUNC_TEST(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NCM_TYPE_SPLINE_FUNC_TEST))
#define NCM_IS_SPLINE_FUNC_TEST_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NCM_TYPE_SPLINE_FUNC_TEST))
#define NCM_SPLINE_FUNC_TEST_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NCM_TYPE_SPLINE_FUNC_TEST, NcmSplineFuncTestClass))

typedef struct _NcmSplineFuncTestClass NcmSplineFuncTestClass;
typedef struct _NcmSplineFuncTest NcmSplineFuncTest;
typedef struct _NcmSplineFuncTestPrivate NcmSplineFuncTestPrivate;

struct _NcmSplineFuncTestClass
{
  /*< private > */
  GObjectClass parent_class;
};

struct _NcmSplineFuncTest
{
  /*< private >*/
  GObject parent_instance;
  NcmSplineFuncTestPrivate *priv;
};

/**
 * NcmSplineFuncTestType:
 * @NCM_SPLINE_FUNC_TEST_TYPE_POLYNOMIAL: polynomial interpolation.
 * @NCM_SPLINE_FUNC_TEST_TYPE_COSINE: cosine series.
 * @NCM_SPLINE_FUNC_TEST_TYPE_RBF: RBF interpolation.
 * @NCM_SPLINE_FUNC_TEST_TYPE_USER: user supplied function.
 *
 * Enum to choose which base function to be used in the test suite.
 */
typedef enum _NcmSplineFuncTestType
{
  NCM_SPLINE_FUNC_TEST_TYPE_POLYNOMIAL,
  NCM_SPLINE_FUNC_TEST_TYPE_COSINE,
  NCM_SPLINE_FUNC_TEST_TYPE_RBF,
  NCM_SPLINE_FUNC_TEST_TYPE_USER,
} NcmSplineFuncTestType;

typedef void (*NcmSplineFuncTestPrepare) (gpointer p);

/**
 * NcmSplineFuncTestTypePDF:
 * @NCM_SPLINE_FUNC_TEST_TYPE_PDF_FLAT: a flat PDF (see #ncm_rng_uniform_gen).
 * @NCM_SPLINE_FUNC_TEST_TYPE_PDF_NORMAL: a gaussian PDF (see #ncm_rng_gaussian_gen).
 *
 * Enum to choose which PDF type.
 */
typedef enum _NcmSplineFuncTestTypePDF
{
  NCM_SPLINE_FUNC_TEST_TYPE_PDF_FLAT,
  NCM_SPLINE_FUNC_TEST_TYPE_PDF_NORMAL,
} NcmSplineFuncTestTypePDF;

GType ncm_spline_func_test_get_type (void) G_GNUC_CONST;

NcmSplineFuncTest *ncm_spline_func_test_new (void);
NcmSplineFuncTest *ncm_spline_func_test_ref (NcmSplineFuncTest *sft);

void ncm_spline_func_test_unref (NcmSplineFuncTest *sft);
void ncm_spline_func_test_clear (NcmSplineFuncTest **sft);

void ncm_spline_func_test_set_type (NcmSplineFuncTest *sft, NcmSplineFuncTestType type);

void ncm_spline_func_test_set_ngrid (NcmSplineFuncTest *sft, const guint ngrid);
guint ncm_spline_func_test_get_ngrid (NcmSplineFuncTest *sft);

void ncm_spline_func_test_set_seed (NcmSplineFuncTest *sft, const gulong seed);
gulong ncm_spline_func_test_get_seed (NcmSplineFuncTest *sft);

void ncm_spline_func_test_set_params_info (NcmSplineFuncTest *sft, NcmMatrix *par_info);
void ncm_spline_func_test_set_params_info_all (NcmSplineFuncTest *sft, const guint npar, const gdouble p1, const gdouble p2);
NcmMatrix *ncm_spline_func_test_get_params_info (NcmSplineFuncTest *sft);
NcmVector *ncm_spline_func_test_peek_current_params (NcmSplineFuncTest *sft);


void ncm_spline_func_test_set_xi (NcmSplineFuncTest *sft, const gdouble xi);
gdouble ncm_spline_func_test_get_xi (NcmSplineFuncTest *sft);

void ncm_spline_func_test_set_xf (NcmSplineFuncTest *sft, const gdouble xf);
gdouble ncm_spline_func_test_get_xf (NcmSplineFuncTest *sft);

void ncm_spline_func_test_set_rel_error (NcmSplineFuncTest *sft, const gdouble rel_error);
gdouble ncm_spline_func_test_get_rel_error (NcmSplineFuncTest *sft);

void ncm_spline_func_test_set_scale (NcmSplineFuncTest *sft, const gdouble scale);
gdouble ncm_spline_func_test_get_scale (NcmSplineFuncTest *sft);

void ncm_spline_func_test_set_out_threshold (NcmSplineFuncTest *sft, const gdouble out_threshold);
gdouble ncm_spline_func_test_get_out_threshold (NcmSplineFuncTest *sft);

void ncm_spline_func_test_set_user_gsl_function (NcmSplineFuncTest *sft, gsl_function *F);
void ncm_spline_func_test_set_prepare_user_function (NcmSplineFuncTest *sft, NcmSplineFuncTestPrepare F_prepare);

void ncm_spline_func_test_prepare (NcmSplineFuncTest *sft, NcmSplineFuncType ftype, NcmSplineFuncTestTypePDF pdftype);

void ncm_spline_func_test_set_one_grid_stats (NcmSplineFuncTest *sft);

void ncm_spline_func_test_monte_carlo (NcmSplineFuncTest *sft, guint nsim);

void ncm_spline_func_test_log_vals_one_grid_stats (NcmSplineFuncTest *sft);

void ncm_spline_func_test_save_grid_functions_to_txt (NcmSplineFuncTest *sft, gchar *fname);

void ncm_spline_func_test_save_knots_to_txt (NcmSplineFuncTest *sft, gchar *fname);

void ncm_spline_func_test_monte_carlo_and_save_to_txt (NcmSplineFuncTest *sft, guint nsim, gchar *fname);

void ncm_spline_func_test_log_vals_mc_stats (NcmSplineFuncTest *sft);

G_END_DECLS

#endif /* _NCM_SPLINE_FUNC_TEST_H_ */

