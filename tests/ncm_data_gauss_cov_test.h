/***************************************************************************
 *            ncm_data_gauss_cov_test.h
 *
 *  Wed June 19 14:25:31 2013
 *  Copyright  2013  Sandro Dias Pinto Vitenti
 *  <sandro@isotfware.com.br>
 ****************************************************************************/
/*
 * ncm_data_gauss_cov_test.h
 * Copyright (C) 2013 Sandro Dias Pinto Vitenti <sandro@isotfware.com.br>
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

#ifndef _NCM_DATA_GAUSS_COV_TEST_H_
#define _NCM_DATA_GAUSS_COV_TEST_H_

#include <glib-object.h>

G_BEGIN_DECLS

#define NCM_TYPE_DATA_GAUSS_COV_TEST             (ncm_data_gauss_cov_test_get_type ())
#define NCM_DATA_GAUSS_COV_TEST(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NCM_TYPE_DATA_GAUSS_COV_TEST, NcmDataGaussCovTest))
#define NCM_DATA_GAUSS_COV_TEST_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NCM_TYPE_DATA_GAUSS_COV_TEST, NcmDataGaussCovTestClass))
#define NCM_IS_DATA_GAUSS_COV_TEST(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NCM_TYPE_DATA_GAUSS_COV_TEST))
#define NCM_IS_DATA_GAUSS_COV_TEST_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NCM_TYPE_DATA_GAUSS_COV_TEST))
#define NCM_DATA_GAUSS_COV_TEST_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NCM_TYPE_DATA_GAUSS_COV_TEST, NcmDataGaussCovTestClass))

typedef struct _NcmDataGaussCovTestClass NcmDataGaussCovTestClass;
typedef struct _NcmDataGaussCovTest NcmDataGaussCovTest;


struct _NcmDataGaussCovTestClass
{
  NcmDataGaussCovClass parent_class;
};

struct _NcmDataGaussCovTest
{
  NcmDataGaussCov parent_instance;
  gdouble a, b, c, d;
};

GType ncm_data_gauss_cov_test_get_type (void) G_GNUC_CONST;

NcmData *ncm_data_gauss_cov_test_new (void);
void ncm_data_gauss_cov_test_gen_cov (NcmDataGaussCovTest *gcov_test);
void ncm_data_gauss_cov_test_mean_func (NcmDataGaussCov *gauss, NcmMSet *mset, NcmVector *vp);

G_END_DECLS

#endif /* _NCM_DATA_GAUSS_COV_TEST_H_ */

