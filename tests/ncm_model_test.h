/***************************************************************************
 *            ncm_model_test.h
 *
 *  Fri May 18 16:00:21 2012
 *  Copyright  2012  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2012 <vitenti@uel.br>
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

#ifndef _NCM_MODEL_TEST_H_
#define _NCM_MODEL_TEST_H_

#include <glib-object.h>

G_BEGIN_DECLS

#define NCM_TYPE_MODEL_TEST (ncm_model_test_get_type ())
#define NCM_TYPE_MODEL_TEST_CHILD (ncm_model_test_child_get_type ())
#define NCM_TYPE_MODEL_TEST_CHILD_CHILD (ncm_model_test_child_child_get_type ())

G_DECLARE_DERIVABLE_TYPE (NcmModelTest, ncm_model_test, NCM, MODEL_TEST, NcmModel)
G_DECLARE_DERIVABLE_TYPE (NcmModelTestChild, ncm_model_test_child, NCM, MODEL_TEST_CHILD, NcmModelTest)
G_DECLARE_DERIVABLE_TYPE (NcmModelTestChildChild, ncm_model_test_child_child, NCM, MODEL_TEST_CHILD_CHILD, NcmModelTestChild)

struct _NcmModelTestClass
{
  NcmModelClass parent_class;
};

struct _NcmModelTestChildClass
{
  NcmModelTestClass parent_class;
};

struct _NcmModelTestChildChildClass
{
  NcmModelTestChildClass parent_class;
};

#define SPARAM_LEN1 4
#define SPARAM_LEN2 3
#define SPARAM_LEN3 2

#define VPARAM_LEN1 3
#define VPARAM_LEN2 2
#define VPARAM_LEN3 0

#define SPARAM_LEN_TOT (SPARAM_LEN1 + SPARAM_LEN2 + SPARAM_LEN3)
#define VPARAM_LEN_TOT (VPARAM_LEN1 + VPARAM_LEN2 + VPARAM_LEN3)

extern gchar  *name_tot[3];
extern gchar  *nick_tot[3];
extern gchar  *s_symbol_tot[SPARAM_LEN_TOT],        *v_symbol_tot[VPARAM_LEN_TOT];
extern gchar  *s_name_tot[SPARAM_LEN_TOT],          *v_name_tot[VPARAM_LEN_TOT];
extern gdouble s_lb_tot[SPARAM_LEN_TOT],             v_lb_tot[VPARAM_LEN_TOT];
extern gdouble s_ub_tot[SPARAM_LEN_TOT],             v_ub_tot[VPARAM_LEN_TOT];
extern gdouble s_scale_tot[SPARAM_LEN_TOT],          v_scale_tot[VPARAM_LEN_TOT];
extern gdouble s_abstol_tot[SPARAM_LEN_TOT],         v_abstol_tot[VPARAM_LEN_TOT];
extern gdouble s_defval_tot[SPARAM_LEN_TOT],         v_defval_tot[VPARAM_LEN_TOT];
extern NcmParamType s_ftype_tot[SPARAM_LEN_TOT], v_ftype_tot[VPARAM_LEN_TOT];
extern guint v_len_tot[VPARAM_LEN_TOT];

G_END_DECLS

#endif /* _NCM_MODEL_TEST_H_ */

