/***************************************************************************
 *            ncm_calc.h
 *
 *  Mon March 21 15:32:16 2016
 *  Copyright  2016  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * ncm_calc.h
 * Copyright (C) 2016 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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

#ifndef _NCM_CALC_H_
#define _NCM_CALC_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_model.h>

G_BEGIN_DECLS

#define NCM_TYPE_CALC             (ncm_calc_get_type ())
#define NCM_CALC(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NCM_TYPE_CALC, NcmCalc))
#define NCM_CALC_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NCM_TYPE_CALC, NcmCalcClass))
#define NCM_IS_CALC(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NCM_TYPE_CALC))
#define NCM_IS_CALC_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NCM_TYPE_CALC))
#define NCM_CALC_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NCM_TYPE_CALC, NcmCalcClass))

typedef struct _NcmCalcClass NcmCalcClass;
typedef struct _NcmCalc NcmCalc;

typedef void (*NcmCalcPrepare0) (NcmCalc *calc);
typedef void (*NcmCalcPrepare1) (NcmCalc *calc, NcmModel *m1);
typedef void (*NcmCalcPrepare2) (NcmCalc *calc, NcmModel *m1, NcmModel *m2);
typedef void (*NcmCalcPrepare3) (NcmCalc *calc, NcmModel *m1, NcmModel *m2, NcmModel *m3);
typedef void (*NcmCalcPrepare4) (NcmCalc *calc, NcmModel *m1, NcmModel *m2, NcmModel *m3, NcmModel *m4);
typedef void (*NcmCalcPrepare5) (NcmCalc *calc, NcmModel *m1, NcmModel *m2, NcmModel *m3, NcmModel *m4, NcmModel *m5);
typedef void (*NcmCalcPrepare6) (NcmCalc *calc, NcmModel *m1, NcmModel *m2, NcmModel *m3, NcmModel *m4, NcmModel *m5, NcmModel *m6);

#define NCM_CALC_PREPARE(prepare_ptr) ((NcmCalcPrepare0 *) (prepare_ptr))

struct _NcmCalcClass
{
  /*< private > */
  GObjectClass parent_class;
  guint ndep;
  GArray *dep_list;
  NcmCalcPrepare0 prepare;
};

struct _NcmCalc
{
  /*< private > */
  GObject parent_instance;
  gdouble reltol;
  gdouble abstol;
  GPtrArray *ctrl;
};

GType ncm_calc_get_type (void) G_GNUC_CONST;

void ncm_calc_class_set_num_dep (NcmCalcClass *calc_class, guint ndep);
void ncm_calc_class_set_dep (NcmCalcClass *calc_class, guint p, GType dep_model);
void ncm_calc_class_check (NcmCalcClass *calc_class);

void ncm_calc_prepare_array (NcmCalc *calc, NcmModel **ma);
void ncm_calc_prepare_if_needed_array (NcmCalc *calc, NcmModel **ma);
void ncm_calc_prepare_if_needed_vargs (NcmCalc *calc, ...) G_GNUC_NULL_TERMINATED;

void ncm_calc_set_reltol (NcmCalc *calc, const gdouble reltol);
void ncm_calc_set_abstol (NcmCalc *calc, const gdouble abstol);

gdouble ncm_calc_get_reltol (NcmCalc *calc);
gdouble ncm_calc_get_abstol (NcmCalc *calc);

#define NCM_CALC_DEFAULT_RELTOL (1.0e-7)
#define NCM_CALC_DEFAULT_ABSTOL (0.0)
#define NCM_CALC_MAX_DEPS (6)

G_END_DECLS

#endif /* _NCM_CALC_H_ */

