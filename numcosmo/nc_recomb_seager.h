/***************************************************************************
 *            nc_recomb_seager.h
 *
 *  Mon November 05 18:28:36 2012
 *  Copyright  2012  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) 2012 Sandro Dias Pinto Vitenti <sandro@isoftware.com.br>
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

#ifndef _NC_RECOMB_SEAGER_H_
#define _NC_RECOMB_SEAGER_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/nc_recomb.h>

G_BEGIN_DECLS

#define NC_TYPE_RECOMB_SEAGER             (nc_recomb_seager_get_type ())
#define NC_RECOMB_SEAGER(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_RECOMB_SEAGER, NcRecombSeager))
#define NC_RECOMB_SEAGER_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_RECOMB_SEAGER, NcRecombSeagerClass))
#define NC_IS_RECOMB_SEAGER(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_RECOMB_SEAGER))
#define NC_IS_RECOMB_SEAGER_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_RECOMB_SEAGER))
#define NC_RECOMB_SEAGER_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_RECOMB_SEAGER, NcRecombSeagerClass))

typedef struct _NcRecombSeagerClass NcRecombSeagerClass;
typedef struct _NcRecombSeager NcRecombSeager;

struct _NcRecombSeagerClass
{
  /*< private >*/
  NcRecombClass parent_class;
};

/**
 * NcRecombSeagerOpt:
 * @NC_RECOM_SEAGER_OPT_FUDGE: Whether to use fudge factor in the case_B recombination fitting formulas.
 * @NC_RECOM_SEAGER_OPT_FUDGE_GAUSS_COR: Whether to use gaussian correction in the case_B recombination fitting formulas.
 * @NC_RECOM_SEAGER_OPT_ALL: All options.
 * 
 * FIXME
 * 
 */
typedef enum _NcRecombSeagerOpt
{
  NC_RECOM_SEAGER_OPT_FUDGE           = 1 << 0,
  NC_RECOM_SEAGER_OPT_FUDGE_GAUSS_COR = 1 << 1,
  NC_RECOM_SEAGER_OPT_ALL = (1 << 2) - 1, /*< private >*/
  NC_RECOM_SEAGER_OPT_LEN,                /*< skip >*/
} NcRecombSeagerOpt;

struct _NcRecombSeager
{
  /*< private >*/
  NcRecomb parent_instance;
  NcRecombSeagerOpt opts;
  gdouble H_fudge;
  gdouble AGauss1, AGauss2;
  gdouble zGauss1, zGauss2;
  gdouble wGauss1, wGauss2;
  gpointer cvode;
  gboolean init;
  CVRhsFn ion;
  CVDlsDenseJacFn ion_J;
  N_Vector y0;
  N_Vector y;
  N_Vector abstol;
  guint n;
};

GType nc_recomb_seager_get_type (void) G_GNUC_CONST;

NcRecombSeager *nc_recomb_seager_new (void);
NcRecombSeager *nc_recomb_seager_new_full (gdouble init_frac, gdouble zi, gdouble prec);
NcRecombSeager *nc_recomb_seager_ref (NcRecombSeager *recomb_seager);
void nc_recomb_seager_free (NcRecombSeager *recomb_seager);
void nc_recomb_seager_clear (NcRecombSeager **recomb_seager);

void nc_recomb_seager_set_options (NcRecombSeager *recomb_seager, NcRecombSeagerOpt opts);
NcRecombSeagerOpt nc_recomb_seager_get_options (NcRecombSeager *recomb_seager);

gdouble nc_recomb_seager_pequignot_HI_case_B (NcRecombSeager *recomb_seager, NcHICosmo *cosmo, const gdouble Tm);
gdouble nc_recomb_seager_pequignot_HI_case_B_dTm (NcRecombSeager *recomb_seager, NcHICosmo *cosmo, const gdouble Tm);
gdouble nc_recomb_seager_hummer_HeI_case_B (NcRecombSeager *recomb_seager, NcHICosmo *cosmo, const gdouble Tm);
gdouble nc_recomb_seager_hummer_HeI_case_B_dTm (NcRecombSeager *recomb_seager, NcHICosmo *cosmo, const gdouble Tm);

G_END_DECLS

#endif /* _NC_RECOMB_SEAGER_H_ */
