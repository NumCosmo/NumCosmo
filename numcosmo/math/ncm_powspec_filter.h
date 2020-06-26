/***************************************************************************
 *            ncm_powspec_filter.h
 *
 *  Fri June 17 10:11:55 2016
 *  Copyright  2016  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * ncm_powspec_filter.h
 * Copyright (C) 2016 Sandro Dias Pinto Vitenti <sandro@isoftware.com.br>
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

#ifndef _NCM_POWSPEC_FILTER_H_
#define _NCM_POWSPEC_FILTER_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_powspec.h>
#include <numcosmo/math/ncm_fftlog.h>
#include <numcosmo/math/ncm_spline2d.h>

G_BEGIN_DECLS

#define NCM_TYPE_POWSPEC_FILTER             (ncm_powspec_filter_get_type ())
#define NCM_POWSPEC_FILTER(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NCM_TYPE_POWSPEC_FILTER, NcmPowspecFilter))
#define NCM_POWSPEC_FILTER_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NCM_TYPE_POWSPEC_FILTER, NcmPowspecFilterClass))
#define NCM_IS_POWSPEC_FILTER(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NCM_TYPE_POWSPEC_FILTER))
#define NCM_IS_POWSPEC_FILTER_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NCM_TYPE_POWSPEC_FILTER))
#define NCM_POWSPEC_FILTER_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NCM_TYPE_POWSPEC_FILTER, NcmPowspecFilterClass))

typedef struct _NcmPowspecFilterClass NcmPowspecFilterClass;
typedef struct _NcmPowspecFilter NcmPowspecFilter;

struct _NcmPowspecFilterClass
{
  /*< private >*/
  GObjectClass parent_class;
};

/**
 * NcmPowspecFilterType:
 * @NCM_POWSPEC_FILTER_TYPE_TOPHAT: Apply the top-hat filter.
 * @NCM_POWSPEC_FILTER_TYPE_GAUSS: Apply the Gaussian filter.
 *
 * Filter type to apply to the power spectrum.
 * See also #NcmFftlogTophatwin2 and #NcmFftlogGausswin2, for the top-hat and Gaussian filters, respectively.
 *
 */
typedef enum _NcmPowspecFilterType
{
  NCM_POWSPEC_FILTER_TYPE_TOPHAT = 0,
  NCM_POWSPEC_FILTER_TYPE_GAUSS,
  /* < private > */
  NCM_POWSPEC_FILTER_TYPE_LEN, /*< skip >*/
} NcmPowspecFilterType;

struct _NcmPowspecFilter
{
  /*< private >*/
  GObject parent_instance;
  NcmPowspec *ps;
  NcmPowspecFilterType type;
  NcmFftlog *fftlog;
  gdouble lnr0;
  gdouble lnk0;
  gdouble Lk;
  gdouble zi;
  gdouble zf;
  gboolean calibrated;
  gdouble reltol;
  gdouble reltol_z;
  NcmSpline2d *var;
  NcmSpline2d *dvar;
  NcmModelCtrl *ctrl;
  gboolean constructed;
};

GType ncm_powspec_filter_get_type (void) G_GNUC_CONST;

NcmPowspecFilter *ncm_powspec_filter_new (NcmPowspec *ps, NcmPowspecFilterType type);
NcmPowspecFilter *ncm_powspec_filter_ref (NcmPowspecFilter *psf);

void ncm_powspec_filter_free (NcmPowspecFilter *psf);
void ncm_powspec_filter_clear (NcmPowspecFilter **psf);

void ncm_powspec_filter_prepare (NcmPowspecFilter *psf, NcmModel *model);
void ncm_powspec_filter_prepare_if_needed (NcmPowspecFilter *psf, NcmModel *model);

void ncm_powspec_filter_set_type (NcmPowspecFilter *psf, NcmPowspecFilterType type);
void ncm_powspec_filter_set_lnr0 (NcmPowspecFilter *psf, gdouble lnr0);
void ncm_powspec_filter_set_best_lnr0 (NcmPowspecFilter *psf);

void ncm_powspec_filter_set_zi (NcmPowspecFilter *psf, gdouble zi);
void ncm_powspec_filter_set_zf (NcmPowspecFilter *psf, gdouble zf);

void ncm_powspec_filter_require_zi (NcmPowspecFilter *psf, gdouble zi);
void ncm_powspec_filter_require_zf (NcmPowspecFilter *psf, gdouble zf);

gdouble ncm_powspec_filter_get_r_min (NcmPowspecFilter *psf);
gdouble ncm_powspec_filter_get_r_max (NcmPowspecFilter *psf);

gdouble ncm_powspec_filter_eval_lnvar_lnr (NcmPowspecFilter *psf, const gdouble z, const gdouble lnr);
gdouble ncm_powspec_filter_eval_var (NcmPowspecFilter *psf, const gdouble z, const gdouble r);
gdouble ncm_powspec_filter_eval_var_lnr (NcmPowspecFilter *psf, const gdouble z, const gdouble lnr);
gdouble ncm_powspec_filter_eval_sigma_lnr (NcmPowspecFilter *psf, const gdouble z, const gdouble lnr);
gdouble ncm_powspec_filter_eval_sigma (NcmPowspecFilter *psf, const gdouble z, const gdouble r);

gdouble ncm_powspec_filter_eval_dvar_dlnr (NcmPowspecFilter *psf, const gdouble z, const gdouble lnr);
gdouble ncm_powspec_filter_eval_dlnvar_dlnr (NcmPowspecFilter *psf, const gdouble z, const gdouble lnr);
gdouble ncm_powspec_filter_eval_dlnvar_dr (NcmPowspecFilter *psf, const gdouble z, const gdouble lnr);

gdouble ncm_powspec_filter_eval_dnvar_dlnrn (NcmPowspecFilter *psf, const gdouble z, const gdouble lnr, guint n);
gdouble ncm_powspec_filter_eval_dnlnvar_dlnrn (NcmPowspecFilter *psf, const gdouble z, const gdouble lnr, guint n);

gdouble ncm_powspec_filter_volume_rm3 (NcmPowspecFilter *psf);

#define NCM_POWSPEC_FILTER_DEFAULT_SIZE (200)

G_END_DECLS

#endif /* _NCM_POWSPEC_FILTER_H_ */

