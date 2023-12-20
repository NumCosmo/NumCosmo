/***************************************************************************
 *            ncm_powspec_filter.h
 *
 *  Fri June 17 10:11:55 2016
 *  Copyright  2016  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * ncm_powspec_filter.h
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

#ifndef _NCM_POWSPEC_FILTER_H_
#define _NCM_POWSPEC_FILTER_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_powspec.h>
#include <numcosmo/math/ncm_fftlog.h>
#include <numcosmo/math/ncm_spline2d.h>

G_BEGIN_DECLS

#define NCM_TYPE_POWSPEC_FILTER (ncm_powspec_filter_get_type ())

G_DECLARE_FINAL_TYPE (NcmPowspecFilter, ncm_powspec_filter, NCM, POWSPEC_FILTER, GObject)

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

NcmPowspecFilter *ncm_powspec_filter_new (NcmPowspec *ps, NcmPowspecFilterType type);
NcmPowspecFilter *ncm_powspec_filter_ref (NcmPowspecFilter *psf);

void ncm_powspec_filter_free (NcmPowspecFilter *psf);
void ncm_powspec_filter_clear (NcmPowspecFilter **psf);

void ncm_powspec_filter_prepare (NcmPowspecFilter *psf, NcmModel *model);
void ncm_powspec_filter_prepare_if_needed (NcmPowspecFilter *psf, NcmModel *model);

void ncm_powspec_filter_set_type (NcmPowspecFilter *psf, NcmPowspecFilterType type);
void ncm_powspec_filter_set_lnr0 (NcmPowspecFilter *psf, gdouble lnr0);
void ncm_powspec_filter_set_best_lnr0 (NcmPowspecFilter *psf);

void ncm_powspec_filter_set_reltol (NcmPowspecFilter *psf, const gdouble reltol);
void ncm_powspec_filter_set_reltol_z (NcmPowspecFilter *psf, const gdouble reltol_z);

void ncm_powspec_filter_set_zi (NcmPowspecFilter *psf, gdouble zi);
void ncm_powspec_filter_set_zf (NcmPowspecFilter *psf, gdouble zf);

void ncm_powspec_filter_require_zi (NcmPowspecFilter *psf, gdouble zi);
void ncm_powspec_filter_require_zf (NcmPowspecFilter *psf, gdouble zf);

NcmPowspecFilterType ncm_powspec_filter_get_filter_type (NcmPowspecFilter *psf);

gdouble ncm_powspec_filter_get_reltol (NcmPowspecFilter *psf);
gdouble ncm_powspec_filter_get_reltol_z (NcmPowspecFilter *psf);

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

NcmPowspec *ncm_powspec_filter_peek_powspec (NcmPowspecFilter *psf);

#define NCM_POWSPEC_FILTER_DEFAULT_SIZE (200)

G_END_DECLS

#endif /* _NCM_POWSPEC_FILTER_H_ */

