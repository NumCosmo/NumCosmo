/***************************************************************************
 *            ncm_powspec_corr3d.h
 *
 *  Wed March 20 15:22:26 2019
 *  Copyright  2019  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/

/*
 * ncm_powspec_corr3d.h
 * Copyright (C) 2019 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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

#ifndef _NCM_POWSPEC_CORR3D_H_
#define _NCM_POWSPEC_CORR3D_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_powspec.h>
#include <numcosmo/math/ncm_fftlog.h>
#include <numcosmo/math/ncm_spline2d.h>

G_BEGIN_DECLS

#define NCM_TYPE_POWSPEC_CORR3D (ncm_powspec_corr3d_get_type ())

G_DECLARE_FINAL_TYPE (NcmPowspecCorr3d, ncm_powspec_corr3d, NCM, POWSPEC_CORR3D, NcmPowspec)

NcmPowspecCorr3d *ncm_powspec_corr3d_new (NcmPowspec * ps);
NcmPowspecCorr3d *ncm_powspec_corr3d_ref (NcmPowspecCorr3d *psc);

void ncm_powspec_corr3d_free (NcmPowspecCorr3d *psc);
void ncm_powspec_corr3d_clear (NcmPowspecCorr3d **psc);

void ncm_powspec_corr3d_set_reltol (NcmPowspecCorr3d *psc, const gdouble reltol);
void ncm_powspec_corr3d_set_reltol_z (NcmPowspecCorr3d *psc, const gdouble reltol_z);

gdouble ncm_powspec_corr3d_get_reltol (NcmPowspecCorr3d *psc);
gdouble ncm_powspec_corr3d_get_reltol_z (NcmPowspecCorr3d *psc);

void ncm_powspec_corr3d_prepare (NcmPowspecCorr3d *psc, NcmModel *model);
void ncm_powspec_corr3d_prepare_if_needed (NcmPowspecCorr3d *psc, NcmModel *model);

void ncm_powspec_corr3d_set_lnr0 (NcmPowspecCorr3d *psc, gdouble lnr0);
void ncm_powspec_corr3d_set_best_lnr0 (NcmPowspecCorr3d *psc);

void ncm_powspec_corr3d_set_zi (NcmPowspecCorr3d *psc, gdouble zi);
void ncm_powspec_corr3d_set_zf (NcmPowspecCorr3d *psc, gdouble zf);

gdouble ncm_powspec_corr3d_get_r_min (NcmPowspecCorr3d *psc);
gdouble ncm_powspec_corr3d_get_r_max (NcmPowspecCorr3d *psc);

gdouble ncm_powspec_corr3d_eval_xi_lnr (NcmPowspecCorr3d *psc, const gdouble z, const gdouble lnr);
gdouble ncm_powspec_corr3d_eval_xi (NcmPowspecCorr3d *psc, const gdouble z, const gdouble r);

#define NCM_POWSPEC_CORR3D_DEFAULT_SIZE (200)

G_END_DECLS

#endif /* _NCM_POWSPEC_CORR3D_H_ */

