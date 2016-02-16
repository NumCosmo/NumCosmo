/***************************************************************************
 *            nc_power_spectrum.h
 *
 *  Tue February 16 17:01:03 2016
 *  Copyright  2016  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * nc_power_spectrum.h
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

#ifndef _NC_POWER_SPECTRUM_H_
#define _NC_POWER_SPECTRUM_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/nc_hicosmo.h>
#include <numcosmo/nc_hireion.h>
#include <numcosmo/math/ncm_model_ctrl.h>

G_BEGIN_DECLS

#define NC_TYPE_POWER_SPECTRUM             (nc_power_spectrum_get_type ())
#define NC_POWER_SPECTRUM(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_POWER_SPECTRUM, NcPowerSpectrum))
#define NC_POWER_SPECTRUM_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_POWER_SPECTRUM, NcPowerSpectrumClass))
#define NC_IS_POWER_SPECTRUM(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_POWER_SPECTRUM))
#define NC_IS_POWER_SPECTRUM_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_POWER_SPECTRUM))
#define NC_POWER_SPECTRUM_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_POWER_SPECTRUM, NcPowerSpectrumClass))

typedef struct _NcPowerSpectrumClass NcPowerSpectrumClass;
typedef struct _NcPowerSpectrum NcPowerSpectrum;

struct _NcPowerSpectrumClass
{
  /*< private > */
  GObjectClass parent_class;
  void (*prepare) (NcPowerSpectrum *powspec, NcHICosmo *cosmo);
  gdouble (*eval) (NcPowerSpectrum *powspec, NcHICosmo *cosmo, const gdouble z, const gdouble k);
};

struct _NcPowerSpectrum
{
  /*< private > */
  GObject parent_instance;
  gdouble zi;
  gdouble zf;
  gdouble kmin;
  gdouble kmax;
  NcmModelCtrl *ctrl_cosmo;
};

GType nc_power_spectrum_get_type (void) G_GNUC_CONST;

NcPowerSpectrum *nc_power_spectrum_ref (NcPowerSpectrum *powspec);
void nc_power_spectrum_free (NcPowerSpectrum *powspec);
void nc_power_spectrum_clear (NcPowerSpectrum **powspec);

void nc_power_spectrum_set_zi (NcPowerSpectrum *powspec, const gdouble zi);
void nc_power_spectrum_set_zf (NcPowerSpectrum *powspec, const gdouble zf);

void nc_power_spectrum_set_kmin (NcPowerSpectrum *powspec, const gdouble kmin);
void nc_power_spectrum_set_kmax (NcPowerSpectrum *powspec, const gdouble kmax);

G_INLINE_FUNC void nc_power_spectrum_prepare (NcPowerSpectrum *powspec, NcHICosmo *cosmo);
G_INLINE_FUNC void nc_power_spectrum_prepare_if_needed (NcPowerSpectrum *powspec, NcHICosmo *cosmo);
G_INLINE_FUNC gdouble nc_power_spectrum_eval (NcPowerSpectrum *powspec, NcHICosmo *cosmo, const gdouble z, const gdouble k);

G_END_DECLS

#endif /* _NC_POWER_SPECTRUM_H_ */


#ifndef _NC_POWER_SPECTRUM_INLINE_H_
#define _NC_POWER_SPECTRUM_INLINE_H_
#ifdef NUMCOSMO_HAVE_INLINE

G_BEGIN_DECLS

G_INLINE_FUNC void 
nc_power_spectrum_prepare (NcPowerSpectrum *powspec, NcHICosmo *cosmo)
{
  NC_POWER_SPECTRUM_GET_CLASS (powspec)->prepare (powspec, cosmo);
}

G_INLINE_FUNC void
nc_power_spectrum_prepare_if_needed (NcPowerSpectrum *powspec, NcHICosmo *cosmo)
{
  gboolean cosmo_up = ncm_model_ctrl_update (powspec->ctrl_cosmo, NCM_MODEL (cosmo));

  if (cosmo_up)
    nc_power_spectrum_prepare (powspec, cosmo);
}


G_INLINE_FUNC gdouble 
nc_power_spectrum_eval (NcPowerSpectrum *powspec, NcHICosmo *cosmo, const gdouble z, const gdouble k)
{
  return NC_POWER_SPECTRUM_GET_CLASS (powspec)->eval (powspec, cosmo, z, k);
}

G_END_DECLS

#endif /* NUMCOSMO_HAVE_INLINE */
#endif /* _NC_POWER_SPECTRUM_INLINE_H_ */
