/***************************************************************************
 *            nc_hiprim_two_fluids.h
 *
 *  Tue Apr 23 11:14:19 2024
 *  Copyright  2024  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_hiprim_two_fluids.h
 * Copyright (C) 2024 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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

#ifndef _NC_HIPRIM_TWO_FLUIDS_H_
#define _NC_HIPRIM_TWO_FLUIDS_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_c.h>
#include <numcosmo/nc_hiprim.h>

G_BEGIN_DECLS

#define NC_TYPE_HIPRIM_TWO_FLUIDS             (nc_hiprim_two_fluids_get_type ())
#define NC_HIPRIM_TWO_FLUIDS(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_HIPRIM_TWO_FLUIDS, NcHIPrimTwoFluids))
#define NC_HIPRIM_TWO_FLUIDS_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_HIPRIM_TWO_FLUIDS, NcHIPrimTwoFluidsClass))
#define NC_IS_HIPRIM_TWO_FLUIDS(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_HIPRIM_TWO_FLUIDS))
#define NC_IS_HIPRIM_TWO_FLUIDS_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_HIPRIM_TWO_FLUIDS))
#define NC_HIPRIM_TWO_FLUIDS_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_HIPRIM_TWO_FLUIDS, NcHIPrimTwoFluidsClass))

typedef struct _NcHIPrimTwoFluidsClass NcHIPrimTwoFluidsClass;
typedef struct _NcHIPrimTwoFluids NcHIPrimTwoFluids;

/**
 * NcHIPrimTwoFluidsSParams:
 * @NC_HIPRIM_TWO_FLUIDS_LN10E10ASA: Amplitude of the adiabatic scalar mode $\ln(10^{10}\mathcal{A}_\mathrm{s})$
 * @NC_HIPRIM_TWO_FLUIDS_T_SA_RATIO: Tensor-to-scalar ratio $r$
 * @NC_HIPRIM_TWO_FLUIDS_LNK0: Logarithm of the mode $k_0$ in $\mathrm{Mpc}^{-1}$.
 * @NC_HIPRIM_TWO_FLUIDS_LNW: Logarithm of the equation of state parameter $w$.
 * @NC_HIPRIM_TWO_FLUIDS_N_T: Spectral index of the tensor power spectrum.
 *
 * Parameters for the two fluids primordial power spectrum.
 *
 */
typedef enum /*< enum,underscore_name=NC_HIPRIM_TWO_FLUIDS_SPARAMS >*/
{
  NC_HIPRIM_TWO_FLUIDS_LN10E10ASA,
  NC_HIPRIM_TWO_FLUIDS_T_SA_RATIO,
  NC_HIPRIM_TWO_FLUIDS_LNK0,
  NC_HIPRIM_TWO_FLUIDS_LNW,
  NC_HIPRIM_TWO_FLUIDS_N_T,
  /* < private > */
  NC_HIPRIM_TWO_FLUIDS_SPARAM_LEN, /*< skip >*/
} NcHIPrimTwoFluidsSParams;

struct _NcHIPrimTwoFluidsClass
{
  /*< private >*/
  NcHIPrimClass parent_class;
};

struct _NcHIPrimTwoFluids
{
  /*< private >*/
  NcHIPrim parent_instance;
};

GType nc_hiprim_two_fluids_get_type (void) G_GNUC_CONST;

#define NC_HIPRIM_TWO_FLUIDS_DEFAULT_LN10E10ASA (3.179)
#define NC_HIPRIM_TWO_FLUIDS_DEFAULT_T_SA_RATIO (0.2)
#define NC_HIPRIM_TWO_FLUIDS_DEFAULT_LNK0 (-5.0 * M_LN10)
#define NC_HIPRIM_TWO_FLUIDS_DEFAULT_LNW (-4.0 * M_LN10)
#define NC_HIPRIM_TWO_FLUIDS_DEFAULT_N_T (0.0)

NcHIPrimTwoFluids *nc_hiprim_two_fluids_new (void);

void nc_hiprim_two_fluids_set_use_default_calib (NcHIPrimTwoFluids *two_fluids, gboolean use_default_calib);
gboolean nc_hiprim_two_fluids_get_use_default_calib (NcHIPrimTwoFluids *two_fluids);

void nc_hiprim_two_fluids_set_lnk_lnw_spline (NcHIPrimTwoFluids *two_fluids, NcmSpline2d *lnSA_powspec_lnk_lnw);
NcmSpline2d *nc_hiprim_two_fluids_peek_lnk_lnw_spline (NcHIPrimTwoFluids *two_fluids);

G_END_DECLS

#endif /* _NC_HIPRIM_TWO_FLUIDS_H_ */

