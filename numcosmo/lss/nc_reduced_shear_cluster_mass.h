/***************************************************************************
 *            nc_reduced_shear_cluster_mass.h
 *
 *  Mon Mar 19 15:42:23 2018
 *  Copyright  2018  Mariana Penna Lima
 *  <pennalima@gmail.com>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Mariana Penna Lima 2018 <pennalima@gmail.com> 
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

#ifndef _NC_REDUCED_SHEAR_CLUSTER_MASS_H_
#define _NC_REDUCED_SHEAR_CLUSTER_MASS_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_model.h>
#include <numcosmo/nc_hicosmo.h>
#include <numcosmo/lss/nc_density_profile.h>
#include <numcosmo/lss/nc_wl_surface_mass_density.h>


#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_rng.h>
#include <gsl/gsl_multifit_nlin.h>
#endif /* NUMCOSMO_GIR_SCAN */

G_BEGIN_DECLS

#define NC_TYPE_REDUCED_SHEAR_CLUSTER_MASS            (nc_reduced_shear_cluster_mass_get_type ())
#define NC_REDUCED_SHEAR_CLUSTER_MASS(obj)            (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_REDUCED_SHEAR_CLUSTER_MASS, NcReducedShearClusterMass))
#define NC_REDUCED_SHEAR_CLUSTER_MASS_CLASS(klass)    (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_REDUCED_SHEAR_CLUSTER_MASS, NcReducedShearClusterMassClass))
#define NC_IS_REDUCED_SHEAR_CLUSTER_MASS(obj)         (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_REDUCED_SHEAR_CLUSTER_MASS))
#define NC_IS_REDUCED_SHEAR_CLUSTER_MASS_CLASS(klass) (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_REDUCED_SHEAR_CLUSTER_MASS))
#define NC_REDUCED_SHEAR_CLUSTER_MASS_GET_CLASS(obj)  (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_REDUCED_SHEAR_CLUSTER_MASS, NcReducedShearClusterMassClass))

typedef struct _NcReducedShearClusterMassClass NcReducedShearClusterMassClass;
typedef struct _NcReducedShearClusterMass NcReducedShearClusterMass;

/**
 * NcReducedShearClusterMassParams:
 * @NC_REDUCED_SHEAR_CLUSTER_MASS_A: shear calibration parameter
 * @NC_REDUCED_SHEAR_CLUSTER_MASS_B: shear calibration parameter
 * @NC_REDUCED_SHEAR_CLUSTER_MASS_C: shear calibration parameter
 * @NC_REDUCED_SHEAR_CLUSTER_MASS_XP: pivot parameter
 * @NC_REDUCED_SHEAR_CLUSTER_MASS_VSIGMA: Voigt profile parameter, $\sigma$ is the standard deviation of the Gaussian distribution
 * @NC_REDUCED_SHEAR_CLUSTER_MASS_VGAMMA: Voigt profile parameter, $\Gamma$ is the width of the  Lorentzian profile
 * 
 * These parameters refers to the shear calibration, and 
 * the Voigt profile. 
 * 
 * In particular, we consider the shear calibration as defined for STEP and STEP2 simulations:
 * $$\hat{g}_{Teo} = (1 + m)g + c,$$
 * where $g$ is the reduced shear, the multiplicative bias $m$ is  
 * $$ m =\left\{
 *   \begin{array}{c l}	
 *        a \frac{r_{gal}}{r_{PSF}} + b & \frac{r_{gal}}{r_{PSF}} < x_p \\
 *        b & \frac{r_{gal}}{r_{PSF}} \geq x_p
 *   \end{array}\right. $$
 * and $c = constant.$
 * $r_{gal}$ and $r_{PSF}$ are the galaxy and PSF (point spread function) sizes, respectively.
 * See [Applegate (2014)][XApplegate2014].  
 * 
 */
typedef enum _NcReducedShearClusterMassParams
{
  NC_REDUCED_SHEAR_CLUSTER_MASS_A = 0, 
  NC_REDUCED_SHEAR_CLUSTER_MASS_B, 
  NC_REDUCED_SHEAR_CLUSTER_MASS_C,
  NC_REDUCED_SHEAR_CLUSTER_MASS_XP, 
  NC_REDUCED_SHEAR_CLUSTER_MASS_VSIGMA,
  NC_REDUCED_SHEAR_CLUSTER_MASS_VGAMMA,
  /* < private > */
  NC_REDUCED_SHEAR_CLUSTER_MASS_SPARAM_LEN, /*< skip >*/
} NcReducedShearClusterMassParams;

#define NC_REDUCED_SHEAR_CLUSTER_MASS_DEFAULT_A  (0.0)
#define NC_REDUCED_SHEAR_CLUSTER_MASS_DEFAULT_B (0.0)
#define NC_REDUCED_SHEAR_CLUSTER_MASS_DEFAULT_C (0.0)
#define NC_REDUCED_SHEAR_CLUSTER_MASS_DEFAULT_XP (0.2)
#define NC_REDUCED_SHEAR_CLUSTER_MASS_DEFAULT_VSIGMA (0.3)
#define NC_REDUCED_SHEAR_CLUSTER_MASS_DEFAULT_VGAMMA (0.05)

#define NC_REDUCED_SHEAR_CLUSTER_MASS_DEFAULT_PARAMS_ABSTOL (0.0)

struct _NcReducedShearClusterMassClass
{
  /*< private >*/
  NcmModelClass parent_class;
};

struct _NcReducedShearClusterMass
{
  /*< private >*/
  NcmModel parent_instance;
  gdouble R_Mpc;
  gdouble nzbins;
  const gsl_multifit_fdfsolver_type *T;
  gsl_multifit_fdfsolver *s;
  gdouble *workz;
};

GType nc_reduced_shear_cluster_mass_get_type (void) G_GNUC_CONST;

NCM_MSET_MODEL_DECLARE_ID (nc_reduced_shear_cluster_mass);

NcReducedShearClusterMass *nc_reduced_shear_cluster_mass_new (void);
NcReducedShearClusterMass *nc_reduced_shear_cluster_mass_ref (NcReducedShearClusterMass *rscm);
void nc_reduced_shear_cluster_mass_free (NcReducedShearClusterMass *rscm);
void nc_reduced_shear_cluster_mass_clear (NcReducedShearClusterMass **rscm);

gdouble nc_reduced_shear_cluster_mass_P_z_gth_gobs (NcReducedShearClusterMass *rscm, NcHICosmo *cosmo, const gdouble z, const gdouble g_th, const gdouble g_obs);
gdouble nc_reduced_shear_cluster_mass_posterior_no_shear_calibration (NcReducedShearClusterMass *rscm, NcHICosmo *cosmo, const gdouble z, const gdouble g_obs);

G_END_DECLS

#endif /* _NC_REDUCED_SHEAR_CLUSTER_MASS_H_ */
