/***************************************************************************
 *            nc_multiplicity_func.h
 *
 *  Mon Jun 28 15:09:13 2010
 *  Copyright  2010  Mariana Penna Lima
 *  <pennalima@gmail.com>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Mariana Penna Lima 2012 <pennalima@gmail.com>
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

#ifndef _NC_MULTIPLICITY_FUNC_H_
#define _NC_MULTIPLICITY_FUNC_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/nc_hicosmo.h>

G_BEGIN_DECLS

#define NC_TYPE_MULTIPLICITY_FUNC             (nc_multiplicity_func_get_type ())
#define NC_MULTIPLICITY_FUNC(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_MULTIPLICITY_FUNC, NcMultiplicityFunc))
#define NC_MULTIPLICITY_FUNC_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_MULTIPLICITY_FUNC, NcMultiplicityFuncClass))
#define NC_IS_MULTIPLICITY_FUNC(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_MULTIPLICITY_FUNC))
#define NC_IS_MULTIPLICITY_FUNC_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_MULTIPLICITY_FUNC))
#define NC_MULTIPLICITY_FUNC_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_MULTIPLICITY_FUNC, NcMultiplicityFuncClass))

typedef struct _NcMultiplicityFuncClass NcMultiplicityFuncClass;
typedef struct _NcMultiplicityFunc NcMultiplicityFunc;

struct _NcMultiplicityFuncClass
{
  /*< private >*/
  GObjectClass parent_class;
  gdouble (*eval) (NcMultiplicityFunc *mulf, NcHICosmo *cosmo, gdouble sigma, gdouble z);
};

/**
 * NcMultiplicityFuncMassDef:
 * @NC_MULTIPLICITY_FUNC_MASS_DEF_MEAN: halo mass defined in terms of the mean density $\rho_\mathrm{bg} = \rho_m(z)$
 * @NC_MULTIPLICITY_FUNC_MASS_DEF_CRITICAL: halo mass defined in terms of the critical density $\rho_\mathrm{bg} = \rho_\mathrm{crit}(z)$
 * @NC_MULTIPLICITY_FUNC_MASS_DEF_VIRIAL: halo mass defined in terms of virial overdensity times the critical density $\rho_\mathrm{bg} = \rho_\mathrm{crit
 * @NC_MULTIPLICITY_FUNC_MASS_DEF_FOF: friends of friends 
 *
 * Spherical overdensity halo mass: $$M_\Delta = \frac{4\pi}{3} \Delta \rho_\mathrm{bg} r_\Delta^3,$$
 * where $\rho_\mathrm{bg}$ is the background density of the universe at redshift z, $\rho_\mathrm{bg} (z)$.
 * For @NC_HALO_DENSITY_PROFILE_MASS_DEF_VIRIAL, the parameter #NcHaloDensityProfile:log10MDelta is ignored and
 * \begin{equation}\label{def:DVir}
 * \Delta_\mathrm{Vir} = 18 \pi^2 + 82 x - 39 x^2, \quad x \equiv \Omega_m(z) - 1.
 * \end{equation}
 *
 */
typedef enum _NcMultiplicityFuncMassDef
{
  NC_MULTIPLICITY_FUNC_MASS_DEF_MEAN = 0,
  NC_MULTIPLICITY_FUNC_MASS_DEF_CRITICAL,
  NC_MULTIPLICITY_FUNC_MASS_DEF_VIRIAL,
  NC_MULTIPLICITY_FUNC_MASS_DEF_FOF,
  /* < private > */
  NC_HALO_DENSITY_PROFILE_MASS_DEF_LEN, /*< skip >*/
} NcMultiplicityFuncMassDef;

struct _NcMultiplicityFunc
{
  /*< private >*/
  GObject parent_instance;
};

GType nc_multiplicity_func_get_type (void) G_GNUC_CONST;

NcMultiplicityFunc *nc_multiplicity_func_new_from_name (gchar *multiplicity_name);
gdouble nc_multiplicity_func_eval (NcMultiplicityFunc *mulf, NcHICosmo *cosmo, gdouble sigma, gdouble z);
void nc_multiplicity_func_free (NcMultiplicityFunc *mulf);
void nc_multiplicity_func_clear (NcMultiplicityFunc **mulf);



G_END_DECLS

#endif /* _NC_MULTIPLICITY_FUNC_H_ */
