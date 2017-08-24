/***************************************************************************
 *            nc_wl_surface_mass_density.c
 *
 *  Tue Aug 15 17:22:45 2017
 *  Copyright  2017  Mariana Penna Lima
 *  <pennalima@gmail.com>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Mariana Penna Lima 2017 <pennalima@gmail.com>
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

/**
 * SECTION:nc_wl_surface_mass_density
 * @title: NcWLSurfaceMassDensity
 * @short_description: Weak lensing surface mass density
 *
 * This object implements FIXME
 * FIXME
 * $
 *  \newcommand{\RH}{{R_H}}
 *  \newcommand{\RHc}{{R^\mathrm{c}_H}}
 * $
 * 
 * The Hubble radius (or scale) is defined as the inverse of the Hubble 
 * function $H(z)$ [nc_hicosmo_H()],
 * \begin{equation}\label{eq:def:RHc}
 * \RH = \frac{c}{H(z)}, \qquad \RH_0 = \frac{c}{H_0},
 * \end{equation}
 * where $c$ is the speed of light [ncm_c_c()], $z$ is the redshift and 
 * $H_0 \equiv H(0)$ is the Hubble parameter [nc_hicosmo_H0()]. Similarly, 
 * we also define the comoving Hubble radius as 
 * \begin{equation}\label{eq:def:DH}
 * \RHc(z) = \frac{c}{aH(z)} = \frac{c(1+z)}{a_0H(z)}, \qquad \RHc_0 = \frac{c}{a_0H_0}
 * \end{equation}
 * where ${}_0$ subscript means that the function is calculated at the 
 * present time and the redshift $z$ is defined by the expression
 * $$1 + z = \frac{a_0}{a}.$$
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "lss/nc_wl_surface_mass_density.h"
#include "nc_distance.h"
#include "lss/nc_density_profile.h"
#include "math/integral.h"
#include "math/ncm_c.h"
#include "math/ncm_cfg.h"
#include "math/ncm_spline_cubic_notaknot.h"

enum
{
  PROP_0,
  PROP_MISC,
  PROP_2H,
  PROP_NW,
  PROP_CG,
  PROP_SIZE,
};

G_DEFINE_TYPE (NcWLSurfaceMassDensity, nc_wl_surface_mass_density, G_TYPE_OBJECT);

static void
nc_wl_surface_mass_density_init (NcWLSurfaceMassDensity *smd)
{
  //smd->comoving_wl_surface_mass_density_spline = NULL;

  smd->ctrl = ncm_model_ctrl_new (NULL);

  smd->misc = FALSE;
  smd->twoh = FALSE;
  smd->cg   = FALSE;
  smd->nw   = FALSE;
}

static void
_nc_wl_surface_mass_density_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcWLSurfaceMassDensity *smd = NC_WL_SURFACE_MASS_DENSITY (object);
  g_return_if_fail (NC_IS_WL_SURFACE_MASS_DENSITY (object));

  switch (prop_id)
  {
    //case PROP_ZF:
      //smd->zf = g_value_get_boolean (value);
      //ncm_model_ctrl_force_update (smd->ctrl);
      //break;
    case PROP_MISC:
      smd->misc = g_value_get_boolean (value);
      break;
    case PROP_2H:
      smd->twoh = g_value_get_boolean (value);
      break;
    case PROP_CG:
      smd->cg = g_value_get_boolean (value);
      break;
    case PROP_NW:
      smd->nw = g_value_get_boolean (value);
      break;  
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_wl_surface_mass_density_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcWLSurfaceMassDensity *smd = NC_WL_SURFACE_MASS_DENSITY (object);
  g_return_if_fail (NC_IS_WL_SURFACE_MASS_DENSITY (object));

  switch (prop_id)
  {
    //case PROP_ZF:
      //g_value_set_double (value, smd->zf);
      //break;
    case PROP_MISC:
      g_value_set_boolean (value, smd->misc);
      break;
    case PROP_2H:
      g_value_set_boolean (value, smd->twoh);
      break;
    case PROP_CG:
      g_value_set_boolean (value, smd->cg);
      break;
    case PROP_NW:
      g_value_set_boolean (value, smd->nw);
      break;  
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_wl_surface_mass_density_constructed (GObject *object)
{
  /* Chain up : start */
  G_OBJECT_CLASS (nc_wl_surface_mass_density_parent_class)->constructed (object);
  {

  }
}

static void
_nc_wl_surface_mass_density_dispose (GObject *object)
{
  NcWLSurfaceMassDensity *smd = NC_WL_SURFACE_MASS_DENSITY (object);

  //ncm_ode_spline_clear (&smd->comoving_distance_spline);

  ncm_model_ctrl_clear (&smd->ctrl);

  /* Chain up : end */
  G_OBJECT_CLASS (nc_wl_surface_mass_density_parent_class)->dispose (object);
}

static void
_nc_wl_surface_mass_density_finalize (GObject *object)
{

  /* Chain up : end */
  G_OBJECT_CLASS (nc_wl_surface_mass_density_parent_class)->finalize (object);
}

static void
nc_wl_surface_mass_density_class_init (NcWLSurfaceMassDensityClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);

  object_class->constructed  = &_nc_wl_surface_mass_density_constructed;
  object_class->set_property = &_nc_wl_surface_mass_density_set_property;
  object_class->get_property = &_nc_wl_surface_mass_density_get_property;
  object_class->dispose      = &_nc_wl_surface_mass_density_dispose;
  object_class->finalize     = &_nc_wl_surface_mass_density_finalize;

  g_object_class_install_property (object_class,
                                   PROP_MISC,
                                   g_param_spec_boolean ("misc",
                                                         NULL,
                                                         "Miscentering term",
                                                         FALSE,
                                                         G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_2H,
                                   g_param_spec_boolean ("twoh",
                                                         NULL,
                                                         "Two-halo term",
                                                         FALSE,
                                                         G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_CG,
                                   g_param_spec_boolean ("cg",
                                                         NULL,
                                                         "Baryonic mass of the central galaxy term",
                                                         FALSE,
                                                         G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_NW,
                                   g_param_spec_boolean ("nw",
                                                         NULL,
                                                         "Non-weak shear term",
                                                         FALSE,
                                                         G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  
}

/**
 * nc_wl_surface_mass_density_new:
 * @zf: final redshift $z_f$
 *
 * Creates a new #NcWLSurfaceMassDensity object optimized to perform distance calculations
 * to redshift up to $z_f$.
 *
 * Returns: a new #NcWLSurfaceMassDensity
 */
NcWLSurfaceMassDensity *
nc_wl_surface_mass_density_new (gdouble zf)
{
  return g_object_new (NC_TYPE_DISTANCE, "zf", zf, NULL);
}

/**
 * nc_wl_surface_mass_density_ref:
 * @smd: a #NcWLSurfaceMassDensity
 *
 * FIXME
 *
 * Returns: (transfer full): FIXME
 */
NcWLSurfaceMassDensity *
nc_wl_surface_mass_density_ref (NcWLSurfaceMassDensity *smd)
{
  return g_object_ref (smd);
}

/**
 * nc_wl_surface_mass_density_free:
 * @smd: a #NcWLSurfaceMassDensity
 *
 * FIXME
 *
 */
void
nc_wl_surface_mass_density_free (NcWLSurfaceMassDensity *smd)
{
  g_object_unref (smd);
}

/**
 * nc_wl_surface_mass_density_clear:
 * @smd: a #NcWLSurfaceMassDensity
 *
 * FIXME
 *
 */
void
nc_wl_surface_mass_density_clear (NcWLSurfaceMassDensity **smd)
{
  g_clear_object (smd);
}


/**
 * nc_wl_surface_mass_density_prepare:
 * @smd: a #NcWLSurfaceMassDensity
 * @cosmo: a #NcHICosmo
 *
 * FIXME
 *
 */
void
nc_wl_surface_mass_density_prepare (NcWLSurfaceMassDensity *smd, NcHICosmo *cosmo)
{
  /*IMPLEMENT FIXME */
  ncm_model_ctrl_update (smd->ctrl, NCM_MODEL (cosmo));

  return;
}

