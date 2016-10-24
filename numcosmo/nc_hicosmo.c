/***************************************************************************
 *            nc_hicosmo.c
 *
 *  Tue May 29 19:23:52 2007
 *  Copyright  2007  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2012 <sandro@isoftware.com.br>
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
 * SECTION:nc_hicosmo
 * @title: NcHICosmo
 * @short_description: Abstract class for implementing homogeneous and isotropic cosmological models.
 *
 * FIXME
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "nc_hicosmo.h"
#include "nc_hiprim.h"
#include "nc_hireion.h"
#include "math/ncm_serialize.h"
#include "math/ncm_cfg.h"
#include "math/ncm_mset_func_list.h"

#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>

G_DEFINE_ABSTRACT_TYPE (NcHICosmo, nc_hicosmo, NCM_TYPE_MODEL);

static void
nc_hicosmo_init (NcHICosmo *cosmo)
{
  cosmo->prim  = NULL;
  cosmo->reion = NULL;
  cosmo->T     = gsl_root_fsolver_brent;
  cosmo->s     = gsl_root_fsolver_alloc (cosmo->T);
}

static void
nc_hicosmo_dispose (GObject *object)
{
  NcHICosmo *cosmo = NC_HICOSMO (object);

  nc_hiprim_clear (&cosmo->prim);
  nc_hireion_clear (&cosmo->reion);

  /* Chain up : end */
  G_OBJECT_CLASS (nc_hicosmo_parent_class)->dispose (object);
}

static void
nc_hicosmo_finalize (GObject *object)
{
  NcHICosmo *cosmo = NC_HICOSMO (object);
  gsl_root_fsolver_free (cosmo->s);
  
  /* Chain up : end */
  G_OBJECT_CLASS (nc_hicosmo_parent_class)->finalize (object);
}

static gboolean _nc_hicosmo_valid (NcmModel *model);
NCM_MSET_MODEL_REGISTER_ID (nc_hicosmo, NC_TYPE_HICOSMO);
static void _nc_hicosmo_add_submodel (NcmModel *model, NcmModel *submodel);

static gdouble _nc_hicosmo_H0 (NcHICosmo *cosmo);
static gdouble _nc_hicosmo_Omega_b0 (NcHICosmo *cosmo);
static gdouble _nc_hicosmo_Omega_g0 (NcHICosmo *cosmo);
static gdouble _nc_hicosmo_Omega_nu0 (NcHICosmo *cosmo);
static gdouble _nc_hicosmo_Omega_r0 (NcHICosmo *cosmo);
static gdouble _nc_hicosmo_Omega_c0 (NcHICosmo *cosmo);
static gdouble _nc_hicosmo_Omega_t0 (NcHICosmo *cosmo);
static gdouble _nc_hicosmo_T_gamma0 (NcHICosmo *cosmo);
static gdouble _nc_hicosmo_Yp_4He (NcHICosmo *cosmo);
static gdouble _nc_hicosmo_z_lss (NcHICosmo *cosmo);
static gdouble _nc_hicosmo_as_drag (NcHICosmo *cosmo);
static gdouble _nc_hicosmo_xb (NcHICosmo *cosmo);

static gdouble _nc_hicosmo_E2 (NcHICosmo *cosmo, gdouble z);
static gdouble _nc_hicosmo_dE2_dz (NcHICosmo *cosmo, gdouble z);
static gdouble _nc_hicosmo_d2E2_dz2 (NcHICosmo *cosmo, gdouble z);
static gdouble _nc_hicosmo_bgp_cs2 (NcHICosmo *cosmo, gdouble z);
static gdouble _nc_hicosmo_Dc (NcHICosmo *cosmo, gdouble z);

static void
nc_hicosmo_class_init (NcHICosmoClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);
  NcmModelClass *model_class = NCM_MODEL_CLASS (klass);

  object_class->dispose  = &nc_hicosmo_dispose;
  object_class->finalize = &nc_hicosmo_finalize;

  ncm_model_class_set_name_nick (model_class, "Abstract class for HI cosmological models.", "NcHICosmo");
  ncm_model_class_add_params (model_class, 0, 0, 1);

  ncm_mset_model_register_id (model_class,
                              "NcHICosmo",
                              "Homogeneous and isotropic cosmological models.",
                              NULL,
                              FALSE,
                              NCM_MSET_MODEL_MAIN);

  ncm_model_class_check_params_info (model_class);

  model_class->valid        = &_nc_hicosmo_valid;
  model_class->add_submodel = &_nc_hicosmo_add_submodel;

  klass->H0        = &_nc_hicosmo_H0;
  klass->Omega_b0  = &_nc_hicosmo_Omega_b0;
  klass->Omega_g0  = &_nc_hicosmo_Omega_g0;
  klass->Omega_nu0 = &_nc_hicosmo_Omega_nu0;
  klass->Omega_r0  = &_nc_hicosmo_Omega_r0;
  klass->Omega_c0  = &_nc_hicosmo_Omega_c0;
  klass->Omega_t0  = &_nc_hicosmo_Omega_t0;
  klass->T_gamma0  = &_nc_hicosmo_T_gamma0;
  klass->Yp_4He    = &_nc_hicosmo_Yp_4He;
  klass->z_lss     = &_nc_hicosmo_z_lss;
  klass->as_drag   = &_nc_hicosmo_as_drag;
  klass->xb        = &_nc_hicosmo_xb;
  klass->E2        = &_nc_hicosmo_E2;
  klass->dE2_dz    = &_nc_hicosmo_dE2_dz;
  klass->d2E2_dz2  = &_nc_hicosmo_d2E2_dz2;
  klass->bgp_cs2   = &_nc_hicosmo_bgp_cs2;
  klass->Dc        = &_nc_hicosmo_Dc;
}

static gdouble _nc_hicosmo_H0 (NcHICosmo *cosmo)        { g_error ("nc_hicosmo_H0: model `%s' does not implement this function.", G_OBJECT_TYPE_NAME (cosmo)); return 0.0; }
static gdouble _nc_hicosmo_Omega_b0 (NcHICosmo *cosmo)  { g_error ("nc_hicosmo_Omega_b0: model `%s' does not implement this function.", G_OBJECT_TYPE_NAME (cosmo)); return 0.0;  }
static gdouble _nc_hicosmo_Omega_g0 (NcHICosmo *cosmo)  { g_error ("nc_hicosmo_Omega_g0: model `%s' does not implement this function.", G_OBJECT_TYPE_NAME (cosmo)); return 0.0;  }
static gdouble _nc_hicosmo_Omega_nu0 (NcHICosmo *cosmo) { g_error ("nc_hicosmo_Omega_nu0: model `%s' does not implement this function.", G_OBJECT_TYPE_NAME (cosmo)); return 0.0;  }
static gdouble _nc_hicosmo_Omega_r0 (NcHICosmo *cosmo)  { g_error ("nc_hicosmo_Omega_r0: model `%s' does not implement this function.", G_OBJECT_TYPE_NAME (cosmo)); return 0.0;  }
static gdouble _nc_hicosmo_Omega_c0 (NcHICosmo *cosmo)  { g_error ("nc_hicosmo_Omega_c0: model `%s' does not implement this function.", G_OBJECT_TYPE_NAME (cosmo)); return 0.0;  }
static gdouble _nc_hicosmo_Omega_t0 (NcHICosmo *cosmo)  { g_error ("nc_hicosmo_Omega_t0: model `%s' does not implement this function.", G_OBJECT_TYPE_NAME (cosmo)); return 0.0;  }
static gdouble _nc_hicosmo_T_gamma0 (NcHICosmo *cosmo)  { g_error ("nc_hicosmo_T_gamma0: model `%s' does not implement this function.", G_OBJECT_TYPE_NAME (cosmo)); return 0.0;  }
static gdouble _nc_hicosmo_Yp_4He (NcHICosmo *cosmo)    { g_error ("nc_hicosmo_Yp_4He: model `%s' does not implement this function.", G_OBJECT_TYPE_NAME (cosmo)); return 0.0;  }
static gdouble _nc_hicosmo_z_lss (NcHICosmo *cosmo)     { g_error ("nc_hicosmo_z_lss: model `%s' does not implement this function.", G_OBJECT_TYPE_NAME (cosmo)); return 0.0;  }
static gdouble _nc_hicosmo_as_drag (NcHICosmo *cosmo)   { g_error ("nc_hicosmo_as_drag: model `%s' does not implement this function.", G_OBJECT_TYPE_NAME (cosmo)); return 0.0;  }
static gdouble _nc_hicosmo_xb (NcHICosmo *cosmo)        { g_error ("nc_hicosmo_xb: model `%s' does not implement this function.", G_OBJECT_TYPE_NAME (cosmo)); return 0.0;  }

static gdouble _nc_hicosmo_E2 (NcHICosmo *cosmo, gdouble z)       { g_error ("nc_hicosmo_E2: model `%s' does not implement this function.", G_OBJECT_TYPE_NAME (cosmo)); return 0.0;  }
static gdouble _nc_hicosmo_dE2_dz (NcHICosmo *cosmo, gdouble z)   { g_error ("nc_hicosmo_dE2_dz: model `%s' does not implement this function.", G_OBJECT_TYPE_NAME (cosmo)); return 0.0;  }
static gdouble _nc_hicosmo_d2E2_dz2 (NcHICosmo *cosmo, gdouble z) { g_error ("nc_hicosmo_d2E2_dz2: model `%s' does not implement this function.", G_OBJECT_TYPE_NAME (cosmo)); return 0.0;  }
static gdouble _nc_hicosmo_bgp_cs2 (NcHICosmo *cosmo, gdouble z)  { g_error ("nc_hicosmo_bgp_cs2: model `%s' does not implement this function.", G_OBJECT_TYPE_NAME (cosmo)); return 0.0;  }
static gdouble _nc_hicosmo_Dc (NcHICosmo *cosmo, gdouble z)       { g_error ("nc_hicosmo_Dc: model `%s' does not implement this function.", G_OBJECT_TYPE_NAME (cosmo)); return 0.0;  }

static gboolean
_nc_hicosmo_valid (NcmModel *model)
{
  if (!NCM_MODEL_CLASS (nc_hicosmo_parent_class)->valid (model))
    return FALSE;
  /* Chain up : start */

  return (nc_hicosmo_E2 (NC_HICOSMO (model), 0.0) >= 0.0);
}

static void 
_nc_hicosmo_add_submodel (NcmModel *model, NcmModel *submodel)
{
  /* Chain up : start */
  NCM_MODEL_CLASS (nc_hicosmo_parent_class)->add_submodel (model, submodel);
  {
    NcHICosmo *cosmo = NC_HICOSMO (model);

    if (ncm_model_id (submodel) == nc_hiprim_id ())
    {
      nc_hiprim_clear (&cosmo->prim);
      cosmo->prim = nc_hiprim_ref (NC_HIPRIM (submodel));
    }
    else if (ncm_model_id (submodel) == nc_hireion_id ())
    {
      nc_hireion_clear (&cosmo->reion);
      cosmo->reion = nc_hireion_ref (NC_HIREION (submodel));
    }
  }
}

/**
 * nc_hicosmo_set_H0_impl: (skip)
 * @model_class: a #NcmModelClass
 * @f: an implementation of H0.
 *
 * Sets the implementation of H0 to @f.
 *
 */
NCM_MODEL_SET_IMPL_FUNC(NC_HICOSMO,NcHICosmo,nc_hicosmo,NcHICosmoFunc0,H0)
/**
 * nc_hicosmo_set_Omega_b0_impl: (skip)
 * @model_class: a #NcmModelClass
 * @f: FIXME
 *
 * FIXME
 *
 */
NCM_MODEL_SET_IMPL_FUNC(NC_HICOSMO,NcHICosmo,nc_hicosmo,NcHICosmoFunc0,Omega_b0)
/**
 * nc_hicosmo_set_Omega_g0_impl: (skip)
 * @model_class: a #NcmModelClass
 * @f: FIXME
 *
 * FIXME
 *
 */
NCM_MODEL_SET_IMPL_FUNC(NC_HICOSMO,NcHICosmo,nc_hicosmo,NcHICosmoFunc0,Omega_g0)
/**
 * nc_hicosmo_set_Omega_nu0_impl: (skip)
 * @model_class: a #NcmModelClass
 * @f: FIXME
 *
 * FIXME
 *
 */
NCM_MODEL_SET_IMPL_FUNC(NC_HICOSMO,NcHICosmo,nc_hicosmo,NcHICosmoFunc0,Omega_nu0)
/**
 * nc_hicosmo_set_Omega_r0_impl: (skip)
 * @model_class: a #NcmModelClass
 * @f: FIXME
 *
 * FIXME
 *
 */
NCM_MODEL_SET_IMPL_FUNC(NC_HICOSMO,NcHICosmo,nc_hicosmo,NcHICosmoFunc0,Omega_r0)
/**
 * nc_hicosmo_set_Omega_c0_impl: (skip)
 * @model_class: a #NcmModelClass
 * @f: FIXME
 *
 * FIXME
 *
 */
NCM_MODEL_SET_IMPL_FUNC(NC_HICOSMO,NcHICosmo,nc_hicosmo,NcHICosmoFunc0,Omega_c0)
/**
 * nc_hicosmo_set_Omega_t0_impl: (skip)
 * @model_class: a #NcmModelClass
 * @f: FIXME
 *
 * FIXME
 *
 */
NCM_MODEL_SET_IMPL_FUNC(NC_HICOSMO,NcHICosmo,nc_hicosmo,NcHICosmoFunc0,Omega_t0)
/**
 * nc_hicosmo_set_T_gamma0_impl: (skip)
 * @model_class: a #NcmModelClass
 * @f: FIXME
 *
 * FIXME
 *
 */
NCM_MODEL_SET_IMPL_FUNC(NC_HICOSMO,NcHICosmo,nc_hicosmo,NcHICosmoFunc0,T_gamma0)
/**
 * nc_hicosmo_set_Yp_4He_impl: (skip)
 * @model_class: a #NcmModelClass
 * @f: FIXME
 *
 * FIXME
 *
 */
NCM_MODEL_SET_IMPL_FUNC(NC_HICOSMO,NcHICosmo,nc_hicosmo,NcHICosmoFunc0,Yp_4He)
/**
 * nc_hicosmo_set_z_lss_impl: (skip)
 * @model_class: a #NcmModelClass
 * @f: FIXME
 *
 * FIXME
 *
 */
NCM_MODEL_SET_IMPL_FUNC(NC_HICOSMO,NcHICosmo,nc_hicosmo,NcHICosmoFunc0,z_lss)

/**
 * nc_hicosmo_set_as_drag_impl: (skip)
 * @model_class: a #NcmModelClass
 * @f: FIXME
 *
 * FIXME
 *
 */
NCM_MODEL_SET_IMPL_FUNC(NC_HICOSMO,NcHICosmo,nc_hicosmo,NcHICosmoFunc0,as_drag)

/**
 * nc_hicosmo_set_xb_impl: (skip)
 * @model_class: a #NcmModelClass
 * @f: FIXME
 *
 * FIXME
 *
 */
NCM_MODEL_SET_IMPL_FUNC(NC_HICOSMO,NcHICosmo,nc_hicosmo,NcHICosmoFunc0,xb)

/**
 * nc_hicosmo_set_E2_impl: (skip)
 * @model_class: a #NcmModelClass
 * @f: FIXME
 *
 * Normalized Hubble function squared, $E^2(z)$.
 *
 */
NCM_MODEL_SET_IMPL_FUNC(NC_HICOSMO,NcHICosmo,nc_hicosmo,NcHICosmoFunc1Z,E2)

/**
 * nc_hicosmo_set_dE2_dz_impl: (skip)
 * @model_class: a #NcmModelClass
 * @f: FIXME
 *
 * First derivative with respect to the redshift of the normalized Hubble function squared, $\frac{dE^2(z)}{dz}$. 
 *
 */
NCM_MODEL_SET_IMPL_FUNC(NC_HICOSMO,NcHICosmo,nc_hicosmo,NcHICosmoFunc1Z,dE2_dz)

/**
 * nc_hicosmo_set_d2E2_dz2_impl: (skip)
 * @model_class: a #NcmModelClass
 * @f: FIXME
 *
 * Second derivative with respect to the redshift of the normalized Hubble function squared, $\frac{d^2E^2(z)}{dz^2}$.
 *
 */
NCM_MODEL_SET_IMPL_FUNC(NC_HICOSMO,NcHICosmo,nc_hicosmo,NcHICosmoFunc1Z,d2E2_dz2)

/**
 * nc_hicosmo_set_bgp_cs2_impl: (skip)
 * @model_class: a #NcmModelClass
 * @f: FIXME
 *
 * FIXME
 *
 */
NCM_MODEL_SET_IMPL_FUNC(NC_HICOSMO,NcHICosmo,nc_hicosmo,NcHICosmoFunc1Z,bgp_cs2)

/**
 * nc_hicosmo_set_Dc_impl: (skip)
 * @model_class: a #NcmModelClass
 * @f: FIXME
 *
 * FIXME
 *
 */
NCM_MODEL_SET_IMPL_FUNC(NC_HICOSMO,NcHICosmo,nc_hicosmo,NcHICosmoFunc1Z,Dc)

/**
 * nc_hicosmo_new_from_name:
 * @parent_type: parent's #GType
 * @cosmo_name: Cosmological model's name
 *
 * Creates a new instance of @cosmo_name,
 * asserting that it descends from @parent_type.
 *
 * Returns: newly created @cosmo_name object.
 */
NcHICosmo *
nc_hicosmo_new_from_name (GType parent_type, gchar *cosmo_name)
{
  GObject *obj = ncm_serialize_global_from_string (cosmo_name);
  GType model_type = G_OBJECT_TYPE (obj);

  if (!g_type_is_a (model_type, parent_type))
    g_error ("nc_hicosmo_new_from_name: NcHICosmo %s do not descend from %s.", cosmo_name, g_type_name (parent_type));
  return NC_HICOSMO (obj);
}

/**
 * nc_hicosmo_ref:
 * @cosmo: a #NcHICosmo
 *
 * Increases the reference count of @cosmo by one.
 *
 * Returns: (transfer full): @cosmo.
 */
NcHICosmo *
nc_hicosmo_ref (NcHICosmo *cosmo)
{
  return g_object_ref (cosmo);
}

/**
 * nc_hicosmo_free:
 * @cosmo: a #NcHICosmo
 *
 * Decreases the reference count of @cosmo by one.
 *
 */
void
nc_hicosmo_free (NcHICosmo *cosmo)
{
  g_object_unref (cosmo);
}

/**
 * nc_hicosmo_clear:
 * @cosmo: a #NcHICosmo
 *
 * The reference count of @cosmo is decreased and the pointer is set to NULL.
 *
 */
void
nc_hicosmo_clear (NcHICosmo **cosmo)
{
  g_clear_object (cosmo);
}

static void
_nc_hicosmo_log_all_models_go (GType model_type, guint n)
{
  guint nc, i, j;
  GType *models = g_type_children (model_type, &nc);
  for (i = 0; i < nc; i++)
  {
    guint ncc;
    GType *modelsc = g_type_children (models[i], &ncc);

    g_message ("#  ");
    for (j = 0; j < n; j++) g_message (" ");
    g_message ("%s\n", g_type_name (models[i]));
    if (ncc)
      _nc_hicosmo_log_all_models_go (models[i], n + 2);

    g_free (modelsc);
  }
  g_free (models);
}

/**
 * nc_hicosmo_log_all_models:
 * @parent: #GType of the parent model
 *
 * Logs all models descending from @parent.
 *
 */
void
nc_hicosmo_log_all_models (GType parent)
{
  g_message ("# Registred NcHICosmo:%s are:\n", g_type_name (parent));
  _nc_hicosmo_log_all_models_go (parent, 0);
}

typedef struct _zt_params
{
  NcHICosmo *cosmo;
}zt_params;

static gdouble
_nc_hicosmo_zt_func (gdouble z, void *params)
{
  zt_params *p = (zt_params *) params;
  return nc_hicosmo_q (p->cosmo, z);
}

/**
 * nc_hicosmo_zt:
 * @cosmo: a #NcHICosmo
 * @z_max: maximum redshift $z_\mathrm{max}$
 *
 * Computes the deceleration-acceleration transition redshift, $z_t$ (find 
 * numerically the first root of $q(z)$ in the interval $[0,z_\mathrm{max}]$).
 * If $z_t$ is not found, i.e., $q(z) \neq 0$ in the entire redshift interval, 
 * the function returns NAN.
 * 
 * Redshift interval: $[0.0, @z_max]$.
 *
 * Returns: the transition redshift $z_t$ or NAN if not found.
 */
gdouble
nc_hicosmo_zt (NcHICosmo *cosmo, gdouble z_max)
{
  gint status;
  gint iter = 0, max_iter = 10000;
  gdouble zt;
  gdouble z_hi = 0.0;
  gdouble z_lo = 0.0;
  const gdouble step   = (1.0e-1 > z_max) ? (z_max / 10.0) : 1.0e-1;
  const gdouble reltol = 1.0e-7;
  const gdouble q0     = nc_hicosmo_q (cosmo, z_lo);
  gsl_function F;
  zt_params params;

  params.cosmo = cosmo;
  
  F.function = &_nc_hicosmo_zt_func;
  F.params   = &params;

  do
  {
    z_hi += step;
  } while (q0 * nc_hicosmo_q (cosmo, z_hi) >= 0.0 && (z_hi < z_max));

  if (z_hi > z_max)
    return GSL_NAN;
  
  gsl_root_fsolver_set (cosmo->s, &F, z_lo, z_hi);

  /*
  printf ("using %s method\n", 
          gsl_root_fsolver_name (cosmo->s));

  printf ("%5s [%9s, %9s] %9s %10s %9s\n",
          "iter", "lower", "upper", "root", 
          "err", "err(est)");

  printf ("zmin = %.5g zmax = %.5g\n", z_lo, z_max);
  */

  do
  {
    iter++;
    status = gsl_root_fsolver_iterate (cosmo->s);

    zt     = gsl_root_fsolver_root (cosmo->s);
    z_lo   = gsl_root_fsolver_x_lower (cosmo->s);
    z_hi   = gsl_root_fsolver_x_upper (cosmo->s);
    status = gsl_root_test_interval (z_lo, z_hi, 0.0, reltol);

    /*
     if (status == GSL_SUCCESS)
     printf ("Converged:\n");

     printf ("%5d [%.7f, %.7f] %.7f %.7f\n",
     iter, z_lo, z_max,
     zt, z_max - z_lo);
     */
  }
  while (status == GSL_CONTINUE && iter < max_iter);

  if (status != GSL_SUCCESS)
    return GSL_NAN;
  else
    return zt;  
}

/*
 * Inlined functions
 */
/**
 * nc_hicosmo_H0: (virtual H0)
 * @cosmo: a #NcHICosmo
 *
 * The value of the Hubble constant in unit of $\mathrm{m}\,\mathrm{s}^{-1}\,\mathrm{kpc}^{-1}$,
 * see ncm_c_kpc().
 *
 * Returns: $H_0 \left[\mathrm{m}\,\mathrm{s}^{-1}\,\mathrm{kpc}^{-1}\right]$
 */
/**
 * nc_hicosmo_RH_Mpc:
 * @cosmo: a #NcHICosmo
 *
 * Calculates the Hubble radius in unit of
 * Mpc, i.e., $R_H = (c / (H_0 \times 1\,\mathrm{Mpc}))$. 
 *
 * Returns: $R_H \left[\mathrm{Mpc}\right]$.
 */
/**
 * nc_hicosmo_RH_planck:
 * @cosmo: a #NcHICosmo
 *
 * Calculates the Hubble radius in unit of
 * Mpc, i.e., $R_H = (c / (H_0 \times l_\mathrm{planck}))$.
 * See ncm_c_planck_length().
 *
 * Returns: $R_H \left[l_\mathrm{planck}\right]$.
 */
/**
 * nc_hicosmo_h:
 * @cosmo: a #NcHICosmo
 *
 * Reduced Hubble constant, $h \equiv H_0 / (1\times\mathrm{m}\mathrm{s}^{-1}\mathrm{kpc}^{-1})$.
 *
 * Returns: $h$.
 */
/**
 * nc_hicosmo_h2:
 * @cosmo: a #NcHICosmo
 *
 * Reduced Hubble constant [nc_hicosmo_h()] squared $h^2$.
 *
 * Returns: $h^2$.
 */
/**
 * nc_hicosmo_Omega_b0: (virtual Omega_b0)
 * @cosmo: a #NcHICosmo
 *
 * Dimensionless baryon density today $\Omega_{b0} = \rho_{b0} / \rho_{\mathrm{crit}0}$,
 * see nc_hicosmo_crit_density().
 *
 * Returns: $\Omega_{b0}$
 */
/**
 * nc_hicosmo_Omega_c0: (virtual Omega_c0)
 * @cosmo: a #NcHICosmo
 *
 * Dimensionless cold dark matter density today $\Omega_{c0} = \rho_{c0} / \rho_{\mathrm{crit}0}$,
 * see nc_hicosmo_crit_density().
 *
 * Returns: $\Omega_{c0}$
 */
/**
 * nc_hicosmo_Omega_g0: (virtual Omega_g0)
 * @cosmo: a #NcHICosmo
 *
 * Dimensionless photon density today $\Omega_{\gamma0} = \rho_{\gamma0} / \rho_{\mathrm{crit}0}$,
 * see nc_hicosmo_crit_density().
 *
 * Returns: $\Omega_{\gamma0}$
 */
/**
 * nc_hicosmo_Omega_nu0: (virtual Omega_nu0)
 * @cosmo: a #NcHICosmo
 *
 * Dimensionless relativistic neutrinos density today $\Omega_{\nu0} = \rho_{\nu0} / \rho_{\mathrm{crit}0}$,
 * see nc_hicosmo_crit_density().
 *
 * Returns: $\Omega_{\nu0}$
 */
/**
 * nc_hicosmo_Omega_m0:
 * @cosmo: a #NcHICosmo
 *
 * Dimensionless total dust density today $\Omega_{m0} = \rho_{m0} / \rho_{\mathrm{crit}0}$,
 * see nc_hicosmo_crit_density().
 *
 * Returns: $\Omega_{m0}$.
 */
/**
 * nc_hicosmo_Omega_r0: (virtual Omega_r0)
 * @cosmo: a #NcHICosmo
 *
 * Dimensionless total radiation density today $\Omega_{r0} = \rho_{r0} / \rho_{\mathrm{crit}0}$,
 * see nc_hicosmo_crit_density().
 *
 * Returns: $\Omega_{r0}$
 */
/**
 * nc_hicosmo_Omega_b0h2:
 * @cosmo: a #NcHICosmo
 *
 * Dimensionless baryon density today [nc_hicosmo_Omega_b0()] times $h^2$.
 *
 * Returns: $\Omega_{b0}h^2$.
 */
/**
 * nc_hicosmo_Omega_c0h2:
 * @cosmo: a #NcHICosmo
 *
 * Dimensionless cold dark matter density today [nc_hicosmo_Omega_c0()] times $h^2$.
 *
 * Returns: $\Omega_{c0}h^2$.
 */
/**
 * nc_hicosmo_Omega_g0h2:
 * @cosmo: a #NcHICosmo
 *
 * Dimensionless photon density today [nc_hicosmo_Omega_g0()] times $h^2$.
 *
 * Returns: $\Omega_{\gamma0}h^2$.
 */
/**
 * nc_hicosmo_Omega_nu0h2:
 * @cosmo: a #NcHICosmo
 *
 * Dimensionless relativistic neutrinos density today [nc_hicosmo_Omega_nu0()] times $h^2$.
 *
 * Returns: $\Omega_{\nu0}h^2$.
 */
/**
 * nc_hicosmo_Omega_m0h2:
 * @cosmo: a #NcHICosmo
 *
 * Dimensionless total dust density today [nc_hicosmo_Omega_m0()] times $h^2$.
 *
 * Returns: $\Omega_{m0}h^2$.
 */
/**
 * nc_hicosmo_Omega_r0h2:
 * @cosmo: a #NcHICosmo
 *
 * Dimensionless total radiation density today [nc_hicosmo_Omega_r0()] times $h^2$.
 *
 * Returns: $\Omega_{r0}h^2$.
 */
/**
 * nc_hicosmo_Omega_t0: (virtual Omega_t0)
 * @cosmo: a #NcHICosmo
 *
 * Dimensionless total matter density today $\Omega_{t0} = \rho_{t0} / \rho_{\mathrm{crit}0}$,
 * see nc_hicosmo_crit_density().
 *
 * Returns: $\Omega_{t0}$
 */
/**
 * nc_hicosmo_Omega_k0:
 * @cosmo: a #NcHICosmo
 *
 * The curvature parameter today, $\Omega_{k0}$.
 *
 * Returns: $\Omega_{k0}$.
 */

/**
 * nc_hicosmo_T_gamma0: (virtual T_gamma0)
 * @cosmo: a #NcHICosmo
 *
 * Gets the cosmic microwave background radiation temperature today.
 *
 * Returns: $T_{\gamma0} \left[\mathrm{K}\right]$.
 */
/**
 * nc_hicosmo_Yp_4He: (virtual Yp_4He)
 * @cosmo: a #NcHICosmo
 *
 * Gets the primordial Helium mass fraction, i.e., 
 * $$Y_p = \frac{m_\mathrm{He}n_\mathrm{He}}
 * {m_\mathrm{He}n_\mathrm{He} + m_\mathrm{H}n_\mathrm{H}},$$ where $m_\mathrm{He}$, 
 * $n_\mathrm{He}$, $m_\mathrm{H}$ and $m_\mathrm{H}$ are respectively Helium-4 mass and 
 * number density and Hydrogen-1 mass and number density.
 *
 * Returns: $Y_p$.
 */
/**
 * nc_hicosmo_Yp_1H:
 * @cosmo: a #NcHICosmo
 *
 * The primordial hydrogen mass fraction $$Y_{\text{1H}p} = 1 - Y_p,$$
 * where $Y_p$ is the helium mass fraction, see nc_hicosmo_Yp_4He().
 *
 * Returns: $Y_{\text{1H}p}$.
 */
/**
 * nc_hicosmo_XHe:
 * @cosmo: a #NcHICosmo
 *
 * The primordial Helium to Hydrogen ratio $$X_\text{He} = 
 * \frac{n_\text{He}}{n_\text{H}} = \frac{m_\text{1H}}{m_\text{4He}}
 * \frac{Y_p}{Y_{\text{1H}p}},$$ see nc_hicosmo_Yp_1H() and nc_hicosmo_Yp_4He().
 * 
 * Returns: The primordial Helium to Hydrogen ratio $X_\text{He}$.
 */
/**
 * nc_hicosmo_crit_density:
 * @cosmo: a #NcHICosmo
 *
 * Calculares the critical density $\rho_\mathrm{crit}$ using 
 * ncm_c_crit_density_h2() $\times$ nc_hicosmo_h2().
 * 
 * Returns: The critical density $\rho_{\mathrm{crit}0}$.
 */
/**
 * nc_hicosmo_baryon_density:
 * @cosmo: a #NcHICosmo
 * 
 * Calculares the baryon density $\rho_{b0} = \rho_{\mathrm{crit}0} \Omega_{b0}$ 
 * using nc_hicosmo_crit_density() $\times$ nc_hicosmo_Omega_b0().
 * 
 * Returns: The baryon density $\rho_{b0}$.
 */
/**
 * nc_hicosmo_He_number_density:
 * @cosmo: a #NcHICosmo
 *
 * Calculares the Helium-4 number density $n_\mathrm{4He} = Y_p n_{b0} / m_\mathrm{4He}$
 * using nc_hicosmo_Yp_4He() $\times$ nc_hicosmo_baryon_density() / ncm_c_rest_energy_4He().
 * 
 * Returns: The baryon density $n_\mathrm{4He}$.
 */
/**
 * nc_hicosmo_H_number_density:
 * @cosmo: a #NcHICosmo
 *
 * Calculares the Hydrogen-1 number density $n_\mathrm{1H} = Y_{\mathrm{1H}p} n_{b0} / m_\mathrm{1H}$
 * using nc_hicosmo_Yp_1H() $\times$ nc_hicosmo_baryon_density() / ncm_c_rest_energy_1H().
 * 
 * Returns: The baryon density $n_\mathrm{1H}$.
 */

/**
 * nc_hicosmo_z_lss: (virtual z_lss)
 * @cosmo: a #NcHICosmo
 *
 * FIXME
 *
 * Returns: FIXME
 */
/**
 * nc_hicosmo_as_drag: (virtual as_drag)
 * @cosmo: a #NcHICosmo
 *
 * FIXME
 *
 * Returns: FIXME
 */
/**
 * nc_hicosmo_xb: (virtual xb)
 * @cosmo: a #NcHICosmo
 * 
 * FIXME
 * 
 * Returns: FIXME
 */

/**
 * nc_hicosmo_E2Omega_t:
 * @cosmo: a #NcHICosmo
 * @z: redshift $z$
 *
 * The value of the normalized Hubble function square times
 * the total dimensionless density $\Omega_t$.
 *
 * Returns: $E^2\Omega_t$.
 */
/**
 * nc_hicosmo_H:
 * @cosmo: a #NcHICosmo
 * @z: redshift $z$
 *
 * The value of the Hubble function in unity of $\mathrm{m}\,\mathrm{s}^{-1}\,\mathrm{kpc}^{-1}$,
 * see ncm_c_kpc().
 *
 * Returns: $H(z) \left[\mathrm{m}\,\mathrm{s}^{-1}\,\mathrm{kpc}^{-1}\right]$
 */
/**
 * nc_hicosmo_dH_dz:
 * @cosmo: a #NcHICosmo
 * @z: redshift $z$
 *
 * FIXME
 *
 * Returns: FIXME
 */
/**
 * nc_hicosmo_E:
 * @cosmo: a #NcHICosmo
 * @z: redshift $z$
 *
 * This function computes the normalized Hubble function $E(z)$.
 *
 * Returns: $E(z)$.
 */
/**
 * nc_hicosmo_E2: (virtual E2)
 * @cosmo: a #NcHICosmo
 * @z: redshift $z$
 *
 * Normalized Hubble function squared.
 *
 * Returns: $H^2 / H_0^2$.
 */
/**
 * nc_hicosmo_Em2:
 * @cosmo: a #NcHICosmo
 * @z: redshift $z$
 *
 * This function computes the inverse of the square normalized Hubble function.
 *
 * Returns: $E(z)^{-2}$.
 */
/**
 * nc_hicosmo_dE2_dz: (virtual dE2_dz)
 * @cosmo: a #NcHICosmo
 * @z: redshift $z$
 *
 * FIXME
 *
 * Returns: FIXME
 */
/**
 * nc_hicosmo_d2E2_dz2: (virtual d2E2_dz2)
 * @cosmo: a #NcHICosmo
 * @z: redshift $z$
 *
 * FIXME
 *
 * Returns: FIXME
 */

/**
 * nc_hicosmo_bgp_cs2: (virtual bgp_cs2)
 * @cosmo: a #NcHICosmo
 * @z: redshift $z$
 *
 * Baryon-photon plasma speed of sound squared, 
 * $$c_s^{b\gamma2} = (\dot{\rho}_b + \dot{\rho}_\gamma) / (p_b + p_\gamma).$$
 *
 * Returns: $c_s^{b\gamma2}$.
 */
/**
 * nc_hicosmo_Dc: (virtual Dc)
 * @cosmo: a #NcHICosmo
 * @z: redshift $z$
 *
 * FIXME
 *
 * Returns: FIXME
 */
/**
 * nc_hicosmo_powspec: (virtual powspec)
 * @cosmo: a #NcHICosmo
 * @k: wavenumber $k$
 *
 * FIXME
 *
 * Returns: FIXME
 */

/**
 * nc_hicosmo_q:
 * @cosmo: a #NcHICosmo
 * @z: redshift $z$
 *
 * FIXME
 *
 * Returns: FIXME
 */
/**
 * nc_hicosmo_nec:
 * @cosmo: a #NcHICosmo
 * @z: redshift $z$
 *
 * FIXME
 *
 * Returns: FIXME
 */
/**
 * nc_hicosmo_dec:
 * @cosmo: a #NcHICosmo
 * @z: redshift $z$
 *
 * FIXME
 *
 * Returns: FIXME
 */
/**
 * nc_hicosmo_wec:
 * @cosmo: a #NcHICosmo
 * @z: redshift $z$
 *
 * FIXME
 *
 * Returns: FIXME
 */
/**
 * nc_hicosmo_qp:
 * @cosmo: a #NcHICosmo
 * @z: redshift $z$
 *
 * FIXME
 *
 * Returns: FIXME
 */
/**
 * nc_hicosmo_j:
 * @cosmo: a #NcHICosmo
 * @z: redshift $z$
 *
 * FIXME
 *
 * Returns: FIXME
 */
/**
 * nc_hicosmo_mqE2:
 * @cosmo: a #NcHICosmo
 * @z: redshift $z$
 *
 * Calculates $-q(z)E^2(z)$.
 *
 * Returns: $-q(z)E^2(z)$.
 */

/**
 * nc_hicosmo_abs_alpha:
 * @cosmo: a #NcHICosmo
 * @x: redshift variable $x = 1 + z$
 *
 * FIXME
 *
 * Returns: FIXME
 */
/**
 * nc_hicosmo_x_alpha:
 * @cosmo: a #NcHICosmo
 * @alpha: redshift $\alpha$
 *
 * FIXME
 *
 * Returns: FIXME
 */

/**
 * nc_hicosmo_peek_prim:
 * @cosmo: a #NcHICosmo
 *
 * FIXME
 *
 * Returns: (transfer none): the #NcHIPrim submodel.
 */
/**
 * nc_hicosmo_peek_reion: 
 * @cosmo: a #NcHICosmo
 *
 * FIXME
 *
 * Returns: (transfer none): the #NcHIReion submodel.
 */

/*
 * nc_hicosmo_sigma8:
 * @cosmo: a #NcHICosmo
 * @psf: a #NcmPowspecFilter
 *
 * Computes the variance of the power spectrum at $R = 8 / h$ Mpc at redshift 0.
 *
 * Returns: $\sigma_8$
 *
 */
gdouble
nc_hicosmo_sigma8 (NcHICosmo *cosmo, NcmPowspecFilter *psf)
{
  if (psf->type != NCM_POWSPEC_FILTER_TYPE_TOPHAT)
    g_error ("nc_hicosmo_sigma8: sigma_8 is defined with a tophat filter, but psf is another type of filter.");

  return ncm_powspec_filter_eval_sigma (psf, 0.0, 8.0 / nc_hicosmo_h (cosmo));
}

#define _NC_HICOSMO_FUNC0_TO_FLIST(fname) \
static void _nc_hicosmo_flist_##fname (NcmMSetFuncList *flist, NcmMSet *mset, const gdouble *x, gdouble *res) \
{ \
  NcHICosmo *cosmo = NC_HICOSMO (ncm_mset_peek (mset, nc_hicosmo_id ())); \
  res[0] = nc_hicosmo_##fname (cosmo); \
}

_NC_HICOSMO_FUNC0_TO_FLIST (H0)
_NC_HICOSMO_FUNC0_TO_FLIST (RH_Mpc)
_NC_HICOSMO_FUNC0_TO_FLIST (RH_planck)
_NC_HICOSMO_FUNC0_TO_FLIST (h)
_NC_HICOSMO_FUNC0_TO_FLIST (h2)
_NC_HICOSMO_FUNC0_TO_FLIST (Omega_b0)
_NC_HICOSMO_FUNC0_TO_FLIST (Omega_c0)
_NC_HICOSMO_FUNC0_TO_FLIST (Omega_g0)
_NC_HICOSMO_FUNC0_TO_FLIST (Omega_nu0)
_NC_HICOSMO_FUNC0_TO_FLIST (Omega_m0)
_NC_HICOSMO_FUNC0_TO_FLIST (Omega_r0)
_NC_HICOSMO_FUNC0_TO_FLIST (Omega_t0)
_NC_HICOSMO_FUNC0_TO_FLIST (Omega_k0)
_NC_HICOSMO_FUNC0_TO_FLIST (Omega_b0h2)
_NC_HICOSMO_FUNC0_TO_FLIST (Omega_c0h2)
_NC_HICOSMO_FUNC0_TO_FLIST (Omega_g0h2)
_NC_HICOSMO_FUNC0_TO_FLIST (Omega_nu0h2)
_NC_HICOSMO_FUNC0_TO_FLIST (Omega_m0h2)
_NC_HICOSMO_FUNC0_TO_FLIST (Omega_r0h2)
_NC_HICOSMO_FUNC0_TO_FLIST (T_gamma0)
_NC_HICOSMO_FUNC0_TO_FLIST (Yp_4He)
_NC_HICOSMO_FUNC0_TO_FLIST (Yp_1H)
_NC_HICOSMO_FUNC0_TO_FLIST (XHe)
_NC_HICOSMO_FUNC0_TO_FLIST (z_lss)
_NC_HICOSMO_FUNC0_TO_FLIST (as_drag)
_NC_HICOSMO_FUNC0_TO_FLIST (xb)

static void 
_nc_hicosmo_flist_sigma8 (NcmMSetFuncList *flist, NcmMSet *mset, const gdouble *x, gdouble *res) \
{
  NcHICosmo *cosmo = NC_HICOSMO (ncm_mset_peek (mset, nc_hicosmo_id ()));

  ncm_powspec_filter_prepare_if_needed (NCM_POWSPEC_FILTER (flist->obj), NCM_MODEL (cosmo));
  
  res[0] = nc_hicosmo_sigma8 (cosmo, NCM_POWSPEC_FILTER (flist->obj));
}

#define _NC_HICOSMO_FUNC1_TO_FLIST(fname) \
static void _nc_hicosmo_flist_##fname (NcmMSetFuncList *flist, NcmMSet *mset, const gdouble *x, gdouble *res) \
{ \
 NcHICosmo *cosmo = NC_HICOSMO (ncm_mset_peek (mset, nc_hicosmo_id ())); \
 res[0] = nc_hicosmo_##fname (cosmo, x[0]); \
}

_NC_HICOSMO_FUNC1_TO_FLIST (H)
_NC_HICOSMO_FUNC1_TO_FLIST (dH_dz)
_NC_HICOSMO_FUNC1_TO_FLIST (E)
_NC_HICOSMO_FUNC1_TO_FLIST (E2)
_NC_HICOSMO_FUNC1_TO_FLIST (Em2)
_NC_HICOSMO_FUNC1_TO_FLIST (dE2_dz)
_NC_HICOSMO_FUNC1_TO_FLIST (d2E2_dz2)
_NC_HICOSMO_FUNC1_TO_FLIST (q)
_NC_HICOSMO_FUNC1_TO_FLIST (nec)
_NC_HICOSMO_FUNC1_TO_FLIST (dec)
_NC_HICOSMO_FUNC1_TO_FLIST (wec)
_NC_HICOSMO_FUNC1_TO_FLIST (qp)
_NC_HICOSMO_FUNC1_TO_FLIST (j)
_NC_HICOSMO_FUNC1_TO_FLIST (mqE2)
_NC_HICOSMO_FUNC1_TO_FLIST (zt)

void
_nc_hicosmo_register_functions (void)
{
  ncm_mset_func_list_register ("H0",          "H_0",                        "NcHICosmo", "Hubble constant",                           G_TYPE_NONE, _nc_hicosmo_flist_H0,          0, 1);
  ncm_mset_func_list_register ("RH_Mpc",      "R_{H} [\\mathrm{Mpc}]",      "NcHICosmo", "Hubble radius today in Mpc",                G_TYPE_NONE, _nc_hicosmo_flist_RH_Mpc,      0, 1);
  ncm_mset_func_list_register ("RH_planck",   "R_{H} [l_\\mathrm{planck}]", "NcHICosmo", "Hubble radius today in l_planck",           G_TYPE_NONE, _nc_hicosmo_flist_RH_planck,   0, 1);
  ncm_mset_func_list_register ("h",           "h",                          "NcHICosmo", "Dimensionless Hubble constant",             G_TYPE_NONE, _nc_hicosmo_flist_h,           0, 1);
  ncm_mset_func_list_register ("h2",          "h^2",                        "NcHICosmo", "Dimensionless Hubble constant squared",     G_TYPE_NONE, _nc_hicosmo_flist_h2,          0, 1);
  ncm_mset_func_list_register ("Omega_b0",    "\\Omega_{b0}",               "NcHICosmo", "Baryons density today",                     G_TYPE_NONE, _nc_hicosmo_flist_Omega_b0,    0, 1);
  ncm_mset_func_list_register ("Omega_c0",    "\\Omega_{c0}",               "NcHICosmo", "CDM density today",                         G_TYPE_NONE, _nc_hicosmo_flist_Omega_c0,    0, 1);
  ncm_mset_func_list_register ("Omega_g0",    "\\Omega_{\\gamma0}",         "NcHICosmo", "Photons density today",                     G_TYPE_NONE, _nc_hicosmo_flist_Omega_g0,    0, 1);
  ncm_mset_func_list_register ("Omega_nu0",   "\\Omega_{\\nu0}",            "NcHICosmo", "Ultra-relativistic neutrino density today", G_TYPE_NONE, _nc_hicosmo_flist_Omega_nu0,   0, 1);
  ncm_mset_func_list_register ("Omega_m0",    "\\Omega_{m0}",               "NcHICosmo", "Total dust matter density today",           G_TYPE_NONE, _nc_hicosmo_flist_Omega_m0,    0, 1);
  ncm_mset_func_list_register ("Omega_r0",    "\\Omega_{r0}",               "NcHICosmo", "Total radiation density today",             G_TYPE_NONE, _nc_hicosmo_flist_Omega_r0,    0, 1);
  ncm_mset_func_list_register ("Omega_t0",    "\\Omega_{t0}",               "NcHICosmo", "Total energy density today",                G_TYPE_NONE, _nc_hicosmo_flist_Omega_t0,    0, 1);
  ncm_mset_func_list_register ("Omega_k0",    "\\Omega_{k0}",               "NcHICosmo", "Curvature scale today",                     G_TYPE_NONE, _nc_hicosmo_flist_Omega_k0,    0, 1);
  ncm_mset_func_list_register ("Omega_b0h2",  "\\Omega_{b0}h^2",            "NcHICosmo", "Baryons density today times h^2",           G_TYPE_NONE, _nc_hicosmo_flist_Omega_b0h2,  0, 1);
  ncm_mset_func_list_register ("Omega_c0h2",  "\\Omega_{c0}h^2",            "NcHICosmo", "CDM density today times h^2",               G_TYPE_NONE, _nc_hicosmo_flist_Omega_c0h2,  0, 1);
  ncm_mset_func_list_register ("Omega_g0h2",  "\\Omega_{\\gamma0}h^2",      "NcHICosmo", "Photons density today times h^2",           G_TYPE_NONE, _nc_hicosmo_flist_Omega_g0h2,  0, 1);
  ncm_mset_func_list_register ("Omega_nu0h2", "\\Omega_{\\mu0}h^2",         "NcHICosmo", "UR Neutrinos density today times h^2",      G_TYPE_NONE, _nc_hicosmo_flist_Omega_nu0h2, 0, 1);
  ncm_mset_func_list_register ("Omega_m0h2",  "\\Omega_{m0}h^2",            "NcHICosmo", "Total dust matter density today times h^2", G_TYPE_NONE, _nc_hicosmo_flist_Omega_m0h2,  0, 1);
  ncm_mset_func_list_register ("Omega_r0h2",  "\\Omega_{r0}h^2",            "NcHICosmo", "Total radiation density today times h^2",   G_TYPE_NONE, _nc_hicosmo_flist_Omega_r0h2,  0, 1);
  ncm_mset_func_list_register ("T_gamma0",    "T_{\\gamma0}",               "NcHICosmo", "Photons temperature today",                 G_TYPE_NONE, _nc_hicosmo_flist_T_gamma0,    0, 1);
  ncm_mset_func_list_register ("Yp_4He",      "Y_\\mathrm{p}",              "NcHICosmo", "Primordial Helium mass fraction",           G_TYPE_NONE, _nc_hicosmo_flist_Yp_4He,      0, 1);
  ncm_mset_func_list_register ("Yp_1H",       "Y_{\\mathrm{1H}p}",          "NcHICosmo", "Primordial Hydrogen mass fraction",         G_TYPE_NONE, _nc_hicosmo_flist_Yp_1H,       0, 1);
  ncm_mset_func_list_register ("XHe",         "X_\\mathrm{HeI}",            "NcHICosmo", "Primordial Helium abundance",               G_TYPE_NONE, _nc_hicosmo_flist_XHe,         0, 1);
  ncm_mset_func_list_register ("z_lss",       "z_\\mathrm{lss}",            "NcHICosmo", "Redshift at lss",                           G_TYPE_NONE, _nc_hicosmo_flist_z_lss,       0, 1);
  ncm_mset_func_list_register ("as_drag",     "r_\\mathrm{adrag}",          "NcHICosmo", "As_drag",                                   G_TYPE_NONE, _nc_hicosmo_flist_as_drag,     0, 1);
  ncm_mset_func_list_register ("xb",          "x_b",                        "NcHICosmo", "Bounce scale",                              G_TYPE_NONE, _nc_hicosmo_flist_xb,          0, 1);

  ncm_mset_func_list_register ("sigma8",      "\\sigma_8",                  "NcHICosmo", "sigma8",                                    NCM_TYPE_POWSPEC_FILTER, _nc_hicosmo_flist_sigma8,  0, 1);
  
  ncm_mset_func_list_register ("H",           "H",                               "NcHICosmo", "Hubble function",                           G_TYPE_NONE, _nc_hicosmo_flist_H,        1, 1);
  ncm_mset_func_list_register ("dH_dz",       "\\mathrm{d}H/\\mathrm{d}z",       "NcHICosmo", "Derivative of the Hubble function",         G_TYPE_NONE, _nc_hicosmo_flist_dH_dz,    1, 1);
  ncm_mset_func_list_register ("E",           "E",                               "NcHICosmo", "Hubble function over H_0",                  G_TYPE_NONE, _nc_hicosmo_flist_E,        1, 1);
  ncm_mset_func_list_register ("E2",          "E2",                              "NcHICosmo", "Hubble function over H_0 squared",          G_TYPE_NONE, _nc_hicosmo_flist_E2,       1, 1);
  ncm_mset_func_list_register ("Em2",         "E^{-2}",                          "NcHICosmo", "One over Hubble function over H_0 squared", G_TYPE_NONE, _nc_hicosmo_flist_Em2,      1, 1);
  ncm_mset_func_list_register ("dE2_dz",      "\\mathrm{d}E^2/\\mathrm{d}z",     "NcHICosmo", "First derivative of E2",                    G_TYPE_NONE, _nc_hicosmo_flist_dE2_dz,   1, 1);
  ncm_mset_func_list_register ("d2E2_dz2",    "\\mathrm{d}^2E^2/\\mathrm{d}z^2", "NcHICosmo", "Second derivative of E2",                   G_TYPE_NONE, _nc_hicosmo_flist_d2E2_dz2, 1, 1);
  ncm_mset_func_list_register ("q",           "q",                               "NcHICosmo", "Deceleration function (also SEC >0)",       G_TYPE_NONE, _nc_hicosmo_flist_q,        1, 1);
  ncm_mset_func_list_register ("nec",         "\\mathrm{nec}",                   "NcHICosmo", "NEC violation function (>0)",               G_TYPE_NONE, _nc_hicosmo_flist_nec,      1, 1);
  ncm_mset_func_list_register ("dec",         "\\mathrm{dec}",                   "NcHICosmo", "DEC violation function (>0)",               G_TYPE_NONE, _nc_hicosmo_flist_dec,      1, 1);
  ncm_mset_func_list_register ("wec",         "\\mathrm{wec}",                   "NcHICosmo", "WEC violation function (>0)",               G_TYPE_NONE, _nc_hicosmo_flist_wec,      1, 1);
  ncm_mset_func_list_register ("qp",          "\\mathrm{d}q/\\mathrm{d}z",       "NcHICosmo", "Derivative of the deceleration function",   G_TYPE_NONE, _nc_hicosmo_flist_qp,       1, 1);
  ncm_mset_func_list_register ("j",           "j",                               "NcHICosmo", "Jerk function",                             G_TYPE_NONE, _nc_hicosmo_flist_j,        1, 1);
  ncm_mset_func_list_register ("mqE2",        "-qE^2",                           "NcHICosmo", "Effective geometric OmegaL",                G_TYPE_NONE, _nc_hicosmo_flist_mqE2,     1, 1);
  ncm_mset_func_list_register ("zt",          "z_t",                             "NcHICosmo", "Transition redshift",                       G_TYPE_NONE, _nc_hicosmo_flist_zt,       1, 1);
}
