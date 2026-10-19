/***************************************************************************
 *            nc_xcor_kernel_weak_lensing.c
 *
 *  Tue July 14 12:00:00 2015
 *  Copyright  2015  Cyrille Doux
 *  <cdoux@apc.in2p3.fr>
 *  Sat December 27 20:21:01 2025
 *  Copyright  2025  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) 2015 Cyrille Doux <cdoux@apc.in2p3.fr>
 * Copyright (C) 2025 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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

/**
 * NcXcorKernelWeakLensing:
 *
 * Implementation of #NcXcorKernel for galaxy weak lensing.
 *
 * The kernel is given by:
 * \begin{equation}
 *    W^{\kappa_{\mathrm{gal}}} = \frac{3}{2} \frac{\Omega_m H_0^2}{c} \frac{(1+z)}{H(z)} \chi(z) \int dz^\prime \frac{dn}{dz^\prime} \frac{\chi(z^\prime) - \chi(z)}{\chi(z^\prime)}
 * \end{equation}
 * where $\frac{dn}{dz}$ is the redshift distribution of galaxies.
 *
 * With #NcXcorKernelWeakLensing:doia set, the kernel adds the nonlinear alignment (NLA)
 * intrinsic shear, in the convention of CCL's WeakLensingTracer with use_A_ia,
 * \begin{equation}
 *    W^{\mathrm{IA}} = -A_\mathrm{IA} \left(\frac{1+z}{1+z_\mathrm{piv}}\right)^{\eta_\mathrm{IA}}
 *    \frac{C_1 \rho_\mathrm{crit} \Omega_m}{D(k, z)} \frac{H(z)}{c} \frac{dn}{dz},
 * \end{equation}
 * with $C_1 \rho_\mathrm{crit} = 5 \times 10^{-14}\, h^{-2} M_\odot^{-1} \mathrm{Mpc}^3
 * \times \rho_\mathrm{crit}$, weighted by $j_\ell(k\chi) / (k\chi)^2$ and the spin-2
 * factor $\sqrt{(\ell+2)!/(\ell-2)!}$, as the lensing term. The growth factor is read
 * from a power spectrum, $D(k, z) = \sqrt{P_D(k, z) / P_D(k, 0)}$, as the redshift-space
 * distortion term of #NcXcorKernelGal reads the growth rate. $P_D$ is
 * #NcXcorKernelWeakLensing:growth-powspec when set and the kernel's power spectrum
 * otherwise. The NLA model divides by the linear growth factor: with a nonlinear kernel
 * power spectrum, set growth-powspec to the linear one, since the ratio of nonlinear
 * spectra includes the nonlinear evolution in $D$.
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "ncm/core/ncm_cfg.h"
#include "nc/xcor/nc_xcor.h"
#include "nc/xcor/nc_xcor_kernel_weak_lensing.h"
#include "nc/xcor/nc_xcor_kernel_component.h"
#include "nc/xcor/nc_xcor_lensing_efficiency.h"

#include "ncm/integration/ncm_integrate.h"
#include "ncm/spline/ncm_spline_gsl.h"
#include "ncm/spline/ncm_spline_func.h"
#include "ncm/core/ncm_util.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_math.h>
#include <gsl/gsl_randist.h>
#endif /* NUMCOSMO_GIR_SCAN */

struct _NcXcorKernelWeakLensing
{
  /*< private >*/
  NcXcorKernel parent_instance;

  /* Lensing efficiency */
  NcXcorLensingEfficiency *lens_eff;

  /* Redshift */
  NcmSpline *dn_dz_user;
  NcmSpline *dn_dz;
  gdouble dn_dz_zmin;
  gdouble dn_dz_zmax;
  gdouble dn_dz_min;
  gdouble dn_dz_max;

  gdouble nbar;
  gdouble intr_shear;

  gdouble noise;
  gboolean doia;
  NcmPowspec *growth_ps;
  NcXcorKernelComponent *wl_comp;
  NcXcorKernelComponent *ia_comp;
};

enum
{
  PROP_0,
  PROP_DN_DZ,
  PROP_NBAR,
  PROP_INTR_SHEAR,
  PROP_DOIA,
  PROP_GROWTH_POWSPEC,
  PROP_SIZE,
};

G_DEFINE_TYPE (NcXcorKernelWeakLensing, nc_xcor_kernel_weak_lensing, NC_TYPE_XCOR_KERNEL)

static gdouble _nc_xcor_kernel_weak_lensing_lens_eff_eval_source (NcXcorLensingEfficiency *lens_eff, gdouble z);
static void _nc_xcor_kernel_weak_lensing_lens_eff_get_z_range (NcXcorLensingEfficiency *lens_eff, gdouble *zmin, gdouble *zmax);
static gdouble _nc_xcor_kernel_weak_lensing_dndz (NcXcorKernelWeakLensing *xclkg, gdouble z);
static gdouble _nc_xcor_kernel_weak_lensing_ia_amplitude (NcXcorKernelWeakLensing *xclkg, gdouble z);
static NcmPowspec *_nc_xcor_kernel_weak_lensing_peek_growth_ps (NcXcorKernelWeakLensing *xclkg);

NC_XCOR_LENSING_EFFICIENCY_DEFINE_TYPE (NC, XCOR_KERNEL_WEAK_LENSING_LENS_EFF,
                                        NcXcorKernelWeakLensingLensEff,
                                        nc_xcor_kernel_weak_lensing_lens_eff,
                                        _nc_xcor_kernel_weak_lensing_lens_eff_eval_source,
                                        _nc_xcor_kernel_weak_lensing_lens_eff_get_z_range,
                                        NcXcorKernelWeakLensing *)

/*
 * Weak Lensing Component Definition
 * Handles the weak lensing convergence kernel
 */

typedef struct _WeakLensingComponentData
{
  NcDistance *dist;
  NcmPowspec *ps;
  NcXcorLensingEfficiency *lens_eff;
  gdouble dn_dz_zmax;
} WeakLensingComponentData;

#define _NC_XCOR_KERNEL_COMPONENT_WEAK_LENSING_GET_DATA(comp) \
        ((WeakLensingComponentData *) ((guint8 *) (comp) + sizeof (NcXcorKernelComponent)))

static gdouble _wl_component_eval_kernel (NcXcorKernelComponent *comp, NcHICosmo *cosmo, gdouble chi, gdouble k);
static gdouble _wl_component_eval_window (NcXcorKernelComponent *comp, NcHICosmo *cosmo, const NcXcorKinetic *xck);
static gdouble _wl_component_eval_kfactor (NcXcorKernelComponent *comp, NcHICosmo *cosmo, const NcXcorKinetic *xck, gdouble k);
static void _wl_component_eval_kfactor_vec (NcXcorKernelComponent *comp, NcHICosmo *cosmo, const NcXcorKinetic *xck, NcmVector *k, NcmVector *out);
static gdouble _wl_component_eval_prefactor (NcXcorKernelComponent *comp, NcHICosmo *cosmo, gdouble k, gint l);
static void _wl_component_get_limits (NcXcorKernelComponent *comp, NcHICosmo *cosmo, gdouble *chi_min, gdouble *chi_max, gdouble *k_min, gdouble *k_max);
static void _wl_component_data_clear (WeakLensingComponentData *data);
static NcXcorKernelComponent *_nc_xcor_kernel_component_weak_lensing_new (NcDistance *dist, NcmPowspec *ps, NcXcorLensingEfficiency *lens_eff, gdouble dn_dz_zmax);

NC_XCOR_KERNEL_COMPONENT_DEFINE_TYPE (NC, XCOR_KERNEL_COMPONENT_WEAK_LENSING,
                                      NcXcorKernelComponentWeakLensing,
                                      nc_xcor_kernel_component_weak_lensing,
                                      _wl_component_eval_kernel,
                                      _wl_component_eval_window,
                                      _wl_component_eval_kfactor,
                                      _wl_component_eval_kfactor_vec,
                                      _wl_component_eval_prefactor,
                                      _wl_component_get_limits,
                                      WeakLensingComponentData,
                                      _wl_component_data_clear)

/*
 * Intrinsic Alignment Component Definition
 * Handles the NLA intrinsic shear: -A_IA(z) C1 rho_crit Omega_m / D(k, z) * dn_dz(z) * E(z)
 */

typedef struct _IAComponentData
{
  NcXcorKernelWeakLensing *xclkg;
  NcDistance *dist;
  NcmPowspec *ps;
  NcmPowspec *growth_ps;
} IAComponentData;

#define _NC_XCOR_KERNEL_COMPONENT_IA_GET_DATA(comp) \
        ((IAComponentData *) ((guint8 *) (comp) + sizeof (NcXcorKernelComponent)))

static gdouble _ia_component_eval_kernel (NcXcorKernelComponent *comp, NcHICosmo *cosmo, gdouble chi, gdouble k);
static gdouble _ia_component_eval_window (NcXcorKernelComponent *comp, NcHICosmo *cosmo, const NcXcorKinetic *xck);
static gdouble _ia_component_eval_kfactor (NcXcorKernelComponent *comp, NcHICosmo *cosmo, const NcXcorKinetic *xck, gdouble k);
static void _ia_component_eval_kfactor_vec (NcXcorKernelComponent *comp, NcHICosmo *cosmo, const NcXcorKinetic *xck, NcmVector *k, NcmVector *out);
static gdouble _ia_component_eval_prefactor (NcXcorKernelComponent *comp, NcHICosmo *cosmo, gdouble k, gint l);
static void _ia_component_get_limits (NcXcorKernelComponent *comp, NcHICosmo *cosmo, gdouble *chi_min, gdouble *chi_max, gdouble *k_min, gdouble *k_max);
static void _ia_component_data_clear (IAComponentData *data);
static NcXcorKernelComponent *_nc_xcor_kernel_component_ia_new (NcXcorKernelWeakLensing *xclkg, NcDistance *dist, NcmPowspec *ps, NcmPowspec *growth_ps);

NC_XCOR_KERNEL_COMPONENT_DEFINE_TYPE (NC, XCOR_KERNEL_COMPONENT_IA,
                                      NcXcorKernelComponentIA,
                                      nc_xcor_kernel_component_ia,
                                      _ia_component_eval_kernel,
                                      _ia_component_eval_window,
                                      _ia_component_eval_kfactor,
                                      _ia_component_eval_kfactor_vec,
                                      _ia_component_eval_prefactor,
                                      _ia_component_get_limits,
                                      IAComponentData,
                                      _ia_component_data_clear)

/* Kernel Parameters */

#define VECTOR   (NCM_MODEL (xclkg))
#define IA_A     (ncm_model_orig_param_get (VECTOR, NC_XCOR_KERNEL_WEAK_LENSING_IA_A))
#define IA_ETA   (ncm_model_orig_param_get (VECTOR, NC_XCOR_KERNEL_WEAK_LENSING_IA_ETA))
#define IA_Z_PIV (ncm_model_orig_param_get (VECTOR, NC_XCOR_KERNEL_WEAK_LENSING_IA_Z_PIV))

static void
nc_xcor_kernel_weak_lensing_init (NcXcorKernelWeakLensing *xclkg)
{
  xclkg->lens_eff   = NULL;
  xclkg->dn_dz_user = NULL;
  xclkg->dn_dz      = NULL;
  xclkg->dn_dz_zmin = 0.0;
  xclkg->dn_dz_zmax = 0.0;
  xclkg->dn_dz_min  = 0.0;
  xclkg->dn_dz_max  = 0.0;
  xclkg->nbar       = 0.0;
  xclkg->intr_shear = 0.0;
  xclkg->noise      = 0.0;
  xclkg->doia       = FALSE;
  xclkg->growth_ps  = NULL;
  xclkg->wl_comp    = NULL;
  xclkg->ia_comp    = NULL;
}

static void
_nc_xcor_kernel_weak_lensing_take_dndz (NcXcorKernelWeakLensing *xclkg, NcmSpline *dn_dz)
{
  NcmVector *z_vec      = ncm_spline_peek_xv (dn_dz);
  const guint z_vec_len = ncm_vector_len (z_vec);

  /* The caller's spline is kept as given, for the property and serialization; the
   * kernel normalizes a copy of its own, so setting the spline leaves the caller's
   * unchanged and a kernel rebuilt from its serialization normalizes the same values. */
  ncm_spline_clear (&xclkg->dn_dz_user);
  ncm_spline_clear (&xclkg->dn_dz);

  xclkg->dn_dz_user = dn_dz;
  xclkg->dn_dz      = ncm_spline_copy (dn_dz);
  xclkg->dn_dz_zmin = ncm_vector_get (z_vec, 0);
  xclkg->dn_dz_zmax = ncm_vector_get (z_vec, z_vec_len - 1);

  {
    const gdouble ngal = ncm_spline_eval_integ (xclkg->dn_dz, xclkg->dn_dz_zmin, xclkg->dn_dz_zmax);

    ncm_vector_scale (ncm_spline_peek_yv (xclkg->dn_dz), 1.0 / ngal);
    ncm_spline_prepare (xclkg->dn_dz);

    xclkg->dn_dz_min = ncm_spline_eval (xclkg->dn_dz, xclkg->dn_dz_zmin);
    xclkg->dn_dz_max = ncm_spline_eval (xclkg->dn_dz, xclkg->dn_dz_zmax);
  }

  nc_xcor_kernel_mark_outdated (NC_XCOR_KERNEL (xclkg));
}

static void
_nc_xcor_kernel_weak_lensing_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcXcorKernelWeakLensing *xclkg = NC_XCOR_KERNEL_WEAK_LENSING (object);

  g_return_if_fail (NC_IS_XCOR_KERNEL_WEAK_LENSING (object));

  switch (prop_id)
  {
    case PROP_DN_DZ:
      _nc_xcor_kernel_weak_lensing_take_dndz (xclkg, g_value_dup_object (value));
      break;
    case PROP_NBAR:
      xclkg->nbar = g_value_get_double (value);
      break;
    case PROP_INTR_SHEAR:
      xclkg->intr_shear = g_value_get_double (value);
      break;
    case PROP_DOIA:
      xclkg->doia = g_value_get_boolean (value);
      break;
    case PROP_GROWTH_POWSPEC:
      ncm_powspec_clear (&xclkg->growth_ps);
      xclkg->growth_ps = g_value_dup_object (value);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_nc_xcor_kernel_weak_lensing_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcXcorKernelWeakLensing *xclkg = NC_XCOR_KERNEL_WEAK_LENSING (object);

  g_return_if_fail (NC_IS_XCOR_KERNEL_WEAK_LENSING (object));

  switch (prop_id)
  {
    case PROP_DN_DZ:
      g_value_set_object (value, xclkg->dn_dz_user);
      break;
    case PROP_NBAR:
      g_value_set_double (value, xclkg->nbar);
      break;
    case PROP_INTR_SHEAR:
      g_value_set_double (value, xclkg->intr_shear);
      break;
    case PROP_DOIA:
      g_value_set_boolean (value, xclkg->doia);
      break;
    case PROP_GROWTH_POWSPEC:
      g_value_set_object (value, xclkg->growth_ps);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static gdouble _nc_xcor_kernel_weak_lensing_lens_eff_eval_source (NcXcorLensingEfficiency *lens_eff, gdouble z);
static void _nc_xcor_kernel_weak_lensing_lens_eff_get_z_range (NcXcorLensingEfficiency *lens_eff, gdouble *zmin, gdouble *zmax);

static void
_nc_xcor_kernel_weak_lensing_constructed (GObject *object)
{
  NcXcorKernelWeakLensing *xclkg = NC_XCOR_KERNEL_WEAK_LENSING (object);
  NcXcorKernel *xclk             = NC_XCOR_KERNEL (xclkg);

  /* Chain up : middle */
  G_OBJECT_CLASS (nc_xcor_kernel_weak_lensing_parent_class)->constructed (object);
  {
    {
      NcDistance *dist                             = nc_xcor_kernel_peek_dist (xclk);
      NcXcorKernelWeakLensingLensEff *lens_eff_obj = g_object_new (
        nc_xcor_kernel_weak_lensing_lens_eff_get_type (),
        "distance", dist,
        NULL);

      lens_eff_obj->data = xclkg;
      xclkg->lens_eff    = NC_XCOR_LENSING_EFFICIENCY (lens_eff_obj);
    }

    /* Noise level */
    xclkg->noise = gsl_pow_2 (xclkg->intr_shear) / xclkg->nbar;

    /* Create weak lensing component */
    {
      NcDistance *dist = nc_xcor_kernel_peek_dist (xclk);
      NcmPowspec *ps   = nc_xcor_kernel_peek_powspec (xclk);

      g_assert_null (xclkg->wl_comp);
      xclkg->wl_comp = _nc_xcor_kernel_component_weak_lensing_new (dist, ps, xclkg->lens_eff, xclkg->dn_dz_zmax);

      /* Component for the intrinsic alignment */
      if (xclkg->doia)
      {
        g_assert_null (xclkg->ia_comp);
        xclkg->ia_comp = _nc_xcor_kernel_component_ia_new (xclkg, dist, ps, _nc_xcor_kernel_weak_lensing_peek_growth_ps (xclkg));
      }
    }
  }
}

static void
_nc_xcor_kernel_weak_lensing_dispose (GObject *object)
{
  NcXcorKernelWeakLensing *xclkg = NC_XCOR_KERNEL_WEAK_LENSING (object);

  ncm_spline_clear (&xclkg->dn_dz_user);
  ncm_spline_clear (&xclkg->dn_dz);
  nc_xcor_lensing_efficiency_clear (&xclkg->lens_eff);
  nc_xcor_kernel_component_clear (&xclkg->wl_comp);
  nc_xcor_kernel_component_clear (&xclkg->ia_comp);
  ncm_powspec_clear (&xclkg->growth_ps);

  /* Chain up : end */
  G_OBJECT_CLASS (nc_xcor_kernel_weak_lensing_parent_class)->dispose (object);
}

static void
_nc_xcor_kernel_weak_lensing_finalize (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (nc_xcor_kernel_weak_lensing_parent_class)->finalize (object);
}

static gdouble _nc_xcor_kernel_weak_lensing_eval_limber_z (NcXcorKernel *xclk, NcHICosmo *cosmo, gdouble z, const NcXcorKinetic *xck, gint l);
static gdouble _nc_xcor_kernel_weak_lensing_eval_limber_z_prefactor (NcXcorKernel *xclk, NcHICosmo *cosmo, gint l);
static void _nc_xcor_kernel_weak_lensing_prepare (NcXcorKernel *xclk, NcHICosmo *cosmo);
static void _nc_xcor_kernel_weak_lensing_add_noise (NcXcorKernel *xclk, NcmVector *vp1, NcmVector *vp2, guint lmin);
guint _nc_xcor_kernel_weak_lensing_obs_len (NcXcorKernel *xclk);
guint _nc_xcor_kernel_weak_lensing_obs_params_len (NcXcorKernel *xclk);
static void _nc_xcor_kernel_weak_lensing_get_z_range (NcXcorKernel *xclk, gdouble *zmin, gdouble *zmax, gdouble *zmid);
static GPtrArray *_nc_xcor_kernel_weak_lensing_get_component_list (NcXcorKernel *xclk);

static void
nc_xcor_kernel_weak_lensing_class_init (NcXcorKernelWeakLensingClass *klass)
{
  GObjectClass *object_class      = G_OBJECT_CLASS (klass);
  NcXcorKernelClass *parent_class = NC_XCOR_KERNEL_CLASS (klass);
  NcmModelClass *model_class      = NCM_MODEL_CLASS (klass);

  model_class->set_property = &_nc_xcor_kernel_weak_lensing_set_property;
  model_class->get_property = &_nc_xcor_kernel_weak_lensing_get_property;
  object_class->constructed = &_nc_xcor_kernel_weak_lensing_constructed;
  object_class->dispose     = &_nc_xcor_kernel_weak_lensing_dispose;
  object_class->finalize    = &_nc_xcor_kernel_weak_lensing_finalize;

  ncm_model_class_set_name_nick (model_class, "Xcor weak lensing", "Xcor-WL");
  ncm_model_class_add_params (model_class, NC_XCOR_KERNEL_WEAK_LENSING_SPARAM_LEN, NC_XCOR_KERNEL_WEAK_LENSING_VPARAM_LEN, PROP_SIZE);

  g_object_class_install_property (object_class,
                                   PROP_DN_DZ,
                                   g_param_spec_object ("dndz",
                                                        NULL,
                                                        "Source redshift distribution",
                                                        NCM_TYPE_SPLINE,
                                                        G_PARAM_CONSTRUCT | G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_NBAR,
                                   g_param_spec_double ("nbar",
                                                        NULL,
                                                        "nbar (galaxy angular density)",
                                                        0.0, GSL_POSINF, 0.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  g_object_class_install_property (object_class,
                                   PROP_INTR_SHEAR,
                                   g_param_spec_double ("intr-shear",
                                                        NULL,
                                                        "Intrinsic galaxy shear",
                                                        0.0, 20.0, 0.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcXcorKernelWeakLensing:doia:
   *
   * Whether to include the nonlinear alignment (NLA) intrinsic shear. Adds a
   * component with kernel $-A_\mathrm{IA}(z) C_1 \rho_\mathrm{crit} \Omega_m
   * E(z)\, \mathrm{d}n/\mathrm{d}z / D(k, z)$ weighted by $j_\ell(k\chi) /
   * (k\chi)^2$, where $D(k, z) = \sqrt{P_D(k, z) / P_D(k, 0)}$ is the growth factor
   * read from #NcXcorKernelWeakLensing:growth-powspec and $A_\mathrm{IA}(z) = A_\mathrm{IA}
   * [(1+z)/(1+z_\mathrm{piv})]^{\eta_\mathrm{IA}}$.
   */
  g_object_class_install_property (object_class,
                                   PROP_DOIA,
                                   g_param_spec_boolean ("doia",
                                                         NULL,
                                                         "Do intrinsic alignment",
                                                         FALSE,
                                                         G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcXcorKernelWeakLensing:growth-powspec:
   *
   * The power spectrum the intrinsic alignment term reads its growth factor from,
   * $D(k, z) = \sqrt{P_D(k, z) / P_D(k, 0)}$; %NULL (the default) uses the kernel's
   * #NcXcorKernel:powspec. The NLA model divides by the linear growth factor: with a
   * nonlinear #NcXcorKernel:powspec, set this to the linear spectrum. Used only with
   * #NcXcorKernelWeakLensing:doia.
   */
  g_object_class_install_property (object_class,
                                   PROP_GROWTH_POWSPEC,
                                   g_param_spec_object ("growth-powspec",
                                                        NULL,
                                                        "Power spectrum of the intrinsic alignment growth factor",
                                                        NCM_TYPE_POWSPEC,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  ncm_model_class_set_sparam (model_class, NC_XCOR_KERNEL_WEAK_LENSING_IA_A,
                              "A_{\\mathrm{IA}}", "ia_A",
                              -10.0, 10.0, 0.1,
                              NC_XCOR_KERNEL_WEAK_LENSING_DEFAULT_PARAMS_ABSTOL, NC_XCOR_KERNEL_WEAK_LENSING_DEFAULT_IA_A, NCM_PARAM_TYPE_FIXED);

  ncm_model_class_set_sparam (model_class, NC_XCOR_KERNEL_WEAK_LENSING_IA_ETA,
                              "\\eta_{\\mathrm{IA}}", "ia_eta",
                              -10.0, 10.0, 0.1,
                              NC_XCOR_KERNEL_WEAK_LENSING_DEFAULT_PARAMS_ABSTOL, NC_XCOR_KERNEL_WEAK_LENSING_DEFAULT_IA_ETA, NCM_PARAM_TYPE_FIXED);

  ncm_model_class_set_sparam (model_class, NC_XCOR_KERNEL_WEAK_LENSING_IA_Z_PIV,
                              "z_{\\mathrm{piv}}", "ia_z_piv",
                              0.0, 10.0, 0.1,
                              NC_XCOR_KERNEL_WEAK_LENSING_DEFAULT_PARAMS_ABSTOL, NC_XCOR_KERNEL_WEAK_LENSING_DEFAULT_IA_Z_PIV, NCM_PARAM_TYPE_FIXED);

  /* Check for errors in parameters initialization */
  ncm_model_class_check_params_info (model_class);

  parent_class->eval_limber_z           = &_nc_xcor_kernel_weak_lensing_eval_limber_z;
  parent_class->eval_limber_z_prefactor = &_nc_xcor_kernel_weak_lensing_eval_limber_z_prefactor;
  parent_class->prepare                 = &_nc_xcor_kernel_weak_lensing_prepare;
  parent_class->add_noise               = &_nc_xcor_kernel_weak_lensing_add_noise;

  parent_class->obs_len            = &_nc_xcor_kernel_weak_lensing_obs_len;
  parent_class->obs_params_len     = &_nc_xcor_kernel_weak_lensing_obs_params_len;
  parent_class->get_z_range        = &_nc_xcor_kernel_weak_lensing_get_z_range;
  parent_class->get_component_list = &_nc_xcor_kernel_weak_lensing_get_component_list;

  ncm_model_class_add_impl_flag (model_class, NC_XCOR_KERNEL_IMPL_ALL);
}

static gdouble
_nc_xcor_kernel_weak_lensing_lens_eff_eval_source (NcXcorLensingEfficiency *lens_eff, gdouble z)
{
  NcXcorKernelWeakLensingLensEff *data_obj = NC_XCOR_KERNEL_WEAK_LENSING_LENS_EFF (lens_eff);
  NcXcorKernelWeakLensing *xclkg           = data_obj->data;

  return _nc_xcor_kernel_weak_lensing_dndz (xclkg, z);
}

static void
_nc_xcor_kernel_weak_lensing_lens_eff_get_z_range (NcXcorLensingEfficiency *lens_eff, gdouble *zmin, gdouble *zmax)
{
  NcXcorKernelWeakLensingLensEff *data_obj = NC_XCOR_KERNEL_WEAK_LENSING_LENS_EFF (lens_eff);
  NcXcorKernelWeakLensing *xclkg           = data_obj->data;

  *zmin = xclkg->dn_dz_zmin;
  *zmax = xclkg->dn_dz_zmax;
}

/*
 * Implementation of the weak lensing kernel helper functions
 */

static gdouble
_nc_xcor_kernel_weak_lensing_dndz (NcXcorKernelWeakLensing *xclkg, gdouble z)
{
  const gdouble alpha = 1.0e-2;

  if (z < xclkg->dn_dz_zmin)
    return xclkg->dn_dz_min * exp (-gsl_pow_2 ((z - xclkg->dn_dz_zmin) / alpha));

  if (z > xclkg->dn_dz_zmax)
    return xclkg->dn_dz_max * exp (-gsl_pow_2 ((z - xclkg->dn_dz_zmax) / alpha));

  return ncm_spline_eval (xclkg->dn_dz, z);
}

/* Returns C1 rho_crit A_IA(z), with C1 = 5e-14 h^-2 M_sun^-1 Mpc^3: the h^-2 cancels the
 * h^2 of the critical density, so the product is dimensionless. */
static gdouble
_nc_xcor_kernel_weak_lensing_ia_amplitude (NcXcorKernelWeakLensing *xclkg, gdouble z)
{
  const gdouble C1_rho_crit = 5.0e-14 * ncm_c_crit_mass_density_h2_solar_mass_Mpc3 ();

  return C1_rho_crit * IA_A * pow ((1.0 + z) / (1.0 + IA_Z_PIV), IA_ETA);
}

static NcmPowspec *
_nc_xcor_kernel_weak_lensing_peek_growth_ps (NcXcorKernelWeakLensing *xclkg)
{
  if (xclkg->growth_ps != NULL)
    return xclkg->growth_ps;

  return nc_xcor_kernel_peek_powspec (NC_XCOR_KERNEL (xclkg));
}

static void
_nc_xcor_kernel_weak_lensing_prepare (NcXcorKernel *xclk, NcHICosmo *cosmo)
{
  NcXcorKernelWeakLensing *xclkg = NC_XCOR_KERNEL_WEAK_LENSING (xclk);
  NcDistance *dist               = nc_xcor_kernel_peek_dist (xclk);
  NcmPowspec *ps                 = nc_xcor_kernel_peek_powspec (xclk);

  nc_distance_prepare_if_needed (dist, cosmo);
  ncm_powspec_prepare_if_needed (ps, NCM_MODEL (cosmo));

  if (xclkg->growth_ps != NULL)
    ncm_powspec_prepare_if_needed (xclkg->growth_ps, NCM_MODEL (cosmo));

  g_assert_cmpfloat (0.0, ==, 0.0);
  g_assert (gsl_finite (ncm_spline_eval (xclkg->dn_dz, 0.0)));
  g_assert (gsl_finite (ncm_spline_eval (xclkg->dn_dz, xclkg->dn_dz_zmax)));

  nc_xcor_lensing_efficiency_prepare (xclkg->lens_eff, cosmo);

  g_assert_nonnull (xclkg->wl_comp);
  nc_xcor_kernel_component_prepare (xclkg->wl_comp, cosmo);

  if (xclkg->doia)
  {
    g_assert_nonnull (xclkg->ia_comp);
    nc_xcor_kernel_component_prepare (xclkg->ia_comp, cosmo);
  }
}

/*
 * Implementation of the new Component interface.
 */

static void
_wl_component_data_clear (WeakLensingComponentData *data)
{
  /* No need to clear, these are weak references from parent kernel */
}

static gdouble
_wl_component_eval_kernel (NcXcorKernelComponent *comp, NcHICosmo *cosmo, gdouble chi, gdouble k)
{
  WeakLensingComponentData *data = _NC_XCOR_KERNEL_COMPONENT_WEAK_LENSING_GET_DATA (comp);
  const gdouble z                = nc_distance_inv_comoving (data->dist, cosmo, chi);
  const gdouble powspec          = ncm_powspec_eval (data->ps, NCM_MODEL (cosmo), z, k / nc_hicosmo_RH_Mpc (cosmo));
  const gdouble lens_eff_z       = nc_xcor_lensing_efficiency_eval (data->lens_eff, z);
  const gdouble kernel           = (1.0 + z) / chi * lens_eff_z;
  const gdouble operator_k       = 1.0 / gsl_pow_2 (k);

  return operator_k * kernel * sqrt (powspec);
}

static gdouble
_wl_component_eval_window (NcXcorKernelComponent *comp, NcHICosmo *cosmo, const NcXcorKinetic *xck)
{
  WeakLensingComponentData *data = _NC_XCOR_KERNEL_COMPONENT_WEAK_LENSING_GET_DATA (comp);

  return (1.0 + xck->z) / xck->chi_z * nc_xcor_lensing_efficiency_eval (data->lens_eff, xck->z);
}

static gdouble
_wl_component_eval_kfactor (NcXcorKernelComponent *comp, NcHICosmo *cosmo, const NcXcorKinetic *xck, gdouble k)
{
  WeakLensingComponentData *data = _NC_XCOR_KERNEL_COMPONENT_WEAK_LENSING_GET_DATA (comp);
  const gdouble powspec          = ncm_powspec_eval (data->ps, NCM_MODEL (cosmo), xck->z, k / nc_hicosmo_RH_Mpc (cosmo));

  return sqrt (powspec) / gsl_pow_2 (k);
}

static void
_wl_component_eval_kfactor_vec (NcXcorKernelComponent *comp, NcHICosmo *cosmo, const NcXcorKinetic *xck, NcmVector *k, NcmVector *out)
{
  WeakLensingComponentData *data = _NC_XCOR_KERNEL_COMPONENT_WEAK_LENSING_GET_DATA (comp);
  const gdouble RH_Mpc           = nc_hicosmo_RH_Mpc (cosmo);
  const guint len                = ncm_vector_len (k);
  guint i;

  ncm_vector_scale (k, 1.0 / RH_Mpc);
  ncm_powspec_eval_vec (data->ps, NCM_MODEL (cosmo), xck->z, k, out);
  ncm_vector_scale (k, RH_Mpc);

  for (i = 0; i < len; i++)
    ncm_vector_set (out, i, sqrt (ncm_vector_get (out, i)) / gsl_pow_2 (ncm_vector_get (k, i)));
}

static gdouble
_wl_component_eval_prefactor (NcXcorKernelComponent *comp, NcHICosmo *cosmo, gdouble k, gint l)
{
  const gdouble cosmo_factor = 1.5 * nc_hicosmo_Omega_m0 (cosmo);
  const gdouble lfactor      = sqrt ((l + 2.0) * (l + 1.0) * l * (l - 1.0));

  return cosmo_factor * lfactor;
}

static void
_wl_component_get_limits (NcXcorKernelComponent *comp, NcHICosmo *cosmo, gdouble *chi_min, gdouble *chi_max, gdouble *k_min, gdouble *k_max)
{
  WeakLensingComponentData *data = _NC_XCOR_KERNEL_COMPONENT_WEAK_LENSING_GET_DATA (comp);
  NcDistance *dist               = data->dist;
  NcmPowspec *ps                 = data->ps;

  nc_distance_prepare_if_needed (dist, cosmo);
  ncm_powspec_prepare_if_needed (ps, NCM_MODEL (cosmo));

  *chi_min = nc_distance_comoving (dist, cosmo, 1.0e-6);
  *chi_max = nc_distance_comoving (dist, cosmo, data->dn_dz_zmax);
  *k_min   = ncm_powspec_get_kmin (ps) * nc_hicosmo_RH_Mpc (cosmo);
  *k_max   = ncm_powspec_get_kmax (ps) * nc_hicosmo_RH_Mpc (cosmo);
}

static NcXcorKernelComponent *
_nc_xcor_kernel_component_weak_lensing_new (NcDistance *dist, NcmPowspec *ps, NcXcorLensingEfficiency *lens_eff, gdouble dn_dz_zmax)
{
  NcXcorKernelComponent *comp    = g_object_new (nc_xcor_kernel_component_weak_lensing_get_type (), NULL);
  WeakLensingComponentData *data = _NC_XCOR_KERNEL_COMPONENT_WEAK_LENSING_GET_DATA (comp);

  data->dist       = dist;
  data->ps         = ps;
  data->lens_eff   = lens_eff;
  data->dn_dz_zmax = dn_dz_zmax;

  return comp;
}

/*
 * Implementation of the intrinsic alignment component.
 *
 * The kernel is the density kernel E(z) dn_dz(z) times -C1 rho_crit A_IA(z) Omega_m /
 * D(k, z), with j_l(k chi) / (k chi)^2 as the lensing term: it has the 1 / k^2 of the
 * lensing component and a 1 / chi^2 in its window. The growth factor is read from
 * the growth power spectrum P_D, D(k, z) = sqrt (P_D(k, z) / P_D(k, 0)), so the k factor
 * is sqrt (P(k, z) P_D(k, 0) / P_D(k, z)) / k^2. Omega_m and the spin-2 factor are in
 * the prefactor.
 */

static void
_ia_component_data_clear (IAComponentData *data)
{
  /* No need to clear, these are weak references from parent kernel */
}

static gdouble
_ia_component_eval_kernel (NcXcorKernelComponent *comp, NcHICosmo *cosmo, gdouble chi, gdouble k)
{
  IAComponentData *data = _NC_XCOR_KERNEL_COMPONENT_IA_GET_DATA (comp);
  const gdouble z       = nc_distance_inv_comoving (data->dist, cosmo, chi);
  const gdouble E_z     = nc_hicosmo_E (cosmo, z);
  const gdouble k_RH    = k / nc_hicosmo_RH_Mpc (cosmo);
  const gdouble powspec = ncm_powspec_eval (data->ps, NCM_MODEL (cosmo), z, k_RH);
  const gdouble P_D_z   = ncm_powspec_eval (data->growth_ps, NCM_MODEL (cosmo), z, k_RH);
  const gdouble P_D_0   = ncm_powspec_eval (data->growth_ps, NCM_MODEL (cosmo), 0.0, k_RH);
  const gdouble dn_dz_z = _nc_xcor_kernel_weak_lensing_dndz (data->xclkg, z);
  const gdouble A_z     = _nc_xcor_kernel_weak_lensing_ia_amplitude (data->xclkg, z);

  return -A_z *dn_dz_z *E_z *sqrt (powspec *P_D_0 / P_D_z) / gsl_pow_2 (k * chi);
}

static gdouble
_ia_component_eval_window (NcXcorKernelComponent *comp, NcHICosmo *cosmo, const NcXcorKinetic *xck)
{
  IAComponentData *data = _NC_XCOR_KERNEL_COMPONENT_IA_GET_DATA (comp);
  const gdouble dn_dz_z = _nc_xcor_kernel_weak_lensing_dndz (data->xclkg, xck->z);
  const gdouble A_z     = _nc_xcor_kernel_weak_lensing_ia_amplitude (data->xclkg, xck->z);

  return -A_z * dn_dz_z * xck->E_z / gsl_pow_2 (xck->chi_z);
}

static gdouble
_ia_component_eval_kfactor (NcXcorKernelComponent *comp, NcHICosmo *cosmo, const NcXcorKinetic *xck, gdouble k)
{
  IAComponentData *data = _NC_XCOR_KERNEL_COMPONENT_IA_GET_DATA (comp);
  const gdouble k_RH    = k / nc_hicosmo_RH_Mpc (cosmo);
  const gdouble powspec = ncm_powspec_eval (data->ps, NCM_MODEL (cosmo), xck->z, k_RH);
  const gdouble P_D_z   = ncm_powspec_eval (data->growth_ps, NCM_MODEL (cosmo), xck->z, k_RH);
  const gdouble P_D_0   = ncm_powspec_eval (data->growth_ps, NCM_MODEL (cosmo), 0.0, k_RH);

  return sqrt (powspec * P_D_0 / P_D_z) / gsl_pow_2 (k);
}

static void
_ia_component_eval_kfactor_vec (NcXcorKernelComponent *comp, NcHICosmo *cosmo, const NcXcorKinetic *xck, NcmVector *k, NcmVector *out)
{
  IAComponentData *data = _NC_XCOR_KERNEL_COMPONENT_IA_GET_DATA (comp);
  const gdouble RH_Mpc  = nc_hicosmo_RH_Mpc (cosmo);
  const guint len       = ncm_vector_len (k);
  guint i;

  ncm_vector_scale (k, 1.0 / RH_Mpc);
  ncm_powspec_eval_vec (data->ps, NCM_MODEL (cosmo), xck->z, k, out);

  for (i = 0; i < len; i++)
  {
    const gdouble k_i   = ncm_vector_get (k, i);
    const gdouble P_D_z = ncm_powspec_eval (data->growth_ps, NCM_MODEL (cosmo), xck->z, k_i);
    const gdouble P_D_0 = ncm_powspec_eval (data->growth_ps, NCM_MODEL (cosmo), 0.0, k_i);

    ncm_vector_set (out, i, sqrt (ncm_vector_get (out, i) * P_D_0 / P_D_z) / gsl_pow_2 (k_i * RH_Mpc));
  }

  ncm_vector_scale (k, RH_Mpc);
}

static gdouble
_ia_component_eval_prefactor (NcXcorKernelComponent *comp, NcHICosmo *cosmo, gdouble k, gint l)
{
  const gdouble lfactor = sqrt ((l + 2.0) * (l + 1.0) * l * (l - 1.0));

  return nc_hicosmo_Omega_m0 (cosmo) * lfactor;
}

static void
_ia_component_get_limits (NcXcorKernelComponent *comp, NcHICosmo *cosmo, gdouble *chi_min, gdouble *chi_max, gdouble *k_min, gdouble *k_max)
{
  IAComponentData *data          = _NC_XCOR_KERNEL_COMPONENT_IA_GET_DATA (comp);
  NcDistance *dist               = data->dist;
  NcmPowspec *ps                 = data->ps;
  NcXcorKernelWeakLensing *xclkg = data->xclkg;

  nc_distance_prepare_if_needed (dist, cosmo);
  ncm_powspec_prepare_if_needed (ps, NCM_MODEL (cosmo));
  ncm_powspec_prepare_if_needed (data->growth_ps, NCM_MODEL (cosmo));

  /* Bounded by the dn_dz support, as the clustering component of NcXcorKernelGal. */
  *chi_min = nc_distance_comoving (dist, cosmo, MAX (xclkg->dn_dz_zmin, 1.0e-6));
  *chi_max = nc_distance_comoving (dist, cosmo, xclkg->dn_dz_zmax);
  *k_min   = ncm_powspec_get_kmin (ps) * nc_hicosmo_RH_Mpc (cosmo);
  *k_max   = ncm_powspec_get_kmax (ps) * nc_hicosmo_RH_Mpc (cosmo);
}

static NcXcorKernelComponent *
_nc_xcor_kernel_component_ia_new (NcXcorKernelWeakLensing *xclkg, NcDistance *dist, NcmPowspec *ps, NcmPowspec *growth_ps)
{
  NcXcorKernelComponent *comp = g_object_new (nc_xcor_kernel_component_ia_get_type (), NULL);
  IAComponentData *data       = _NC_XCOR_KERNEL_COMPONENT_IA_GET_DATA (comp);

  data->xclkg     = xclkg;
  data->dist      = dist;
  data->ps        = ps;
  data->growth_ps = growth_ps;

  return comp;
}

/*
 * Kernel implementation.
 */

/* chi^2 times the radial weight (1 + z) g(z) / chi, with chi cancelled: written as the
 * product, chi = 0 at z = 0 gives 0 * inf instead of the limit 0.
 *
 * The intrinsic alignment term is its window times chi^2, divided by the 3/2 Omega_m of
 * the prefactor, over D(k, z) of the growth power spectrum at the Limber k = (l + 1/2) /
 * chi. It is skipped where dn_dz vanishes, which includes chi = 0, where k is infinite. */
static gdouble
_nc_xcor_kernel_weak_lensing_eval_limber_z (NcXcorKernel *xclk, NcHICosmo *cosmo, gdouble z, const NcXcorKinetic *xck, gint l)
{
  NcXcorKernelWeakLensing *xclkg = NC_XCOR_KERNEL_WEAK_LENSING (xclk);
  gdouble res                    = xck->chi_z * (1.0 + z) * nc_xcor_lensing_efficiency_eval (xclkg->lens_eff, z);

  if (xclkg->doia)
  {
    const gdouble dn_dz_z = _nc_xcor_kernel_weak_lensing_dndz (xclkg, z);

    if (dn_dz_z != 0.0)
    {
      NcmPowspec *ps     = _nc_xcor_kernel_weak_lensing_peek_growth_ps (xclkg);
      const gdouble k    = (l + 0.5) / (xck->chi_z * nc_hicosmo_RH_Mpc (cosmo));
      const gdouble D_kz = sqrt (ncm_powspec_eval (ps, NCM_MODEL (cosmo), z, k) /
                                 ncm_powspec_eval (ps, NCM_MODEL (cosmo), 0.0, k));
      const gdouble A_z = _nc_xcor_kernel_weak_lensing_ia_amplitude (xclkg, z);

      res -= A_z * dn_dz_z * xck->E_z / (1.5 * D_kz);
    }
  }

  return res;
}

static gdouble
_nc_xcor_kernel_weak_lensing_eval_limber_z_prefactor (NcXcorKernel *xclk, NcHICosmo *cosmo, gint l)
{
  const gdouble lfactor = sqrt ((l + 2.0) * (l + 1.0) * l * (l - 1.0));
  const gdouble nu      = l + 0.5;

  return 1.5 * nc_hicosmo_Omega_m0 (cosmo) * lfactor / (nu * nu);
}

static void
_nc_xcor_kernel_weak_lensing_add_noise (NcXcorKernel *xclk, NcmVector *vp1, NcmVector *vp2, guint lmin)
{
  NcXcorKernelWeakLensing *xclkg = NC_XCOR_KERNEL_WEAK_LENSING (xclk);

  /* ncm_vector_add_constant (vp1, NOISE_BIAS); / * take noise_bias into account * / */
  ncm_vector_memcpy (vp2, vp1);
  ncm_vector_add_constant (vp2, xclkg->noise);
}

guint
_nc_xcor_kernel_weak_lensing_obs_len (NcXcorKernel *xclk)
{
  return 2;
}

guint
_nc_xcor_kernel_weak_lensing_obs_params_len (NcXcorKernel *xclk)
{
  return 1;
}

static void
_nc_xcor_kernel_weak_lensing_get_z_range (NcXcorKernel *xclk, gdouble *zmin, gdouble *zmax, gdouble *zmid)
{
  NcXcorKernelWeakLensing *xclkg = NC_XCOR_KERNEL_WEAK_LENSING (xclk);

  *zmin = 0.0;
  *zmax = xclkg->dn_dz_zmax;
  *zmid = ncm_vector_get (ncm_spline_peek_xv (xclkg->dn_dz), ncm_vector_get_max_index (ncm_spline_peek_yv (xclkg->dn_dz))) / 2.0;
}

static GPtrArray *
_nc_xcor_kernel_weak_lensing_get_component_list (NcXcorKernel *xclk)
{
  NcXcorKernelWeakLensing *xclkg = NC_XCOR_KERNEL_WEAK_LENSING (xclk);
  GPtrArray *comp_list           = g_ptr_array_new_with_free_func (g_object_unref);

  if (xclkg->wl_comp != NULL)
    g_ptr_array_add (comp_list, g_object_ref (xclkg->wl_comp));

  if (xclkg->doia && (xclkg->ia_comp != NULL))
    g_ptr_array_add (comp_list, g_object_ref (xclkg->ia_comp));

  return comp_list;
}

/**
 * nc_xcor_kernel_weak_lensing_new:
 * @dist: a #NcDistance
 * @ps: a #NcmPowspec
 * @dn_dz: a #NcmSpline
 * @nbar: a gdouble, gal density
 * @intr_shear: a gdouble, intrinsic galaxy shear
 *
 * Returns: a #NcXcorKernelWeakLensing
 */
NcXcorKernelWeakLensing *
nc_xcor_kernel_weak_lensing_new (NcDistance *dist, NcmPowspec *ps, NcmSpline *dn_dz, gdouble nbar, gdouble intr_shear)
{
  NcXcorKernelWeakLensing *xclkg = g_object_new (NC_TYPE_XCOR_KERNEL_WEAK_LENSING,
                                                 "dist", dist,
                                                 "powspec", ps,
                                                 "dndz", dn_dz,
                                                 "nbar", nbar,
                                                 "intr-shear", intr_shear,
                                                 NULL);

  return xclkg;
}

