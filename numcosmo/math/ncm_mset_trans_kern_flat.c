/***************************************************************************
 *            ncm_mset_trans_kern_flat.c
 *
 *  Thu October 02 13:37:11 2014
 *  Copyright  2014  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * ncm_mset_trans_kern_flat.c
 * Copyright (C) 2014 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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
 * SECTION:ncm_mset_trans_kern_flat
 * @title: NcmMSetTransKernFlat
 * @short_description: Multivariate flat sampler
 *
 * This object subclasses NcmMSetTransKern and implements a multivariate flat sampler.
 *
 * Implementation of a multivariate flat sampler, offering a straightforward method for
 * generating random parameter vectors with multivariate parameters. This sampler
 * generates vectors uniformly distributed in the hypercube defined by the specified
 * bounds for each parameter.
 *
 * **Key Functionality:**
 *
 * - Generates random parameter vectors with multivariate parameters.
 * - Utilizes a simple flat sampling method within the hypercube defined by parameter
 *   bounds.
 *
 * This implementation is particularly useful when a basic flat sampling approach is
 * needed for generating random parameter vectors with multivariate parameters,
 * especially in scenarios where limited information is available about the posterior
 * distribution.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_mset_trans_kern_flat.h"
#include "math/ncm_c.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_randist.h>
#endif /* NUMCOSMO_GIR_SCAN */

enum
{
  PROP_0,
};

struct _NcmMSetTransKernFlat
{
  /*< private >*/
  NcmMSetTransKern parent_instance;
  gdouble parea;
};

G_DEFINE_TYPE (NcmMSetTransKernFlat, ncm_mset_trans_kern_flat, NCM_TYPE_MSET_TRANS_KERN)

static void
ncm_mset_trans_kern_flat_init (NcmMSetTransKernFlat *tkernf)
{
  tkernf->parea = 0.0;
}

static void
_ncm_mset_trans_kern_flat_finalize (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (ncm_mset_trans_kern_flat_parent_class)->finalize (object);
}

static void _ncm_mset_trans_kern_flat_set_mset (NcmMSetTransKern *tkern, NcmMSet *mset);
static void _ncm_mset_trans_kern_flat_generate (NcmMSetTransKern *tkern, NcmVector *theta, NcmVector *thetastar, NcmRNG *rng);
static gdouble _ncm_mset_trans_kern_flat_pdf (NcmMSetTransKern *tkern, NcmVector *theta, NcmVector *thetastar);
static const gchar *_ncm_mset_trans_kern_flat_get_name (NcmMSetTransKern *tkern);

static void
ncm_mset_trans_kern_flat_class_init (NcmMSetTransKernFlatClass *klass)
{
  GObjectClass *object_class         = G_OBJECT_CLASS (klass);
  NcmMSetTransKernClass *tkern_class = NCM_MSET_TRANS_KERN_CLASS (klass);

  object_class->finalize = &_ncm_mset_trans_kern_flat_finalize;

  tkern_class->set_mset = &_ncm_mset_trans_kern_flat_set_mset;
  tkern_class->generate = &_ncm_mset_trans_kern_flat_generate;
  tkern_class->pdf      = &_ncm_mset_trans_kern_flat_pdf;
  tkern_class->get_name = &_ncm_mset_trans_kern_flat_get_name;
}

static void
_ncm_mset_trans_kern_flat_set_mset (NcmMSetTransKern *tkern, NcmMSet *mset)
{
  NcmMSetTransKernFlat *tkernf = NCM_MSET_TRANS_KERN_FLAT (tkern);
  NcmMSet *mset0               = ncm_mset_trans_kern_peek_mset (tkern);
  guint fparam_len             = ncm_mset_fparam_len (mset0);
  guint i;

  tkernf->parea = 1.0;

  for (i = 0; i < fparam_len; i++)
  {
    const gdouble lb = ncm_mset_fparam_get_lower_bound (mset0, i);
    const gdouble ub = ncm_mset_fparam_get_upper_bound (mset0, i);

    tkernf->parea *= ub - lb;
    g_assert_cmpfloat (tkernf->parea, >, 0.0);
  }
}

static void
_ncm_mset_trans_kern_flat_generate (NcmMSetTransKern *tkern, NcmVector *theta, NcmVector *thetastar, NcmRNG *rng)
{
  NcmMSet *mset    = ncm_mset_trans_kern_peek_mset (tkern);
  guint fparam_len = ncm_mset_fparam_len (mset);
  guint i;

  for (i = 0; i < fparam_len; i++)
  {
    const gdouble lb  = ncm_mset_fparam_get_lower_bound (mset, i);
    const gdouble ub  = ncm_mset_fparam_get_upper_bound (mset, i);
    const gdouble val = lb + (ub - lb) * ncm_rng_uniform01_pos_gen (rng);

    ncm_vector_set (thetastar, i, val);
  }
}

static gdouble
_ncm_mset_trans_kern_flat_pdf (NcmMSetTransKern *tkern, NcmVector *theta, NcmVector *thetastar)
{
  return 1.0 / NCM_MSET_TRANS_KERN_FLAT (tkern)->parea;
}

static const gchar *
_ncm_mset_trans_kern_flat_get_name (NcmMSetTransKern *tkern)
{
  return "Multivariate Flat Sampler";
}

/**
 * ncm_mset_trans_kern_flat_new:
 *
 * New NcmMSetTransKern flat.
 *
 * Returns: (transfer full): a new #NcmMSetTransKernFlat.
 *
 */
NcmMSetTransKernFlat *
ncm_mset_trans_kern_flat_new (void)
{
  NcmMSetTransKernFlat *tkernf = g_object_new (NCM_TYPE_MSET_TRANS_KERN_FLAT,
                                               NULL);

  return tkernf;
}

