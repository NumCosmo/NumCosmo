/***************************************************************************
 *            nc_powspec_ml_cbe.c
 *
 *  Tue April 05 10:42:14 2016
 *  Copyright  2016  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * nc_powspec_ml_cbe.c
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

/**
 * SECTION:nc_powspec_ml_cbe
 * @title: NcPowspecMLCBE
 * @short_description: linear matter power spectrum from CLASS backend.
 * 
 * Provides the linear matter power spectrum using the CLASS backend #NcCBE.
 * 
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "nc_powspec_ml_cbe.h"
#include "nc_powspec_ml_transfer.h"
#include "lss/nc_transfer_func_eh.h"
#include "nc_hiprim.h"

enum
{
  PROP_0,
  PROP_CBE,
	PROP_SIZE
};

G_DEFINE_TYPE (NcPowspecMLCBE, nc_powspec_ml_cbe, NC_TYPE_POWSPEC_ML);

static void
nc_powspec_ml_cbe_init (NcPowspecMLCBE *ps_cbe)
{
  NcTransferFunc *tf = nc_transfer_func_eh_new ();
  ps_cbe->cbe  = NULL;
  ps_cbe->lnPk = NULL;
  ps_cbe->eh   = NC_POWSPEC_ML (nc_powspec_ml_transfer_new (tf));

  nc_transfer_func_free (tf);
}

static void
_nc_powspec_ml_cbe_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcPowspecMLCBE *ps_cbe = NC_POWSPEC_ML_CBE (object);
  g_return_if_fail (NC_IS_POWSPEC_ML_CBE (object));

  switch (prop_id)
  {
    case PROP_CBE:
      nc_powspec_ml_cbe_set_cbe (ps_cbe, g_value_get_object (value));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_powspec_ml_cbe_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcPowspecMLCBE *ps_cbe = NC_POWSPEC_ML_CBE (object);
  g_return_if_fail (NC_IS_POWSPEC_ML_CBE (object));

  switch (prop_id)
  {
    case PROP_CBE:
      g_value_set_object (value, nc_powspec_ml_cbe_peek_cbe (ps_cbe));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_powspec_ml_cbe_constructed (GObject *object)
{
  /* Chain up : start */
  G_OBJECT_CLASS (nc_powspec_ml_cbe_parent_class)->constructed (object);
  {
    NcPowspecMLCBE *ps_cbe = NC_POWSPEC_ML_CBE (object);
    if (ps_cbe->cbe == NULL)
    {
      ps_cbe->cbe = nc_cbe_new ();
    }
  }
}

static void
_nc_powspec_ml_cbe_dispose (GObject *object)
{
  NcPowspecMLCBE *ps_cbe = NC_POWSPEC_ML_CBE (object);

  nc_cbe_clear (&ps_cbe->cbe);
  ncm_spline2d_clear (&ps_cbe->lnPk);
  nc_powspec_ml_clear (&ps_cbe->eh);
  
  /* Chain up : end */
  G_OBJECT_CLASS (nc_powspec_ml_cbe_parent_class)->dispose (object);
}

static void
_nc_powspec_ml_cbe_finalize (GObject *object)
{

  /* Chain up : end */
  G_OBJECT_CLASS (nc_powspec_ml_cbe_parent_class)->finalize (object);
}

static void _nc_powspec_ml_cbe_prepare (NcmPowspec *powspec, NcmModel *model);
static gdouble _nc_powspec_ml_cbe_eval (NcmPowspec *powspec, NcmModel *model, const gdouble z, const gdouble k);
static void _nc_powspec_ml_cbe_get_nknots (NcmPowspec *powspec, guint *Nz, guint *Nk);

static void
nc_powspec_ml_cbe_class_init (NcPowspecMLCBEClass *klass)
{
  GObjectClass *object_class     = G_OBJECT_CLASS (klass);
  NcmPowspecClass *powspec_class = NCM_POWSPEC_CLASS (klass);

  object_class->set_property = &_nc_powspec_ml_cbe_set_property;
  object_class->get_property = &_nc_powspec_ml_cbe_get_property;

  object_class->constructed  = &_nc_powspec_ml_cbe_constructed;
  object_class->dispose      = &_nc_powspec_ml_cbe_dispose;
  object_class->finalize     = &_nc_powspec_ml_cbe_finalize;

  g_object_class_install_property (object_class,
                                   PROP_CBE,
                                   g_param_spec_object ("cbe",
                                                        NULL,
                                                        "Class backend object",
                                                        NC_TYPE_CBE,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  
  powspec_class->prepare    = &_nc_powspec_ml_cbe_prepare;
  powspec_class->eval       = &_nc_powspec_ml_cbe_eval;
  powspec_class->get_nknots = &_nc_powspec_ml_cbe_get_nknots;
}

static void 
_nc_powspec_ml_cbe_prepare (NcmPowspec *powspec, NcmModel *model)
{
  NcHICosmo *cosmo       = NC_HICOSMO (model);
  NcPowspecMLCBE *ps_cbe = NC_POWSPEC_ML_CBE (powspec);

  g_assert (NC_IS_HICOSMO (model));
  g_assert (ncm_model_peek_submodel_by_mid (model, nc_hiprim_id ()) != NULL);

  nc_cbe_set_calc_transfer (ps_cbe->cbe, TRUE);
  nc_cbe_set_max_matter_pk_z (ps_cbe->cbe, powspec->zf);
  nc_cbe_set_max_matter_pk_k (ps_cbe->cbe, NC_POWSPEC_ML_CBE_INTERN_KMAX);

  nc_cbe_prepare_if_needed (ps_cbe->cbe, cosmo);

  ncm_spline2d_clear (&ps_cbe->lnPk);

  ps_cbe->lnPk = nc_cbe_get_matter_ps (ps_cbe->cbe);

  ncm_powspec_prepare_if_needed (NCM_POWSPEC (ps_cbe->eh), model);
}

static gdouble 
_nc_powspec_ml_cbe_eval (NcmPowspec *powspec, NcmModel *model, const gdouble z, const gdouble k)
{
  NcPowspecMLCBE *ps_cbe = NC_POWSPEC_ML_CBE (powspec);
  if (k < NC_POWSPEC_ML_CBE_INTERN_KMIN)
  {
    const gdouble lnkmin = log (NC_POWSPEC_ML_CBE_INTERN_KMIN);
    const gdouble match = exp (ncm_spline2d_eval (ps_cbe->lnPk, lnkmin, z)) / ncm_powspec_eval (NCM_POWSPEC (ps_cbe->eh), model, z, NC_POWSPEC_ML_CBE_INTERN_KMIN);

    return match * ncm_powspec_eval (NCM_POWSPEC (ps_cbe->eh), model, z, k);
  }
  else if (k > NC_POWSPEC_ML_CBE_INTERN_KMAX)
  {
    const gdouble lnkmax = log (NC_POWSPEC_ML_CBE_INTERN_KMAX);
    const gdouble match = exp (ncm_spline2d_eval (ps_cbe->lnPk, lnkmax, z)) / ncm_powspec_eval (NCM_POWSPEC (ps_cbe->eh), model, z, NC_POWSPEC_ML_CBE_INTERN_KMAX);

    return match * ncm_powspec_eval (NCM_POWSPEC (ps_cbe->eh), model, z, k);
  }
  else
  {
    return exp (ncm_spline2d_eval (ps_cbe->lnPk, log (k), z));
  }
}

static void 
_nc_powspec_ml_cbe_get_nknots (NcmPowspec *powspec, guint *Nz, guint *Nk)
{
  NcPowspecMLCBE *ps_cbe = NC_POWSPEC_ML_CBE (powspec);

  Nz[0] = ncm_vector_len (ps_cbe->lnPk->yv);
  Nk[0] = ncm_vector_len (ps_cbe->lnPk->xv);
}

/**
 * nc_powspec_ml_cbe_new:
 * 
 * Creates a new #NcPowspecMLCBE from a new #NcCBE.
 * 
 * Returns: (transfer full): the newly created #NcPowspecMLCBE.
 */
NcPowspecMLCBE *
nc_powspec_ml_cbe_new (void)
{
  NcPowspecMLCBE *ps_cbe = g_object_new (NC_TYPE_POWSPEC_ML_CBE,
                                         NULL);

  return ps_cbe;
}

/**
 * nc_powspec_ml_cbe_new_full:
 * @cbe: a #NcCBE
 * 
 * Creates a new #NcPowspecMLCBE from @cbe.
 * 
 * Returns: (transfer full): the newly created #NcPowspecMLCBE.
 */
NcPowspecMLCBE *
nc_powspec_ml_cbe_new_full (NcCBE *cbe)
{
  NcPowspecMLCBE *ps_cbe = g_object_new (NC_TYPE_POWSPEC_ML_CBE,
                                         "cbe", cbe, 
                                         NULL);

  return ps_cbe;
}


/**
 * nc_powspec_ml_cbe_set_tf:
 * @ps_cbe: a #NcPowspecMLCBE
 * @cbe: a #NcCBE
 * 
 * Sets the #NcCBE to @cbe.
 * 
 */
void 
nc_powspec_ml_cbe_set_cbe (NcPowspecMLCBE *ps_cbe, NcCBE *cbe)
{
  g_clear_object (&ps_cbe->cbe);
  ps_cbe->cbe = nc_cbe_ref (cbe);
}

/**
 * nc_powspec_ml_cbe_peek_cbe:
 * @ps_cbe: a #NcPowspecMLCBE
 * 
 * Peeks the #NcCBE inside @ps_cbe.
 * 
 * Returns: (transfer none): the #NcCBE inside @ps_cbe.
 */
NcCBE *
nc_powspec_ml_cbe_peek_cbe (NcPowspecMLCBE *ps_cbe)
{
  return ps_cbe->cbe;
}

