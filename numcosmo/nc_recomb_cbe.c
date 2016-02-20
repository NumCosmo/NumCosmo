/***************************************************************************
 *            nc_recomb_cbe.c
 *
 *  Wed November 18 15:26:55 2015
 *  Copyright  2015  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>>
 ****************************************************************************/
/*
 * nc_recomb_cbe.c
 * Copyright (C) 2015 Sandro Dias Pinto Vitenti <sandro@isoftware.com.br>
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
 * SECTION:nc_recomb_cbe
 * @title: NcRecombCBE
 * @short_description: Cosmic recombination by Class.
 * @include: numcosmo/nc_recomb_cbe.h
 *
 * Cosmic recobination as implemeted by Class.
 * For more details see: #NcCBE.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "nc_recomb_cbe.h"
#include "math/ncm_spline_cubic_notaknot.h"

enum
{
  PROP_0,
  PROP_CBE
};

G_DEFINE_TYPE (NcRecombCBE, nc_recomb_cbe, NC_TYPE_RECOMB);

static void
nc_recomb_cbe_init (NcRecombCBE *recomb_cbe)
{
  NcRecomb *recomb = NC_RECOMB (recomb_cbe);
  
  recomb_cbe->cbe = NULL;

  recomb->Xe_s = ncm_spline_cubic_notaknot_new ();
}

static void
nc_recomb_cbe_dispose (GObject *object)
{
  NcRecombCBE *recomb_cbe = NC_RECOMB_CBE (object);

  nc_cbe_clear (&recomb_cbe->cbe);
  
	/* Chain up : end */  
  G_OBJECT_CLASS (nc_recomb_cbe_parent_class)->dispose (object);
}

static void
nc_recomb_cbe_finalize (GObject *object)
{

	/* Chain up : end */  
  G_OBJECT_CLASS (nc_recomb_cbe_parent_class)->finalize (object);
}

static void
nc_recomb_cbe_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcRecombCBE *recomb_cbe = NC_RECOMB_CBE (object);
  g_return_if_fail (NC_IS_RECOMB_CBE (object));

  switch (prop_id)
  {
    case PROP_CBE:
      nc_recomb_cbe_set_cbe (recomb_cbe, g_value_get_object (value));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
nc_recomb_cbe_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcRecombCBE *recomb_cbe = NC_RECOMB_CBE (object);
  g_return_if_fail (NC_IS_RECOMB_CBE (object));

  switch (prop_id)
  {
    case PROP_CBE:
      g_value_set_object (value, nc_recomb_cbe_peek_cbe (recomb_cbe));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void nc_recomb_cbe_prepare (NcRecomb *recomb, NcHICosmo *cosmo);

static void
nc_recomb_cbe_class_init (NcRecombCBEClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);
	NcRecombClass *recomb_class = NC_RECOMB_CLASS (klass);

  object_class->set_property = nc_recomb_cbe_set_property;
  object_class->get_property = nc_recomb_cbe_get_property;
  object_class->dispose      = nc_recomb_cbe_dispose;
  object_class->finalize     = nc_recomb_cbe_finalize;

  g_object_class_install_property (object_class,
                                   PROP_CBE,
                                   g_param_spec_object ("cbe",
                                                        NULL,
                                                        "Class backend",
                                                        NC_TYPE_CBE,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

	recomb_class->prepare = &nc_recomb_cbe_prepare;
}

static void
nc_recomb_cbe_prepare (NcRecomb *recomb, NcHICosmo *cosmo)
{
  NcRecombCBE *recomb_cbe = NC_RECOMB_CBE (recomb);

  nc_cbe_thermodyn_prepare (recomb_cbe->cbe, cosmo);

  ncm_spline_clear (&recomb->Xe_s);
  ncm_spline_clear (&recomb->tau_s);
  ncm_spline_clear (&recomb->dtau_dlambda_s);

  recomb->Xe_s           = nc_cbe_thermodyn_get_Xe (recomb_cbe->cbe);
  recomb->tau_s          = ncm_spline_copy_empty (recomb->Xe_s);
	recomb->dtau_dlambda_s = ncm_spline_copy_empty (recomb->Xe_s);

	_nc_recomb_prepare_tau_splines (recomb, cosmo);
}

/**
 * nc_recomb_cbe_new:
 * 
 * Creates a new #NcRecombCBE using default properties.
 * 
 * Returns: (transfer full): a new #NcRecombCBE.
 */
NcRecombCBE *
nc_recomb_cbe_new (void)
{
  NcCBE *cbe = nc_cbe_new ();
  return g_object_new (NC_TYPE_RECOMB_CBE,
                       "cbe", cbe, 
                       NULL);
}

/**
 * nc_recomb_cbe_full_new:
 * @cbe: a #NcCBE object
 * 
 * Creates a new #NcRecombCBE using default properties
 * and @cbe as the Class backend object #NcCBE.
 * 
 * Returns: (transfer full): a new #NcRecombCBE.
 */
NcRecombCBE *
nc_recomb_cbe_full_new (NcCBE *cbe)
{
  return g_object_new (NC_TYPE_RECOMB_CBE,
                       "cbe", cbe, 
                       NULL);
}

/**
 * nc_recomb_cbe_ref:
 * @recomb_cbe: a #NcRecombCBE.
 *
 * Increases the reference count of @recomb_cbe.
 *
 * Returns: (transfer full): @recomb_cbe.
 */
NcRecombCBE *
nc_recomb_cbe_ref (NcRecombCBE *recomb_cbe)
{
  return NC_RECOMB_CBE (g_object_ref (recomb_cbe));
}

/**
 * nc_recomb_cbe_free:
 * @recomb_cbe: a #NcRecombCBE.
 *
 * Decreases the reference count of @recomb_cbe.
 *
 */
void
nc_recomb_cbe_free (NcRecombCBE *recomb_cbe)
{
  g_object_unref (recomb_cbe);
}

/**
 * nc_recomb_cbe_clear:
 * @recomb_cbe: a #NcRecombCBE.
 *
 * Decreases the reference count of *@recomb_cbe if
 * *@recomb_cbe is not NULL, then sets *@recomb_cbe to NULL.
 *
 */
void
nc_recomb_cbe_clear (NcRecombCBE **recomb_cbe)
{
  g_clear_object (recomb_cbe);
}

/**
 * nc_recomb_cbe_set_cbe:
 * @recomb_cbe: a #NcRecombCBE
 * @cbe: a #NcCBE
 * 
 * Sets @cbe as the Class backend to be used.
 *
 */
void 
nc_recomb_cbe_set_cbe (NcRecombCBE *recomb_cbe, NcCBE *cbe)
{
  nc_cbe_ref (cbe);
  nc_cbe_clear (&recomb_cbe->cbe);
  recomb_cbe->cbe = cbe;
  nc_cbe_set_thermodyn (cbe, TRUE);
}

/**
 * nc_recomb_cbe_peek_cbe:
 * @recomb_cbe: a #NcRecombCBE.
 *
 * Peeks the currently used #NcCBE.
 * 
 * Returns: (transfer none): the used #NcCBE.
 */
NcCBE *
nc_recomb_cbe_peek_cbe (NcRecombCBE *recomb_cbe)
{
  return recomb_cbe->cbe;
}
