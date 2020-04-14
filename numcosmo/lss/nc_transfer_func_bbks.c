/***************************************************************************
 *            nc_transfer_func_bbks.c
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

/**
 * SECTION:nc_transfer_func_bbks
 * @title: NcTransferFuncBBKS
 * @short_description: Bardeen, Bond, Kaiser and Szalay (BBKS) transfer function.
 * 
 * Bardeen, Bond, Kaiser and Szalay (BBKS) transfer function, see [Sugiyama (1995)][XSugiyama1995].
 * 
 * 
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "lss/nc_transfer_func_bbks.h"
#include "nc_enum_types.h"

struct _NcTransferFuncBBKSPrivate
{
  NcTransferFuncBBKSType type;
  gdouble c1;
  gdouble c2;
  gdouble c3;
  gdouble c4;
  gdouble c5_wm; /* c5_wm = c5/wm */
  gdouble h;
};

enum
{
  PROP_0,
  PROP_TYPE,
};

G_DEFINE_TYPE_WITH_PRIVATE (NcTransferFuncBBKS, nc_transfer_func_bbks, NC_TYPE_TRANSFER_FUNC);

static void
nc_transfer_func_bbks_init (NcTransferFuncBBKS *tf_bbks)
{
  NcTransferFuncBBKSPrivate * const self = tf_bbks->priv = nc_transfer_func_bbks_get_instance_private (tf_bbks);

  self->c1    = 0.0;
  self->c2    = 0.0;
  self->c3    = 0.0;
  self->c4    = 0.0;
  self->c5_wm = 0.0;
  self->h     = 0.0;
}

static void
_nc_transfer_func_bbks_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcTransferFuncBBKS *tf_bbks = NC_TRANSFER_FUNC_BBKS (object);
  g_return_if_fail (NC_IS_TRANSFER_FUNC_BBKS (object));

  switch (prop_id)
  {
    case PROP_TYPE:
      nc_transfer_func_bbks_set_type (tf_bbks, g_value_get_enum (value));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_transfer_func_bbks_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcTransferFuncBBKS *tf_bbks = NC_TRANSFER_FUNC_BBKS (object);
  NcTransferFuncBBKSPrivate * const self = tf_bbks->priv;
  g_return_if_fail (NC_IS_TRANSFER_FUNC_BBKS (object));

  switch (prop_id)
  {
    case PROP_TYPE:
      g_value_set_enum (value, self->type);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_transfer_func_bbks_finalize (GObject *object)
{

  /* Chain up : end */
  G_OBJECT_CLASS (nc_transfer_func_bbks_parent_class)->finalize (object);
}

static void _nc_transfer_func_bbks_prepare (NcTransferFunc *tf, NcHICosmo *cosmo);
static gdouble _nc_transfer_func_bbks_calc (NcTransferFunc *tf, gdouble kh);

static void
nc_transfer_func_bbks_class_init (NcTransferFuncBBKSClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);
  NcTransferFuncClass* parent_class = NC_TRANSFER_FUNC_CLASS (klass);

  object_class->set_property = &_nc_transfer_func_bbks_set_property;
  object_class->get_property = &_nc_transfer_func_bbks_get_property;
  object_class->finalize     = &_nc_transfer_func_bbks_finalize;

  g_object_class_install_property (object_class,
                                   PROP_TYPE,
                                   g_param_spec_enum ("type",
                                                      NULL,
                                                      "BBKS variant type",
                                                      NC_TYPE_TRANSFER_FUNC_BBKS_TYPE, NC_TRANSFER_FUNC_BBKS_TYPE_NOBARYONS,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  
  parent_class->prepare = &_nc_transfer_func_bbks_prepare;
  parent_class->calc    = &_nc_transfer_func_bbks_calc;
}

/**
 * nc_transfer_func_bbks_new:
 *
 * FIXME
 *
 * Returns: A new #NcTransferFunc.
 */
NcTransferFunc *
nc_transfer_func_bbks_new ()
{
  return g_object_new (NC_TYPE_TRANSFER_FUNC_BBKS, NULL);
}

static void
_nc_transfer_func_bbks_prepare (NcTransferFunc *tf, NcHICosmo *cosmo)
{
  NcTransferFuncBBKS *tf_bbks = NC_TRANSFER_FUNC_BBKS (tf);
  NcTransferFuncBBKSPrivate * const self = tf_bbks->priv;
  
  const gdouble T_0 = nc_hicosmo_T_gamma0 (cosmo);
  const gdouble c1  = 3.89;
  const gdouble c2  = gsl_pow_2 (16.1);
  const gdouble c3  = gsl_pow_3 (5.46);
  const gdouble c4  = gsl_pow_4 (6.71);
  const gdouble c5  = gsl_pow_2 (T_0 / 2.7);   /* CMB: (T_0/2.7)^2 = (2.725/2.7)^2 */
  const gdouble h   = nc_hicosmo_h (cosmo);
  const gdouble h2  = h * h;
  const gdouble Ob  = nc_hicosmo_Omega_b0 (cosmo);
  const gdouble Om  = nc_hicosmo_Omega_m0 (cosmo);
  const gdouble wm  = Om * h2;

  self->c1    = c1;
  self->c2    = c2;
  self->c3    = c3;
  self->c4    = c4;
  self->h     = h;

  switch (self->type)
  {
    case NC_TRANSFER_FUNC_BBKS_TYPE_NOBARYONS:
      self->c5_wm = c5 / wm;
      break;
    case NC_TRANSFER_FUNC_BBKS_TYPE_BARYONS:
      self->c5_wm = (c5 / wm) / exp (- Ob - sqrt (2.0 * h) * Ob / Om);
      break;
    case NC_TRANSFER_FUNC_BBKS_TYPE_CCL:
      self->c5_wm = (1.0 / wm) / exp (- Ob - sqrt (2.0 * h) * Ob / Om); /* Check why they modify it like this, is it an typo? */
      break;
    default:
      g_assert_not_reached ();
      break;
  }

  if (c5 == 0.0)
    g_warning ("_nc_transfer_func_bbks_prepare: no radiation universe, BBKS is not defined transfer function will be exaclty = 1.");
}

static gdouble
_nc_transfer_func_bbks_calc (NcTransferFunc *tf, gdouble kh)
{
  NcTransferFuncBBKS *tf_bbks = NC_TRANSFER_FUNC_BBKS (tf);
  NcTransferFuncBBKSPrivate * const self = tf_bbks->priv;
  
  const gdouble k  = kh * self->h;
  const gdouble q  = k * self->c5_wm;
  const gdouble q1 = 2.34 * q;
  const gdouble q2 = q * q;
  const gdouble q3 = q2 * q;
  const gdouble q4 = q3 * q;

  return (q1 == 0.0 ? 1.0 : (log1p (q1) / q1)) * pow (1.0 + self->c1 * q + self->c2 * q2 + self->c3 * q3 + self->c4 * q4, -1.0 / 4.0);
}

/**
 * nc_transfer_func_bbks_set_type:
 * @tf_bbks: a #NcTransferFuncBBKS
 * @bbks_type: a #NcTransferFuncBBKSType
 * 
 * Sets BBKS variant type.
 * 
 */
void 
nc_transfer_func_bbks_set_type (NcTransferFuncBBKS *tf_bbks, NcTransferFuncBBKSType bbks_type)
{
  NcTransferFuncBBKSPrivate * const self = tf_bbks->priv;
  self->type = bbks_type;
}
