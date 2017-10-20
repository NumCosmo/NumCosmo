/***************************************************************************
 *            nc_hipert_first_order.c
 *
 *  Mon October 09 16:58:16 2017
 *  Copyright  2017  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * nc_hipert_first_order.c
 * Copyright (C) 2017 Sandro Dias Pinto Vitenti <sandro@isoftware.com.br>
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
 * SECTION:nc_hipert_first_order
 * @title: NcHIPertFirstOrder
 * @short_description: Base object for implementing first order perturbation in a Friedmann background.
 *
 * FIXME
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "perturbations/nc_hipert_first_order.h"
#include "nc_recomb_seager.h"
#include "nc_enum_types.h"

#ifndef NUMCOSMO_GIR_SCAN
#include "math/rcm.h"
#endif /* NUMCOSMO_GIR_SCAN */

struct _NcHIPertFirstOrderPrivate
{
  NcHIPertGrav *grav;
  GPtrArray *comps;
  NcHIPertBGVar *bg_var;
  NcHIPertCompGauge gauge;
};

enum
{
  PROP_0,
  PROP_GAUGE,
  PROP_GRAV,
  PROP_CARRAY,
  PROP_DIST,
  PROP_RECOMB,
  PROP_SCALEFACTOR,
  PROP_LEN,
};

G_DEFINE_TYPE (NcHIPertFirstOrder, nc_hipert_first_order, NC_TYPE_HIPERT_BOLTZMANN);

static void
nc_hipert_first_order_init (NcHIPertFirstOrder *fo)
{
  fo->priv         = G_TYPE_INSTANCE_GET_PRIVATE (fo, NC_TYPE_HIPERT_FIRST_ORDER, NcHIPertFirstOrderPrivate);
  fo->priv->grav   = NULL;
  fo->priv->comps  = g_ptr_array_new ();
  fo->priv->bg_var = nc_hipert_bg_var_new ();
  fo->priv->gauge  = NC_HIPERT_GRAV_GAUGE_LEN;

  if (FALSE)
  {
    gint *adj;
    gint node_num = 20;
    gint adj_max  = 20 * 19; 
    gint adj_num  = 0;
    gint *adj_row;
    gint bandwidth;
    gint i;
    gint *perm;
    gint *perm_inv;

    adj_row  = g_new0 (gint, node_num + 1);
    adj      = g_new0 (gint, adj_max);
    perm     = g_new0 (gint, node_num);
    perm_inv = g_new0 (gint, node_num);

    //graph_01_adj (, adj_num, adj_row, adj);

    adj_set (node_num, adj_max, &adj_num, adj_row, adj, -1, -1 );
    adj_set (node_num, adj_max, &adj_num, adj_row, adj, 1, 10 );
    adj_set (node_num, adj_max, &adj_num, adj_row, adj, 10, 5 );
    adj_set (node_num, adj_max, &adj_num, adj_row, adj, 10, 4 );
    adj_set (node_num, adj_max, &adj_num, adj_row, adj, 10, 15 );
    adj_set (node_num, adj_max, &adj_num, adj_row, adj, 5, 4 );
    adj_set (node_num, adj_max, &adj_num, adj_row, adj, 5, 7 );
    adj_set (node_num, adj_max, &adj_num, adj_row, adj, 1, 20 );
    adj_set (node_num, adj_max, &adj_num, adj_row, adj, 1, 19 );
    adj_set (node_num, adj_max, &adj_num, adj_row, adj, 1, 18 );

    g_message ("#\n");
    adj_print ( node_num, adj_num, adj_row, adj, "  Adjacency matrix:" );

    adj_show ( node_num, adj_num, adj_row, adj );

    bandwidth = adj_bandwidth ( node_num, adj_num, adj_row, adj );

    g_message ("#    ADJ bandwidth = %d\n#\n", bandwidth);

    genrcm ( node_num, adj_num, adj_row, adj, perm );

    perm_inverse3 ( node_num, perm, perm_inv );

    g_message ("#\n#    The RCM permutation and inverse:\n#\n");

    for ( i = 0; i < node_num; i++ )
    {
      g_message ("#    %8d  %8d  %8d\n", i + 1, perm[i], perm_inv[i]);
    }

    g_message ("#\n#    Permuted adjacency matrix:\n#\n");

    adj_perm_show ( node_num, adj_num, adj_row, adj, perm, perm_inv );

    bandwidth = adj_perm_bandwidth ( node_num, adj_num, adj_row, adj,
                                    perm, perm_inv );

    g_message ("#\n#    ADJ (permuted) bandwidth = %d\n#\n", bandwidth);

    g_free (adj);
    g_free (adj_row);
    g_free (perm);
    g_free (perm_inv);
  }
}

static void
_nc_hipert_first_order_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcHIPertFirstOrder *fo = NC_HIPERT_FIRST_ORDER (object);
  g_return_if_fail (NC_IS_HIPERT_FIRST_ORDER (object));

  switch (prop_id)
  {
    case PROP_GAUGE:
      nc_hipert_first_order_set_gauge (fo, g_value_get_enum (value));
      break;    
    case PROP_GRAV:
      nc_hipert_first_order_set_grav (fo, g_value_get_object (value));
      break;    
    case PROP_CARRAY:
    {
      NcmObjArray *oa = (NcmObjArray *) g_value_get_boxed (value);
      if (oa != NULL)
      {
        guint i;
        for (i = 0; i < oa->len; i++)
        {
          NcHIPertComp *comp = NC_HIPERT_COMP (ncm_obj_array_peek (oa, i));
          nc_hipert_first_order_add_comp (fo, comp);
        }
      }
      break;
    }
    case PROP_DIST:
      nc_hipert_bg_var_set_dist (fo->priv->bg_var, g_value_get_object (value));
      break;    
    case PROP_RECOMB:
      nc_hipert_bg_var_set_recomb (fo->priv->bg_var, g_value_get_object (value));
      break;    
    case PROP_SCALEFACTOR:
      nc_hipert_bg_var_set_scalefactor (fo->priv->bg_var, g_value_get_object (value));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_hipert_first_order_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcHIPertFirstOrder *fo = NC_HIPERT_FIRST_ORDER (object);
  g_return_if_fail (NC_IS_HIPERT_FIRST_ORDER (object));

  switch (prop_id)
  {
    case PROP_GAUGE:
      g_value_set_enum (value, nc_hipert_first_order_get_gauge (fo));
      break;    
    case PROP_GRAV:
      g_value_take_object (value, nc_hipert_first_order_get_grav (fo));
      break;    
    case PROP_CARRAY:
    {
      NcmObjArray *oa = ncm_obj_array_new ();
      guint i;

      for (i = 0; i < fo->priv->comps->len; i++)
      {
        NcHIPertComp *comp = g_ptr_array_index (fo->priv->comps, i);
        ncm_obj_array_add (oa, G_OBJECT (comp));
      }

      g_value_take_boxed (value, oa);
      break;
    }    
    case PROP_DIST:
      g_value_take_object (value, nc_hipert_bg_var_get_dist (fo->priv->bg_var));
      break;    
    case PROP_RECOMB:
      g_value_take_object (value, nc_hipert_bg_var_get_recomb (fo->priv->bg_var));
      break;    
    case PROP_SCALEFACTOR:
      g_value_take_object (value, nc_hipert_bg_var_get_scalefactor (fo->priv->bg_var));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_hipert_first_order_dispose (GObject *object)
{
  NcHIPertFirstOrder *fo = NC_HIPERT_FIRST_ORDER (object);

  nc_hipert_grav_clear (&fo->priv->grav);
  
  if (fo->priv->comps != NULL)
  {
    guint i;
    for (i = 0; i < fo->priv->comps->len; i++)
    {
      NcHIPertComp *comp = NC_HIPERT_COMP (g_ptr_array_index (fo->priv->comps, i));

      nc_hipert_comp_clear (&comp);
      g_ptr_array_index (fo->priv->comps, i) = comp; 
    }
  }  
  g_clear_pointer (&fo->priv->comps, (GDestroyNotify) g_ptr_array_unref);

  nc_hipert_bg_var_clear (&fo->priv->bg_var);
  
  /* Chain up : end */
  G_OBJECT_CLASS (nc_hipert_first_order_parent_class)->dispose (object);
}

static void
_nc_hipert_first_order_finalize (GObject *object)
{

  /* Chain up : end */
  G_OBJECT_CLASS (nc_hipert_first_order_parent_class)->finalize (object);
}

static void
nc_hipert_first_order_class_init (NcHIPertFirstOrderClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);

  g_type_class_add_private (klass, sizeof (NcHIPertFirstOrderPrivate));

  object_class->set_property = &_nc_hipert_first_order_set_property;
  object_class->get_property = &_nc_hipert_first_order_get_property;
  object_class->dispose      = &_nc_hipert_first_order_dispose;
  object_class->finalize     = &_nc_hipert_first_order_finalize;

  g_object_class_install_property (object_class,
                                   PROP_GAUGE,
                                   g_param_spec_enum ("gauge",
                                                      NULL,
                                                      "Gauge",
                                                      NC_TYPE_HIPERT_GRAV_GAUGE, NC_HIPERT_GRAV_GAUGE_SYNCHRONOUS, 
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_GRAV,
                                   g_param_spec_object ("grav",
                                                        NULL,
                                                        "Gravitation object",
                                                        NC_TYPE_HIPERT_GRAV,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_CARRAY,
                                   g_param_spec_boxed ("comp-array",
                                                       NULL,
                                                       "Components array",
                                                       NCM_TYPE_OBJ_ARRAY,
                                                       G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_DIST,
                                   g_param_spec_object ("distance",
                                                        NULL,
                                                        "Distance object",
                                                        NC_TYPE_DISTANCE,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_RECOMB,
                                   g_param_spec_object ("recomb",
                                                        NULL,
                                                        "Recombination object",
                                                        NC_TYPE_RECOMB,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_SCALEFACTOR,
                                   g_param_spec_object ("scalefactor",
                                                        NULL,
                                                        "Scale factor object",
                                                        NC_TYPE_SCALEFACTOR,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
}

/**
 * nc_hipert_first_order_new:
 * 
 * Creates a new #NcHIPertFirstOrder.
 * 
 * Returns: (transfer full): the newly instantiated #NcHIPertFirstOrder.
 */
NcHIPertFirstOrder *
nc_hipert_first_order_new (void)
{
  NcDistance *dist = nc_distance_new (1.0);
  NcRecomb *recomb = NC_RECOMB (nc_recomb_seager_new ());
  NcScalefactor *a = nc_scalefactor_new (0, 1.0, dist);

  NcHIPertFirstOrder *fo = nc_hipert_first_order_new_full (dist, recomb, a);

  nc_distance_free (dist);
  nc_recomb_free (recomb);
  nc_scalefactor_free (a);

  return fo;
}

/**
 * nc_hipert_first_order_new_full:
 * @dist: a #NcDistance
 * @recomb: a #NcRecomb
 * @a: a #NcScalefactor
 * 
 * Creates a new #NcHIPertFirstOrder.
 * 
 * Returns: (transfer full): the newly instantiated #NcHIPertFirstOrder.
 */
NcHIPertFirstOrder *
nc_hipert_first_order_new_full (NcDistance *dist, NcRecomb *recomb, NcScalefactor *a)
{
  NcHIPertFirstOrder *fo = g_object_new (NC_TYPE_HIPERT_FIRST_ORDER,
                                         "distance",    dist,
                                         "recomb",      recomb,
                                         "scalefactor", a,
                                         NULL);
  return fo;
}


/**
 * nc_hipert_first_order_ref:
 * @fo: a #NcHIPertFirstOrder
 *
 * Increases the reference count of @fo.
 *
 * Returns: (transfer full): @fo.
 */
NcHIPertFirstOrder *
nc_hipert_first_order_ref (NcHIPertFirstOrder *fo)
{
  return g_object_ref (fo);
}

/**
 * nc_hipert_first_order_free:
 * @fo: a #NcHIPertFirstOrder
 *
 * Decreases the reference count of @fo.
 *
 */
void 
nc_hipert_first_order_free (NcHIPertFirstOrder *fo)
{
  g_object_unref (fo);
}

/**
 * nc_hipert_first_order_clear:
 * @fo: a #NcHIPertFirstOrder
 *
 * Decreases the reference count of *@fo and sets the pointer *@fo to NULL.
 *
 */
void 
nc_hipert_first_order_clear (NcHIPertFirstOrder **fo)
{
  g_clear_object (fo);
}

static void
_nc_hipert_first_order_prepare_internal (NcHIPertFirstOrder *fo)
{
  if (fo->priv->grav != NULL)
  {
    GArray *grav_drho = NULL, *grav_v = NULL, *grav_dp = NULL, *grav_Pi = NULL;
    guint i;
    
    nc_hipert_comp_get_Tscalar_coupling (NC_HIPERT_COMP (fo->priv->grav), &grav_drho, &grav_v, &grav_dp, &grav_Pi);

    for (i = 0; i < fo->priv->comps->len; i++)
    {
      NcHIPertComp *comp = NC_HIPERT_COMP (g_ptr_array_index (fo->priv->comps, i));
      if (comp != NULL)
      {
        GArray *drho_i = NULL, *v_i = NULL, *dp_i = NULL, *Pi_i = NULL;
        nc_hipert_comp_get_Tscalar_coupling (comp, &drho_i, &v_i, &dp_i, &Pi_i);

        g_clear_pointer (&drho_i, g_array_unref);
        g_clear_pointer (&v_i,    g_array_unref);
        g_clear_pointer (&dp_i,   g_array_unref);
        g_clear_pointer (&Pi_i,   g_array_unref);
      }
    }


    
    g_clear_pointer (&grav_drho, g_array_unref);
    g_clear_pointer (&grav_v,    g_array_unref);
    g_clear_pointer (&grav_dp,   g_array_unref);
    g_clear_pointer (&grav_Pi,   g_array_unref);
  }
}

/**
 * nc_hipert_first_order_set_gauge:
 * @fo: a #NcHIPertFirstOrder
 * @gauge: a #NcHIPertCompGauge
 *
 * Sets the gauge to be used in the first order system.
 *
 */
void 
nc_hipert_first_order_set_gauge (NcHIPertFirstOrder *fo, NcHIPertCompGauge gauge)
{
  if (gauge != fo->priv->gauge)
  {
    guint i;
    if (fo->priv->grav != NULL)
      nc_hipert_comp_set_gauge (NC_HIPERT_COMP (fo->priv->grav), gauge);

    for (i = 0; i < fo->priv->comps->len; i++)
    {
      NcHIPertComp *comp = NC_HIPERT_COMP (g_ptr_array_index (fo->priv->comps, i));
      if (comp != NULL)
        nc_hipert_comp_set_gauge (comp, gauge);
    }

    fo->priv->gauge = gauge;
    _nc_hipert_first_order_prepare_internal (fo);
  }
}

/**
 * nc_hipert_first_order_get_gauge:
 * @fo: a #NcHIPertFirstOrder
 * 
 * Gets the gauge used by @fo.
 * 
 * Returns: the gauge used by @fo
 */
NcHIPertCompGauge 
nc_hipert_first_order_get_gauge (NcHIPertFirstOrder *fo)
{
  return fo->priv->gauge;
}

/**
 * nc_hipert_first_order_set_grav:
 * @fo: a #NcHIPertFirstOrder
 * @grav: a #NcHIPertGrav
 *
 * Sets the gravitation object.
 *
 */
void 
nc_hipert_first_order_set_grav (NcHIPertFirstOrder *fo, NcHIPertGrav *grav)
{
  nc_hipert_grav_clear (&fo->priv->grav);
  if (grav != NULL)
  {
    fo->priv->grav = nc_hipert_grav_ref (grav);
    nc_hipert_comp_set_gauge (NC_HIPERT_COMP (fo->priv->grav), fo->priv->gauge);
    _nc_hipert_first_order_prepare_internal (fo);
  }
}

/**
 * nc_hipert_first_order_get_grav:
 * @fo: a #NcHIPertFirstOrder
 * 
 * Gets the gravitation #NcHIPertGrav object.
 * 
 * Returns: (transfer full) (nullable): the #NcHIPertGrav object used by @fo.
 */
NcHIPertGrav *
nc_hipert_first_order_get_grav (NcHIPertFirstOrder *fo)
{
  return (fo->priv->grav != NULL) ? nc_hipert_grav_ref (fo->priv->grav) : fo->priv->grav;
}

/**
 * nc_hipert_first_order_peek_grav:
 * @fo: a #NcHIPertFirstOrder
 * 
 * Peeks the #NcHIPertGrav object.
 * 
 * Returns: (transfer none) (nullable): the #NcHIPertGrav object used by @fo.
 */
NcHIPertGrav *
nc_hipert_first_order_peek_grav (NcHIPertFirstOrder *fo)
{
  return fo->priv->grav;
}

/**
 * nc_hipert_first_order_add_comp:
 * @fo: a #NcHIPertFirstOrder
 * @comp: a #NcHIPertComp
 *
 * Adds a new component @comp to the system.
 *
 */
void 
nc_hipert_first_order_add_comp (NcHIPertFirstOrder *fo, NcHIPertComp *comp)
{
  const guint len    = nc_hipert_bg_var_len (fo->priv->bg_var);
  NcHIPertBGVarID id = nc_hipert_comp_get_id (comp);

  if (NC_IS_HIPERT_GRAV (comp))
    g_error ("nc_hipert_first_order_add_comp: the gravitation component should be added using `nc_hipert_first_order_set_grav'.");

  g_assert_cmpint (id, >=, 0);
  g_assert_cmpint (id, <, len);
  g_ptr_array_set_size (fo->priv->comps, len);

  if (g_ptr_array_index (fo->priv->comps, id) != NULL)
    g_warning ("nc_hipert_first_order_add_comp: component with `%d' (%s) already included, ignoring...", id, G_OBJECT_TYPE_NAME (comp));
  else
  {
    g_ptr_array_index (fo->priv->comps, id) = nc_hipert_comp_ref (comp);
    nc_hipert_comp_set_gauge (comp, fo->priv->gauge);
    _nc_hipert_first_order_prepare_internal (fo);
  }
}
