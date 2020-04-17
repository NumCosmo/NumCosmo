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

#include <nvector/nvector_serial.h>

#include <cvode/cvode.h>
#include <cvode/cvode_ls.h>
#include <arkode/arkode.h>
#include <arkode/arkode_arkstep.h>
#include <arkode/arkode_ls.h>

#include <sundials/sundials_matrix.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <sunlinsol/sunlinsol_dense.h>
#include <sunmatrix/sunmatrix_band.h>
#include <sunlinsol/sunlinsol_band.h>

#include <sundials/sundials_types.h> 
#endif /* NUMCOSMO_GIR_SCAN */

typedef struct _NcHIPertFirstOrderVar
{
  gint src;
  gint index;
  GArray *deps;
} NcHIPertFirstOrderVar;

struct _NcHIPertFirstOrderPrivate
{
  NcHIPertGrav *grav;
  GPtrArray *comps;
  GPtrArray *active_comps;
  GArray *vars;
  GArray *perm;
  GArray *perm_inv;
  NcHIPertBGVar *bg_var;
  NcHIPertGravGauge gauge;
  gpointer cvode;
  gboolean cvode_init;
  gpointer arkode;
  gboolean arkode_init;
  guint cur_sys_size;
  N_Vector y;
  N_Vector abstol_v;
  SUNMatrix A;
  SUNLinearSolver LS;
  gdouble reltol;
  gdouble abstol;
  gint mupper;
  gint mlower;
  NcHIPertFirstOrderInteg integ;
  NcHIPertGravTScalar *T_scalar_i;
  NcHIPertGravTScalar *T_scalar_tot;
  NcHIPertGravScalar *G_scalar;
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
  PROP_RELTOL,
  PROP_ABSTOL,
  PROP_INTEG,
  PROP_LEN,
};

G_DEFINE_TYPE_WITH_PRIVATE (NcHIPertFirstOrder, nc_hipert_first_order, NC_TYPE_HIPERT_BOLTZMANN);

void 
_nc_hipert_first_order_clear_var (NcHIPertFirstOrderVar *var)
{
  g_clear_pointer (&var->deps, g_array_unref);
}

static void
nc_hipert_first_order_init (NcHIPertFirstOrder *fo)
{
  NcHIPertFirstOrderPrivate * const self = fo->priv = nc_hipert_first_order_get_instance_private (fo); 
  
  self->grav         = NULL;
  self->comps        = g_ptr_array_new ();
  self->active_comps = g_ptr_array_new ();
  self->vars         = g_array_new (TRUE, TRUE, sizeof (NcHIPertFirstOrderVar));
  self->perm         = g_array_new (FALSE, FALSE, sizeof (gint));
  self->perm_inv     = g_array_new (FALSE, FALSE, sizeof (gint));
  self->bg_var       = nc_hipert_bg_var_new ();
  self->gauge        = NC_HIPERT_GRAV_GAUGE_LEN;

  self->cvode        = NULL;
  self->cvode_init   = FALSE;
  self->arkode       = NULL;
  self->arkode_init  = FALSE;

  self->reltol       = 0.0;
  self->abstol       = 0.0;

  self->mupper       = 0;
  self->mlower       = 0;

  self->integ        = NC_HIPERT_FIRST_ORDER_INTEG_LEN;
  self->cur_sys_size = 0;
  self->y            = NULL;
  self->abstol_v     = NULL;

  self->A            = NULL;
  self->LS           = NULL;

  self->T_scalar_i   = nc_hipert_grav_T_scalar_new ();
  self->T_scalar_tot = nc_hipert_grav_T_scalar_new ();
  self->G_scalar     = nc_hipert_grav_scalar_new ();

  g_array_set_clear_func (self->vars, (GDestroyNotify)_nc_hipert_first_order_clear_var);  
  g_ptr_array_set_free_func (self->active_comps, (GDestroyNotify)nc_hipert_comp_free);
}

static void
_nc_hipert_first_order_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcHIPertFirstOrder *fo = NC_HIPERT_FIRST_ORDER (object);
  NcHIPertFirstOrderPrivate * const self = fo->priv;
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
      nc_hipert_bg_var_set_dist (self->bg_var, g_value_get_object (value));
      break;    
    case PROP_RECOMB:
      nc_hipert_bg_var_set_recomb (self->bg_var, g_value_get_object (value));
      break;    
    case PROP_SCALEFACTOR:
      nc_hipert_bg_var_set_scalefactor (self->bg_var, g_value_get_object (value));
      break;
    case PROP_RELTOL:
      nc_hipert_first_order_set_reltol (fo, g_value_get_double (value));
      break;
    case PROP_ABSTOL:
      nc_hipert_first_order_set_abstol (fo, g_value_get_double (value));
      break;
    case PROP_INTEG:
      nc_hipert_first_order_set_integ (fo, g_value_get_enum (value));
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
  NcHIPertFirstOrderPrivate * const self = fo->priv;
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

      for (i = 0; i < self->comps->len; i++)
      {
        NcHIPertComp *comp = g_ptr_array_index (self->comps, i);
        if (comp != NULL)
          ncm_obj_array_add (oa, G_OBJECT (comp));
      }

      g_value_take_boxed (value, oa);
      break;
    }    
    case PROP_DIST:
      g_value_take_object (value, nc_hipert_bg_var_get_dist (self->bg_var));
      break;    
    case PROP_RECOMB:
      g_value_take_object (value, nc_hipert_bg_var_get_recomb (self->bg_var));
      break;    
    case PROP_SCALEFACTOR:
      g_value_take_object (value, nc_hipert_bg_var_get_scalefactor (self->bg_var));
      break;
    case PROP_RELTOL:
      g_value_set_double (value, nc_hipert_first_order_get_reltol (fo));
      break;
    case PROP_ABSTOL:
      g_value_set_double (value, nc_hipert_first_order_get_abstol (fo));
      break;
    case PROP_INTEG:
      g_value_set_enum (value, nc_hipert_first_order_get_integ (fo));
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
  NcHIPertFirstOrderPrivate * const self = fo->priv;
  
  nc_hipert_grav_clear (&self->grav);
  
  if (self->comps != NULL)
  {
    guint i;
    for (i = 0; i < self->comps->len; i++)
    {
      NcHIPertComp *comp = NC_HIPERT_COMP (g_ptr_array_index (self->comps, i));

      nc_hipert_comp_clear (&comp);
      g_ptr_array_index (self->comps, i) = comp; 
    }
  }

  g_clear_pointer (&self->comps,        g_ptr_array_unref);
  g_clear_pointer (&self->active_comps, g_ptr_array_unref);

  g_clear_pointer (&self->vars,         g_array_unref);
  g_clear_pointer (&self->perm,         g_array_unref);
  g_clear_pointer (&self->perm_inv,     g_array_unref);

  nc_hipert_bg_var_clear (&self->bg_var);
  
  /* Chain up : end */
  G_OBJECT_CLASS (nc_hipert_first_order_parent_class)->dispose (object);
}

static void
_nc_hipert_first_order_finalize (GObject *object)
{
  NcHIPertFirstOrder *fo = NC_HIPERT_FIRST_ORDER (object);
  NcHIPertFirstOrderPrivate * const self = fo->priv;

  if (self->cvode != NULL)
  {
    CVodeFree (&self->cvode);
    self->cvode      = NULL;
    self->cvode_init = FALSE;
  }
  if (self->arkode != NULL)
  {
    ARKStepFree (&self->arkode);
    self->arkode      = NULL;
    self->arkode_init = FALSE;
  }

  if (self->A != NULL)
    SUNMatDestroy (self->A);

  if (self->LS != NULL)
  {
    gint flag = SUNLinSolFree (self->LS);
    NCM_CVODE_CHECK (&flag, "SUNLinSolFree", 1, );
  }

  g_clear_pointer (&self->y, N_VDestroy);
  g_clear_pointer (&self->abstol_v, N_VDestroy);

  self->cur_sys_size = 0;

  g_clear_pointer (&self->T_scalar_i,   nc_hipert_grav_T_scalar_free);
  g_clear_pointer (&self->T_scalar_tot, nc_hipert_grav_T_scalar_free);

  g_clear_pointer (&self->G_scalar,     nc_hipert_grav_scalar_free);
  
  /* Chain up : end */
  G_OBJECT_CLASS (nc_hipert_first_order_parent_class)->finalize (object);
}

static void
nc_hipert_first_order_class_init (NcHIPertFirstOrderClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);

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
  g_object_class_install_property (object_class,
                                   PROP_RELTOL,
                                   g_param_spec_double ("reltol",
                                                        NULL,
                                                        "Relative tolerance",
                                                        0.0, 1.0, NCM_DEFAULT_PRECISION,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_ABSTOL,
                                   g_param_spec_double ("abstol",
                                                        NULL,
                                                        "Absolute tolerance tolerance",
                                                        0.0, G_MAXDOUBLE, 0.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_INTEG,
                                   g_param_spec_enum ("integ",
                                                        NULL,
                                                        "ODE integrator",
                                                        NC_TYPE_HIPERT_FIRST_ORDER_INTEG,
                                                        NC_HIPERT_FIRST_ORDER_INTEG_ARKODE,
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
  NcScalefactor *a = nc_scalefactor_new (1.0, dist);

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
_nc_hipert_first_order_add_pad (GArray *a, gint pad)
{
  guint i;
  if (pad == 0)
    return;
  
  for (i = 0; i < a->len; i++)
  {
    if (g_array_index (a, gint, i) >= 0)
    {
      g_array_index (a, gint, i) += pad;
    }
  }
}

#define APPEND(a,b) (g_array_append_vals ((a), (b)->data, (b)->len))

gint __cmp_gint (gconstpointer a, gconstpointer b) { const gint *u = a; const gint *v = b; return (u[0] < v[0]) ? -1 : ((u[0] > v[0]) ? 1 : 0); }

static void
_nc_hipert_first_order_solve_deps (NcHIPertFirstOrder *fo, NcHIPertGravInfo *ginfo, NcHIPertGravTScalarInfo *Tsinfo, GArray *deps, guint r)
{
  gboolean subs = FALSE;
  guint i;

  if (r > 9)
    g_error ("_nc_hipert_first_order_solve_deps: too many recursion levels.");
  
  for (i = 0; i < deps->len; )
  {
    const gint v = g_array_index (deps, gint, i);
    /*printf ("%d %d %d\n", i, deps->len, v);*/
    
    if (v < 0)
    {
      g_array_remove_index (deps, i);
      subs = TRUE;
      switch (v)
      {
        case NC_HIPERT_GRAV_SELEM_PHI:
          /*g_message ("Appending phi     %d!\n", ginfo->phi_deps->len);*/
          APPEND (deps, ginfo->phi_deps);
          break;
        case NC_HIPERT_GRAV_SELEM_DSIGMA:
/*          g_message ("Appending dsigma  %d!\n", ginfo->dsigma_deps->len);*/
          APPEND (deps, ginfo->dsigma_deps);
          break;
        case NC_HIPERT_GRAV_SELEM_PSI:
          /*g_message ("Appending psi     %d!\n", ginfo->psi_deps->len);*/
          APPEND (deps, ginfo->psi_deps);
          break;
        case NC_HIPERT_GRAV_SELEM_DOTPSI:
          /*g_message ("Appending dotphi  %d!\n", ginfo->dotpsi_deps->len);*/
          APPEND (deps, ginfo->dotpsi_deps);
          break;
        case NC_HIPERT_GRAV_SELEM_DRHO:
          /*g_message ("Appending drho    %d!\n", Tsinfo->drho_deps->len);*/
          APPEND (deps, Tsinfo->drho_deps);
          break;
        case NC_HIPERT_GRAV_SELEM_RHOPPV:
          /*g_message ("Appending rhoppv  %d!\n", Tsinfo->rhoppv_deps->len);*/
          APPEND (deps, Tsinfo->rhoppv_deps);
          break;
        case NC_HIPERT_GRAV_SELEM_DP:
          /*g_message ("Appending dp      %d!\n", Tsinfo->dp_deps->len);*/
          APPEND (deps, Tsinfo->dp_deps);
          break;
        case NC_HIPERT_GRAV_SELEM_DPI:
          /*g_message ("Appending dPi     %d!\n", Tsinfo->dPi_deps->len);*/
          APPEND (deps, Tsinfo->dPi_deps);
          break;
        default:
          g_assert_not_reached ();
          break;
      }
    }
    else
      i++;
  }

  if (subs)
  {
    _nc_hipert_first_order_solve_deps (fo, ginfo, Tsinfo, deps, r + 1);
  }
  else
  {
    g_array_sort (deps, &__cmp_gint);
    if (deps->len > 1)
    {
      gint last = g_array_index (deps, gint, 0);
      for (i = 1; i < deps->len; )
      {
        gint v = g_array_index (deps, gint, i);
        if (v == last)
        {
          g_array_remove_index (deps, i);
        }
        else
        {
          last = v;
          i++;
        }
      }
    }
  }
}

static void
_nc_hipert_first_order_arrange_vars (NcHIPertFirstOrder *fo)
{
  NcHIPertFirstOrderPrivate * const self = fo->priv;
  gint *adj;
  gint node_num = self->vars->len;
  gchar *Jrow   = g_new0 (gchar, node_num + 1);
  gint adj_max  = node_num * (node_num - 1); 
  gint adj_num  = 0;
  gint *adj_row;
  gint *perm;
  gint *perm_inv;
  gint orig_mupper;
  gint orig_mlower;
  gint i;

  adj_row  = g_new0 (gint, node_num + 1);
  adj      = g_new0 (gint, adj_max);

  g_array_set_size (self->perm,     node_num);
  g_array_set_size (self->perm_inv, node_num);
  
  perm     = &g_array_index (self->perm,     gint, 0);
  perm_inv = &g_array_index (self->perm_inv, gint, 0);

  if (FALSE)
  {
    gint lll = self->vars->len - 1;
    g_array_append_val (g_array_index (self->vars, NcHIPertFirstOrderVar, 0).deps, lll);
    lll--;
    if (lll >= 0)
      g_array_append_val (g_array_index (self->vars, NcHIPertFirstOrderVar, 1).deps, lll);
  }
  
  adj_set (node_num, adj_max, &adj_num, adj_row, adj, -1, -1 );

  ncm_message ("#\n# Original jacobian:\n#\n");
  for (i = 0; i < node_num; i++)
  {
    NcHIPertFirstOrderVar var = g_array_index (self->vars, NcHIPertFirstOrderVar, i);
    gint j;

    for (j = 0; j < node_num; j++)
    {
      if (i == j)
        Jrow[j] = 'D'; 
      else
        Jrow[j] = '.';
    }
    Jrow[j] = '\0';
    
    for (j = 0; j < var.deps->len; j++)
    {
      gint dep = g_array_index (var.deps, gint, j);
      adj_set (node_num, adj_max, &adj_num, adj_row, adj, var.index + 1, dep + 1);
      if (var.index != dep)
        Jrow[dep] = 'X';
      /*printf ("%d %d %d %d\n", var.index, dep, i, j);*/
    }
    ncm_message ("#  %s\n", Jrow);
  }
/*  
  g_message ("#\n");
  adj_print ( node_num, adj_num, adj_row, adj, "  Adjacency matrix:" );

  adj_show ( node_num, adj_num, adj_row, adj );

  bandwidth = adj_bandwidth ( node_num, adj_num, adj_row, adj );

  g_message ("#    ADJ bandwidth = %d\n#\n", bandwidth);
*/
  genrcm ( node_num, adj_num, adj_row, adj, perm ); /* REMEMBER: perm[new_index]     == old_index !!!! */
  perm_inverse3 ( node_num, perm, perm_inv );       /* REMEMBER: perm_inv[old_index] == new_index !!!! */
  
  /*g_message ("#\n#    The RCM permutation and inverse:\n#\n");*/

  orig_mupper = 0;
  orig_mlower = 0;

  for (i = 0; i < node_num; i++)
  {
    NcHIPertFirstOrderVar *var = &g_array_index (self->vars, NcHIPertFirstOrderVar, i);
    var->index = perm[i] - 1;
    g_assert_cmpint (perm[perm_inv[i]-1], ==, i + 1);
    g_assert_cmpint (perm_inv[perm[i]-1], ==, i + 1);
    /*g_message ("#    %8d  %8d  %8d | %8d  %8d\n", i + 1, perm[i], perm_inv[i], perm[perm_inv[i]-1], perm_inv[perm[i]-1]);*/
    {
      NcHIPertFirstOrderVar var = g_array_index (self->vars, NcHIPertFirstOrderVar, i);
      gint j;
      for (j = 0; j < var.deps->len; j++)
      {
        gint dep = g_array_index (var.deps, gint, j);

        orig_mupper = MAX (orig_mupper, dep - i);
        orig_mlower = MAX (orig_mlower, i - dep);
      }    
    }
  }
  g_message ("#\n#  ADJ (non-permuted) bandwidth = (%d, %d)\n", orig_mupper, orig_mlower);
  
/*
  g_message ("#\n#    Permuted adjacency matrix:\n#\n");

  adj_perm_show ( node_num, adj_num, adj_row, adj, perm, perm_inv );
*/
  
  self->mupper = 0;
  self->mlower = 0;
    
  ncm_message ("#\n# Reordered jacobian:\n#\n");
  for (i = 0; i < node_num; i++)
  {
    NcHIPertFirstOrderVar var = g_array_index (self->vars, NcHIPertFirstOrderVar, perm[i] - 1);
    gint j;

    for (j = 0; j < node_num; j++)
    {
      if (i == j)
        Jrow[j] = 'D'; 
      else
        Jrow[j] = '.';
    }
    Jrow[j] = '\0';

    for (j = 0; j < var.deps->len; j++)
    {
      gint dep = perm_inv[g_array_index (var.deps, gint, j)] - 1;

      self->mupper = MAX (self->mupper, dep - i);
      self->mlower = MAX (self->mlower, i - dep);

      if (i != dep)
        Jrow[dep] = 'X';
    }
    ncm_message ("#  %s\n", Jrow);
  }

  g_free (Jrow);
  g_message ("#\n#  ADJ (permuted) bandwidth = (%d, %d)\n", self->mupper, self->mlower);

  if ((orig_mupper + orig_mlower) <= (self->mupper + self->mlower))
  {
    for (i = 0; i < node_num; i++)
    {
      g_array_index (self->perm,     gint, i) = i;
      g_array_index (self->perm_inv, gint, i) = i;
    }    
  }

  g_free (adj);
  g_free (adj_row);
}

static void
_nc_hipert_first_order_prepare_internal (NcHIPertFirstOrder *fo)
{
  NcHIPertFirstOrderPrivate * const self = fo->priv;
  if (self->grav != NULL)
  {
    NcHIPertGravInfo *ginfo         = nc_hipert_grav_get_G_scalar_info (self->grav);
    NcHIPertGravTScalarInfo *Tsinfo = nc_hipert_grav_T_scalar_info_new ();
    const guint grav_ndyn           = nc_hipert_grav_ndyn_var (self->grav);

    guint i, pad = 0;

    g_array_set_size (self->vars, 0);
      
    /* Adding gravitation potentials to the variables list */
    for (i = 0; i < grav_ndyn; i++)
    {
      GArray *grav_dyn_var_i_deps = nc_hipert_grav_get_deps (self->grav, i);
      NcHIPertFirstOrderVar var = {-1, self->vars->len, grav_dyn_var_i_deps};

      _nc_hipert_first_order_add_pad (grav_dyn_var_i_deps, pad);

      g_array_append_val (self->vars, var);
    }
    pad = self->vars->len;

    for (i = 0; i < self->comps->len; i++)
    {
      NcHIPertComp *comp = NC_HIPERT_COMP (g_ptr_array_index (self->comps, i));
      if (comp != NULL)
      {
        NcHIPertGravTScalarInfo *Tsinfo_i = nc_hipert_comp_get_T_scalar_info (comp);
        guint ndyn = nc_hipert_comp_ndyn_var (comp);
        guint j;

        nc_hipert_grav_T_scalar_info_add_pad (Tsinfo_i, pad);
        nc_hipert_grav_T_scalar_info_append (Tsinfo, Tsinfo_i);
        nc_hipert_grav_T_scalar_info_free (Tsinfo_i);

        for (j = 0; j < ndyn; j++)
        {
          GArray *comp_j_deps       = nc_hipert_comp_get_deps (comp, j);
          NcHIPertFirstOrderVar var = {i, self->vars->len, comp_j_deps};

          _nc_hipert_first_order_add_pad (comp_j_deps, pad);
          
          g_array_append_val (self->vars, var);
        }
        pad = self->vars->len;
      }
    }

    if (TRUE)
    {
      guint i;

      ncm_cfg_msg_sepa ();
      g_message ("# phi deps:    ");
      for (i = 0; i < ginfo->phi_deps->len; i++)
      {
        g_message (" %2d", g_array_index (ginfo->phi_deps, gint, i));
      }
      g_message ("\n");
      g_message ("# dsigma deps: ");
      for (i = 0; i < ginfo->dsigma_deps->len; i++)
      {
        g_message (" %2d", g_array_index (ginfo->dsigma_deps, gint, i));
      }
      g_message ("\n");
      g_message ("# psi deps:    ");
      for (i = 0; i < ginfo->psi_deps->len; i++)
      {
        g_message (" %2d", g_array_index (ginfo->psi_deps, gint, i));
      }
      g_message ("\n");
      g_message ("# dotpsi deps: ");
      for (i = 0; i < ginfo->dotpsi_deps->len; i++)
      {
        g_message (" %2d", g_array_index (ginfo->dotpsi_deps, gint, i));
      }
      g_message ("\n");


      g_message ("# drho deps:   ");
      for (i = 0; i < Tsinfo->drho_deps->len; i++)
      {
        g_message (" %2d", g_array_index (Tsinfo->drho_deps, gint, i));
      }
      g_message ("\n");
      g_message ("# rhoppv deps: ");
      for (i = 0; i < Tsinfo->rhoppv_deps->len; i++)
      {
        g_message (" %2d", g_array_index (Tsinfo->rhoppv_deps, gint, i));
      }
      g_message ("\n");
      g_message ("# dp deps:     ");
      for (i = 0; i < Tsinfo->dp_deps->len; i++)
      {
        g_message (" %2d", g_array_index (Tsinfo->dp_deps, gint, i));
      }
      g_message ("\n");
      g_message ("# dPi deps:    ");
      for (i = 0; i < Tsinfo->dPi_deps->len; i++)
      {
        g_message (" %2d", g_array_index (Tsinfo->dPi_deps, gint, i));
      }
      g_message ("\n");
    }
    
    {
      guint i;
      for (i = 0; i < self->vars->len; i++)
      {
        NcHIPertFirstOrderVar var = g_array_index (self->vars, NcHIPertFirstOrderVar, i);
        _nc_hipert_first_order_solve_deps (fo, ginfo, Tsinfo, var.deps, 0);
      }      
    }
    
    _nc_hipert_first_order_arrange_vars (fo);
      
    nc_hipert_grav_T_scalar_info_free (Tsinfo);
    nc_hipert_grav_info_free (ginfo);
  }
}

/**
 * nc_hipert_first_order_set_gauge:
 * @fo: a #NcHIPertFirstOrder
 * @gauge: a #NcHIPertGravGauge
 *
 * Sets the gauge to be used in the first order system.
 *
 */
void 
nc_hipert_first_order_set_gauge (NcHIPertFirstOrder *fo, NcHIPertGravGauge gauge)
{
  NcHIPertFirstOrderPrivate * const self = fo->priv;
  if (gauge != self->gauge)
  {
    guint i;
    if (self->grav != NULL)
      nc_hipert_grav_set_gauge (self->grav, gauge);

    for (i = 0; i < self->comps->len; i++)
    {
      NcHIPertComp *comp = NC_HIPERT_COMP (g_ptr_array_index (self->comps, i));
      if (comp != NULL)
        nc_hipert_comp_set_gauge (comp, gauge);
    }

    self->gauge = gauge;
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
NcHIPertGravGauge 
nc_hipert_first_order_get_gauge (NcHIPertFirstOrder *fo)
{
  NcHIPertFirstOrderPrivate * const self = fo->priv;
  return self->gauge;
}

/**
 * nc_hipert_first_order_set_reltol:
 * @fo: a #NcHIPertFirstOrder
 * @reltol: relative tolerance used during the integration
 * 
 * Sets the relative tolerance to @reltol.
 * 
 */
void 
nc_hipert_first_order_set_reltol (NcHIPertFirstOrder *fo, const gdouble reltol)
{
  NcHIPertFirstOrderPrivate * const self = fo->priv;
  self->reltol = reltol;
}

/**
 * nc_hipert_first_order_set_abstol:
 * @fo: a #NcHIPertFirstOrder
 * @abstol: absolute tolerance used during the integration
 * 
 * Sets the absolute tolerance to @abstol.
 * 
 */
void 
nc_hipert_first_order_set_abstol (NcHIPertFirstOrder *fo, const gdouble abstol)
{
  NcHIPertFirstOrderPrivate * const self = fo->priv;
  self->abstol = abstol;
}

/**
 * nc_hipert_first_order_get_reltol:
 * @fo: a #NcHIPertFirstOrder
 *
 * Gets the relative tolerance.
 * 
 * Returns: the current relative tolerance.
 */
gdouble 
nc_hipert_first_order_get_reltol (NcHIPertFirstOrder *fo)
{
  NcHIPertFirstOrderPrivate * const self = fo->priv;
  return self->reltol;
}

/**
 * nc_hipert_first_order_get_abstol:
 * @fo: a #NcHIPertFirstOrder
 *
 * Gets the absolute tolerance.
 * 
 * Returns: the current absolute tolerance.
 */
gdouble 
nc_hipert_first_order_get_abstol (NcHIPertFirstOrder *fo)
{
  NcHIPertFirstOrderPrivate * const self = fo->priv;
  return self->abstol;
}

/**
 * nc_hipert_first_order_set_integ:
 * @fo: a #NcHIPertFirstOrder
 * @integ: integrator type #NcHIPertFirstOrderInteg
 *
 * Sets the integrator to be used.
 * 
 */
void
nc_hipert_first_order_set_integ (NcHIPertFirstOrder *fo, NcHIPertFirstOrderInteg integ)
{
  NcHIPertFirstOrderPrivate * const self = fo->priv;
  if (self->integ != integ)
  {
    self->integ = integ;
  }
}

/**
 * nc_hipert_first_order_get_integ:
 * @fo: a #NcHIPertFirstOrder
 *
 * Gets the integrator used.
 * 
 * Returns: the current integrator used by @fo.
 */
NcHIPertFirstOrderInteg 
nc_hipert_first_order_get_integ (NcHIPertFirstOrder *fo)
{
  NcHIPertFirstOrderPrivate * const self = fo->priv;
  return self->integ;
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
  NcHIPertFirstOrderPrivate * const self = fo->priv;
  nc_hipert_grav_clear (&self->grav);
  if (grav != NULL)
  {
    self->grav = nc_hipert_grav_ref (grav);
    nc_hipert_grav_set_gauge (self->grav, self->gauge);
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
  NcHIPertFirstOrderPrivate * const self = fo->priv;
  return (self->grav != NULL) ? nc_hipert_grav_ref (self->grav) : self->grav;
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
  NcHIPertFirstOrderPrivate * const self = fo->priv;
  return self->grav;
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
  NcHIPertFirstOrderPrivate * const self = fo->priv;
  const guint len    = nc_hipert_bg_var_len (self->bg_var);
  NcHIPertBGVarID id = nc_hipert_comp_get_id (comp);

  if (NC_IS_HIPERT_GRAV (comp))
    g_error ("nc_hipert_first_order_add_comp: the gravitation component should be added using `nc_hipert_first_order_set_grav'.");

  g_assert_cmpint (id, >=, 0);
  g_assert_cmpint (id, <, len);
  g_ptr_array_set_size (self->comps, len);

  if (g_ptr_array_index (self->comps, id) != NULL)
    g_warning ("nc_hipert_first_order_add_comp: component with `%d' (%s) already included, ignoring...", id, G_OBJECT_TYPE_NAME (comp));
  else
  {
    g_ptr_array_index (self->comps, id) = nc_hipert_comp_ref (comp);

    g_ptr_array_add (self->active_comps, nc_hipert_comp_ref (comp));
      
    nc_hipert_comp_set_gauge (comp, self->gauge);
    _nc_hipert_first_order_prepare_internal (fo);
  }
}

typedef struct _NcHIPertFirstOrderWS
{
  NcHIPertFirstOrder *fo;
  NcHIPertBGVarYDY *ydy;
  NcHICosmo *cosmo;
} NcHIPertFirstOrderWS;

static gint
_nc_hipert_first_order_f (realtype t, N_Vector y, N_Vector ydot, gpointer f_data)
{
  NcHIPertFirstOrderWS *ws = (NcHIPertFirstOrderWS *) f_data;
  NcHIPertFirstOrder *fo = ws->fo;
  NcHIPertFirstOrderPrivate * const self = fo->priv;
  NcHICosmo *cosmo       = ws->cosmo;
  NcHIPertBGVarYDY *ydy  = ws->ydy;
  NcHIPertBGVar *bg_var  = self->bg_var;
  const guint ncomps     = self->active_comps->len;
  guint i;

  /*printf ("Getting background variables at % 22.15g | %p %p %p\n", t, ydy, y, ydot);fflush (stdout);*/
  nc_hicosmo_get_bg_var (cosmo, t, bg_var);

  ydy->y  = N_VGetArrayPointer (y);
  ydy->dy = N_VGetArrayPointer (ydot);

  N_VConst (0.0, ydot);

  nc_hipert_grav_T_scalar_set_zero (self->T_scalar_tot);
  
  for (i = 0; i < ncomps; i++)
  {
    NcHIPertComp *comp = g_ptr_array_index (self->active_comps, i);

    nc_hipert_grav_T_scalar_set_zero (self->T_scalar_i);

    /*printf ("Getting T_scalar at % 22.15g for comp %d\n", t, i);fflush (stdout);*/
    nc_hipert_comp_get_T_scalar (comp, bg_var, ydy, self->T_scalar_i);
    nc_hipert_grav_T_scalar_add (self->T_scalar_tot, self->T_scalar_tot, self->T_scalar_i);
  }

  nc_hipert_grav_scalar_set_zero (self->G_scalar);
  nc_hipert_grav_get_G_scalar (self->grav, bg_var, ydy, self->T_scalar_tot, self->G_scalar);

  nc_hipert_grav_get_dy_scalar (self->grav, bg_var, ydy, self->T_scalar_tot, self->G_scalar);
  
  for (i = 0; i < ncomps; i++)
  {
    NcHIPertComp *comp = g_ptr_array_index (self->active_comps, i);

    nc_hipert_comp_get_dy_scalar (comp, bg_var, ydy, self->T_scalar_tot, self->G_scalar);
  }

  return 0;
}

static gdouble
_nc_hipert_first_order_set_init_cond (NcHIPertFirstOrder *fo, NcHICosmo *cosmo, const gdouble k)
{
  NcHIPertFirstOrderPrivate * const self = fo->priv;
  gint i;

  for (i = 0; i < self->cur_sys_size; i++)
  {
    NV_Ith_S (self->y, i) = 1.0;
  }

  return 0.0;
}

static void
_nc_hipert_first_order_alloc_integrator (NcHIPertFirstOrder *fo)
{
  NcHIPertFirstOrderPrivate * const self = fo->priv;
  gint flag;

  if (self->cur_sys_size != self->vars->len)
  {
    if (self->cvode != NULL)
    {
      CVodeFree (&self->cvode);
      self->cvode      = NULL;
      self->cvode_init = FALSE;
    }
    if (self->arkode != NULL)
    {
      ARKStepFree (&self->arkode);
      self->arkode      = NULL;
      self->arkode_init = FALSE;
    }

    g_clear_pointer (&self->y, N_VDestroy);
    g_clear_pointer (&self->abstol_v, N_VDestroy);

    self->cur_sys_size = self->vars->len;

    self->y        = N_VNew_Serial (self->cur_sys_size);
    self->abstol_v = N_VNew_Serial (self->cur_sys_size);

    if (self->A != NULL)
      SUNMatDestroy (self->A);

    if (self->LS != NULL)
    {
      flag = SUNLinSolFree (self->LS);
      NCM_CVODE_CHECK (&flag, "SUNLinSolFree", 1, );
    }

    self->A = SUNBandMatrix (self->cur_sys_size, self->mupper, self->mlower);
    NCM_CVODE_CHECK ((gpointer)self->A, "SUNBandMatrix", 0, );

    self->LS = SUNBandLinearSolver (self->y, self->A);
    NCM_CVODE_CHECK ((gpointer)self->LS, "SUNDenseLinearSolver", 0, );
  }
}

static void
_nc_hipert_first_order_prepare_integrator (NcHIPertFirstOrder *fo, const gdouble t0)
{
  NcHIPertFirstOrderPrivate * const self = fo->priv;
  gint flag;

  N_VConst (self->abstol, self->abstol_v);

  switch (self->integ)
  {
    case NC_HIPERT_FIRST_ORDER_INTEG_CVODE:
    {
      if (self->cvode_init)
      {
        self->cvode = CVodeCreate (CV_BDF); /*CVodeCreate (CV_BDF, CV_NEWTON);*/
      
        flag = CVodeInit (self->cvode, &_nc_hipert_first_order_f, t0, self->y);
        NCM_CVODE_CHECK (&flag, "CVodeInit", 1, );

        flag = CVodeSVtolerances (self->cvode, self->reltol, self->abstol_v);
        NCM_CVODE_CHECK (&flag, "CVodeSVtolerances", 1, );

        flag = CVodeSetMaxNumSteps (self->cvode, 0);
        NCM_CVODE_CHECK (&flag, "CVodeSetMaxNumSteps", 1, );

        flag = CVodeSetLinearSolver (self->cvode, self->LS, self->A);
        NCM_CVODE_CHECK (&flag, "CVDlsSetLinearSolver", 1, );
        
        //flag = CVodeSetJacFn (self->cvode, NULL /*J*/);
        //NCM_CVODE_CHECK (&flag, "CVodeSetJacFn", 1, );

        flag = CVodeSetInitStep (self->cvode, fabs (t0) * self->reltol);
        NCM_CVODE_CHECK (&flag, "CVodeSetInitStep", 1, );

        self->cvode_init = TRUE;
      }
      else
      {    
        flag = CVodeReInit (self->cvode, t0, self->y);
        NCM_CVODE_CHECK (&flag, "CVodeInit", 1, );

        flag = CVodeSetInitStep (self->cvode, fabs (t0) * self->reltol);
        NCM_CVODE_CHECK (&flag, "CVodeSetInitStep", 1, );
      }  
      break;
    }
    case NC_HIPERT_FIRST_ORDER_INTEG_ARKODE:
    {
#define INTTYPE _nc_hipert_first_order_f, NULL
      if (!self->arkode_init)
      {
        self->arkode = ARKStepCreate (INTTYPE, t0, self->y);
        NCM_CVODE_CHECK (&self->arkode, "ARKStepCreate", 0, );

        flag = ARKStepSVtolerances (self->arkode, self->reltol, self->abstol_v);
        NCM_CVODE_CHECK (&flag, "ARKStepSVtolerances", 1, );

        flag = ARKStepSetMaxNumSteps (self->arkode, 0);
        NCM_CVODE_CHECK (&flag, "ARKStepSetMaxNumSteps", 1, );

        flag = ARKStepSetLinearSolver (self->arkode, self->LS, self->A);
        NCM_CVODE_CHECK (&flag, "ARKStepSetLinearSolver", 1, );
        
        flag = ARKStepSetLinear (self->arkode, 1);
        NCM_CVODE_CHECK (&flag, "ARKStepSetLinear", 1, );

        //flag = ARKStepSetJacFn (self->cvode, NULL /*J*/);
        //NCM_CVODE_CHECK (&flag, "ARKStepSetJacFn", 1, );

        //flag = ARKStepSetOrder (self->arkode, 7);
        //NCM_CVODE_CHECK (&flag, "ARKStepSetOrder", 1, );

        //flag = ARKStepSetERKTableNum (self->arkode, FEHLBERG_13_7_8);
        //NCM_CVODE_CHECK (&flag, "ARKStepSetERKTableNum", 1, );

        flag = ARKStepSetInitStep (self->arkode, fabs (t0) * self->reltol);
        NCM_CVODE_CHECK (&flag, "ARKStepSetInitStep", 1, );

        self->arkode_init = TRUE;
      }
      else
      {
        flag = ARKStepReInit (self->arkode, INTTYPE, t0, self->y);
        NCM_CVODE_CHECK (&flag, "ARKStepInit", 1, );

        flag = ARKStepSetInitStep (self->arkode, fabs (t0) * self->reltol);
        NCM_CVODE_CHECK (&flag, "ARKStepSetInitStep", 1, );
      }
      break;
    }
    default:
      g_error ("_nc_hipert_first_order_prepare_integrator: integrator %d not supported.", self->integ);
      break;
  }
}

/**
 * nc_hipert_first_order_prepare:
 * @fo: a #NcHIPertFirstOrder
 * @cosmo: a #NcHICosmo
 *
 * Adds a new component @comp to the system.
 *
 */
void 
nc_hipert_first_order_prepare (NcHIPertFirstOrder *fo, NcHICosmo *cosmo)
{
  NcHIPertFirstOrderPrivate * const self = fo->priv;
  NcHIPertBGVarYDY *ydy = nc_hipert_bg_var_ydy_new ();
  NcHIPertFirstOrderWS userdata = {fo, ydy, cosmo};
  gdouble tf = 1.0;
  gdouble t0;
  gint flag;

  if (self->vars->len == 0)
  {
    g_warning ("nc_hipert_first_order_prepare: empty system, nothing to do.");
    return;
  }
  
  nc_hipert_bg_var_prepare_if_needed (self->bg_var, cosmo);

  _nc_hipert_first_order_alloc_integrator (fo);
  t0 = _nc_hipert_first_order_set_init_cond (fo, cosmo, 1.0);
  _nc_hipert_first_order_prepare_integrator (fo, t0);

  switch (self->integ)
  {
    case NC_HIPERT_FIRST_ORDER_INTEG_CVODE:
    {
      flag = CVodeSetStopTime (self->cvode, tf);
      NCM_CVODE_CHECK (&flag, "CVodeSetStopTime", 1, );
    
      flag = CVodeSetUserData (self->cvode, &userdata);
      NCM_CVODE_CHECK (&flag, "CVodeSetUserData", 1, );

      while (TRUE)
      {
        gdouble t;
        
        flag = CVode (self->cvode, tf, self->y, &t, CV_ONE_STEP);
        NCM_CVODE_CHECK (&flag, "CVode", 1, );

        if (t == tf)
          break;
      }
    
      break;
    }
    case NC_HIPERT_FIRST_ORDER_INTEG_ARKODE:
    {
      flag = ARKStepSetStopTime (self->arkode, tf);
      NCM_CVODE_CHECK (&flag, "CVodeSetStopTime", 1, );

      flag = ARKStepSetUserData (self->arkode, &userdata);
      NCM_CVODE_CHECK (&flag, "ARKStepSetUserData", 1, );

      while (TRUE)
      {
        gdouble t;
        gint i;
        
        flag = ARKStepEvolve (self->arkode, tf, self->y, &t, ARK_ONE_STEP);
        NCM_CVODE_CHECK (&flag, "ARKStepEvolve", 1, );

        printf ("% 22.15g ", t);
        for (i = 0; i < self->cur_sys_size; i++)
        {
          printf ("% 22.15g ", NV_Ith_S (self->y, i));
        }
        printf ("\n");

        if (t == tf)
          break;
      }
      
      break;
    }      
    default:
      g_error ("_nc_hipert_first_order_prepare_integrator: integrator %d not supported.", self->integ);
      break;
  }

  nc_hipert_bg_var_ydy_free (ydy);
}


