/* -*- Mode: C; indent-tabs-mode: nil; c-basic-offset: 2; tab-width: 2 -*-  */
/***************************************************************************
 *            nc_hipert_comp_pb.c
 *
 *  Fri October 13 11:10:24 2017
 *  Copyright  2017  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * nc_hipert_comp_pb.c
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
 * SECTION:nc_hipert_comp_pb
 * @title: NcHIPertCompPB
 * @short_description: Photon-Baryon plasma compoment
 *
 * FIXME
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "perturbations/nc_hipert_comp_pb.h"

struct _NcHIPertCompPBPrivate
{
  guint lmax;
};

enum
{
  PROP_0,
  PROP_LMAX
};

G_DEFINE_TYPE_WITH_PRIVATE (NcHIPertCompPB, nc_hipert_comp_pb, NC_TYPE_HIPERT_COMP);

static void
nc_hipert_comp_pb_init (NcHIPertCompPB *pb)
{
  pb->priv = nc_hipert_comp_pb_get_instance_private (pb);

  pb->priv->lmax = 0;
}

static void
_nc_hipert_comp_pb_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcHIPertCompPB *pb = NC_HIPERT_COMP_PB (object);
  g_return_if_fail (NC_IS_HIPERT_COMP_PB (object));

  switch (prop_id)
  {
    case PROP_LMAX:
      nc_hipert_comp_pb_set_lmax (pb, g_value_get_uint (value));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_hipert_comp_pb_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcHIPertCompPB *pb = NC_HIPERT_COMP_PB (object);
  g_return_if_fail (NC_IS_HIPERT_COMP_PB (object));

  switch (prop_id)
  {
    case PROP_LMAX:
      g_value_set_uint (value, nc_hipert_comp_pb_get_lmax (pb));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_hipert_comp_pb_dispose (GObject *object)
{

  /* Chain up : end */
  G_OBJECT_CLASS (nc_hipert_comp_pb_parent_class)->dispose (object);
}

static void
_nc_hipert_comp_pb_finalize (GObject *object)
{
  
  /* Chain up : end */
  G_OBJECT_CLASS (nc_hipert_comp_pb_parent_class)->finalize (object);
}

NC_HIPERT_BG_VAR_ID_FUNC_IMPL (nc_hipert_comp_pb, NcHIPertCompPB);

static guint _nc_hipert_comp_pb_ndyn_var (NcHIPertComp *comp);
static GArray *_nc_hipert_comp_pb_get_deps (NcHIPertComp *comp, guint vindex);
static NcHIPertGravTScalarInfo *_nc_hipert_comp_pb_get_T_scalar_info (NcHIPertComp *comp);
static void _nc_hipert_comp_pb_get_T_scalar (NcHIPertComp *comp, NcHIPertBGVar *bg_var, NcHIPertBGVarYDY *ydy, NcHIPertGravTScalar *T_scalar);
static void _nc_hipert_comp_pb_get_dy_scalar (NcHIPertComp *comp, NcHIPertBGVar *bg_var, NcHIPertBGVarYDY *ydy, NcHIPertGravTScalar *T_scalar, NcHIPertGravScalar *G_scalar);

static void
nc_hipert_comp_pb_class_init (NcHIPertCompPBClass *klass)
{
  GObjectClass *object_class    = G_OBJECT_CLASS (klass);
  NcHIPertCompClass *comp_class = NC_HIPERT_COMP_CLASS (klass);

  object_class->set_property = &_nc_hipert_comp_pb_set_property;
  object_class->get_property = &_nc_hipert_comp_pb_get_property;
  object_class->dispose      = &_nc_hipert_comp_pb_dispose;
  object_class->finalize     = &_nc_hipert_comp_pb_finalize;

  g_object_class_install_property (object_class,
                                   PROP_LMAX,
                                   g_param_spec_uint ("l-max",
                                                     NULL,
                                                     "l_max",
                                                     4, G_MAXUINT, 12,
                                                     G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  nc_hipert_bg_var_class_register_id ("NcHIPertCompPB", 
                                      "First order Photon-Baryons background variables", 
                                      NULL,
                                      0);

  comp_class->ndyn_var          = &_nc_hipert_comp_pb_ndyn_var;
  comp_class->get_deps          = &_nc_hipert_comp_pb_get_deps;
  comp_class->get_T_scalar_info = &_nc_hipert_comp_pb_get_T_scalar_info;
  comp_class->get_T_scalar      = &_nc_hipert_comp_pb_get_T_scalar;
  comp_class->get_dy_scalar     = &_nc_hipert_comp_pb_get_dy_scalar;
}

static guint 
_nc_hipert_comp_pb_ndyn_var (NcHIPertComp *comp)
{
  NcHIPertCompPB *pb = NC_HIPERT_COMP_PB (comp);
  return (pb->priv->lmax + 1) + 2;
}

#define LEN(a) (sizeof (a) / sizeof (*a))
#define APPEND(a,b) (g_array_append_vals ((a), (b), LEN (b)))

static GArray *
_nc_hipert_comp_pb_get_deps (NcHIPertComp *comp, guint vindex)
{
  NcHIPertCompPB *pb = NC_HIPERT_COMP_PB (comp);
  GArray *deps = g_array_new (TRUE, TRUE, sizeof (gint));

  switch (vindex)
  {
    case NC_HIPERT_COMP_PB_VAR_DELTA_B:
    {
      gint deps_a[] = {
        NC_HIPERT_COMP_PB_VAR_V_B,
        NC_HIPERT_GRAV_SELEM_DOTPSI, 
        NC_HIPERT_GRAV_SELEM_DSIGMA};

      APPEND (deps, deps_a);
      break;
    }
    case NC_HIPERT_COMP_PB_VAR_V_B:
    {
      gint deps_a[] = {
        NC_HIPERT_COMP_PB_VAR_DELTA_B,
        NC_HIPERT_COMP_PB_VAR_V_G,
        NC_HIPERT_GRAV_SELEM_PHI};

      APPEND (deps, deps_a);
      break;
    }
    case NC_HIPERT_COMP_PB_VAR_DELTA_G:
    {
      gint deps_a[] = {
        NC_HIPERT_COMP_PB_VAR_V_G, 
        NC_HIPERT_GRAV_SELEM_DOTPSI, 
        NC_HIPERT_GRAV_SELEM_DSIGMA};

      APPEND (deps, deps_a);
      break;
    }
    case NC_HIPERT_COMP_PB_VAR_V_G:
    {
      gint deps_a[] = {
        NC_HIPERT_COMP_PB_VAR_DELTA_G, 
        NC_HIPERT_COMP_PB_VAR_THETA_G, 
        NC_HIPERT_COMP_PB_VAR_V_B, 
        NC_HIPERT_GRAV_SELEM_PHI};
      
      APPEND (deps, deps_a);
      break;
    }
    case NC_HIPERT_COMP_PB_VAR_F_G3:
    {
      gint deps_a[] = { 
        NC_HIPERT_COMP_PB_VAR_THETA_G,
        NC_HIPERT_COMP_PB_VAR_F_G (4)
        };
      
      APPEND (deps, deps_a);
      break;
    }
    
    default:
    {
      if (vindex < NC_HIPERT_COMP_PB_VAR_F_G (pb->priv->lmax))
      {
        gint deps_a[] = { 
          vindex - 1,
          vindex + 1
          };

        APPEND (deps, deps_a);
      }
      else if (vindex == NC_HIPERT_COMP_PB_VAR_F_G (pb->priv->lmax))
      {
        gint deps_a[] = { 
          vindex - 1,
          };

        APPEND (deps, deps_a);
      }
      else
        g_assert_not_reached ();
      break;
    }
  }
  
  return deps;
}

static NcHIPertGravTScalarInfo * 
_nc_hipert_comp_pb_get_T_scalar_info (NcHIPertComp *comp)
{
  NcHIPertGravTScalarInfo *Tsinfo = nc_hipert_grav_T_scalar_info_new ();
  
  gint drho_deps_a[]   = {NC_HIPERT_COMP_PB_VAR_DELTA_B, NC_HIPERT_COMP_PB_VAR_DELTA_G};
  gint rhoppv_deps_a[] = {NC_HIPERT_COMP_PB_VAR_V_B, NC_HIPERT_COMP_PB_VAR_V_G};
  gint dPi_deps_a[]    = {NC_HIPERT_COMP_PB_VAR_THETA_G};

  APPEND (Tsinfo->drho_deps,   drho_deps_a);
  APPEND (Tsinfo->rhoppv_deps, rhoppv_deps_a);
  APPEND (Tsinfo->dPi_deps,    dPi_deps_a);

  return Tsinfo;
}

static void 
_nc_hipert_comp_pb_get_T_scalar (NcHIPertComp *comp, NcHIPertBGVar *bg_var, NcHIPertBGVarYDY *ydy, NcHIPertGravTScalar *T_scalar)
{ 
//  const gdouble delta_b = nc_hipert_bg_var_ydy_get_y_i (ydy, NC_HIPERT_COMP_PB_VAR_DELTA_B);
//  const gdouble delta_g = nc_hipert_bg_var_ydy_get_y_i (ydy, NC_HIPERT_COMP_PB_VAR_DELTA_G);
  
}

static void 
_nc_hipert_comp_pb_get_dy_scalar (NcHIPertComp *comp, NcHIPertBGVar *bg_var, NcHIPertBGVarYDY *ydy, NcHIPertGravTScalar *T_scalar, NcHIPertGravScalar *G_scalar)
{ 

}

#undef APPEND

/**
 * nc_hipert_comp_pb_new:
 *
 * Creates a new #NcHIPertCompPB.
 *
 * Returns: (transfer full): the newly created #NcHIPertCompPB.
 */
NcHIPertCompPB *
nc_hipert_comp_pb_new (void)
{
  NcHIPertCompPB *pb = g_object_new (NC_TYPE_HIPERT_COMP_PB,
                                     NULL);

  return pb;
}

/**
 * nc_hipert_comp_pb_ref:
 * @pb: a #NcHIPertCompPB
 *
 * Increases the reference count of @pb.
 *
 * Returns: (transfer full): @pb.
 */
NcHIPertCompPB *
nc_hipert_comp_pb_ref (NcHIPertCompPB *pb)
{
  return g_object_ref (pb);
}

/**
 * nc_hipert_comp_pb_free:
 * @pb: a #NcHIPertCompPB
 *
 * Decreases the reference count of @pb.
 *
 */
void 
nc_hipert_comp_pb_free (NcHIPertCompPB *pb)
{
  g_object_unref (pb);
}

/**
 * nc_hipert_comp_pb_clear:
 * @pb: a #NcHIPertCompPB
 *
 * Decreases the reference count of *@pb and sets the pointer *@pb to NULL.
 *
 */
void 
nc_hipert_comp_pb_clear (NcHIPertCompPB **pb)
{
  g_clear_object (pb);
}

/**
 * nc_hipert_comp_pb_set_lmax:
 * @pb: a #NcHIPertCompPB
 * @lmax: the maximum momentum of the photon distribution $\ell_\mathrm{max}$
 *
 * Sets the maximum momentum of the photon distribution to $\ell_\mathrm{max}=$@lmax.
 *
 */
void 
nc_hipert_comp_pb_set_lmax (NcHIPertCompPB *pb, guint lmax)
{
  if (pb->priv->lmax != lmax)
  {
    pb->priv->lmax = lmax;
  }
}

/**
 * nc_hipert_comp_pb_get_lmax:
 * @pb: a #NcHIPertCompPB
 *
 * Returns: the maximum momentum of the photon distribution to $\ell_\mathrm{max}$.
 */
guint 
nc_hipert_comp_pb_get_lmax (NcHIPertCompPB *pb)
{
  return pb->priv->lmax;
}
