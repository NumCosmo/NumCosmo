/* -*- Mode: C; indent-tabs-mode: nil; c-basic-offset: 2; tab-width: 2 -*-  */
/***************************************************************************
 *            nc_hipert_grav_einstein.c
 *
 *  Fri October 13 10:37:11 2017
 *  Copyright  2017  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * nc_hipert_grav_einstein.c
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
 * SECTION:nc_hipert_grav_einstein
 * @title: NcHIPertGravEinstein
 * @short_description: First order Einstein equations on a Friedmann background.
 *
 * FIXME
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "perturbations/nc_hipert_grav_einstein.h"

struct _NcHIPertGravEinsteinPrivate
{
  gint a;
};

enum
{
  PROP_0,
  PROP_NHOC
};

G_DEFINE_TYPE_WITH_PRIVATE (NcHIPertGravEinstein, nc_hipert_grav_einstein, NC_TYPE_HIPERT_GRAV);

static void
nc_hipert_grav_einstein_init (NcHIPertGravEinstein *gr)
{
  gr->priv = nc_hipert_grav_einstein_get_instance_private (gr);
}

static void
_nc_hipert_grav_einstein_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  g_return_if_fail (NC_IS_HIPERT_GRAV_EINSTEIN (object));

  switch (prop_id)
  {
    case PROP_NHOC:
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_hipert_grav_einstein_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  g_return_if_fail (NC_IS_HIPERT_GRAV_EINSTEIN (object));

  switch (prop_id)
  {
    case PROP_NHOC:
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_hipert_grav_einstein_dispose (GObject *object)
{

  /* Chain up : end */
  G_OBJECT_CLASS (nc_hipert_grav_einstein_parent_class)->dispose (object);
}

static void
_nc_hipert_grav_einstein_finalize (GObject *object)
{

  /* Chain up : end */
  G_OBJECT_CLASS (nc_hipert_grav_einstein_parent_class)->finalize (object);
}

NC_HIPERT_BG_VAR_ID_FUNC_IMPL (nc_hipert_grav_einstein, NcHIPertGravEinstein)

static guint _nc_hipert_grav_einstein_ndyn_var (NcHIPertGrav *grav);
static GArray *_nc_hipert_grav_einstein_get_deps (NcHIPertGrav *grav, guint vindex);

static NcHIPertGravInfo *_nc_hipert_grav_einstein_get_G_scalar_info (NcHIPertGrav *grav);
static void _nc_hipert_grav_einstein_get_G_scalar (NcHIPertGrav *grav, NcHIPertBGVar *bg_var, NcHIPertBGVarYDY *ydy, NcHIPertGravTScalar *T_scalar, NcHIPertGravScalar *G_scalar);
static void _nc_hipert_grav_einstein_get_dy_scalar (NcHIPertGrav *grav, NcHIPertBGVar *bg_var, NcHIPertBGVarYDY *ydy, NcHIPertGravTScalar *T_scalar, NcHIPertGravScalar *G_scalar);

static void
nc_hipert_grav_einstein_class_init (NcHIPertGravEinsteinClass *klass)
{
  GObjectClass *object_class    = G_OBJECT_CLASS (klass);
  NcHIPertGravClass *grav_class = NC_HIPERT_GRAV_CLASS (klass);

  object_class->set_property = &_nc_hipert_grav_einstein_set_property;
  object_class->get_property = &_nc_hipert_grav_einstein_get_property;
  object_class->dispose      = &_nc_hipert_grav_einstein_dispose;
  object_class->finalize     = &_nc_hipert_grav_einstein_finalize;

  g_object_class_install_property (object_class,
                                   PROP_NHOC,
                                   g_param_spec_int ("nhoc",
                                                     NULL,
                                                     "nhoc",
                                                     G_MININT, G_MAXINT, 0,
                                                     G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  nc_hipert_bg_var_class_register_id ("NcHIPertGravEinstein", 
                                      "First order Einstein equations background variables", 
                                      NULL,
                                      0);

  
  grav_class->ndyn_var          = &_nc_hipert_grav_einstein_ndyn_var;
  grav_class->get_deps          = &_nc_hipert_grav_einstein_get_deps;
  grav_class->get_G_scalar_info = &_nc_hipert_grav_einstein_get_G_scalar_info;
  grav_class->get_G_scalar      = &_nc_hipert_grav_einstein_get_G_scalar;
  grav_class->get_dy_scalar     = &_nc_hipert_grav_einstein_get_dy_scalar;
}

#define LEN(a) (sizeof (a) / sizeof (*(a)))
#define APPEND(a,b) (g_array_append_vals ((a), (b), LEN (b)))

static guint 
_nc_hipert_grav_einstein_ndyn_var (NcHIPertGrav *grav)
{
  guint ndyn_var = 0;
  switch (nc_hipert_grav_get_gauge (grav))
  {
    case NC_HIPERT_GRAV_GAUGE_SYNCHRONOUS:
    case NC_HIPERT_GRAV_GAUGE_NEWTONIAN:
      ndyn_var = 1;
      break;
    case NC_HIPERT_GRAV_GAUGE_CONST_CURV:
    case NC_HIPERT_GRAV_GAUGE_CONST_EXP:
      ndyn_var = 0;
      break;
    default:
      g_assert_not_reached ();
      break;
  }

  return ndyn_var;
}

static GArray *
_nc_hipert_grav_einstein_get_deps (NcHIPertGrav *grav, guint vindex)
{
  GArray *deps = g_array_new (TRUE, TRUE, sizeof (gint));
  
  switch (nc_hipert_grav_get_gauge (grav))
  {
    case NC_HIPERT_GRAV_GAUGE_SYNCHRONOUS:
    {
      switch (vindex)
      {
        case NC_HIPERT_GRAV_DYN_VAR (0):
        {
          NcHIPertGravSElem deps_a[] = {NC_HIPERT_GRAV_SELEM_DSIGMA, NC_HIPERT_GRAV_SELEM_RHOPPV};
          APPEND (deps, deps_a);
          break;
        }
        default:
          g_assert_not_reached ();
          break;
      }
      break;
    }
    case NC_HIPERT_GRAV_GAUGE_NEWTONIAN:
    {
      switch (vindex)
      {
        case NC_HIPERT_GRAV_DYN_VAR (0):
        {
          NcHIPertGravSElem deps_a[] = {NC_HIPERT_GRAV_SELEM_PHI, NC_HIPERT_GRAV_SELEM_RHOPPV};
          APPEND (deps, deps_a);
          break;
        }
        default:
          g_assert_not_reached ();
          break;
      }
      break;
    }
    case NC_HIPERT_GRAV_GAUGE_CONST_CURV:
    {
      break;
    }
    case NC_HIPERT_GRAV_GAUGE_CONST_EXP:
    {
      break;
    }
    default:
      g_assert_not_reached ();
      break;
  }

  return deps;
}

static NcHIPertGravInfo * 
_nc_hipert_grav_einstein_get_G_scalar_info (NcHIPertGrav *grav)
{
  NcHIPertGravInfo *ginfo = nc_hipert_grav_info_new ();
    
  switch (nc_hipert_grav_get_gauge (grav))
  {
    case NC_HIPERT_GRAV_GAUGE_SYNCHRONOUS:
    {
      NcHIPertGravSElem dsigma_deps_a[] = {NC_HIPERT_GRAV_SELEM_PSI, NC_HIPERT_GRAV_SELEM_DRHO, NC_HIPERT_GRAV_SELEM_RHOPPV};
      NcHIPertGravSElem psi_deps_a[]    = {NC_HIPERT_GRAV_DYN_VAR (0)};
      NcHIPertGravSElem dotpsi_deps_a[] = {NC_HIPERT_GRAV_SELEM_DSIGMA, NC_HIPERT_GRAV_SELEM_RHOPPV};

      APPEND (ginfo->dsigma_deps, dsigma_deps_a);
      APPEND (ginfo->psi_deps,    psi_deps_a);
      APPEND (ginfo->dotpsi_deps, dotpsi_deps_a);

      break;
    }
    case NC_HIPERT_GRAV_GAUGE_NEWTONIAN:
    {
      NcHIPertGravSElem phi_deps_a[]    = {NC_HIPERT_GRAV_SELEM_PSI, NC_HIPERT_GRAV_SELEM_DPI};
      NcHIPertGravSElem psi_deps_a[]    = {NC_HIPERT_GRAV_DYN_VAR (0)};
      NcHIPertGravSElem dotpsi_deps_a[] = {NC_HIPERT_GRAV_SELEM_PHI, NC_HIPERT_GRAV_SELEM_RHOPPV};

      APPEND (ginfo->phi_deps,    phi_deps_a);
      APPEND (ginfo->psi_deps,    psi_deps_a);
      APPEND (ginfo->dotpsi_deps, dotpsi_deps_a);

      break;
    }
    case NC_HIPERT_GRAV_GAUGE_CONST_CURV:
      break;
    default:
      g_assert_not_reached ();
      break;
  }

  return ginfo;
}

static void 
_nc_hipert_grav_einstein_get_G_scalar (NcHIPertGrav *grav, NcHIPertBGVar *bg_var, NcHIPertBGVarYDY *ydy, NcHIPertGravTScalar *T_scalar, NcHIPertGravScalar *G_scalar)
{ 
 
}

static void 
_nc_hipert_grav_einstein_get_dy_scalar (NcHIPertGrav *grav, NcHIPertBGVar *bg_var, NcHIPertBGVarYDY *ydy, NcHIPertGravTScalar *T_scalar, NcHIPertGravScalar *G_scalar)
{
  
}

/**
 * nc_hipert_grav_einstein_new:
 *
 * Creates a new #NcHIPertGravEinstein.
 *
 * Returns: (transfer full): the newly created #NcHIPertGravEinstein.
 */
NcHIPertGravEinstein *
nc_hipert_grav_einstein_new (void)
{
  NcHIPertGravEinstein *gr = g_object_new (NC_TYPE_HIPERT_GRAV_EINSTEIN,
                                           NULL);

  return gr;
}

/**
 * nc_hipert_grav_einstein_ref:
 * @gr: a #NcHIPertGravEinstein
 *
 * Increases the reference count of @gr.
 *
 * Returns: (transfer full): @gr.
 */
NcHIPertGravEinstein *
nc_hipert_grav_einstein_ref (NcHIPertGravEinstein *gr)
{
  return g_object_ref (gr);
}

/**
 * nc_hipert_grav_einstein_free:
 * @gr: a #NcHIPertGravEinstein
 *
 * Decreases the reference count of @gr.
 *
 */
void 
nc_hipert_grav_einstein_free (NcHIPertGravEinstein *gr)
{
  g_object_unref (gr);
}

/**
 * nc_hipert_grav_einstein_clear:
 * @gr: a #NcHIPertGravEinstein
 *
 * Decreases the reference count of *@gr and sets the pointer *@gr to NULL.
 *
 */
void 
nc_hipert_grav_einstein_clear (NcHIPertGravEinstein **gr)
{
  g_clear_object (gr);
}
