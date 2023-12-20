/***************************************************************************
 *            nc_transfer_func_camb.c
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
 * SECTION:nc_transfer_func_camb
 * @title: NcTransferFuncCAMB
 * @short_description: Transfer function using CAMB as backend.
 *
 * FIXME
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "lss/nc_transfer_func_camb.h"
#include "math/ncm_spline_cubic_notaknot.h"
#include "math/ncm_cfg.h"

G_DEFINE_TYPE (NcTransferFuncCAMB, nc_transfer_func_camb, NC_TYPE_TRANSFER_FUNC)

gchar *camb_filename = NULL;

/**
 * nc_transfer_func_camb_new:
 *
 * FIXME
 *
 * Returns: A new #NcTransferFunc.
 */
NcTransferFunc *
nc_transfer_func_camb_new ()
{
  return g_object_new (NC_TYPE_TRANSFER_FUNC_CAMB, NULL);
}

static void
_nc_transfer_func_camb_prepare (NcTransferFunc *tf, NcHICosmo *cosmo)
{
  NcTransferFuncCAMB *tf_CAMB = NC_TRANSFER_FUNC_CAMB (tf);
  FILE *camb_tf;
  gint c, nlines = 0, i = 0, ret;
  GArray *x;
  GArray *y;

  if (tf_CAMB->init)
  {
    NcmVector *T_spline_xv = ncm_spline_peek_xv (tf_CAMB->T_spline);
    NcmVector *T_spline_yv = ncm_spline_peek_yv (tf_CAMB->T_spline);

    x = ncm_vector_get_array (T_spline_xv);
    y = ncm_vector_get_array (T_spline_yv);
    ncm_spline_free (tf_CAMB->T_spline);
  }
  else
  {
    x = g_array_new (FALSE, FALSE, sizeof (gdouble));
    y = g_array_new (FALSE, FALSE, sizeof (gdouble));
  }

  if (camb_filename == NULL)
    g_error ("To use camb transfer function first set the camb_filename extern variable.");

  camb_tf = fopen (camb_filename, "r");

  while ((c = fgetc (camb_tf)) != EOF)
    if (c == '\n')
      nlines++;

  rewind (camb_tf);

  g_array_set_size (x, nlines);
  g_array_set_size (y, nlines);

  i = 0;

  while ((ret = fscanf (camb_tf, " %lg %lg \n", &g_array_index (x, gdouble, i), &g_array_index (y, gdouble, i)) != EOF))
  {
    /*printf ("AQUI %d %d % 20.15g % 20.15g\n", ret, i, g_array_index (x, gdouble, i), g_array_index (y, gdouble, i)); */
    g_array_index (x, gdouble, i) = log (g_array_index (x, gdouble, i));
    g_array_index (y, gdouble, i) = log (g_array_index (y, gdouble, i));
    /*printf ("AQU1 %d %d % 20.15g % 20.15g\n", ret, i, g_array_index (x, gdouble, i), g_array_index (y, gdouble, i)); */
    i++;
  }

  tf_CAMB->T_spline = ncm_spline_cubic_notaknot_new ();
  ncm_spline_new_array (tf_CAMB->T_spline, x, y, TRUE);

  g_array_unref (x);
  g_array_unref (y);

  tf_CAMB->init = TRUE;
}

static gdouble
_nc_transfer_func_camb_calc (NcTransferFunc *tf, gdouble kh)
{
  NCM_UNUSED (tf);
  NCM_UNUSED (kh);
  g_assert_not_reached ();
}

static void
nc_transfer_func_camb_init (NcTransferFuncCAMB *tf_camb)
{
  /* TODO: Add initialization code here */
  tf_camb->T_spline = NULL;
  tf_camb->init     = FALSE;
}

static void
_nc_transfer_func_camb_dispose (GObject *object)
{
  NcTransferFuncCAMB *tf_camb = NC_TRANSFER_FUNC_CAMB (object);

  ncm_spline_clear (&tf_camb->T_spline);

  /* Chain up : end */
  G_OBJECT_CLASS (nc_transfer_func_camb_parent_class)->dispose (object);
}

static void
_nc_transfer_func_camb_finalize (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (nc_transfer_func_camb_parent_class)->finalize (object);
}

static void
nc_transfer_func_camb_class_init (NcTransferFuncCAMBClass *klass)
{
  GObjectClass *object_class        = G_OBJECT_CLASS (klass);
  NcTransferFuncClass *parent_class = NC_TRANSFER_FUNC_CLASS (klass);

  parent_class->prepare = &_nc_transfer_func_camb_prepare;
  parent_class->calc    = &_nc_transfer_func_camb_calc;

  object_class->dispose  = _nc_transfer_func_camb_dispose;
  object_class->finalize = _nc_transfer_func_camb_finalize;
}

