/***************************************************************************
 *            confidence_region.c
 *
 *  Mon Jun 11 13:28:00 2007
 *  Copyright  2007  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2012 <sandro@lapsandro>
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
 * SECTION:confidence_region
 * @title: Confidence Region
 * @short_description: FIXME
 *
 * FIXME
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "likelihood/confidence_region.h"
#include "math/ncm_priors.h"
#include "math/util.h"

#include <gsl/gsl_eigen.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_blas.h>

static gboolean ncm_fit_cr_step (NcConfidenceRegion *cr, gdouble x);
static gdouble nc_confidence_region_f (gdouble x, gpointer p);
static gdouble nc_confidence_region_df (gdouble x, gpointer p);
static void nc_confidence_region_fdf (gdouble x, gpointer p, gdouble *y, gdouble *dy);
static gdouble nc_confidence_region_numdiff_df (gdouble x, gpointer p);
static void nc_confidence_region_numdiff_fdf (gdouble x, gpointer p, gdouble *y, gdouble *dy);

static gdouble ncm_fit_cr_root_steffenson (NcConfidenceRegion *cr, gdouble x);
static gdouble ncm_fit_cr_root_brent (NcConfidenceRegion *cr, gdouble x0, gdouble x);

/**
 * nc_confidence_region_new_1d: (skip)
 * @fit: a #NcmFit.
 * @gmid: a #NcmModelID.
 * @pid: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
NcConfidenceRegion *
nc_confidence_region_new_1d (NcmFit *fit, NcmModelID gmid, guint pid)
{
  NcConfidenceRegion *cr = g_slice_new (NcConfidenceRegion);
  NcmMSet *mset = ncm_mset_copy_all (fit->mset);

  if (ncm_mset_param_get_ftype (fit->mset, gmid, pid) != NCM_PARAM_TYPE_FREE)
    g_error ("Cannot find confidence region for a non fitted parameter[%d:%u].", gmid, pid);

  ncm_mset_param_set_ftype (mset, gmid, pid, NCM_PARAM_TYPE_FIXED);

  cr->total_func_eval = 0;
  cr->n           = 1;
  cr->bestfit     = ncm_fit_ref (fit);
  cr->constrained = ncm_fit_copy_new (fit, fit->lh, mset, fit->grad.gtype);
  ncm_mset_free (mset);

  cr->shift[0] = 0.0;
  cr->pi[0].gmid = gmid;
  cr->pi[0].pid = pid;

  cr->search_type = NC_CONFIDENCE_REGION_SEARCH_1D;
  cr->chi2 = -1.0;
  cr->minimize = TRUE;

  cr->prior = g_slice_alloc (2 * sizeof (NcmPriorGauss));
  cr->points = NULL;

  return cr;
}

/**
 * nc_confidence_region_new: (skip)
 * @fit: a #NcmFit.
 * @gmid1: a #NcmModelID.
 * @pid1: FIXME
 * @gmid2: a #NcmModelID.
 * @pid2: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
NcConfidenceRegion *
nc_confidence_region_new (NcmFit *fit, NcmModelID gmid1, guint pid1, NcmModelID gmid2, guint pid2)
{
  NcConfidenceRegion *cr = g_slice_new (NcConfidenceRegion);
  NcmMSet *mset = ncm_mset_copy_all (fit->mset);
  gdouble local_cov_data[4];
  NcmMatrix *local_cov = ncm_matrix_new_data_static (local_cov_data, 2, 2);
  gsl_eigen_symmv_workspace *vw = gsl_eigen_symmv_alloc (2);
  gint fpi1 = ncm_mset_fparam_get_fpi (fit->mset, gmid1, pid1);
  gint fpi2 = ncm_mset_fparam_get_fpi (fit->mset, gmid2, pid2);

  if (fpi1 < 0)
    g_error ("Cannot find confidence region for a non fitted parameter[%02d:%02u].", gmid1, pid1);
  if (fpi2 < 0)
    g_error ("Cannot find confidence region for a non fitted parameter[%02d:%02u].", gmid2, pid2);

  g_assert (!((gmid1 == gmid2) && (pid1 == pid2)));

  if (fpi1 > fpi2)
    cr->inv_order = TRUE;
  else
    cr->inv_order = FALSE;

  ncm_mset_param_set_ftype (mset, gmid1, pid1, NCM_PARAM_TYPE_FIXED);
  ncm_mset_param_set_ftype (mset, gmid2, pid2, NCM_PARAM_TYPE_FIXED);

  cr->total_func_eval = 0;
  cr->n           = 2;
  cr->bestfit     = fit;
  cr->constrained = ncm_fit_copy_new (fit, fit->lh, mset, fit->grad.gtype);
  ncm_mset_free (mset);
  cr->covar_orto = ncm_matrix_new_sunk (2, 2);
  cr->covar_ev   = ncm_vector_new_sunk (2);

  cr->shift[0] = 0.0;
  cr->shift[1] = 0.0;

  cr->pi[0].gmid = gmid1;
  cr->pi[0].pid = pid1;
  cr->pi[1].gmid = gmid2;
  cr->pi[1].pid = pid2;

  ncm_matrix_set (local_cov, 0, 0, ncm_fit_covar_fparam_cov (fit, fpi1, fpi1));
  ncm_matrix_set (local_cov, 0, 1, ncm_fit_covar_fparam_cov (fit, fpi1, fpi2));
  ncm_matrix_set (local_cov, 1, 0, ncm_fit_covar_fparam_cov (fit, fpi2, fpi1));
  ncm_matrix_set (local_cov, 1, 1, ncm_fit_covar_fparam_cov (fit, fpi2, fpi2));

  gsl_eigen_symmv (NCM_MATRIX_GSL (local_cov), ncm_vector_gsl (cr->covar_ev), NCM_MATRIX_GSL (cr->covar_orto), vw);
  gsl_eigen_symmv_free (vw);
  ncm_matrix_transpose (cr->covar_orto);

  cr->search_type = NC_CONFIDENCE_REGION_SEARCH_ELLIPTIC_RADIUS;
  cr->chi2 = -1.0;
  cr->minimize = TRUE;

  cr->prior = g_slice_alloc (2 * sizeof (NcmPriorGauss));
  cr->points = NULL;

  ncm_matrix_free (local_cov);

  return cr;
}

/**
 * nc_confidence_region_free:
 * @cr: a #NcConfidenceRegion
 *
 * FIXME
 *
 * Returns: FIXME
 */
gboolean
nc_confidence_region_free (NcConfidenceRegion *cr)
{
  ncm_fit_free (cr->bestfit);
  ncm_fit_free (cr->constrained);
  if (cr->covar_ev != NULL)
    ncm_vector_free (cr->covar_ev);
  if (cr->covar_orto != NULL)
    ncm_matrix_free (cr->covar_orto);
  if (cr->points != NULL)
    g_list_free (cr->points);

  g_slice_free1 (2 * sizeof (NcmPriorGauss), cr->prior);
  g_slice_free (NcConfidenceRegion, cr);
  return TRUE;
}

/**
 * ncm_fit_cr_get_pi:
 * @cr: FIXME
 * @n: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
NcmMSetPIndex *
ncm_fit_cr_get_pi (NcConfidenceRegion *cr, guint n)
{
  return &cr->pi[n];
}

/**
 * ncm_fit_cr_get_pi_array:
 * @cr: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
NcmMSetPIndex *
ncm_fit_cr_get_pi_array (NcConfidenceRegion *cr)
{
  return cr->pi;
}

/**
 * ncm_fit_cr_get_param: (skip)
 * @cr: FIXME
 * @n: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
gdouble
ncm_fit_cr_get_param (NcConfidenceRegion *cr, guint n)
{
  return ncm_mset_param_get (cr->constrained->mset, cr->pi[n].gmid, cr->pi[n].pid);
}

/**
 * ncm_fit_cr_get_bf_param: (skip)
 * @cr: FIXME
 * @n: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
gdouble
ncm_fit_cr_get_bf_param (NcConfidenceRegion *cr, guint n)
{
  return ncm_mset_param_get (cr->bestfit->mset, cr->pi[n].gmid, cr->pi[n].pid);
}

/**
 * ncm_fit_cr_get_param_shift: (skip)
 * @cr: FIXME
 * @n: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
gdouble
ncm_fit_cr_get_param_shift (NcConfidenceRegion *cr, guint n)
{
  return cr->shift[n];
}

/**
 * ncm_fit_cr_set_param_shift: (skip)
 * @cr: FIXME
 * @n: FIXME
 * @s: FIXME
 *
 * FIXME
 *
 */
void
ncm_fit_cr_set_param_shift (NcConfidenceRegion *cr, guint n, gdouble s)
{
  cr->shift[n] = s;
}

/**
 * ncm_fit_cr_points_add: (skip)
 * @points: FIXME
 * @x: FIXME
 * @y: FIXME
 * @theta: FIXME
 * @p1: FIXME
 * @p2: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
gboolean
ncm_fit_cr_points_add (GList **points, gdouble x, gdouble y, gdouble theta, gdouble p1, gdouble p2)
{
  NcConfidenceRegion2dPoint *crp_new = g_slice_new (NcConfidenceRegion2dPoint);
  crp_new->x = x;
  crp_new->y = y;
  crp_new->theta = theta;
  crp_new->p1 = p1;
  crp_new->p2 = p2;
  *points = g_list_append (*points, crp_new);
  return TRUE;
}

/**
 * ncm_fit_cr_points_exists: (skip)
 * @points: FIXME
 * @x: FIXME
 * @y: FIXME
 * @n: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
gboolean
ncm_fit_cr_points_exists (GList **points, gdouble x, gdouble y, guint n)
{
  GList *pos = g_list_last (*points);
  while (n != 0)
  {
    NcConfidenceRegion2dPoint *crp1, *crp2;
    gdouble t1, t2;
    if (pos == NULL) break;
    crp1 = (NcConfidenceRegion2dPoint *) pos->data;
    pos = g_list_previous (pos);
    if (pos == NULL) break;
    crp2 = (NcConfidenceRegion2dPoint *) pos->data;
    t1 = (x - crp1->x) / (crp2->x - crp1->x);
    t2 = (y - crp1->y) / (crp2->y - crp1->y);
    if (fabs (t1 - t2) < 1e-1 && t1 >= 0.0 && t1 <= 1.0)
      return TRUE;
    n--;
  }
  return FALSE;
}

/**
 * ncm_fit_cr_points_print: (skip)
 * @points: FIXME
 * @out: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
gboolean
ncm_fit_cr_points_print (GList *points, FILE *out)
{
  GList *spoints = points;
  points = g_list_first (spoints);
  while (points)
  {
    NcConfidenceRegion2dPoint *crp = (NcConfidenceRegion2dPoint *) points->data;
    fprintf (out, "% 20.15g % 20.15g % 20.15g % 20.15g % 20.15g\n", crp->p1, crp->p2, crp->x, crp->y, crp->theta);fflush(out);
    points = g_list_next (points);
  }
  points = g_list_first (spoints);
  {
    NcConfidenceRegion2dPoint *crp = (NcConfidenceRegion2dPoint *) points->data;
    fprintf (out, "% 20.15g % 20.15g % 20.15g % 20.15g % 20.15g\n", crp->p1, crp->p2, crp->x, crp->y, crp->theta);fflush(out);
  }
  return TRUE;
}

/**
 * ncm_fit_cr_points_free: (skip)
 * @points: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
gboolean
ncm_fit_cr_points_free (GList *points)
{
  points = g_list_first (points);
  while (points)
  {
    NcConfidenceRegion2dPoint *crp = (NcConfidenceRegion2dPoint *) points->data;
    g_slice_free (NcConfidenceRegion2dPoint, crp);
    points = g_list_next (points);
  }
  g_list_free (points);
  return TRUE;
}

/**
 * ncm_fit_cr:
 * @fit: a #NcmFit
 * @gmid1: FIXME
 * @pid1: FIXME
 * @gmid2: FIXME
 * @pid2: FIXME
 * @p: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
gboolean
ncm_fit_cr (NcmFit *fit, NcmModelID gmid1, guint pid1, NcmModelID gmid2, guint pid2, gdouble p)
{
  gdouble r, theta;
  gdouble r0 = 0.0;
  gdouble step;
  gdouble new_x, new_y;
  gsl_rng *rand = ncm_get_rng();
  NcConfidenceRegion *cr = nc_confidence_region_new (fit, gmid1, pid1, gmid2, pid2);

  cr->chi2 = gsl_cdf_chisq_Qinv (1.0 - p, 2);
  r = sqrt (cr->chi2);
  step = 0.001;

  cr->search_type = NC_CONFIDENCE_REGION_SEARCH_ELLIPTIC_RADIUS;

  for (theta = 0.0; theta < 2.0 * M_PI; theta += step)
  {
    cr->theta = theta;
    r = ncm_fit_cr_root_steffenson (cr, r);
    while (!gsl_finite(r) || r < 0.0)
      r = ncm_fit_cr_root_steffenson (cr, gsl_rng_uniform (rand));

    new_x = ncm_fit_cr_get_param (cr, 0);
    new_y = ncm_fit_cr_get_param (cr, 1);

    if ((theta != 0.0) && fabs((r-r0)/r0) > 0.10)
    {
      theta -= step;
      step /= 2.0;
      if (step > 1e-4)
        continue;
      else
        g_error ("Found a discontinuity");
    }
    step = atan( 0.01 / fabs(r) );

    printf ("\t%g %g %g %g %g\n", new_x, new_y, r, theta, step);
    fflush (stdout);
    r0 = r;
  }

  nc_confidence_region_free (cr);

  return TRUE;
}

#define BASE_SCALE(a) (2.0 * M_PI * sqrt (a) / 100.0)
#define RESCALE (0.5)
#define NMAXTRIES 40
#define TIMEOUT 90.0

/**
 * ncm_fit_cr2: (skip)
 * @fit: a @NcmFit
 * @gmid1: a #NcmModelID.
 * @pid1: FIXME
 * @gmid2: a #NcmModelID.
 * @pid2: FIXME
 * @p: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
GList *
ncm_fit_cr2 (NcmFit *fit, NcmModelID gmid1, guint pid1, NcmModelID gmid2, guint pid2, gdouble p)
{
  NcConfidenceRegion *cr;
  GTimer *iter_timer = g_timer_new ();
  gdouble total_time = 0.0;
  gdouble old_x, old_y;
  gdouble init_x, init_y;
  gsl_rng *rand = ncm_get_rng();
  gint i, counter = -1;
  GList *points = NULL, *final_points = NULL;
  gdouble theta = gsl_rng_uniform (rand) * 2.0 * M_PI;
  gboolean completed = FALSE;
  gdouble second_try = FALSE;

  cr = nc_confidence_region_new (fit, gmid1, pid1, gmid2, pid2);
  cr->chi2 = gsl_cdf_chisq_Qinv (1.0 - p, 2);

  while (!completed)
  {
    gboolean start = TRUE;
    gboolean end = FALSE;
    gboolean count_error = TRUE;
    guint change_scale = 0;
    guint tries = 0;
    guint tries_d = 0;
    guint repet_error = 0;
    gdouble vval, r0, r, scale;

    cr->search_type = NC_CONFIDENCE_REGION_SEARCH_ELLIPTIC_RADIUS;
    scale = sqrt (cr->chi2);
    r0 = 0.0;
    r = scale;
    cr->theta = (NC_CR_COVAR_EV (cr, 0) < NC_CR_COVAR_EV (cr, 0)) ? 0.0 : M_PI / 2.0;
    while ((vval = nc_confidence_region_f (r, cr)) < 0)
    {
      r0 = r;
      r += scale;
    }
    g_message ("#  Jumping to %.15g : %.15g =~ 0.0 [chi2 %.8f %2.4f%%]\n", cr->r, vval, cr->chi2, p * 100.0);

    {
      //	  cr->r = ncm_fit_cr_root_brent (cr, r0, r);
      cr->r = ncm_fit_cr_root_steffenson (cr, r);
      g_message ("#  Corrected to %.15g\n", cr->r);

      completed = TRUE;
      if (second_try)
        completed = TRUE;

      while (!gsl_finite (cr->r) || cr->r < 0.0)
      {
        cr->r = (gsl_rng_uniform (rand) - 0.5) * r0 / 5.0 + r0;
        g_message ("#  Trying[%d] r = %g, theta = %f\n", tries, cr->r, cr->theta * 180.0 / M_PI);
        ncm_mset_copyto (cr->bestfit->mset, cr->constrained->mset);
        cr->r = ncm_fit_cr_root_steffenson (cr, cr->r);
        tries++;
        if (tries == 1000)
          g_error ("Cannot find the border of (%2.4f) C.L.", p);
      }
    }

    ncm_fit_cr_step (cr, cr->r);

    cr->shift[0] = cr->r * cos (cr->theta);
    cr->shift[1] = cr->r * sin (cr->theta);
    ncm_fit_cr_points_add (&points, cr->shift[0], cr->shift[1], GSL_NAN, ncm_fit_cr_get_param (cr, 0), ncm_fit_cr_get_param (cr, 1));
    init_x = cr->shift[0];
    init_y = cr->shift[1];

    cr->search_type = NC_CONFIDENCE_REGION_SEARCH_ELLIPTIC_ANGLE;
    cr->r = BASE_SCALE(cr->chi2);
    cr->theta = 2.0 * M_PI;

    {
      gdouble theta0 = GSL_NAN;
      gdouble first_theta = GSL_NAN;
      gdouble first_shift0 = cr->shift[0];
      gdouble first_shift1 = cr->shift[1];

      old_x = 0.0;
      old_y = 0.0;
      i = 0;
      total_time += g_timer_elapsed (iter_timer, NULL);
      g_timer_start (iter_timer);

      while (TRUE)
      {
        theta0 = M_PI * 0.5, theta = 1.5 * M_PI;
        printf ("# THETA % 20.15g % 20.15g | % 20.15g\n", theta0, theta, cr->theta);
        //cr->theta = ncm_fit_cr_root_brent (cr, theta0, theta);
        cr->theta = ncm_fit_cr_root_steffenson (cr, cr->theta);
        theta0 = cr->theta - 0.5 * M_PI;
        theta  = cr->theta + 0.5 * M_PI;
        printf ("# THETA % 20.15g % 20.15g | % 20.15g\n", theta0, theta, cr->theta);
        tries = 0;
        while (!gsl_finite (cr->theta))
        {
          if (g_timer_elapsed (iter_timer, NULL) > TIMEOUT && !start)
          {
            g_message ("#  Timeout finding next point (%f)\n", g_timer_elapsed (iter_timer, NULL));
            cr->theta = GSL_NAN;
            break;
          }
          if (repet_error > 5)
          {
            g_message ("#  Max consecutive error count reached (%f)\n", g_timer_elapsed (iter_timer, NULL));
            break;
          }
          cr->theta = gsl_rng_uniform (rand) * 2.0 * M_PI;
          g_message ("#  Trying[%d,%d] theta = %g =>%f<=\n", tries, repet_error, cr->theta, g_timer_elapsed (iter_timer, NULL));
          cr->theta = ncm_fit_cr_root_steffenson (cr, cr->theta);
          if (tries > 1)
          {
            g_message ("#  Changing scale...[%g] -> [%g] =>%f<=\n", cr->r, cr->r*RESCALE, g_timer_elapsed (iter_timer, NULL));
            cr->r *= RESCALE;
            change_scale += 5;
          }
          tries++;
          if (count_error) {repet_error++;count_error=FALSE;}
          if (tries == 100)
            break;
        }
        cr->theta = NC_RADIAN_0_2PI (cr->theta);
        first_theta = gsl_finite (cr->theta) ? cr->theta : GSL_NAN;

        if (!start)
        {
          gdouble diff = fabs(cr->theta - theta0);
          gboolean going_back;
          ncm_fit_cr_step (cr, cr->theta);
          going_back = ncm_fit_cr_points_exists (&points, NC_CR_X (cr, cr->r, cr->theta), NC_CR_Y (cr, cr->r, cr->theta), -1);
          diff = NC_RADIAN_0_2PI(diff);

          if ((fabs(1.0 - diff/M_PI) < 0.010) || going_back)
          {
            cr->theta = 2.0 * M_PI * gsl_rng_uniform (rand);
            g_message ("#  Inverse path found [%d], retrying [%d, %d] =>%f<=...\n", going_back, tries_d, repet_error, g_timer_elapsed (iter_timer, NULL));
            if (tries_d > 1)
            {
              cr->theta = theta0;
              g_message ("#  Changing scale...[%g] -> [%g] =>%f<=\n", cr->r, cr->r * RESCALE, g_timer_elapsed (iter_timer, NULL));
              cr->r *= RESCALE;
              change_scale += 5;
            }
            if (count_error) {repet_error++;count_error=FALSE;}
            tries_d++;
            if (repet_error >= 5)
            {
              g_message ("#  Max consecutive error count reached (%f)\n", g_timer_elapsed (iter_timer, NULL));
              cr->theta = GSL_NAN;
            }
            else if (tries_d == NMAXTRIES)
              cr->theta = GSL_NAN;
            else if (g_timer_elapsed (iter_timer, NULL) > TIMEOUT && !start)
            {
              g_message ("#  Timeout finding next point (%f)\n", g_timer_elapsed (iter_timer, NULL));
              cr->theta = GSL_NAN;
            }
            else
              continue;
          }
        }
        start = FALSE;
        tries_d = 0;
        if (!count_error) {count_error=TRUE;}
        else {repet_error = 0;}

        if (change_scale)
        {
          if (change_scale > 0)
            change_scale--;
          if (!change_scale)
          {
            g_message ("#  Returning to default scale\n");
            cr->r = BASE_SCALE(cr->chi2);
            change_scale = 0;
          }
        }

        if (!gsl_finite (cr->theta) || i > 5000)
        {
          g_message ("#  Cannot continue[%d], reverting direction [%d]...\n", i, end);
          if (end) break;
          points = g_list_reverse (points);
          ncm_mset_copyto (cr->bestfit->mset, cr->constrained->mset);
          cr->shift[0] = first_shift0;
          cr->shift[1] = first_shift1;
          theta0 = first_theta + M_PI;
          init_x = old_x;
          init_y = old_y;
          cr->r = BASE_SCALE(cr->chi2);
          change_scale = 0;
          cr->theta = 0.0;
          i = 0;
          end = TRUE;
          total_time += g_timer_elapsed (iter_timer, NULL);
          g_timer_start (iter_timer);
          count_error = TRUE;
          repet_error = 0;
          continue;
        }

        g_message ("# Found point [%lu] (% 12.10f % 12.10f ) (% 12.10f % 12.10f ) theta: %.5f theta0: %.5f\n# Profile params ",
                   cr->total_func_eval,
                   ncm_fit_cr_get_param (cr, 0), ncm_fit_cr_get_param (cr, 1),
                   NC_CR_X (cr, cr->r, cr->theta), NC_CR_Y (cr, cr->r, cr->theta),
                   cr->theta * 180.0 / M_PI, theta0 * 180.0 / M_PI);
        ncm_mset_params_log_vals (cr->constrained->mset);

        ncm_fit_cr_points_add (&points, NC_CR_X (cr, cr->r, cr->theta), NC_CR_Y (cr, cr->r, cr->theta), cr->theta, ncm_fit_cr_get_param (cr, 0), ncm_fit_cr_get_param (cr, 1));

        old_x = cr->shift[0];
        old_y = cr->shift[1];
        ncm_fit_cr_step (cr, cr->theta);
        cr->shift[0] += cr->r * cos (cr->theta);
        cr->shift[1] += cr->r * sin (cr->theta);

        counter--;
        if (!counter)
          break;

        {
          gboolean near_x = fabs(init_x - cr->shift[0]) < cr->r;
          gboolean near_y = fabs(init_y - cr->shift[1]) < cr->r;
          if (i > 10 && near_x && near_y && counter < 0)
          {
            g_message ("#  Start found at [%d], ending...\n",i);
            completed = TRUE;
            break;
            counter = 5;
          }
        }
        total_time += g_timer_elapsed (iter_timer, NULL);
        g_timer_start (iter_timer);
        i++;
        theta0 = cr->theta;
      }
    }
    ncm_fit_cr_points_add (&points, NC_CR_X (cr, cr->r, cr->theta), NC_CR_Y (cr, cr->r, cr->theta), GSL_NAN, ncm_fit_cr_get_param (cr, 0), ncm_fit_cr_get_param (cr, 1));
    if (final_points == NULL)
    {
      final_points = points;
      points = NULL;
    }
    else
    {
      final_points = g_list_concat (final_points, g_list_reverse (points));
      points = NULL;
    }
    if (!completed)
    {
      second_try = TRUE;
      g_message ("#  Trying in another direction %g => %g\n", theta, NC_RADIAN_0_2PI(theta + M_PI));
      theta = NC_RADIAN_0_2PI(theta + M_PI);
      ncm_mset_copyto (cr->bestfit->mset, cr->constrained->mset);
      cr->search_type = NC_CONFIDENCE_REGION_SEARCH_ELLIPTIC_RADIUS;
    }
  }

  total_time += g_timer_elapsed (iter_timer, NULL);
  g_timer_destroy (iter_timer);
  g_message ("#  C.R. (%f), took %fs = %fmin to complete\n", p, total_time, total_time / 60.0);

  return final_points;
}

/**
 * ncm_fit_cr2_fisher: (skip)
 * @fit: a #NcmFit
 * @gmid1: a #NcmModelID.
 * @pid1: FIXME
 * @gmid2: a #NcmModelID.
 * @pid2: FIXME
 * @p: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
GList *
ncm_fit_cr2_fisher (NcmFit *fit, NcmModelID gmid1, guint pid1, NcmModelID gmid2, guint pid2, gdouble p)
{
  NcConfidenceRegion *cr;
  GList *points = NULL;
  gdouble theta;

  cr = nc_confidence_region_new (fit, gmid1, pid1, gmid2, pid2);
  cr->chi2 = gsl_cdf_chisq_Qinv (1.0 - p, 2);
  cr->minimize = FALSE;
  cr->search_type = NC_CONFIDENCE_REGION_SEARCH_ELLIPTIC_ANGLE;
  cr->r = sqrt (cr->chi2);

  cr->shift[0] = 0.0;
  cr->shift[1] = 0.0;

  for (theta = 0.0; theta <= 2.0 * M_PI; theta += 0.01)
  {
    ncm_fit_cr_step (cr, theta);
    ncm_fit_cr_points_add (&points, NC_CR_X (cr, cr->r, cr->theta), NC_CR_Y (cr, cr->r, cr->theta), theta, ncm_fit_cr_get_param (cr, 0), ncm_fit_cr_get_param (cr, 1));
  }

  return points;
}

gboolean
ncm_fit_cr_1dim_diff (NcmFit *fit, NcmModelID gmid, guint pid, gdouble p, gint nu, gdouble *err_inf, gdouble *err_sup)
{
  NcConfidenceRegion *cr = nc_confidence_region_new_1d (fit, gmid, pid);
  gdouble r_min, r_max, r, scale;
  gsl_rng *rand = ncm_get_rng();

  scale = ncm_fit_covar_sd (cr->bestfit, gmid, pid);
  cr->chi2 = gsl_cdf_chisq_Qinv (1.0 - p, nu);
  cr->search_type = NC_CONFIDENCE_REGION_SEARCH_1D;

  r = -scale;
  r = ncm_fit_cr_root_steffenson (cr, r);

  while (!gsl_finite (r) || r > 0.0)
  {
    printf ("# Found r = %g,", r);
    r = -scale * gsl_rng_uniform (rand);
    printf (" trying r = %g\n", r);
    fflush (stdout);
    ncm_mset_copyto (cr->bestfit->mset, cr->constrained->mset);
    r = ncm_fit_cr_root_steffenson (cr, r);
  }
  r_min = r;
  //  printf ("# Lower end found (%g)\n", r);
  r = scale;
  r = ncm_fit_cr_root_steffenson (cr, r);
  while (!gsl_finite(r) || r < 0.0)
  {
    printf ("# Found r = %g,", r);
    r = scale * gsl_rng_uniform (rand);
    printf (" trying r = %g\n", r);
    fflush (stdout);
    ncm_mset_copyto (cr->bestfit->mset, cr->constrained->mset);
    r = ncm_fit_cr_root_steffenson (cr, r);
  }
  r_max = r;
  //  printf ("# Upper end found (%g)\n", r);
  *err_inf = r_min;
  *err_sup = r_max;
  return TRUE;
}

gboolean
ncm_fit_cr_1dim (NcmFit *fit, NcmModelID gmid, guint pid, gdouble p, gint nu, gdouble *err_inf, gdouble *err_sup)
{
  NcConfidenceRegion *cr = nc_confidence_region_new_1d (fit, gmid, pid);
  gdouble r_min, r_max, r, scale;

#define NC_CR_1DIM_SCALE_INCR (1.1)

  scale = ncm_fit_covar_sd (cr->bestfit, gmid, pid);
  cr->chi2 = gsl_cdf_chisq_Qinv (1.0 - p, nu);
  cr->search_type = NC_CONFIDENCE_REGION_SEARCH_1D;

  r = 0.0;
  r_min = -scale;
  while (nc_confidence_region_f (r_min, cr) < 0.0)
  {
    r = r_min;
    r_min *= NC_CR_1DIM_SCALE_INCR;
  }
  r_min = ncm_fit_cr_root_brent (cr, r_min, r);

  r = 0.0;
  r_max = scale;
  while (nc_confidence_region_f (r_max, cr) < 0.0)
  {
    r = r_max;
    r_max *= NC_CR_1DIM_SCALE_INCR;
  }
  r_max = ncm_fit_cr_root_brent (cr, r, r_max);

  *err_inf = r_min;
  *err_sup = r_max;
  return TRUE;
}

gdouble
ncm_fit_type_constrain_error (NcmFit *fit, gdouble p, gint nu, gdouble dir, NcmMSetFunc *func, gdouble z, gboolean walk)
{
  g_assert_not_reached ();
  return 0.0;
  /*
   NcConfidenceRegion *cr = g_slice_new (NcConfidenceRegion);
   NcLikelihood *clh;
   gdouble r;
   gsl_rng *rand = ncm_get_rng();
   gdouble scale = ncm_fit_function_error (fit, func, z, FALSE);

   g_assert (dir != 0.0);
   scale *= dir;

   cr->prior[0].func = func;
   cr->prior[0].z = z;
   cr->prior[0].mean = ncm_mset_func_eval1 (cr->prior[0].func, fit->mset, z);
   cr->prior[0].sigma = cr->prior[0].mean * 1e-4;

   cr->n = 1;

   profile = ncm_fit_params_copy (fit->pt);
   cr_fit = ncm_fit_params_copy (fit->pt);

   clh = nc_likelihood_copy (fit->lh);

   nc_prior_add_gaussian (clh, &cr->prior[0]);

   cr->best = fit;
   cr->constrained = ncm_fit_new_with_types (clh, mset, fit->type, fit->grad.gtype, profile);
   cr->fit = ncm_fit_new_with_types (fit->lh, mset, fit->type, fit->grad.gtype, cr_fit);
   cr->chi2 = gsl_cdf_chisq_Qinv (1.0 - p, nu);

   cr->search_type = NC_CONFIDENCE_REGION_SEARCH_CONSTRAIN_1D;
   cr->shift[0] = cr->prior[0].mean;

   r = scale;
   if (walk)
   r = ncm_fit_cr_root_walkto (cr, GSL_SIGN(scale));
   r = ncm_fit_cr_root_steffenson (cr, r);

   if (FALSE)
   {
     ncm_fit_cr_root_brent (cr, 0, 0);
     r = ncm_fit_cr_root_walkto (cr, GSL_SIGN(scale));
}

while (!gsl_finite(r) || (GSL_SIGN(r) != GSL_SIGN(scale)))
{
  printf ("# r = %g, scale = %g ", r, scale);
  r = scale * gsl_pow_2 (gsl_rng_uniform (rand));
  printf ("| Trying r = %g\n", r);
  fflush (stdout);
  ncm_mset_copyto (cr->bestfit->mset, cr->constrained->mset);
  r = ncm_fit_cr_root_steffenson (cr, r);
}
nc_likelihood_free (clh);
g_slice_free (NcConfidenceRegion, cr);

return r;
*/
}

static gdouble
ncm_fit_cr_root_steffenson (NcConfidenceRegion *cr, gdouble x)
{
  gint status;
  gint iter = 0, max_iter = 1000000;
  const gsl_root_fdfsolver_type *T;
  gsl_root_fdfsolver *s;
  gsl_function_fdf F;
  gdouble x0, prec = 1e-2;

  F.f = &nc_confidence_region_f;
  if (cr->constrained->grad.gtype != NCM_FIT_GRAD_ANALYTICAL)
  {
    F.df = &nc_confidence_region_numdiff_df;
    F.fdf = &nc_confidence_region_numdiff_fdf;
  }
  else
  {
    F.df = &nc_confidence_region_df;
    F.fdf = &nc_confidence_region_fdf;
  }
  F.params = cr;

  T = gsl_root_fdfsolver_steffenson;
  //T = gsl_root_fdfsolver_newton;
  s = gsl_root_fdfsolver_alloc (T);
  gsl_root_fdfsolver_set (s, &F, x);

  do
  {
    iter++;
    status = gsl_root_fdfsolver_iterate (s);

    if (status)
    {
      ncm_mset_params_print_vals (cr->constrained->mset, stdout);
      g_warning ("%s", gsl_strerror (status));
      gsl_root_fdfsolver_free (s);
      return GSL_NAN;
    }
    if (iter > 100)
      prec *= 10;

    x0 = x;
    x = gsl_root_fdfsolver_root (s);
    status = gsl_root_test_delta (x, x0, 0, prec);

    if (!gsl_finite (nc_confidence_region_f (x, cr)))
    {
      x = GSL_NAN;
      break;
    }

    if ( (status == GSL_SUCCESS) && TRUE)
    {
      printf ("# Converged: %d, with prec = %g\n", iter, prec);
    }{
      printf ("#\t[%d]: %g %g\n", iter, x, nc_confidence_region_f (x, cr));
      fflush (stdout);
    }
  }
  while (status == GSL_CONTINUE && iter < max_iter);

  if (fabs(prec) > 1.0)
    x = GSL_NAN;
  gsl_root_fdfsolver_free (s);
  return x;
}

static gdouble
ncm_fit_cr_root_brent (NcConfidenceRegion *cr, gdouble x0, gdouble x)
{
  gint status;
  gint iter = 0, max_iter = 1000000;
  const gsl_root_fsolver_type *T;
  gsl_root_fsolver *s;
  gsl_function F;
  gdouble prec = 1e-5, x1 = x;

  F.function = &nc_confidence_region_f;
  F.params = cr;

  T = gsl_root_fsolver_brent;
  s = gsl_root_fsolver_alloc (T);
  gsl_root_fsolver_set (s, &F, x0, x1);

  //printf ("# Starting at % 20.15g % 20.15g\n", x0, x1);
  do
  {
    iter++;
    status = gsl_root_fsolver_iterate (s);
    if (status)
    {
      //      nc_params_print_all (cr->fit->mset, stdout);
      g_warning ("%s", gsl_strerror (status));
      gsl_root_fsolver_free (s);
      return GSL_NAN;
    }
    if (iter > 100)
      prec *= 10;

    x = gsl_root_fsolver_root (s);
    x0 = gsl_root_fsolver_x_lower (s);
    x1 = gsl_root_fsolver_x_upper (s);
    status = gsl_root_test_interval (x0, x1, 0, prec);

    if (!gsl_finite (nc_confidence_region_f (x, cr)))
    {
      g_debug ("Ops");
      x = GSL_NAN;
      break;
    }

    if ( (status == GSL_SUCCESS) && FALSE)
    {
      printf ("# Converged: %d, with prec = %g\n", iter, prec);
      //}{
      printf ("#\t[%d]: (%g, %g):%g %g\n", iter, x0, x1, x, nc_confidence_region_f (x, cr));
      fflush (stdout);
    }
  }
  while (status == GSL_CONTINUE && iter < max_iter);
  if (fabs(prec) > 1e-4)
    x = GSL_NAN;
  gsl_root_fsolver_free (s);
  return x;
}

static gboolean
ncm_fit_cr_step (NcConfidenceRegion *cr, gdouble x)
{
  guint i = 1;

  do {
    if (i < 1)
      ncm_mset_copyto (cr->bestfit->mset, cr->constrained->mset);
    switch (cr->search_type)
    {
      case NC_CONFIDENCE_REGION_SEARCH_1D:
      {
        NcmMSetPIndex *pi = ncm_fit_cr_get_pi (cr, 0);
        gdouble p_val = ncm_fit_cr_get_bf_param (cr, 0) + cr->shift[0] + x;
        const gdouble p_lb = ncm_mset_param_get_lower_bound (cr->constrained->mset, pi->gmid, pi->pid);
        const gdouble p_ub = ncm_mset_param_get_upper_bound (cr->constrained->mset, pi->gmid, pi->pid);
        p_val = GSL_MAX (p_val, p_lb);
        p_val = GSL_MIN (p_val, p_ub);
        ncm_mset_param_set (cr->constrained->mset, pi->gmid, pi->pid, p_val);
      }
        break;
      case NC_CONFIDENCE_REGION_SEARCH_ELLIPTIC_RADIUS:
      {
        gdouble xp_data[2], x_data[2];
        gsl_vector_view xp_view = gsl_vector_view_array (xp_data, 2);
        gsl_vector_view x_view = gsl_vector_view_array (x_data, 2);

        xp_data[0] = sqrt (NC_CR_COVAR_EV (cr, 0)) * NC_CR_X (cr, x, cr->theta);
        xp_data[1] = sqrt (NC_CR_COVAR_EV (cr, 1)) * NC_CR_Y (cr, x, cr->theta);

        x_data[0] = ncm_fit_cr_get_bf_param (cr, 0);
        x_data[1] = ncm_fit_cr_get_bf_param (cr, 1);

        gsl_blas_dgemv (CblasTrans, 1.0, NCM_MATRIX_GSL(cr->covar_orto), &xp_view.vector, 1.0, &x_view.vector);
        ncm_mset_param_set_pi (cr->constrained->mset, ncm_fit_cr_get_pi_array (cr), x_data, 2);
        break;
      }
      case NC_CONFIDENCE_REGION_SEARCH_ELLIPTIC_ANGLE:
      {
        gdouble xp_data[2], x_data[2];
        gsl_vector_view xp_view = gsl_vector_view_array (xp_data, 2);
        gsl_vector_view x_view = gsl_vector_view_array (x_data, 2);

        xp_data[0] = sqrt (NC_CR_COVAR_EV (cr, 0)) * NC_CR_X (cr, cr->r, x);
        xp_data[1] = sqrt (NC_CR_COVAR_EV (cr, 1)) * NC_CR_Y (cr, cr->r, x);

        x_data[0] = ncm_fit_cr_get_bf_param(cr, 0);
        x_data[1] = ncm_fit_cr_get_bf_param(cr, 1);

        gsl_blas_dgemv (CblasTrans, 1.0, NCM_MATRIX_GSL (cr->covar_orto), &xp_view.vector, 1.0, &x_view.vector);
        ncm_mset_param_set_pi (cr->constrained->mset, ncm_fit_cr_get_pi_array (cr), x_data, 2);
        break;
      }
      case NC_CONFIDENCE_REGION_SEARCH_CONSTRAIN_1D:
        cr->prior[0].mean = ncm_fit_cr_get_bf_param (cr, 0) + cr->shift[0] + x;
        break;
    }
  } while (cr->minimize && !ncm_fit_run (cr->constrained, NCM_FIT_RUN_MSGS_NONE) && i--);
  cr->total_func_eval += cr->constrained->func_eval;
  //  g_debug ("ok[%d][%g]: %g", i, cr->min[0] + x, gsl_pow_2((cr->prior[0].func.f (cr->constrained->cp, cr->prior[0].z) - cr->prior[0].mean) / (cr->prior[0].sigma)));
  if (cr->search_type == NC_CONFIDENCE_REGION_SEARCH_CONSTRAIN_1D &&
      (gsl_pow_2 ((ncm_mset_func_eval1 (cr->prior[0].func, cr->constrained->mset, cr->prior[0].z) - cr->prior[0].mean) / (cr->prior[0].sigma))) > 1e-4)
    return FALSE;

  if (i == -1)
    return FALSE;
  return TRUE;
}

static gdouble
nc_confidence_region_f (gdouble x, gpointer p)
{
  NcConfidenceRegion *cr = (NcConfidenceRegion *)p;
  gdouble f;
  if(!ncm_fit_cr_step (cr, x))
  {
    //    printf("\n# Failed to fit (f)\n#");
    return GSL_POSINF;
  }
  f = cr->constrained->m2lnL - (cr->bestfit->m2lnL + cr->chi2);
  //g_debug ("f %.15g | %.15g %.15g %.15g\n", f, cr->constrained->m2lnL, cr->bestfit->m2lnL, cr->chi2);
  if (!gsl_finite(f))
    f = GSL_POSINF;
  return f;
}

static gdouble
nc_confidence_region_df (gdouble x, gpointer p)
{
  NcConfidenceRegion *cr = (NcConfidenceRegion *)p;
  gdouble df_data[2] = {0.0, 0.0};
  gdouble res = 0.0;
  NcmVector *df = ncm_vector_new_data_static (df_data, cr->n, 1);
  guint d1, d2;
  if (cr->inv_order) { d1 = 1; d2 = 0; }
  else { d1 = 0; d2 = 1; }
  if(!ncm_fit_cr_step (cr, x))
  {
    //    printf("\n# Failed to fit (df)\n#");
    return GSL_POSINF;
  }

  switch (cr->search_type)
  {
    case NC_CONFIDENCE_REGION_SEARCH_1D:
      cr->constrained->grad.m2lnL_grad (cr->constrained, df);
      res = df_data[0];
      break;
    case NC_CONFIDENCE_REGION_SEARCH_ELLIPTIC_RADIUS:
    {
      gdouble dxp_data[2], dx_data[2];
      gsl_vector_view dxp_view = gsl_vector_view_array (dxp_data, 2);
      gsl_vector_view dx_view = gsl_vector_view_array (dx_data, 2);
      cr->constrained->grad.m2lnL_grad (cr->constrained, df);
      dxp_data[0] = sqrt(NC_CR_COVAR_EV (cr, 0)) * cos (cr->theta);
      dxp_data[1] = sqrt(NC_CR_COVAR_EV (cr, 1)) * sin (cr->theta);
      gsl_blas_dgemv (CblasTrans, 1.0, NCM_MATRIX_GSL (cr->covar_orto), &dxp_view.vector, 0.0, &dx_view.vector);
      res = df_data[d1] * dx_data[0] + df_data[d2] * dx_data[1];
    }
      break;
    case NC_CONFIDENCE_REGION_SEARCH_ELLIPTIC_ANGLE:
    {
      gdouble dxp_data[2], dx_data[2];
      gsl_vector_view dxp_view = gsl_vector_view_array (dxp_data, 2);
      gsl_vector_view dx_view = gsl_vector_view_array (dx_data, 2);
      cr->constrained->grad.m2lnL_grad (cr->constrained, df);
      dxp_data[0] = -sqrt(NC_CR_COVAR_EV (cr, 0)) * cr->r * sin (x);
      dxp_data[1] =  sqrt(NC_CR_COVAR_EV (cr, 1)) * cr->r * cos (x);
      gsl_blas_dgemv (CblasTrans, 1.0, NCM_MATRIX_GSL (cr->covar_orto), &dxp_view.vector, 0.0, &dx_view.vector);
      res = df_data[d1] * dx_data[0] + df_data[d2] * dx_data[1];
    }
      break;
    case NC_CONFIDENCE_REGION_SEARCH_CONSTRAIN_1D:
      res = -2.0 * (ncm_mset_func_eval1 (cr->prior[0].func, cr->constrained->mset, cr->prior[0].z) - cr->prior[0].mean) / gsl_pow_2 (cr->prior[0].sigma);
      break;
  }
  //g_debug ("df %g", res);
  ncm_vector_free (df);
  if (!gsl_finite(res))
    res = GSL_POSINF;
  return res;
}

static gdouble
nc_confidence_region_numdiff_df (gdouble x, gpointer p)
{
  gsl_function F;
  gdouble res, err;
  F.function = &nc_confidence_region_f;
  F.params = p;
  res = ncm_numdiff_1 (&F, x, x * 1e-5, &err);
  //gsl_deriv_central (&F, x, x * 1e-5, &res, &err);
  //gsl_deriv_forward (&F, x, x * 1e-5, &res, &err);

  return res;
}

static void
nc_confidence_region_numdiff_fdf (gdouble x, gpointer p, gdouble *y, gdouble *dy)
{
  *dy = nc_confidence_region_numdiff_df (x, p);
  *y = nc_confidence_region_f (x, p);
  return;
}

static void
nc_confidence_region_fdf (gdouble x, gpointer p, gdouble *y, gdouble *dy)
{
  NcConfidenceRegion *cr = (NcConfidenceRegion *)p;
  gdouble df_data[2] = {0.0, 0.0};
  NcmVector *df = ncm_vector_new_data_static (df_data, cr->n, 1);
  guint d1, d2;
  if (cr->inv_order) { d1 = 1; d2 = 0; }
  else { d1 = 0; d2 = 1; }
  if(!ncm_fit_cr_step (cr, x))
  {
    //    printf("\n# Failed to fit (fdf)\n#");
    *y = GSL_POSINF;
    *dy = GSL_POSINF;
    return;
  }

  *y = cr->constrained->m2lnL - (cr->bestfit->m2lnL + cr->chi2);

  switch (cr->search_type)
  {
    case NC_CONFIDENCE_REGION_SEARCH_1D:
      cr->constrained->grad.m2lnL_grad (cr->constrained, df);
      *dy = df_data[0];
      break;
    case NC_CONFIDENCE_REGION_SEARCH_ELLIPTIC_RADIUS:
    {
      gdouble dxp_data[2], dx_data[2];
      gsl_vector_view dxp_view = gsl_vector_view_array (dxp_data, 2);
      gsl_vector_view dx_view = gsl_vector_view_array (dx_data, 2);
      cr->constrained->grad.m2lnL_grad (cr->constrained, df);
      dxp_data[0] = sqrt(NC_CR_COVAR_EV (cr, 0)) * cos (cr->theta);
      dxp_data[1] = sqrt(NC_CR_COVAR_EV (cr, 1)) * sin (cr->theta);
      gsl_blas_dgemv (CblasTrans, 1.0, NCM_MATRIX_GSL (cr->covar_orto), &dxp_view.vector, 0.0, &dx_view.vector);
      *dy = df_data[d1] * dx_data[0] + df_data[d2] * dx_data[1];
    }
      break;
    case NC_CONFIDENCE_REGION_SEARCH_ELLIPTIC_ANGLE:
    {
      gdouble dxp_data[2], dx_data[2];
      gsl_vector_view dxp_view = gsl_vector_view_array (dxp_data, 2);
      gsl_vector_view dx_view = gsl_vector_view_array (dx_data, 2);
      cr->constrained->grad.m2lnL_grad (cr->constrained, df);
      dxp_data[0] = -sqrt (NC_CR_COVAR_EV (cr, 0)) * cr->r * sin (x);
      dxp_data[1] =  sqrt (NC_CR_COVAR_EV (cr, 1)) * cr->r * cos (x);
      gsl_blas_dgemv (CblasTrans, 1.0, NCM_MATRIX_GSL (cr->covar_orto), &dxp_view.vector, 0.0, &dx_view.vector);
      *dy = df_data[d1] * dx_data[0] + df_data[d2] * dx_data[1];
    }
      break;
    case NC_CONFIDENCE_REGION_SEARCH_CONSTRAIN_1D:
      *dy = -2.0 * (ncm_mset_func_eval1 (cr->prior[0].func, cr->constrained->mset, cr->prior[0].z) - cr->prior[0].mean) / gsl_pow_2 (cr->prior[0].sigma);
      break;
  }
  //g_debug ("f %g", *y);
  //g_debug ("df %g", *dy);
  ncm_vector_free (df);
  if (!gsl_finite(*y))
    *y = GSL_POSINF;
  if (!gsl_finite(*dy))
    *dy = GSL_POSINF;

  return;
}
