/***************************************************************************
 *            print_data.c
 *
 *  Sat Sep 10 18:55:42 2011
 *  Copyright  2011  Mariana Penna Lima
 *  <pennalima@gmail.com>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Mariana Penna Lima 2012 <pennalima@gmail.com>
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
 * SECTION:print_data
 * @title: Mass Function Data Printing
 * @short_description: FIXME
 *
 * FIXME
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "lss/print_data.h"
#include "lss/nc_cluster_abundance.h"
#include "nc_data_cluster_ncount.h"

/**
 * nc_mass_funtion_print:
 * @ca_unbinned: FIXME
 * @model: a NcHICosmo
 * @out: FIXME
 * @header: FIXME
 *
 * FIXME
 *
 */
void
nc_mass_function_print (NcmData *data, NcHICosmo *model, FILE *out, gchar *header)
{
  NcDataClusterNCount *ncount = NC_DATA_CLUSTER_NCOUNT (data);
  NcClusterAbundance *cad = ncount->cad;
  gint i, j;
  gint nbins_M = 100;
  gint nbins_z = 20;
  gsl_vector *lnM_nodes = gsl_vector_alloc (nbins_M);
  gsl_vector *z_nodes = gsl_vector_alloc (nbins_z);

  for (i = 0; i < nbins_M; i++)
  {
    gdouble lnM = cad->lnMi + (cad->lnMf - cad->lnMi) / (nbins_M - 1.0) * i;
    gsl_vector_set (lnM_nodes, i, lnM);

    printf ("lnM = %5.5g Me = %5.5g M10 = %5.5g\n", gsl_vector_get (lnM_nodes, i), exp(gsl_vector_get (lnM_nodes, i)), pow(10, gsl_vector_get (lnM_nodes, i)));
  }

  for (j = 0; j < nbins_z; j++)
  {
    gdouble z = cad->zi + (cad->zf - cad->zi) / (nbins_z - 1.0) * j;
    gsl_vector_set (z_nodes, j, z);

    printf ("z = %5.5g\n", gsl_vector_get (z_nodes, j));
  }

  if (header != NULL)
    fprintf (out, "# %s\n# ", header);
  else
    fprintf (out, "# ");

  gsl_histogram2d *h = nc_data_cluster_ncount_hist_lnM_z (data, lnM_nodes, z_nodes);


  fprintf (out, "# z M N/(logM * V) (catalog) dn/dlog10M (theory) Nmi(catalog)(abundance in bins of redshift and mass) Nmi(theory)\n");
  for (j = 0; j < nbins_z - 1; j++)
  {
    gdouble zm, dz, zl, zu, V;
    gsl_histogram2d_get_yrange (h, j, &zl, &zu);
    zm = (zu + zl) / 2.0;
    dz = (zu - zl);
    V = nc_mass_function_dv_dzdomega (cad->mfp, model, zm) * cad->mfp->area_survey * dz;
    for (i = 0; i < nbins_M; i++)
    {
      gdouble ln_ml, ln_mu, Mm, lnMm, log_mu, log_ml;
      gdouble Nmi = gsl_histogram2d_get (h, i, j);
      gdouble dndlog10M, ca_M;
      gsl_histogram2d_get_xrange (h, i, &ln_ml, &ln_mu);
      lnMm = (ln_ml + ln_mu) / 2.0;
      Mm = exp (lnMm);
      log_mu = log10 (exp(ln_mu));
      log_ml = log10 (exp(ln_ml));
      dndlog10M = M_LN10 * nc_mass_function_dn_dlnm (cad->mfp, model, lnMm, zm);
      ca_M = (log_mu - log_ml) * V * dndlog10M;

      //printf ("log-mu = %5.5g log-ml = %5.5g\n", log_mu, log_ml);
      fprintf (out, "% 6.6g % 6.6e % 6.6g % 6.6g % 6.6g % 6.6g\n", zm, Mm, Nmi / ((log_mu - log_ml) * V), dndlog10M, Nmi, ca_M);
    }
    fprintf (out, "\n\n");
  }

  gsl_vector_free (lnM_nodes);
  gsl_vector_free (z_nodes);
  gsl_histogram2d_free (h);
}
