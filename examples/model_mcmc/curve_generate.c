/***************************************************************************
 *            curve_generate.c
 *
 *  Wed July 12 17:01:40 2017
 *  Copyright  2017  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * curve_generate.c
 *
 * Copyright (C) 2017 - Sandro Dias Pinto Vitenti
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */

/*
 * We need to include our new object headers after the 
 * main NumCosmo header.
 * 
 */
#include <numcosmo/numcosmo.h>
#include "nc_curve.h"
#include "nc_curve_linear.h"
#include "nc_data_curve.h"

gint
main (gint argc, gchar *argv[])
{
  /*
   * Here we declare the variables that contain the options
   * and set their default values.
   * 
   */ 
  gint data_len           = 10;
  gint64 seed             = -1;
  gdouble sigma           = 1.0;
  gdouble xl              = 0.0;
  gdouble xu              = 1.0;
  gchar *out              = NULL;
  gchar *curve_name       = NULL;
  NcDataCurve *data_curve = NULL;

  /*
   * The code below uses the GOption framework to deal with the
   * program options.
   * 
   */
  GError *error = NULL;
  GOptionContext *context;
  GOptionEntry entries[] =
  {
    { "length",         'l', 0, G_OPTION_ARG_INT,    &data_len,  "Number of data points to generate.",          NULL         },
    { "seed",           'S', 0, G_OPTION_ARG_INT64,  &seed,      "Pseudo-random number generator seed `val'.",  "val"        },
    { "xl",             'L', 0, G_OPTION_ARG_DOUBLE, &xl,        "Sets the lower bound to `val'.",              "val"        },
    { "xu",             'U', 0, G_OPTION_ARG_DOUBLE, &xu,        "Sets the upper bound to `val'.",              "val"        },
    { "sigma",          's', 0, G_OPTION_ARG_DOUBLE, &sigma,     "Sets all standard deviations to `val'.",      "val"        },
    { "curve",          'c', 0, G_OPTION_ARG_STRING, curve_name, "Uses curve object named `curve_name'.",       "curve_name" },
    { "out",            'o', 0, G_OPTION_ARG_STRING, &out,       "Save data object to `filename'.",             "filename"   },
    { NULL }
  };

  context = g_option_context_new ("- generate a data set sampling.");
  g_option_context_set_summary (context, "NcDataCurve sample generator");
  g_option_context_set_description (context, "NcDataCurve sample generator");

  g_option_context_add_main_entries (context, entries, NULL);

  if (!g_option_context_parse (context, &argc, &argv, &error))
  {
    g_print ("option parsing failed: %s\n", error->message);
    exit (1);
  }

  /*
   * Checking if the options have sane values.
   * 
   */
  if (out == NULL)
  {
    g_message ("Choose the output filename, using --out/-o.\n");
    exit (1);
  }
  if (xu <= xl)
  {
    g_message ("Invalid interval [% 22.15g, % 22.15g], use --xl/-L and --xu/-U to set the interval lower/upper bounds.\n",
               xl, xu);
    exit (1);    
  }
  if (sigma <= 0.0)
  {
    g_message ("Invalid standard deviation `% 22.15g', use --sigma/-s to set its value.\n",
               sigma);
    exit (1);    
  }
  if (data_len <= 0)
  {
    g_message ("Invalid number of points `%d', use --length/-l to set its value.\n",
               data_len);
    exit (1);
  }

  /*
   * After initializing the library we need to register our new objects.
   * If these objects were incorporated in the main NumCosmo project
   * these registry would take place inside ncm_cfg_init().
   * 
   */
  ncm_cfg_init ();
  ncm_cfg_register_obj (NC_TYPE_CURVE);  
  ncm_cfg_register_obj (NC_TYPE_CURVE_LINEAR);
  ncm_cfg_register_obj (NC_TYPE_DATA_CURVE);

  /*
   * First of all we initialize the NcDataCurve object. For that we 
   * use the sigma options to set a equal value for all $\sigma_i$
   * and create an equally spaced grid between xl and xu with data_len
   * knots.
   * 
   */
  {
    NcmVector *xv     = ncm_vector_new (data_len);
    NcmVector *sigmav = ncm_vector_new (data_len);
    guint i;

    data_curve = nc_data_curve_new (data_len);
    
    ncm_vector_set_all (sigmav, sigma);

    for (i = 0; i < data_len; i++)
    {
      const gdouble x_i = xl + (xu - xl) / (data_len - 1.0) * i;
      ncm_vector_set (xv, i, x_i);
    }

    nc_data_curve_init_data (data_curve, xv, NULL, sigmav);

    ncm_vector_free (xv);
    ncm_vector_free (sigmav);
  }  

  /*
   * Now that our data object has everything but the actual data points
   * we can use the _resample method of NcmData to generate a sample.
   * 
   * 
   */
  {
    /* If the user has passed a NcCurve model we use it, otherwise we used NcCurveLinear as default. */
    NcCurve *curve    = (curve_name != NULL) ? nc_curve_new_from_name (curve_name) : NC_CURVE (nc_curve_linear_new (xl, xu));

    /* Recasting it to NcmModel to use its methods */
    NcmModel *model   = NCM_MODEL (curve);

    /* Creating a NcmMSet model to contain our model */
    NcmMSet *mset     = ncm_mset_new (model, NULL);

    /* 
     * We need a random number generator, here we create one with a random seed or with the one given by the user
     * 
     */
    NcmRNG *rng       = (seed == -1) ? ncm_rng_new (NULL) : ncm_rng_seeded_new (NULL, seed);

    /* We need a NcmSerialize object to transform our objects to strings and save it to disk */
    NcmSerialize *ser = ncm_serialize_new (NCM_SERIALIZE_OPT_NONE);

    /*
     * This method overwrite the data in data_curve (currently uninitalized)
     * with a random set obtained by the likelihood and the random number 
     * generator `rng'.
     * 
     */
    ncm_data_resample (NCM_DATA (data_curve), mset, rng);

    /*
     * The method ncm_serialize_to_file serializes our data object and saves
     * it to disk using the filename in the variable `out'. One can also 
     * use ncm_serialize_to_binfile to save it in a binary format, which is
     * useful for large data objects.
     * 
     */
    ncm_serialize_to_file (ser, G_OBJECT (data_curve), out);

    ncm_mset_free (mset);
    ncm_rng_free (rng);
  }
  
  nc_data_curve_free (data_curve);
}
