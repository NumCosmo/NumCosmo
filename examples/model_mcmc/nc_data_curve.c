/***************************************************************************
 *            nc_data_curve.c
 *
 *  Wed July 12 13:45:57 2017
 *  Copyright  2017  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * nc_data_curve.c
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
 * SECTION:nc_data_curve
 * @title: NcDataCurve
 * @short_description: Diagonal Gaussian data for #NcCurve models
 * 
 * This is a implementation of #NcmDataGaussDiag desgined to 
 * contain Gaussian distributed data connected to #NcCurve models.
 * 
 * (This is a gtk-doc comment, which always start with two **)
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "nc_data_curve.h"
#include "nc_curve.h"

/* 
 * The #NcmDataGaussDiag object already have a vector for the function
 * values and standard deviations, which are general for any diagonal
 * gaussian data. Here our data points are just the value of the function
 * at a point $x_i$ plus a Gaussian noise $\delta_i$ with variance $\sigma_i^2$.
 * For this reason, we need an additional vector to store our $x_i$
 * values.
 * 
 * For this reason we create a property and an entry in our private structure.
 * 
 */
enum
{
  PROP_0,
  PROP_X,
  PROP_SIZE,
};

struct _NcDataCurvePrivate
{
  NcmVector *xv;
};

/*
 * Here we define the object GType using the macro G_DEFINE_TYPE.
 * This macro basically defines the function nc_data_curve_get_type (void) and
 * everything necessary to define an object in the GLib type system. The
 * since it is not an ABSTRACT object, it creates a GType that can be instantiated.
 * 
 * The last argument provides the GType of the parent object.
 * 
 */
G_DEFINE_TYPE (NcDataCurve, nc_data_curve, NCM_TYPE_DATA_GAUSS_DIAG);

static void
nc_data_curve_init (NcDataCurve *data_curve)
{
  /* 
   * The first step is the creation of the private structure, all allocation and 
   * de-allocation is automatically performed by the GObject framework.
   */
  data_curve->priv     = G_TYPE_INSTANCE_GET_PRIVATE (data_curve, NC_TYPE_DATA_CURVE, NcDataCurvePrivate);

  /*
   * Here we initialize all structure member to null/zero.
   * 
   */
  data_curve->priv->xv = NULL;
  
}

static void
_nc_data_curve_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcDataCurve *data_curve = NC_DATA_CURVE (object);
  g_return_if_fail (NC_IS_DATA_CURVE (object));

  switch (prop_id)
  {
    /*
     * Here we use the ncm_vector_substitute to substitute the vector in the 
     * struct, the last argument impose that both vectors must have the same size.
     * 
     */
    case PROP_X:
      ncm_vector_substitute (&data_curve->priv->xv, g_value_get_object (value), TRUE);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_data_curve_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcDataCurve *data_curve = NC_DATA_CURVE (object);
  g_return_if_fail (NC_IS_DATA_CURVE (object));

  switch (prop_id)
  {
    /*
     * Give back a reference of our vector.
     * 
     */
    case PROP_X:
      g_value_set_object (value, data_curve->priv->xv);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

/*
 * Here we must de-allocate any memory allocated *inside* gobject framework,
 * i.e., we must unref any outside object contained in our object.
 * 
 */
static void
_nc_data_curve_dispose (GObject *object)
{
  NcDataCurve *data_curve = NC_DATA_CURVE (object);

  /*
   * Here we are using ncm_vector_clear to decrease the reference count of
   * xv, after this takes place in the first time the pointer 
   * data_curve->priv->xv will be NULL and any subsequent call to dispose
   * will do nothing.
   * 
   */
  ncm_vector_clear (&data_curve->priv->xv);
  
  /* Chain up : end */
  G_OBJECT_CLASS (nc_data_curve_parent_class)->dispose (object);
}

/*
 * Here we must de-allocate any memory allocated *outside* gobject framework.
 * Nothing to do!
 * 
 */
static void
_nc_data_curve_finalize (GObject *object)
{
  
  /* Chain up : end */
  G_OBJECT_CLASS (nc_data_curve_parent_class)->finalize (object);
}

/*
 * Prototype of the implementation.
 * The #NcmDataGaussDiag object has four virtual functions
 * that must always be implemented by their children:
 * 
 */

/*
 * This function is responsible for getting the model in the model set `mset'
 * and calculating the theoretical values (or the mean values) of each Gaussian 
 * point.
 * 
 */
void _nc_data_curve_mean_func (NcmDataGaussDiag *diag, NcmMSet *mset, NcmVector *vp);

/*
 * This function is responsible for getting the model in the model set `mset'
 * and calculating the standard deviation of each Gaussian point. It must return
 * TRUE if the values were updated or FALSE it they are the same as in the last 
 * call.
 * 
 * This function is optional, it must be implemented only if the variance is not fixed.
 * 
 */
/* gboolean _nc_data_curve_sigma_func (NcmDataGaussDiag *diag, NcmMSet *mset, NcmVector *var); */

/*
 * This function is called anytime the data object must change size, i.e.,
 * when the number of Gaussian points change. Since the parent has some
 * of the data structures we must always chain up at the beggining of
 * the function. See the implementation below.
 * 
 */
void _nc_data_curve_set_size (NcmDataGaussDiag *diag, guint np);

/*
 * This function must return the number of points in the object, 
 * it has a default implementation which just return the size.
 * We do not need to reimplement it here.
 * 
 */
/* guint _nc_data_curve_get_size (NcmDataGaussDiag *diag); */

/*
 * At _class_init we will define all properties and parameters we should
 * also include a default implementation for our virtual function `f'.
 * 
 */
static void
nc_data_curve_class_init (NcDataCurveClass *klass)
{
  GObjectClass* object_class        = G_OBJECT_CLASS (klass);
  NcmDataGaussDiagClass *diag_class = NCM_DATA_GAUSS_DIAG_CLASS (klass); 

  /* 
   * Tells the type system that we have a private struct
   * 
   */
  g_type_class_add_private (klass, sizeof (NcDataCurvePrivate));

  /*
   * The usual GObject hooks must be assigned to GObjectClass.
   * 
   * Note that differently from the NcmModel here we put the
   * set/get functions directly in the GObject structure,
   * the NcmModel is a exception!
   * 
   */    
  object_class->set_property = &_nc_data_curve_set_property;
  object_class->get_property = &_nc_data_curve_get_property;
  object_class->dispose      = &_nc_data_curve_dispose;
  object_class->finalize     = &_nc_data_curve_finalize;

  /*
   * This register our x vector as a property called `x'.
   * 
   */
  g_object_class_install_property (object_class,
                                   PROP_X,
                                   g_param_spec_object ("x",
                                                        NULL,
                                                        "x vector",
                                                        NCM_TYPE_VECTOR,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  
  /*
   * Here we assign our implementation of the NcmDataGaussDiag methods.
   * 
   */
  diag_class->mean_func  = &_nc_data_curve_mean_func;
  diag_class->set_size   = &_nc_data_curve_set_size;
}

void 
_nc_data_curve_mean_func (NcmDataGaussDiag *diag, NcmMSet *mset, NcmVector *vp)
{
  /* Recasting to our local object type */
  NcDataCurve *data_curve = NC_DATA_CURVE (diag);
  /*
   * The mset object contains all models used in a analysis, in each data
   * object the functions should get from the mset the models necessary for
   * their calculations.
   * 
   * 
   */
  NcCurve *curve          = NC_CURVE (ncm_mset_peek (mset, nc_curve_id ()));

  /*
   * Asserting that mset contain a model of the type nc_curve_id().
   * 
   */
  g_assert (curve != NULL);

  {
    /* Getting the number of points */
    const guint len = ncm_data_gauss_diag_get_size (diag);
    guint i;

    for (i = 0 ; i < len; i++)
    {
      /* Getting the value of x_i */
      const gdouble x_i    = ncm_vector_get (data_curve->priv->xv, i);
      /*
       * Here we have the important abstraction layer, we call the function
       * nc_curve_f() to get the mean values, does not matter which is
       * the actual implementation of NcCurve or which are its parameters.
       * 
       */
      const gdouble f_th_i = nc_curve_f (curve, x_i);

      ncm_vector_set (vp, i, f_th_i);
    }
  }  
}

void 
_nc_data_curve_set_size (NcmDataGaussDiag *diag, guint np)
{
  /* Recasting to our local object type */
  NcDataCurve *data_curve = NC_DATA_CURVE (diag);

  /*
   * If the new size is zero or if it is different of the actual value
   * we clear the xv vector.
   * 
   */
  if ((np == 0) || (np != diag->np))
    ncm_vector_clear (&data_curve->priv->xv);

  /*
   * If the requested size is different from zero and it is different 
   * from the current, allocate a new vector.
   * 
   */
  if ((np != 0) && (np != diag->np))
    data_curve->priv->xv = ncm_vector_new (np);

  /*
   * Here we must chain up so the parent can update its internal structures
   * too.
   * 
   */
  /* Chain up : end */
  NCM_DATA_GAUSS_DIAG_CLASS (nc_data_curve_parent_class)->set_size (diag, np);
}

/*
 * This defines the constructor for #NcDataCurve 
 * 
 */ 
/**
 * nc_data_curve_new:
 * @np: number of data points
 * 
 * Creates a new #NcDataCurve for a set of @np points.
 * 
 * Returns: (transfer full): a new #NcDataCurve
 */ 
NcDataCurve *
nc_data_curve_new (const guint np)
{
  /*
   * Here we use the actual GObject constructor to instantiate
   * an object of GType NC_TYPE_DATA_CURVE. Note that it 
   * inherits the #NcmDataGaussDiag properties xl and xu, and they must
   * be assigned during construction.
   * 
   */
  NcDataCurve *data_curve = g_object_new (NC_TYPE_DATA_CURVE,
                                          "n-points", np,
                                          NULL);
  return data_curve;
}

/*
 * The reference increasing function
 * 
 */
/**
 * nc_data_curve_ref:
 * @data_curve: a #NcDataCurve
 * 
 * Increase reference count by one.
 * 
 * Returns: (transfer full): @data_curve.
 */
NcDataCurve *
nc_data_curve_ref (NcDataCurve *data_curve)
{
  return g_object_ref (data_curve);
}

/*
 * The reference decreasing function
 * 
 */
/**
 * nc_data_curve_free:
 * @data_curve: a #NcDataCurve
 * 
 * Decrease reference count by one.
 * 
 */
void 
nc_data_curve_free (NcDataCurve *data_curve)
{
  g_object_unref (data_curve);
}

/*
 * This function decreases the reference count
 * by one only if *data_curve != NULL, and in that case
 * it sets *data_curve to NULL after decreasing the 
 * reference count. It is useful to use these 
 * functions in the dispose hooks.
 * 
 * 
 */
/**
 * nc_data_curve_clear:
 * @data_curve: a #NcDataCurve
 * 
 * Decrease reference count by one if *@data_curve != NULL
 * and sets @data_curve to NULL. 
 * 
 */
void 
nc_data_curve_clear (NcDataCurve **data_curve)
{
  g_clear_object (data_curve);
}

/* 
 * Here we create a simple method that can be used to initialize 
 * our data object.
 * 
 */
/**
 * nc_data_curve_init_data:
 * @data_curve: a #NcDataCurve
 * @xv: a #NcmVector
 * @yv: (nullable): a #NcmVector
 * @sigmav: a #NcmVector
 * 
 * Initialize the @data_curve object using the data in the vectors
 * @xv, @yv, and @sigmav. The @yv vector is optional. 
 * 
 */
void 
nc_data_curve_init_data (NcDataCurve *data_curve, NcmVector *xv, NcmVector *yv, NcmVector *sigmav)
{
  /* Getting a recast version of the same object to use its methods */
  NcmDataGaussDiag *diag = NCM_DATA_GAUSS_DIAG (data_curve);
  
  /*
   * Checking the consistency of the input considering yv as optional.
   * 
   */
  const guint xv_len     = ncm_vector_len (xv);
  const guint yv_len     = (yv != NULL) ? ncm_vector_len (yv) : 0;
  const guint sigmav_len = ncm_vector_len (sigmav);

  g_assert_cmpuint (xv_len, ==, sigmav_len);
  if (yv != NULL)
    g_assert_cmpuint (xv_len, ==, yv_len);

  /*
   * Setting the required size.
   * 
   */
  ncm_data_gauss_diag_set_size (diag, xv_len);

  /*
   * Call the GObject function used to set properties.
   * 
   */
  if (yv != NULL)
  {
    g_object_set (data_curve,
                  "x",     xv,
                  "mean",  yv, 
                  "sigma", sigmav,
                  NULL);

    /*
     * Since we set all necessary properties we can change the
     * object status to `initialized'.
     * 
     */
    ncm_data_set_init (NCM_DATA (data_curve), TRUE);
  }
  else
  {
    g_object_set (data_curve,
                  "x",     xv,
                  "sigma", sigmav,
                  NULL);
  }
}
