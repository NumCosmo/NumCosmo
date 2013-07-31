/***************************************************************************
 *            ncm_bootstrap.h
 *
 *  Fri August 16 11:09:19 2013
 *  Copyright  2013  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * ncm_bootstrap.h
 * Copyright (C) 2013 Sandro Dias Pinto Vitenti <sandro@isoftware.com.br>
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

#ifndef _NCM_BOOTSTRAP_H_
#define _NCM_BOOTSTRAP_H_


#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_rng.h>
#include <gsl/gsl_randist.h>

G_BEGIN_DECLS

#define NCM_TYPE_BOOTSTRAP             (ncm_bootstrap_get_type ())
#define NCM_BOOTSTRAP(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NCM_TYPE_BOOTSTRAP, NcmBootstrap))
#define NCM_BOOTSTRAP_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NCM_TYPE_BOOTSTRAP, NcmBootstrapClass))
#define NCM_IS_BOOTSTRAP(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NCM_TYPE_BOOTSTRAP))
#define NCM_IS_BOOTSTRAP_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NCM_TYPE_BOOTSTRAP))
#define NCM_BOOTSTRAP_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NCM_TYPE_BOOTSTRAP, NcmBootstrapClass))

typedef struct _NcmBootstrapClass NcmBootstrapClass;
typedef struct _NcmBootstrap NcmBootstrap;

struct _NcmBootstrap
{
  /*< private >*/
  GObject parent_instance;
  guint fsize;
  guint bsize;
  GArray *bootstrap_index;
  GArray *increasing_index;
};

struct _NcmBootstrapClass
{
  /*< private >*/
  GObjectClass parent_class;
};

GType ncm_bootstrap_get_type (void) G_GNUC_CONST;

NcmBootstrap *ncm_bootstrap_new (void);
NcmBootstrap *ncm_bootstrap_sized_new (guint fsize);
NcmBootstrap *ncm_bootstrap_full_new (guint fsize, guint bsize);
NcmBootstrap *ncm_bootstrap_ref (NcmBootstrap *bstrap);
void ncm_bootstrap_free (NcmBootstrap *bstrap);
void ncm_bootstrap_clear (NcmBootstrap **bstrap);

void ncm_bootstrap_set_fsize (NcmBootstrap *bstrap, guint fsize);
guint ncm_bootstrap_get_fsize (NcmBootstrap *bstrap);
void ncm_bootstrap_set_bsize (NcmBootstrap *bstrap, guint bsize);
guint ncm_bootstrap_get_bsize (NcmBootstrap *bstrap);

G_INLINE_FUNC void ncm_bootstrap_resample (NcmBootstrap *bstrap);
G_INLINE_FUNC guint ncm_bootstrap_get (NcmBootstrap *bstrap, guint i);

G_END_DECLS

#endif /* _NCM_BOOTSTRAP_H_ */

#ifndef _NCM_BOOTSTRAP_INLINE_H_
#define _NCM_BOOTSTRAP_INLINE_H_
#ifdef NUMCOSMO_HAVE_INLINE

G_BEGIN_DECLS

#define NCM_BOOTSTRAP_RNG_NAME "bootstrap"

G_INLINE_FUNC void 
ncm_bootstrap_resample (NcmBootstrap *bstrap)
{
  NcmRNG *rng = ncm_rng_pool_get (NCM_BOOTSTRAP_RNG_NAME);
  gpointer bdata = bstrap->bootstrap_index->data;
  gpointer idata = bstrap->increasing_index->data;;
  const gsize fsize = bstrap->fsize;
  const gsize bsize = bstrap->bsize;
  const gsize element_size = g_array_get_element_size (bstrap->bootstrap_index);
  ncm_rng_lock (rng);
  gsl_ran_sample (rng->r, bdata, bsize, idata, fsize, element_size);
  ncm_rng_unlock (rng);
  ncm_rng_free (rng);
}

G_INLINE_FUNC guint 
ncm_bootstrap_get (NcmBootstrap *bstrap, guint i)
{
  return g_array_index (bstrap->bootstrap_index, guint, i);
}

G_END_DECLS

#endif /* NUMCOSMO_HAVE_INLINE */
#endif /* _NCM_BOOTSTRAP_INLINE_H_ */
