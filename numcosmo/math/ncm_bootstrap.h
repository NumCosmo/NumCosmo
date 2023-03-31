/* -*- Mode: C; indent-tabs-mode: nil; c-basic-offset: 2; tab-width: 2 -*-  */
/***************************************************************************
 *            ncm_bootstrap.h
 *
 *  Fri August 16 11:09:19 2013
 *  Copyright  2013  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * ncm_bootstrap.h
 * Copyright (C) 2013 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_randist.h>
#endif /* NUMCOSMO_GIR_SCAN */

G_BEGIN_DECLS

#define NCM_TYPE_BOOTSTRAP (ncm_bootstrap_get_type ())

G_DECLARE_FINAL_TYPE (NcmBootstrap, ncm_bootstrap, NCM, BOOTSTRAP, GObject)

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

void ncm_bootstrap_resample (NcmBootstrap *bstrap, NcmRNG *rng);
void ncm_bootstrap_remix (NcmBootstrap *bstrap, NcmRNG *rng);
guint ncm_bootstrap_get (NcmBootstrap *bstrap, guint i);
GArray *ncm_bootstrap_get_sortncomp (NcmBootstrap *bstrap);
gboolean ncm_bootstrap_is_init (NcmBootstrap *bstrap);

#define NCM_BOOTSTRAP_RNG_NAME "bootstrap"

G_END_DECLS

#endif /* _NCM_BOOTSTRAP_H_ */
