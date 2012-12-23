/***************************************************************************
 *            ncm_mset.c
 *
 *  Fri May 25 09:37:54 2012
 *  Copyright  2012  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2012 <sandro@isoftware.com.br>
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
 * SECTION:ncm_mset
 * @title: A Set of NcmModels
 * @short_description: Object representing a set of different NcmModel objects
 *
 * FIXME
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_mset.h"
#include "math/ncm_cfg.h"

G_DEFINE_TYPE (NcmMSet, ncm_mset, G_TYPE_OBJECT);
G_DEFINE_BOXED_TYPE (NcmMSetPIndex, ncm_mset_pindex, ncm_mset_pindex_dup, ncm_mset_pindex_free);

/**
 * ncm_mset_pindex_new:
 * @mid: FIXME
 * @pid: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
NcmMSetPIndex *
ncm_mset_pindex_new (NcmModelID mid, guint pid)
{
  NcmMSetPIndex *pi = g_slice_new0 (NcmMSetPIndex);
  pi->mid = mid;
  pi->pid = pid;
    
  return pi;
}

/**
 * ncm_mset_pindex_dup:
 * @pi: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
NcmMSetPIndex *
ncm_mset_pindex_dup (NcmMSetPIndex *pi)
{
  NcmMSetPIndex *new_pi = ncm_mset_pindex_new (pi->mid, pi->pid);
  return new_pi;
}

/**
 * ncm_mset_pindex_free:
 * @pi: FIXME
 *
 * FIXME
 *
 */
void
ncm_mset_pindex_free (NcmMSetPIndex *pi)
{
  g_slice_free (NcmMSetPIndex, pi);
}

/**
 * ncm_mset_empty_new:
 *
 * FIXME
 *
 * Returns: FIXME
 */
NcmMSet *
ncm_mset_empty_new (void)
{
  return g_object_new (NCM_TYPE_MSET, NULL);
}

/**
 * ncm_mset_new:
 * @model0: FIXME
 * @...: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
NcmMSet *
ncm_mset_new (NcmModel *model0, ...)
{
  NcmMSet *mset;
  va_list ap;

  va_start (ap, model0);
  mset = ncm_mset_newv (model0, ap);
  va_end (ap);

  return mset;
}

/**
 * ncm_mset_newv:
 * @model0: FIXME
 * @ap: FIXME
 *
 * FIXME
 *
 * Returns: (transfer full): FIXME
 */
NcmMSet *
ncm_mset_newv (NcmModel *model0, va_list ap)
{
  NcmMSet *mset = ncm_mset_empty_new ();
  NcmModel *model = NULL;

  ncm_mset_set (mset, model0);

  while ((model = va_arg (ap, NcmModel *)) != NULL)
    ncm_mset_set (mset, model);

  return mset;
}

/**
 * ncm_mset_ref:
 * @mset: a #NcmMSet.
 *
 * FIXME
 *
 * Returns: (transfer full): a new #NcmMSet.
 */
NcmMSet *
ncm_mset_ref (NcmMSet *mset)
{
  return g_object_ref (mset);
}

/**
 * ncm_mset_dup:
 * @mset: a #NcmMSet.
 *
 * FIXME
 *
 * Returns: (transfer full): a new #NcmMSet.
 */
NcmMSet *
ncm_mset_dup (NcmMSet *mset)
{
  NcmMSet *mset_copy = ncm_mset_empty_new ();
  NcmModelID mid;
  for (mid = 0; mid < NCM_MODEL_MAX_ID; mid++)
  {
    if (mset->model[mid] != NULL)
    {
      NcmModel *m_copy = ncm_model_copy (ncm_mset_peek (mset, mid));
      ncm_mset_set (mset_copy, m_copy);
      ncm_model_free (m_copy);
    }
  }
  return mset_copy;
}

/**
 * ncm_mset_copy:
 * @mset: a #NcmMSet.
 *
 * FIXME
 *
 * Returns: (transfer full): a new #NcmMSet.
 */
NcmMSet *
ncm_mset_copy (NcmMSet *mset)
{
  NcmMSet *mset_copy = ncm_mset_empty_new ();
  NcmModelID mid;
  for (mid = 0; mid < NCM_MODEL_MAX_ID; mid++)
  {
    if (mset->model[mid] != NULL)
      ncm_mset_set (mset_copy, ncm_mset_peek (mset, mid));
  }
  return mset_copy;
}

/**
 * ncm_mset_copyto:
 * @mset_src: a #NcmMSet.
 * @mset_dest: a #NcmMSet.
 *
 * FIXME
 *
 */
void
ncm_mset_copyto (NcmMSet *mset_src, NcmMSet *mset_dest)
{
  guint i;
  for (i = 0; i < NCM_MODEL_MAX_ID; i++)
  {
    if (mset_src->model[i] != NULL)
    {
      g_assert (mset_dest->model[i] != NULL);
      ncm_model_copyto (mset_src->model[i], mset_dest->model[i]);
    }
  }
}

/**
 * ncm_mset_free:
 * @mset: FIXME
 *
 * FIXME
 *
 */
void
ncm_mset_free (NcmMSet *mset)
{
  g_object_unref (mset);
}

/**
 * ncm_mset_clear:
 * @mset: FIXME
 *
 * FIXME
 *
 */
void
ncm_mset_clear (NcmMSet **mset)
{
  g_clear_object (mset);
}

/**
 * ncm_mset_peek:
 * @mset: a #NcmMSet.
 * @mid: a #NcmModelID.
 *
 * FIXME
 *
 * Returns: (transfer none): FIXME
 */
NcmModel *
ncm_mset_peek (NcmMSet *mset, NcmModelID mid)
{
  return mset->model[mid];
}

/**
 * ncm_mset_get:
 * @mset: a #NcmMSet.
 * @mid: a #NcmModelID.
 *
 * FIXME
 *
 * Returns: (transfer full): FIXME
 */
NcmModel *
ncm_mset_get (NcmMSet *mset, NcmModelID mid)
{
  g_assert (mset->model[mid]);
  return ncm_model_ref (mset->model[mid]);
}

/**
 * ncm_mset_remove:
 * @mset: a #NcmMSet.
 * @mid: a #NcmModelID.
 *
 * FIXME
 *
 */
void
ncm_mset_remove (NcmMSet *mset, NcmModelID mid)
{
  if (mset->model[mid] != NULL)
  {
    mset->total_len -= ncm_model_len (mset->model[mid]);
    ncm_model_free (mset->model[mid]);
    mset->model[mid] = NULL;
    mset->valid_map = FALSE;
  }
}

/**
 * ncm_mset_set:
 * @mset: FIXME
 * @model: FIXME
 *
 * FIXME
 *
 */
void
ncm_mset_set (NcmMSet *mset, NcmModel *model)
{
  NcmModelID mid = ncm_model_id (model);
  guint model_len = ncm_model_len (model);
  guint i;

  if (mset->model[mid] != NULL)
    ncm_mset_remove (mset, mid);

  mset->model[mid] = ncm_model_ref (model);

  g_array_set_size (mset->fpi_array[mid], model_len);
  for (i = 0; i < model_len; i++)
    g_array_index (mset->fpi_array[mid], gint, i) = -1;

  mset->total_len += model_len;

  mset->valid_map = FALSE;
}

/**
 * ncm_mset_prepare_fparam_map:
 * @mset: FIXME
 *
 * FIXME
 *
 */
void
ncm_mset_prepare_fparam_map (NcmMSet *mset)
{
  NcmModelID mid;
  guint pid;

  mset->fparam_len = 0;
  g_array_set_size (mset->pi_array, 0);

  for (mid = 0; mid < NCM_MODEL_MAX_ID; mid++)
  {
    NcmModel *model = ncm_mset_peek (mset, mid);
    if (model == NULL)
      continue;
    else
    {
      gint n_params = ncm_model_len (model);
      for (pid = 0; pid < n_params; pid++)
      {
        if (ncm_model_param_get_ftype (model, pid) == NCM_PARAM_TYPE_FREE)
        {
          NcmMSetPIndex pi = {mid, pid};
          g_array_append_val (mset->pi_array, pi);
          g_array_index (mset->fpi_array[mid], gint, pid) = mset->fparam_len;
          mset->fparam_len++;
        }
        else
          g_array_index (mset->fpi_array[mid], gint, pid) = -1;
      }
    }
  }
  mset->valid_map = TRUE;
}

/**
 * ncm_mset_total_len:
 * @mset: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
guint
ncm_mset_total_len (NcmMSet *mset)
{
  return mset->total_len;
}

/**
 * ncm_mset_fparam_len:
 * @mset: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
guint
ncm_mset_fparam_len (NcmMSet *mset)
{
  return mset->fparam_len;
}

/**
 * ncm_mset_max_param_name:
 * @mset: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
guint
ncm_mset_max_param_name (NcmMSet *mset)
{
  NcmModelID mid;
  gint i;
  guint name_size = 0;
  for (mid = 0; mid < NCM_MODEL_MAX_ID; mid++)
  {
    NcmModel *model = ncm_mset_peek (mset, mid);
    if (model == NULL)
      continue;
    else
    {
      gint n_params = ncm_model_len (model);
      for (i = 0; i < n_params; i++)
      {
        const gchar *pname = ncm_model_param_name (model, i);
        name_size = GSL_MAX (name_size, strlen (pname));
      }
    }
  }
  return name_size;
}

/**
 * ncm_mset_param_len:
 * @mset: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
guint 
ncm_mset_param_len (NcmMSet *mset)
{
  NcmModelID mid;
  guint len = 0;
  for (mid = 0; mid < NCM_MODEL_MAX_ID; mid++)
  {
    NcmModel *model = ncm_mset_peek (mset, mid);
    if (model == NULL)
      continue;
    else
      len += ncm_model_len (model);
  }
  return len;
}

/**
 * ncm_mset_max_fparam_name:
 * @mset: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
guint
ncm_mset_max_fparam_name (NcmMSet *mset)
{
  guint name_size = 0, i;
  guint free_params_len = ncm_mset_fparams_len (mset);
  for (i = 0; i < free_params_len; i++)
  {
    const gchar *pname = ncm_mset_fparam_name (mset, i);
    name_size = GSL_MAX (name_size, strlen (pname));
  }
  return name_size;
}

/**
 * ncm_mset_max_model_nick:
 * @mset: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
guint
ncm_mset_max_model_nick (NcmMSet *mset)
{
  NcmModelID mid;
  guint nick_size = 0;
  for (mid = 0; mid < NCM_MODEL_MAX_ID; mid++)
  {
    NcmModel *model = ncm_mset_peek (mset, mid);
    if (model == NULL)
      continue;
    else
    {
      const gchar *nick = ncm_model_nick (model);
      nick_size = GSL_MAX (nick_size, strlen (nick));
    }
  }
  return nick_size;
}

/**
 * ncm_mset_pretty_log:
 * @mset: FIXME
 *
 * FIXME
 *
 */
void
ncm_mset_pretty_log (NcmMSet *mset)
{
  NcmModelID mid;
  gint i;
  guint name_size = ncm_mset_max_param_name (mset);
  for (mid = 0; mid < NCM_MODEL_MAX_ID; mid++)
  {
    NcmModel *model = ncm_mset_peek (mset, mid);
    if (model == NULL)
      continue;
    else
    {
      gint n_params = ncm_model_len (model);
      ncm_cfg_msg_sepa ();
      g_message ("# Model[%02d]:\n#   - %s\n", mid, ncm_model_name (model));
      ncm_cfg_msg_sepa ();
      g_message ("# Model parameters\n");
      for (i = 0; i < n_params; i++)
      {
        const gchar *pname = ncm_model_param_name (model, i);
        g_message ("#   - %*s[%02d]: % -20.15g [%s]\n", name_size, pname, i,
                   ncm_model_param_get (model, i),
                   (ncm_model_param_get_ftype (model, i) == NCM_PARAM_TYPE_FREE) ? "FREE" : "FIXED");
      }
    }
  }
}

/**
 * ncm_mset_params_pretty_print:
 * @mset: a #NcmMSet
 * @out: name of the file
 * @header: pointer to the command line
 *
 * This function print the command line (first line, commented), the model nick and parameters' names (second line, commented)
 * and their values indicating if they are fixed or free.
 *
 */
void
ncm_mset_params_pretty_print (NcmMSet *mset, FILE *out, gchar *header)
{
  NcmModelID mid;
  gint i;
  guint name_size = ncm_mset_max_param_name (mset);
  guint nick_size = ncm_mset_max_model_nick (mset);
  name_size = GSL_MAX (name_size, 5);

  if (header != NULL)
    fprintf (out, "# %s\n ", header);
  else
    fprintf (out, "#\n");

  for (mid = 0; mid < NCM_MODEL_MAX_ID; mid++)
  {
    NcmModel *model = ncm_mset_peek (mset, mid);
    if (model == NULL)
      continue;
    else
    {
      const gchar *nick = ncm_model_nick (model);
      gint n_params = ncm_model_len (model);
      for (i = 0; i < n_params; i++)
      {
        const gchar *pname = ncm_model_param_name (model, i);
        fprintf (out, "%*s:%*s %*s % 5.5g\n", nick_size, nick, name_size, pname, name_size,
                 (ncm_model_param_get_ftype (model, i) == NCM_PARAM_TYPE_FREE) ? "FREE" : "FIXED",
                 ncm_model_param_get (model, i));
      }
    }
  }

  return;
}

/**
 * ncm_mset_params_log_vals:
 * @mset: FIXME
 *
 * FIXME
 *
 */
void
ncm_mset_params_log_vals (NcmMSet *mset)
{
  NcmModelID mid;
  gint i;
  g_message ("#");
  for (mid = 0; mid < NCM_MODEL_MAX_ID; mid++)
  {
    NcmModel *model = ncm_mset_peek (mset, mid);
    if (model == NULL)
      continue;
    else
    {
      gint n_params = ncm_model_len (model);
      for (i = 0; i < n_params; i++)
      {
        g_message (" % -20.15g", ncm_model_param_get (model, i));
      }
    }
  }
  g_message ("\n");
}

/**
 * ncm_mset_params_print_vals:
 * @mset: FIXME
 * @out: FIXME
 *
 * FIXME
 *
 */
void
ncm_mset_params_print_vals (NcmMSet *mset, FILE *out)
{
  NcmModelID mid;
  gint i;
  for (mid = 0; mid < NCM_MODEL_MAX_ID; mid++)
  {
    NcmModel *model = ncm_mset_peek (mset, mid);
    if (model == NULL)
      continue;
    else
    {
      gint n_params = ncm_model_len (model);
      for (i = 0; i < n_params; i++)
      {
        fprintf (out, " % -20.15g", ncm_model_param_get (model, i));
      }
    }
  }
  fprintf (out, "\n");
}

/**
 * ncm_mset_params_valid:
 * @mset: FIXME
 *
 * FIXME
 * 
 * Returns: FIXME
 */
gboolean 
ncm_mset_params_valid (NcmMSet *mset)
{
  NcmModelID mid;

  for (mid = 0; mid < NCM_MODEL_MAX_ID; mid++)
  {
    NcmModel *model = ncm_mset_peek (mset, mid);
    if (model == NULL)
      continue;
    else if (!ncm_model_params_valid (model))
      return FALSE;
  }
  return TRUE;
}

/**
 * ncm_mset_cmp:
 * @mset0: FIXME
 * @mset1: FIXME
 * @cmp_model: FIXME
 *
 * Compares @mset0 and @mset1 and returns TRUE if both coitains the same models.
 * If @cmp_model is TRUE compare also if the models are the same.
 * 
 * Returns: FIXME
 */
gboolean 
ncm_mset_cmp (NcmMSet *mset0, NcmMSet *mset1, gboolean cmp_model)
{
  NcmModelID mid;
  for (mid = 0; mid < NCM_MODEL_MAX_ID; mid++)
  {
    NcmModel *model0 = ncm_mset_peek (mset0, mid);
    NcmModel *model1 = ncm_mset_peek (mset1, mid);
    if (model0 == NULL && model1 == NULL)
      continue;
    else if (model0 == NULL && model1 != NULL)
      return FALSE;
    else if (model1 == NULL && model0 != NULL)
      return FALSE;
    else if (cmp_model && !ncm_model_is_equal (model0, model1))
      return FALSE;
  }
  return TRUE;
}

/**
 * ncm_mset_param_set:
 * @mset: a #NcmMSet
 * @mid: FIXME
 * @pid: FIXME
 * @x: FIXME
 *
 * FIXME
 *
 */
void
ncm_mset_param_set (NcmMSet *mset, NcmModelID mid, guint pid, const gdouble x)
{
  ncm_model_param_set (ncm_mset_peek (mset, mid), pid, x);
}

/**
 * ncm_mset_param_get:
 * @mset: a #NcmMSet
 * @mid: FIXME
 * @pid: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
gdouble
ncm_mset_param_get (NcmMSet *mset, NcmModelID mid, guint pid)
{
  return ncm_model_param_get (ncm_mset_peek (mset, mid), pid);
}

/**
 * ncm_mset_orig_param_get:
 * @mset: a #NcmMSet
 * @mid: FIXME
 * @pid: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
gdouble
ncm_mset_orig_param_get (NcmMSet *mset, NcmModelID mid, guint pid)
{
  return ncm_model_orig_param_get (ncm_mset_peek (mset, mid), pid);
}

/**
 * ncm_mset_param_name:
 * @mset: a #NcmMSet
 * @mid: FIXME
 * @pid: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
const gchar *
ncm_mset_param_name (NcmMSet *mset, NcmModelID mid, guint pid)
{
  return ncm_model_param_name (ncm_mset_peek (mset, mid), pid);
}

/**
 * ncm_mset_param_get_scale:
 * @mset: a #NcmMSet
 * @mid: FIXME
 * @pid: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
gdouble
ncm_mset_param_get_scale (NcmMSet *mset, NcmModelID mid, guint pid)
{
  return ncm_model_param_get_scale (ncm_mset_peek (mset, mid), pid);
}

/**
 * ncm_mset_param_get_lower_bound:
 * @mset: a #NcmMSet
 * @mid: FIXME
 * @pid: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
gdouble
ncm_mset_param_get_lower_bound (NcmMSet *mset, NcmModelID mid, guint pid)
{
  return ncm_model_param_get_lower_bound (ncm_mset_peek (mset, mid), pid);
}

/**
 * ncm_mset_param_get_upper_bound:
 * @mset: a #NcmMSet
 * @mid: FIXME
 * @pid: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
gdouble
ncm_mset_param_get_upper_bound (NcmMSet *mset, NcmModelID mid, guint pid)
{
  return ncm_model_param_get_upper_bound (ncm_mset_peek (mset, mid), pid);
}

/**
 * ncm_mset_param_get_abstol:
 * @mset: a #NcmMSet
 * @mid: FIXME
 * @pid: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
gdouble
ncm_mset_param_get_abstol (NcmMSet *mset, NcmModelID mid, guint pid)
{
  return ncm_model_param_get_abstol (ncm_mset_peek (mset, mid), pid);
}

/**
 * ncm_mset_param_set_ftype:
 * @mset: a #NcmMSet
 * @mid: FIXME
 * @pid: FIXME
 * @ftype: FIXME
 *
 * FIXME
 *
 */
void
ncm_mset_param_set_ftype (NcmMSet *mset, NcmModelID mid, guint pid, NcmParamType ftype)
{
  NcmModel *model = ncm_mset_peek (mset, mid);
  g_assert (model != NULL);
  ncm_model_param_set_ftype (model, pid, ftype);
}

/**
 * ncm_mset_param_set_all_ftype:
 * @mset: a #NcmMSet
 * @ftype: FIXME
 *
 * FIXME
 *
 */
void
ncm_mset_param_set_all_ftype (NcmMSet *mset, NcmParamType ftype)
{
  NcmModelID mid;
  gint pid;

  for (mid = 0; mid < NCM_MODEL_MAX_ID; mid++)
  {
    NcmModel *model = ncm_mset_peek (mset, mid);
    if (model == NULL)
      continue;
    else
    {
      gint n_params = ncm_model_len (model);
      for (pid = 0; pid < n_params; pid++)
      {
        ncm_model_param_set_ftype (ncm_mset_peek (mset, mid), pid, ftype);
      }
    }
  }
}

/**
 * ncm_mset_param_set_vector:
 * @mset: a #NcmMSet
 * @params: FIXME
 *
 * Set the models parameters using values from @params.
 *
 */
void 
ncm_mset_param_set_vector (NcmMSet *mset, NcmVector *params)
{
  NcmModelID mid;
  guint j = 0;

  g_assert (ncm_vector_len (params) == ncm_mset_param_len (mset));

  for (mid = 0; mid < NCM_MODEL_MAX_ID; mid++)
  {
    NcmModel *model = ncm_mset_peek (mset, mid);
    if (model == NULL)
      continue;
    else
    {
      gint n_params = ncm_model_len (model);
      NcmVector *mvec = ncm_vector_get_subvector (params, j, n_params);
      ncm_model_params_set_vector (model, mvec);
      ncm_vector_free (mvec);
      j += n_params;
    }
  }
}

/**
 * ncm_mset_param_get_vector:
 * @mset: a #NcmMSet
 * @params: FIXME
 *
 * Set the compontents of @params using the models parameters.
 *
 */
void 
ncm_mset_param_get_vector (NcmMSet *mset, NcmVector *params)
{
  NcmModelID mid;
  gint pid;
  guint j = 0;

  g_assert (ncm_vector_len (params) == ncm_mset_param_len (mset));

  for (mid = 0; mid < NCM_MODEL_MAX_ID; mid++)
  {
    NcmModel *model = ncm_mset_peek (mset, mid);
    if (model == NULL)
      continue;
    else
    {
      gint n_params = ncm_model_len (model);
      for (pid = 0; pid < n_params; pid++)
      {
        ncm_vector_set (params, j++, ncm_model_param_get (model, pid));
      }
    }
  }
}

/**
 * ncm_mset_param_get_ftype:
 * @mset: a #NcmMSet.
 * @mid: a #NcmModelID.
 * @pid: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
NcmParamType
ncm_mset_param_get_ftype (NcmMSet *mset, NcmModelID mid, guint pid)
{
  NcmModel *model = ncm_mset_peek (mset, mid);
  if (model == NULL)
    g_error ("ncm_mset_param_get_ftype: cannot get ftype of mode %d, model not set.", mid);
  return ncm_model_param_get_ftype (model, pid);
}

/**
 * ncm_mset_param_set_pi:
 * @mset: a #NcmMSet
 * @pi: FIXME
 * @x: FIXME
 * @n: FIXME
 *
 * FIXME
 *
 */
void
ncm_mset_param_set_pi (NcmMSet *mset, NcmMSetPIndex *pi, const gdouble *x, guint n)
{
  guint fpi;
  for (fpi = 0; fpi < n; fpi++)
    ncm_mset_param_set (mset, pi[fpi].mid, pi[fpi].pid, x[fpi]);
}

/**
 * ncm_mset_param_get_pi:
 * @mset: a #NcmMSet
 * @pi: FIXME
 * @x: FIXME
 * @n: FIXME
 *
 * FIXME
 *
 */
void
ncm_mset_param_get_pi (NcmMSet *mset, NcmMSetPIndex *pi, gdouble *x, guint n)
{
  guint fpi;
  for (fpi = 0; fpi < n; fpi++)
    x[fpi] = ncm_mset_param_get (mset, pi[fpi].mid, pi[fpi].pid);
}

/**
 * ncm_mset_fparams_get_vector:
 * @mset: a #NcmMSet
 * @x: FIXME
 *
 * FIXME
 *
 */
void
ncm_mset_fparams_get_vector (NcmMSet *mset, NcmVector *x)
{
  guint fpi;

  for (fpi = 0; fpi < mset->fparam_len; fpi++)
  {
    const NcmMSetPIndex pi = g_array_index (mset->pi_array, NcmMSetPIndex, fpi);
    ncm_vector_set (x, fpi, ncm_mset_param_get (mset, pi.mid, pi.pid));
  }
}

/**
 * ncm_mset_fparams_set_vector:
 * @mset: a #NcmMSet
 * @x: FIXME
 *
 * FIXME
 *
 */
void
ncm_mset_fparams_set_vector (NcmMSet *mset, const NcmVector *x)
{
  guint fpi;

  for (fpi = 0; fpi < mset->fparam_len; fpi++)
  {
    const NcmMSetPIndex pi = g_array_index (mset->pi_array, NcmMSetPIndex, fpi);
    ncm_mset_param_set (mset, pi.mid, pi.pid, ncm_vector_get (x, fpi));
  }
}

/**
 * ncm_mset_fparams_set_array:
 * @mset: a #NcmMSet
 * @x: FIXME
 *
 * FIXME
 *
 */
void
ncm_mset_fparams_set_array (NcmMSet *mset, const gdouble *x)
{
  guint fpi;
  for (fpi = 0; fpi < mset->fparam_len; fpi++)
  {
    const NcmMSetPIndex pi = g_array_index (mset->pi_array, NcmMSetPIndex, fpi);
    ncm_mset_param_set (mset, pi.mid, pi.pid, x[fpi]);
  }
}

/**
 * ncm_mset_fparams_set_gsl_vector: (skip)
 * @mset: a #NcmMSet.
 * @x: FIXME
 *
 * FIXME
 *
 */
void
ncm_mset_fparams_set_gsl_vector (NcmMSet *mset, const gsl_vector *x)
{
  guint fpi;
  for (fpi = 0; fpi < mset->fparam_len; fpi++)
  {
    const NcmMSetPIndex pi = g_array_index (mset->pi_array, NcmMSetPIndex, fpi);
    ncm_mset_param_set (mset, pi.mid, pi.pid, gsl_vector_get (x, fpi));
  }
}

/**
 * ncm_mset_fparams_len:
 * @mset: a #NcmMSet
 *
 * FIXME
 *
 * Returns: FIXME
 */
guint
ncm_mset_fparams_len (NcmMSet *mset)
{
  g_assert (mset->valid_map);
  return mset->fparam_len;
}

/**
 * ncm_mset_fparam_name:
 * @mset: a #NcmMSet
 * @n: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
const gchar *
ncm_mset_fparam_name (NcmMSet *mset, guint n)
{
  g_assert (mset->valid_map && n < mset->fparam_len);
  {
    const NcmMSetPIndex pi = g_array_index (mset->pi_array, NcmMSetPIndex, n);
    return ncm_mset_param_name (mset, pi.mid, pi.pid);
  }
}

/**
 * ncm_mset_fparam_get_scale:
 * @mset: a #NcmMSet
 * @n: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
gdouble
ncm_mset_fparam_get_scale (NcmMSet *mset, guint n)
{
  g_assert (mset->valid_map && n < mset->fparam_len);
  {
    const NcmMSetPIndex pi = g_array_index (mset->pi_array, NcmMSetPIndex, n);
    return ncm_mset_param_get_scale (mset, pi.mid, pi.pid);
  }
}

/**
 * ncm_mset_fparam_get_lower_bound:
 * @mset: a #NcmMSet
 * @n: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
gdouble
ncm_mset_fparam_get_lower_bound (NcmMSet *mset, guint n)
{
  g_assert (mset->valid_map && n < mset->fparam_len);
  {
    const NcmMSetPIndex pi = g_array_index (mset->pi_array, NcmMSetPIndex, n);
    return ncm_mset_param_get_lower_bound (mset, pi.mid, pi.pid);
  }
}

/**
 * ncm_mset_fparam_get_upper_bound:
 * @mset: a #NcmMSet
 * @n: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
gdouble
ncm_mset_fparam_get_upper_bound (NcmMSet *mset, guint n)
{
  g_assert (mset->valid_map && n < mset->fparam_len);
  {
    const NcmMSetPIndex pi = g_array_index (mset->pi_array, NcmMSetPIndex, n);
    return ncm_mset_param_get_upper_bound (mset, pi.mid, pi.pid);
  }
}

/**
 * ncm_mset_fparam_get_abstol:
 * @mset: a #NcmMSet
 * @n: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
gdouble
ncm_mset_fparam_get_abstol (NcmMSet *mset, guint n)
{
  g_assert (mset->valid_map && n < mset->fparam_len);
  {
    const NcmMSetPIndex pi = g_array_index (mset->pi_array, NcmMSetPIndex, n);
    return ncm_mset_param_get_abstol (mset, pi.mid, pi.pid);
  }
}

/**
 * ncm_mset_fparam_get:
 * @mset: a #NcmMSet
 * @n: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
gdouble
ncm_mset_fparam_get (NcmMSet *mset, guint n)
{
  g_assert (mset->valid_map && n < mset->fparam_len);
  {
    const NcmMSetPIndex pi = g_array_index (mset->pi_array, NcmMSetPIndex, n);
    return ncm_mset_param_get (mset, pi.mid, pi.pid);
  }
}

/**
 * ncm_mset_fparam_set:
 * @mset: a #NcmMSet
 * @n: FIXME
 * @x: FIXME
 *
 * FIXME
 *
 */
void
ncm_mset_fparam_set (NcmMSet *mset, guint n, const gdouble x)
{
  g_assert (mset->valid_map && n < mset->fparam_len);
  {
    const NcmMSetPIndex pi = g_array_index (mset->pi_array, NcmMSetPIndex, n);
    return ncm_mset_param_set (mset, pi.mid, pi.pid, x);
  }
}

/**
 * ncm_mset_fparam_get_pi:
 * @mset: a #NcmMSet
 * @n: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
NcmMSetPIndex *
ncm_mset_fparam_get_pi (NcmMSet *mset, guint n)
{
  g_assert (mset->valid_map && n < mset->fparam_len);
  return &(g_array_index (mset->pi_array, NcmMSetPIndex, n));
}

/**
 * ncm_mset_fparam_get_fpi:
 * @mset: a #NcmMSet.
 * @mid: a #NcmModelID.
 * @pid: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
gint
ncm_mset_fparam_get_fpi (NcmMSet *mset, NcmModelID mid, guint pid)
{
  g_assert (mset->valid_map);
  return g_array_index (mset->fpi_array[mid], gint, pid);
}

static void
ncm_mset_init (NcmMSet *mset)
{
  guint i;
  memset (mset->model, 0, sizeof (NcmModel *) * NCM_MODEL_MAX_ID);
  mset->pi_array = g_array_sized_new (FALSE, TRUE, sizeof (NcmMSetPIndex), 20);
  for (i = 0; i < NCM_MODEL_MAX_ID; i++)
    mset->fpi_array[i] = g_array_sized_new (FALSE, TRUE, sizeof (gint), 10);
  mset->valid_map = FALSE;
  mset->total_len = 0;
}

static void
ncm_mset_dispose (GObject *object)
{
  NcmMSet *mset = NCM_MSET (object);
  NcmModelID mid;

  for (mid = 0; mid < NCM_MODEL_MAX_ID; mid++)
  {
    if (mset->model[mid] != NULL)
      ncm_model_clear (&mset->model[mid]);
  }

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_mset_parent_class)->dispose (object);
}

static void
ncm_mset_finalize (GObject *object)
{
  NcmMSet *mset = NCM_MSET (object);
  gint i;

  for (i = 0; i < NCM_MODEL_MAX_ID; i++)
    g_array_unref (mset->fpi_array[i]);

  g_array_unref (mset->pi_array);

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_mset_parent_class)->finalize (object);
}

static void
ncm_mset_class_init (NcmMSetClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);

  object_class->dispose = ncm_mset_dispose;
  object_class->finalize = ncm_mset_finalize;
}

