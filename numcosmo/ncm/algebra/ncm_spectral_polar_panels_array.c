/***************************************************************************
 *            ncm_spectral_polar_panels_array.c
 *
 *  Wed Nov 11 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 *  Copyright  2026  Caio Lima de Oliveira
 *  <caiolimadeoliveira@pm.me>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2026 <vitenti@uel.br>
 * Copyright (C) 2026 Caio Lima de Oliveira <caiolimadeoliveira@pm.me>
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
 * NcmSpectralPolarPanelsArray:
 *
 * An ordered set of #NcmSpectralPolarPanels kept and serialized as one object, such as the
 * terms of an expansion of a function over parameters, one #NcmSpectralPolarPanels per
 * term (see ncm_spectral_polar_panels_pack_eval_on_array_combine()).
 *
 * The sets are added with ncm_spectral_polar_panels_array_add(), and
 * ncm_spectral_polar_panels_array_pack() then keeps them in one allocation only, releasing
 * them: a packed array no longer grows, and its sets are read through its pack
 * (ncm_spectral_polar_panels_array_peek_pack()). The array is serialized packed
 * (#NcmSpectralPolarPanelsArray:packed), so that reading a serialization makes no
 * #NcmSpectralPolarPanels.
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "ncm/algebra/ncm_spectral_polar_panels_array.h"
#include "ncm/algebra/ncm_spectral_polar_panels_private.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <math.h>
#include <string.h>
#endif /* NUMCOSMO_GIR_SCAN */

/*
 * The panels combined so far for one set of weights: index[i] is panel i with the offset of
 * its combined coefficients in coeffs, or NCM_SPECTRAL_POLAR_PANELS_COMB_NONE.
 */
struct _NcmSpectralPolarPanelsComb
{
  GArray *index;  /* NcmSpectralPolarPanelsPanel */
  GArray *coeffs; /* gdouble */
};

#define NCM_SPECTRAL_POLAR_PANELS_COMB_NONE G_MAXSIZE

/**
 * ncm_spectral_polar_panels_comb_new: (skip)
 *
 * Creates an empty cache of combined panels for
 * ncm_spectral_polar_panels_pack_eval_on_array_combine().
 *
 * Returns: (transfer full): a new #NcmSpectralPolarPanelsComb.
 */
NcmSpectralPolarPanelsComb *
ncm_spectral_polar_panels_comb_new (void)
{
  NcmSpectralPolarPanelsComb *comb = g_new (NcmSpectralPolarPanelsComb, 1);

  comb->index  = g_array_new (FALSE, FALSE, sizeof (NcmSpectralPolarPanelsPanel));
  comb->coeffs = g_array_new (FALSE, FALSE, sizeof (gdouble));

  return comb;
}

/**
 * ncm_spectral_polar_panels_comb_free: (skip)
 * @comb: a #NcmSpectralPolarPanelsComb
 *
 * Frees @comb.
 */
void
ncm_spectral_polar_panels_comb_free (NcmSpectralPolarPanelsComb *comb)
{
  g_array_unref (comb->index);
  g_array_unref (comb->coeffs);
  g_free (comb);
}

/**
 * ncm_spectral_polar_panels_comb_reset: (skip)
 * @comb: a #NcmSpectralPolarPanelsComb
 *
 * Forgets every panel @comb holds. It must be called whenever the pack or the weights
 * passed to ncm_spectral_polar_panels_pack_eval_on_array_combine() change.
 */
void
ncm_spectral_polar_panels_comb_reset (NcmSpectralPolarPanelsComb *comb)
{
  g_array_set_size (comb->index, 0);
  g_array_set_size (comb->coeffs, 0);
}

/*
 * The sets of an NcmSpectralPolarPanelsArray in one allocation, for
 * ncm_spectral_polar_panels_pack_eval_on_array_combine(): the edges, the initial panels and
 * the panels of the first set, each panel with the largest numbers of coefficients of the
 * sets there and the offset in data of its terms. These are a count, n_kept, in the first
 * slot, followed by the terms of the sets whose block there is not a single zero, in the
 * order of the sets: a slot with the index of the set, a slot with its numbers of
 * coefficients, and the coefficients. Every array lies in the allocation of the pack, after
 * the structure, so that one span holds everything an evaluation reads.
 */
struct _NcmSpectralPolarPanelsPack
{
  gsize size; /* bytes of the allocation */
  guint n_sets;
  guint n_panels;
  guint n_r_edges;
  guint n_theta_edges;
  const NcmSpectralPolarPanelsPanel *panels;
  const gdouble *r_edges;
  const gdouble *theta_edges;
  const guint *first;
  const gdouble *data;
};

/* Bytes rounded up to whole doubles */
#define NCM_SPECTRAL_POLAR_PANELS_PACK_ROUND(b) ((((b) + sizeof (gdouble) - 1) / sizeof (gdouble)) * sizeof (gdouble))

/* A slot of the terms of a pack: two guint in the place of a double */
typedef union _NcmSpectralPolarPanelsPackSlot
{
  gdouble d;
  guint u[2];
} NcmSpectralPolarPanelsPackSlot;

G_STATIC_ASSERT (sizeof (NcmSpectralPolarPanelsPackSlot) == sizeof (gdouble));

/*
 * Whether the block of a term is a single +0.0, the block of a term trim_modes() drops. A -0.0
 * is kept as any other coefficient.
 */
static gboolean
_ncm_spectral_polar_panels_pack_dropped (const NcmSpectralPolarPanelsPanel *pq, const gdouble *c)
{
  return (pq->n_r == 1) && (pq->n_theta == 1) && (c[0] == 0.0) && !signbit (c[0]);
}

/* Aborts unless the n sets spp have the same edges and the same panels, as a pack keeps them once */
static void
_ncm_spectral_polar_panels_pack_check_layout (guint n, NcmSpectralPolarPanels **spp)
{
  NcmSpectralPolarPanels *s0 = spp[0];
  guint q, i;

  for (q = 1; q < n; q++)
  {
    NcmSpectralPolarPanels *sq = spp[q];
    gboolean same              = (sq->panels->len == s0->panels->len) && (sq->r_edges->len == s0->r_edges->len) &&
                                 (sq->theta_edges->len == s0->theta_edges->len) &&
                                 (memcmp (sq->r_edges->data, s0->r_edges->data, sizeof (gdouble) * s0->r_edges->len) == 0) &&
                                 (memcmp (sq->theta_edges->data, s0->theta_edges->data, sizeof (gdouble) * s0->theta_edges->len) == 0);

    for (i = 0; same && (i < s0->panels->len); i++)
    {
      const NcmSpectralPolarPanelsPanel *p0 = &g_array_index (s0->panels, NcmSpectralPolarPanelsPanel, i);
      const NcmSpectralPolarPanelsPanel *pq = &g_array_index (sq->panels, NcmSpectralPolarPanelsPanel, i);

      same = (p0->ra == pq->ra) && (p0->rb == pq->rb) && (p0->ta == pq->ta) && (p0->tb == pq->tb);
    }

    if (!same)
      g_error ("ncm_spectral_polar_panels_array: set %u does not have the edges and panels of set 0.", q);
  }
}

static NcmSpectralPolarPanelsPack *
_ncm_spectral_polar_panels_pack_new (guint n, NcmSpectralPolarPanels **spp)
{
  NcmSpectralPolarPanels *s0 = spp[0];
  const guint n_panels       = s0->panels->len;
  gsize n_data               = 0;
  gsize b_panels, b_r, b_theta, b_first, size;
  NcmSpectralPolarPanelsPack *pack;
  NcmSpectralPolarPanelsPanel *panels;
  gdouble *data;
  gchar *base;
  guint i, q;

  _ncm_spectral_polar_panels_pack_check_layout (n, spp);

  for (i = 0; i < n_panels; i++)
  {
    n_data++;

    for (q = 0; q < n; q++)
    {
      const NcmSpectralPolarPanelsPanel *pq = &g_array_index (spp[q]->panels, NcmSpectralPolarPanelsPanel, i);

      if (!_ncm_spectral_polar_panels_pack_dropped (pq, &g_array_index (spp[q]->coeffs, gdouble, pq->offset)))
        n_data += 2 + (gsize) pq->n_r * pq->n_theta;
    }
  }

  b_panels = NCM_SPECTRAL_POLAR_PANELS_PACK_ROUND (sizeof (NcmSpectralPolarPanelsPanel) * n_panels);
  b_r      = sizeof (gdouble) * s0->r_edges->len;
  b_theta  = sizeof (gdouble) * s0->theta_edges->len;
  b_first  = NCM_SPECTRAL_POLAR_PANELS_PACK_ROUND (sizeof (guint) * s0->first->len);
  size     = NCM_SPECTRAL_POLAR_PANELS_PACK_ROUND (sizeof (NcmSpectralPolarPanelsPack)) + b_panels + b_r + b_theta + b_first +
             sizeof (gdouble) * n_data;

  base   = g_malloc (size);
  pack   = (NcmSpectralPolarPanelsPack *) base;
  base  += NCM_SPECTRAL_POLAR_PANELS_PACK_ROUND (sizeof (NcmSpectralPolarPanelsPack));
  panels = (NcmSpectralPolarPanelsPanel *) base;
  base  += b_panels;

  pack->size          = size;
  pack->n_sets        = n;
  pack->n_panels      = n_panels;
  pack->n_r_edges     = s0->r_edges->len;
  pack->n_theta_edges = s0->theta_edges->len;
  pack->panels        = panels;

  memcpy (base, s0->r_edges->data, b_r);
  pack->r_edges = (const gdouble *) base;
  base         += b_r;

  memcpy (base, s0->theta_edges->data, b_theta);
  pack->theta_edges = (const gdouble *) base;
  base             += b_theta;

  memcpy (base, s0->first->data, sizeof (guint) * s0->first->len);
  pack->first = (const guint *) base;
  base       += b_first;

  data       = (gdouble *) base;
  pack->data = data;
  n_data     = 0;

  for (i = 0; i < n_panels; i++)
  {
    NcmSpectralPolarPanelsPanel *panel = &panels[i];
    NcmSpectralPolarPanelsPackSlot count;

    *panel         = g_array_index (s0->panels, NcmSpectralPolarPanelsPanel, i);
    panel->n_r     = 0;
    panel->n_theta = 0;
    panel->offset  = n_data;
    count.u[0]     = 0;
    count.u[1]     = 0;
    n_data++;

    for (q = 0; q < n; q++)
    {
      const NcmSpectralPolarPanelsPanel *pq = &g_array_index (spp[q]->panels, NcmSpectralPolarPanelsPanel, i);
      const gdouble *cq                     = &g_array_index (spp[q]->coeffs, gdouble, pq->offset);
      const gsize len                       = (gsize) pq->n_r * pq->n_theta;
      NcmSpectralPolarPanelsPackSlot head;

      /* The largest block of every set, as ncm_spectral_polar_panels_eval_on_array_combine() takes it */
      panel->n_r     = MAX (panel->n_r, pq->n_r);
      panel->n_theta = MAX (panel->n_theta, pq->n_theta);

      if (_ncm_spectral_polar_panels_pack_dropped (pq, cq))
        continue;

      head.u[0] = q;
      head.u[1] = 0;
      memcpy (&data[n_data], &head, sizeof (gdouble));
      head.u[0] = pq->n_r;
      head.u[1] = pq->n_theta;
      memcpy (&data[n_data + 1], &head, sizeof (gdouble));
      memcpy (&data[n_data + 2], cq, sizeof (gdouble) * len);
      n_data += 2 + len;
      count.u[0]++;
    }

    memcpy (&data[panel->offset], &count, sizeof (gdouble));
  }

  return pack;
}

/**
 * ncm_spectral_polar_panels_pack_get_size: (skip)
 * @pack: a #NcmSpectralPolarPanelsPack
 *
 * Returns: the number of bytes of @pack, which start at @pack.
 */
gsize
ncm_spectral_polar_panels_pack_get_size (const NcmSpectralPolarPanelsPack *pack)
{
  return pack->size;
}

/**
 * ncm_spectral_polar_panels_pack_find_panel_near: (skip)
 * @pack: a #NcmSpectralPolarPanelsPack
 * @hint: a panel index, usually the result of the previous search
 * @r: a radius
 * @theta: an angle
 *
 * ncm_spectral_polar_panels_find_panel_near() on the first set of @pack.
 *
 * Returns: the index of the panel.
 */
guint
ncm_spectral_polar_panels_pack_find_panel_near (const NcmSpectralPolarPanelsPack *pack, guint hint, gdouble r, gdouble theta)
{
  if (hint < pack->n_panels)
  {
    const NcmSpectralPolarPanelsPanel *panel = &pack->panels[hint];

    if ((r >= panel->ra) && (r < panel->rb) && (theta >= panel->ta) && (theta < panel->tb))
      return hint;
  }

  return _ncm_spectral_polar_panels_find_panel_in (pack->r_edges, pack->n_r_edges, pack->theta_edges, pack->n_theta_edges,
                                                   pack->first, pack->panels, r, theta);
}

/**
 * ncm_spectral_polar_panels_pack_eval_on_array_combine: (skip)
 * @pack: a #NcmSpectralPolarPanelsPack
 * @w: (array): the weights, one per set of @pack
 * @comb: the panels combined so far for this @pack and @w
 * @i: panel index
 * @n: number of points
 * @r: (array length=n): the radii
 * @theta: (array length=n): the angles
 * @out: (array length=n): the values, set
 *
 * Evaluates at the @n points of panel @i the series of $\sum_q w_q F_q$, set $q$ of @pack
 * the expansion of $F_q$, as ncm_spectral_polar_panels_eval_on_array() would on the set
 * ncm_spectral_polar_panels_combine() makes, combining only the coefficients of panel @i
 * instead of those of every panel: the coefficients are added in the order of the sets, and
 * the sets whose block on panel @i is a single zero, which add nothing, are skipped. The
 * combined panel is kept in @comb, so that later calls on panel @i do not combine it again:
 * @comb must be reset (ncm_spectral_polar_panels_comb_reset()) when @pack or @w change. A
 * pack of a single set with weight one is evaluated on its coefficients as they are, the
 * values of ncm_spectral_polar_panels_eval_on_array() on the set, and @comb is not used.
 */
void
ncm_spectral_polar_panels_pack_eval_on_array_combine (const NcmSpectralPolarPanelsPack *pack, const gdouble *w, NcmSpectralPolarPanelsComb *comb,
                                                      guint i, guint n, const gdouble *r, const gdouble *theta, gdouble *out)
{
  NcmSpectralPolarPanelsPanel panel;
  const gdouble *c;

  /* A single set with weight one is evaluated where its terms lie, without @comb */
  if ((pack->n_sets == 1) && (w[0] == 1.0))
  {
    NcmSpectralPolarPanelsPackSlot slot;
    const gdouble *t;

    g_assert_cmpuint (i, <, pack->n_panels);

    panel = pack->panels[i];
    t     = pack->data + panel.offset;
    memcpy (&slot, t, sizeof (gdouble));

    if (slot.u[0] == 0)
    {
      /* The block of the set is a single zero */
      guint j;

      for (j = 0; j < n; j++)
        out[j] = 0.0;

      return;
    }

    memcpy (&slot, t + 2, sizeof (gdouble));
    panel.n_r     = slot.u[0];
    panel.n_theta = slot.u[1];
    panel.offset  = 0;

    _ncm_spectral_polar_panels_eval_points (&panel, t + 3, n, r, theta, out);

    return;
  }

  if (comb->index->len <= i)
  {
    const NcmSpectralPolarPanelsPanel none = {0.0, 0.0, 0.0, 0.0, 0, 0, NCM_SPECTRAL_POLAR_PANELS_COMB_NONE};

    g_assert_cmpuint (i, <, pack->n_panels);

    while (comb->index->len < pack->n_panels)
      g_array_append_val (comb->index, none);
  }

  panel = g_array_index (comb->index, NcmSpectralPolarPanelsPanel, i);

  if (panel.offset == NCM_SPECTRAL_POLAR_PANELS_COMB_NONE)
  {
    const gdouble *t;
    NcmSpectralPolarPanelsPackSlot slot;
    gdouble *cc;
    guint m, n_kept;

    panel = pack->panels[i];
    t     = pack->data + panel.offset;
    memcpy (&slot, t, sizeof (gdouble));
    n_kept = slot.u[0];
    t++;

    panel.offset = comb->coeffs->len;
    g_array_set_size (comb->coeffs, panel.offset + (gsize) panel.n_r * panel.n_theta);
    cc = (gdouble *) comb->coeffs->data + panel.offset;
    memset (cc, 0, sizeof (gdouble) * panel.n_r * panel.n_theta);

    for (m = 0; m < n_kept; m++)
    {
      guint q, n_r, n_theta, l, k;
      gdouble wq;

      memcpy (&slot, t, sizeof (gdouble));
      q = slot.u[0];
      memcpy (&slot, t + 1, sizeof (gdouble));
      n_r     = slot.u[0];
      n_theta = slot.u[1];
      wq      = w[q];

      if (wq != 0.0)
      {
        const gdouble *cq = t + 2;

        for (l = 0; l < n_theta; l++)
        {
          gdouble *row_o       = cc + (gsize) l * panel.n_r;
          const gdouble *row_q = cq + (gsize) l * n_r;

          #pragma omp simd

          for (k = 0; k < n_r; k++)
            row_o[k] += wq * row_q[k];
        }
      }

      t += 2 + (gsize) n_r * n_theta;
    }

    g_array_index (comb->index, NcmSpectralPolarPanelsPanel, i) = panel;
  }

  c            = (const gdouble *) comb->coeffs->data + panel.offset;
  panel.offset = 0;

  _ncm_spectral_polar_panels_eval_points (&panel, c, n, r, theta, out);
}

/**
 * ncm_spectral_polar_panels_pack_get_range: (skip)
 * @pack: a #NcmSpectralPolarPanelsPack
 * @ra: (out): first radial edge
 * @rb: (out): last radial edge
 * @ta: (out): first angular edge
 * @tb: (out): last angular edge
 *
 * ncm_spectral_polar_panels_get_range() of the sets of @pack.
 */
void
ncm_spectral_polar_panels_pack_get_range (const NcmSpectralPolarPanelsPack *pack, gdouble *ra, gdouble *rb, gdouble *ta, gdouble *tb)
{
  *ra = pack->r_edges[0];
  *rb = pack->r_edges[pack->n_r_edges - 1];
  *ta = pack->theta_edges[0];
  *tb = pack->theta_edges[pack->n_theta_edges - 1];
}

/* The pack as NcmSpectralPolarPanelsArray:packed holds it */
static GVariant *
_ncm_spectral_polar_panels_pack_to_variant (const NcmSpectralPolarPanelsPack *pack)
{
  GVariantBuilder panels;
  guint i, m;

  g_variant_builder_init (&panels, G_VARIANT_TYPE ("a(dddduua(uuuad))"));

  for (i = 0; i < pack->n_panels; i++)
  {
    const NcmSpectralPolarPanelsPanel *panel = &pack->panels[i];
    const gdouble *t                         = pack->data + panel->offset;
    NcmSpectralPolarPanelsPackSlot slot;
    GVariantBuilder terms;
    guint n_kept;

    memcpy (&slot, t, sizeof (gdouble));
    n_kept = slot.u[0];
    t++;

    g_variant_builder_init (&terms, G_VARIANT_TYPE ("a(uuuad)"));

    for (m = 0; m < n_kept; m++)
    {
      guint q, n_r, n_theta;

      memcpy (&slot, t, sizeof (gdouble));
      q = slot.u[0];
      memcpy (&slot, t + 1, sizeof (gdouble));
      n_r     = slot.u[0];
      n_theta = slot.u[1];

      g_variant_builder_add (&terms, "(uuu@ad)", q, n_r, n_theta,
                             g_variant_new_fixed_array (G_VARIANT_TYPE ("d"), t + 2, (gsize) n_r * n_theta, sizeof (gdouble)));
      t += 2 + (gsize) n_r * n_theta;
    }

    g_variant_builder_add (&panels, "(dddduu@a(uuuad))", panel->ra, panel->rb, panel->ta, panel->tb, panel->n_r, panel->n_theta,
                           g_variant_builder_end (&terms));
  }

  return g_variant_ref_sink (g_variant_new ("(u@ad@ad@au@a(dddduua(uuuad)))", pack->n_sets,
                                            g_variant_new_fixed_array (G_VARIANT_TYPE ("d"), pack->r_edges, pack->n_r_edges, sizeof (gdouble)),
                                            g_variant_new_fixed_array (G_VARIANT_TYPE ("d"), pack->theta_edges, pack->n_theta_edges, sizeof (gdouble)),
                                            g_variant_new_fixed_array (G_VARIANT_TYPE ("u"), pack->first,
                                                                       (gsize) (pack->n_r_edges - 1) * (pack->n_theta_edges - 1) + 1, sizeof (guint)),
                                            g_variant_builder_end (&panels)));
}

/*
 * The pack of NcmSpectralPolarPanelsArray:packed, checked as ncm_spectral_polar_panels_append()
 * checks a panel: edges ascending, initial panels holding their panels in order, the terms of
 * a panel in ascending sets, inside its largest block.
 */
static NcmSpectralPolarPanelsPack *
_ncm_spectral_polar_panels_pack_from_variant (GVariant *var)
{
  GVariant *r_var      = g_variant_get_child_value (var, 1);
  GVariant *theta_var  = g_variant_get_child_value (var, 2);
  GVariant *first_var  = g_variant_get_child_value (var, 3);
  GVariant *panels_var = g_variant_get_child_value (var, 4);
  const guint n_panels = g_variant_n_children (panels_var);
  gsize n_r_edges = 0, n_theta_edges = 0, n_first = 0, n_data = 0;
  const gdouble *r_edges     = g_variant_get_fixed_array (r_var, &n_r_edges, sizeof (gdouble));
  const gdouble *theta_edges = g_variant_get_fixed_array (theta_var, &n_theta_edges, sizeof (gdouble));
  const guint *first         = g_variant_get_fixed_array (first_var, &n_first, sizeof (guint));
  gsize b_panels, b_r, b_theta, b_first, size;
  NcmSpectralPolarPanelsPack *pack;
  NcmSpectralPolarPanelsPanel *panels;
  gdouble *data;
  gchar *base;
  guint n_sets, i, c;

  g_variant_get_child (var, 0, "u", &n_sets);

  if ((n_sets == 0) || (n_r_edges < 2) || (n_theta_edges < 2) || (n_panels == 0) ||
      (n_first != (n_r_edges - 1) * (n_theta_edges - 1) + 1) || (first[0] != 0) || (first[n_first - 1] != n_panels))
    g_error ("ncm_spectral_polar_panels_array: the packed sets have %u sets, %zu and %zu edges, %u panels and %zu "
             "initial panel offsets, which do not agree.", n_sets, n_r_edges, n_theta_edges, n_panels, n_first);

  for (i = 0; i + 1 < n_r_edges; i++)
    if (!(r_edges[i] < r_edges[i + 1]))
      g_error ("ncm_spectral_polar_panels_array: the packed radial edges are not ascending at index %u.", i);

  for (i = 0; i + 1 < n_theta_edges; i++)
    if (!(theta_edges[i] < theta_edges[i + 1]))
      g_error ("ncm_spectral_polar_panels_array: the packed angular edges are not ascending at index %u.", i);

  for (c = 0; c + 1 < n_first; c++)
    if (first[c] > first[c + 1])
      g_error ("ncm_spectral_polar_panels_array: the packed initial panel offsets decrease at index %u.", c);

  for (i = 0; i < n_panels; i++)
  {
    GVariant *panel_var = g_variant_get_child_value (panels_var, i);
    GVariant *terms_var = g_variant_get_child_value (panel_var, 6);
    const guint n_terms = g_variant_n_children (terms_var);
    guint n_r_max, n_theta_max, m, last_q = 0;

    g_variant_get_child (panel_var, 4, "u", &n_r_max);
    g_variant_get_child (panel_var, 5, "u", &n_theta_max);

    if ((n_r_max == 0) || (n_theta_max == 0))
      g_error ("ncm_spectral_polar_panels_array: packed panel %u has an empty largest block.", i);

    n_data++;

    for (m = 0; m < n_terms; m++)
    {
      GVariant *coeffs_var;
      guint q, n_r, n_theta;
      gsize len = 0;

      g_variant_get_child (terms_var, m, "(uuu@ad)", &q, &n_r, &n_theta, &coeffs_var);
      g_variant_get_fixed_array (coeffs_var, &len, sizeof (gdouble));

      if ((q >= n_sets) || ((m > 0) && (q <= last_q)) || (n_r == 0) || (n_theta == 0) || (n_r > n_r_max) ||
          (n_theta > n_theta_max) || (len != (gsize) n_r * n_theta))
        g_error ("ncm_spectral_polar_panels_array: term %u of packed panel %u (set %u, %u by %u, %zu coefficients) "
                 "is not a block of a set inside the panel's %u by %u.", m, i, q, n_theta, n_r, len, n_theta_max, n_r_max);

      last_q  = q;
      n_data += 2 + len;
      g_variant_unref (coeffs_var);
    }

    g_variant_unref (terms_var);
    g_variant_unref (panel_var);
  }

  b_panels = NCM_SPECTRAL_POLAR_PANELS_PACK_ROUND (sizeof (NcmSpectralPolarPanelsPanel) * n_panels);
  b_r      = sizeof (gdouble) * n_r_edges;
  b_theta  = sizeof (gdouble) * n_theta_edges;
  b_first  = NCM_SPECTRAL_POLAR_PANELS_PACK_ROUND (sizeof (guint) * n_first);
  size     = NCM_SPECTRAL_POLAR_PANELS_PACK_ROUND (sizeof (NcmSpectralPolarPanelsPack)) + b_panels + b_r + b_theta + b_first +
             sizeof (gdouble) * n_data;

  base   = g_malloc (size);
  pack   = (NcmSpectralPolarPanelsPack *) base;
  base  += NCM_SPECTRAL_POLAR_PANELS_PACK_ROUND (sizeof (NcmSpectralPolarPanelsPack));
  panels = (NcmSpectralPolarPanelsPanel *) base;
  base  += b_panels;

  pack->size          = size;
  pack->n_sets        = n_sets;
  pack->n_panels      = n_panels;
  pack->n_r_edges     = n_r_edges;
  pack->n_theta_edges = n_theta_edges;
  pack->panels        = panels;

  memcpy (base, r_edges, b_r);
  pack->r_edges = (const gdouble *) base;
  base         += b_r;

  memcpy (base, theta_edges, b_theta);
  pack->theta_edges = (const gdouble *) base;
  base             += b_theta;

  memcpy (base, first, sizeof (guint) * n_first);
  pack->first = (const guint *) base;
  base       += b_first;

  data       = (gdouble *) base;
  pack->data = data;
  n_data     = 0;

  for (c = 0; c + 1 < n_first; c++)
  {
    const guint ci        = c / (n_theta_edges - 1);
    const guint cj        = c % (n_theta_edges - 1);
    const gdouble cell_ra = r_edges[ci];
    const gdouble cell_rb = r_edges[ci + 1];
    const gdouble cell_ta = theta_edges[cj];
    const gdouble cell_tb = theta_edges[cj + 1];

    for (i = first[c]; i < first[c + 1]; i++)
    {
      NcmSpectralPolarPanelsPanel *panel = &panels[i];
      GVariant *panel_var                = g_variant_get_child_value (panels_var, i);
      GVariant *terms_var                = g_variant_get_child_value (panel_var, 6);
      const guint n_terms                = g_variant_n_children (terms_var);
      NcmSpectralPolarPanelsPackSlot slot;
      guint m;

      g_variant_get_child (panel_var, 0, "d", &panel->ra);
      g_variant_get_child (panel_var, 1, "d", &panel->rb);
      g_variant_get_child (panel_var, 2, "d", &panel->ta);
      g_variant_get_child (panel_var, 3, "d", &panel->tb);
      g_variant_get_child (panel_var, 4, "u", &panel->n_r);
      g_variant_get_child (panel_var, 5, "u", &panel->n_theta);

      if (!(panel->ra < panel->rb) || !(panel->ta < panel->tb) || (panel->ra < cell_ra) || (panel->rb > cell_rb) ||
          (panel->ta < cell_ta) || (panel->tb > cell_tb))
        g_error ("ncm_spectral_polar_panels_array: packed panel %u [%.17g, %.17g] x [%.17g, %.17g] is empty or not "
                 "inside its initial panel %u.", i, panel->ra, panel->rb, panel->ta, panel->tb, c);

      panel->offset = n_data;
      slot.u[0]     = n_terms;
      slot.u[1]     = 0;
      memcpy (&data[n_data], &slot, sizeof (gdouble));
      n_data++;

      for (m = 0; m < n_terms; m++)
      {
        GVariant *coeffs_var;
        const gdouble *cq;
        guint q, n_r, n_theta;
        gsize len = 0;

        g_variant_get_child (terms_var, m, "(uuu@ad)", &q, &n_r, &n_theta, &coeffs_var);
        cq = g_variant_get_fixed_array (coeffs_var, &len, sizeof (gdouble));

        slot.u[0] = q;
        slot.u[1] = 0;
        memcpy (&data[n_data], &slot, sizeof (gdouble));
        slot.u[0] = n_r;
        slot.u[1] = n_theta;
        memcpy (&data[n_data + 1], &slot, sizeof (gdouble));
        memcpy (&data[n_data + 2], cq, sizeof (gdouble) * len);
        n_data += 2 + len;

        g_variant_unref (coeffs_var);
      }

      g_variant_unref (terms_var);
      g_variant_unref (panel_var);
    }
  }

  g_variant_unref (r_var);
  g_variant_unref (theta_var);
  g_variant_unref (first_var);
  g_variant_unref (panels_var);

  return pack;
}

enum
{
  PROP_ARRAY_0,
  PROP_ARRAY_PACKED,
};

/* The type of NcmSpectralPolarPanelsArray:packed */
#define NCM_SPECTRAL_POLAR_PANELS_ARRAY_PACKED_TYPE "(uadadaua(dddduua(uuuad)))"

struct _NcmSpectralPolarPanelsArray
{
  /*< private >*/
  GObject parent_instance;
  NcmObjArray *panels;              /* the sets as added, NULL once packed */
  guint n_sets;                     /* the number of sets */
  NcmSpectralPolarPanelsPack *pack; /* made by ncm_spectral_polar_panels_array_pack() */
};

G_DEFINE_TYPE (NcmSpectralPolarPanelsArray, ncm_spectral_polar_panels_array, G_TYPE_OBJECT)

static void
ncm_spectral_polar_panels_array_init (NcmSpectralPolarPanelsArray *spa)
{
  spa->panels = ncm_obj_array_new ();
  spa->n_sets = 0;
  spa->pack   = NULL;
}

static void
_ncm_spectral_polar_panels_array_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmSpectralPolarPanelsArray *spa = NCM_SPECTRAL_POLAR_PANELS_ARRAY (object);

  g_return_if_fail (NCM_IS_SPECTRAL_POLAR_PANELS_ARRAY (object));

  switch (prop_id)
  {
    case PROP_ARRAY_PACKED:
    {
      GVariant *var = g_value_get_variant (value);

      ncm_obj_array_clear (&spa->panels);
      g_clear_pointer (&spa->pack, g_free);

      if (var != NULL)
      {
        spa->pack   = _ncm_spectral_polar_panels_pack_from_variant (var);
        spa->n_sets = spa->pack->n_sets;
      }
      else
      {
        spa->panels = ncm_obj_array_new ();
        spa->n_sets = 0;
      }

      break;
    }
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_ncm_spectral_polar_panels_array_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmSpectralPolarPanelsArray *spa = NCM_SPECTRAL_POLAR_PANELS_ARRAY (object);

  g_return_if_fail (NCM_IS_SPECTRAL_POLAR_PANELS_ARRAY (object));

  switch (prop_id)
  {
    case PROP_ARRAY_PACKED:

      if (spa->n_sets == 0)
      {
        g_value_set_variant (value, NULL);
      }
      else if (spa->pack != NULL)
      {
        g_value_take_variant (value, _ncm_spectral_polar_panels_pack_to_variant (spa->pack));
      }
      else
      {
        /* A pack made for this only, so that the array still takes sets */
        NcmSpectralPolarPanelsPack *pack = _ncm_spectral_polar_panels_pack_new (spa->n_sets, (NcmSpectralPolarPanels **) spa->panels->pdata);

        g_value_take_variant (value, _ncm_spectral_polar_panels_pack_to_variant (pack));
        g_free (pack);
      }

      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_ncm_spectral_polar_panels_array_dispose (GObject *object)
{
  NcmSpectralPolarPanelsArray *spa = NCM_SPECTRAL_POLAR_PANELS_ARRAY (object);

  ncm_obj_array_clear (&spa->panels);
  g_clear_pointer (&spa->pack, g_free);

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_spectral_polar_panels_array_parent_class)->dispose (object);
}

static void
ncm_spectral_polar_panels_array_class_init (NcmSpectralPolarPanelsArrayClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);

  object_class->set_property = &_ncm_spectral_polar_panels_array_set_property;
  object_class->get_property = &_ncm_spectral_polar_panels_array_get_property;
  object_class->dispose      = &_ncm_spectral_polar_panels_array_dispose;

  /**
   * NcmSpectralPolarPanelsArray:packed:
   *
   * The sets packed, as a #GVariant of type `(uadadaua(dddduua(uuuad)))`: the number of
   * sets, the radial and the angular edges, the offsets of the panels of the initial panels
   * (as ncm_spectral_polar_panels_find_panel() reads them), and for each panel $r_a$, $r_b$,
   * $\theta_a$, $\theta_b$, the largest $K_r$ and $K_\theta$ of the sets there and the terms
   * of the sets whose block is not a single zero, in ascending sets: the set, its $K_r$ and
   * $K_\theta$ and its $K_\theta K_r$ coefficients by angular index. %NULL for no sets. An
   * array set from it is packed.
   */
  g_object_class_install_property (object_class,
                                   PROP_ARRAY_PACKED,
                                   g_param_spec_variant ("packed",
                                                         NULL,
                                                         "The packed sets",
                                                         G_VARIANT_TYPE (NCM_SPECTRAL_POLAR_PANELS_ARRAY_PACKED_TYPE), NULL,
                                                         G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
}

/**
 * ncm_spectral_polar_panels_array_new:
 *
 * Creates an empty set of #NcmSpectralPolarPanels.
 *
 * Returns: (transfer full): a new #NcmSpectralPolarPanelsArray.
 */
NcmSpectralPolarPanelsArray *
ncm_spectral_polar_panels_array_new (void)
{
  return g_object_new (NCM_TYPE_SPECTRAL_POLAR_PANELS_ARRAY, NULL);
}

/**
 * ncm_spectral_polar_panels_array_ref:
 * @spa: a #NcmSpectralPolarPanelsArray
 *
 * Increases the reference count of @spa by one.
 *
 * Returns: (transfer full): @spa.
 */
NcmSpectralPolarPanelsArray *
ncm_spectral_polar_panels_array_ref (NcmSpectralPolarPanelsArray *spa)
{
  return g_object_ref (spa);
}

/**
 * ncm_spectral_polar_panels_array_free:
 * @spa: a #NcmSpectralPolarPanelsArray
 *
 * Decreases the reference count of @spa by one.
 */
void
ncm_spectral_polar_panels_array_free (NcmSpectralPolarPanelsArray *spa)
{
  g_object_unref (spa);
}

/**
 * ncm_spectral_polar_panels_array_clear:
 * @spa: a #NcmSpectralPolarPanelsArray
 *
 * Decreases the reference count of *@spa by one and sets *@spa to %NULL.
 */
void
ncm_spectral_polar_panels_array_clear (NcmSpectralPolarPanelsArray **spa)
{
  g_clear_object (spa);
}

/**
 * ncm_spectral_polar_panels_array_add:
 * @spa: a #NcmSpectralPolarPanelsArray
 * @spp: a #NcmSpectralPolarPanels
 *
 * Appends @spp to @spa, referencing it. @spa must not be packed
 * (ncm_spectral_polar_panels_array_pack()).
 */
void
ncm_spectral_polar_panels_array_add (NcmSpectralPolarPanelsArray *spa, NcmSpectralPolarPanels *spp)
{
  g_assert (NCM_IS_SPECTRAL_POLAR_PANELS (spp));
  g_assert_nonnull (spa->panels);
  ncm_obj_array_add (spa->panels, G_OBJECT (spp));
  spa->n_sets++;
}

/**
 * ncm_spectral_polar_panels_array_len:
 * @spa: a #NcmSpectralPolarPanelsArray
 *
 * Returns: the number of #NcmSpectralPolarPanels in @spa.
 */
guint
ncm_spectral_polar_panels_array_len (NcmSpectralPolarPanelsArray *spa)
{
  return spa->n_sets;
}

/**
 * ncm_spectral_polar_panels_array_pack:
 * @spa: a #NcmSpectralPolarPanelsArray
 *
 * Packs the sets of @spa in one allocation (ncm_spectral_polar_panels_array_peek_pack())
 * and releases them, so that @spa keeps its coefficients once. The sets must have the same
 * edges and panels, and there must be at least one. This changes @spa: it must be called
 * before @spa is shared with other threads. Packing a packed array does nothing.
 */
void
ncm_spectral_polar_panels_array_pack (NcmSpectralPolarPanelsArray *spa)
{
  if (spa->pack != NULL)
    return;

  g_assert_cmpuint (spa->n_sets, >, 0);

  spa->pack = _ncm_spectral_polar_panels_pack_new (spa->n_sets, (NcmSpectralPolarPanels **) spa->panels->pdata);
  ncm_obj_array_clear (&spa->panels);
}

/**
 * ncm_spectral_polar_panels_array_peek_pack: (skip)
 * @spa: a packed #NcmSpectralPolarPanelsArray
 *
 * Gets the sets of @spa packed in one allocation, for
 * ncm_spectral_polar_panels_pack_eval_on_array_combine(). @spa must be packed
 * (ncm_spectral_polar_panels_array_pack()); the pack is freed with @spa.
 *
 * Returns: (transfer none): the pack of @spa.
 */
const NcmSpectralPolarPanelsPack *
ncm_spectral_polar_panels_array_peek_pack (NcmSpectralPolarPanelsArray *spa)
{
  g_assert_nonnull (spa->pack);

  return spa->pack;
}

