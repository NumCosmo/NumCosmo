/***************************************************************************
 *            ncm_powspec_analytic.h
 *
 *  Tue August 26 12:00:00 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) 2026 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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

#ifndef _NCM_POWSPEC_ANALYTIC_H_
#define _NCM_POWSPEC_ANALYTIC_H_

#include <glib-object.h>
#include <glib.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/ncm/powspec/ncm_powspec.h>

G_BEGIN_DECLS

#define NCM_TYPE_POWSPEC_ANALYTIC (ncm_powspec_analytic_get_type ())

G_DECLARE_FINAL_TYPE (NcmPowspecAnalytic, ncm_powspec_analytic, NCM, POWSPEC_ANALYTIC, NcmPowspec)

/**
 * NcmPowspecAnalyticShape:
 * @NCM_POWSPEC_ANALYTIC_SHAPE_POWER_LAW: $T(k) = 1$, so $P \propto k^{n_s}$
 * @NCM_POWSPEC_ANALYTIC_SHAPE_BBKS: the BBKS transfer function
 * @NCM_POWSPEC_ANALYTIC_SHAPE_RATIONAL: $T(k) = 1/(1 + (k/k_\mathrm{eq})^2)$
 *
 * Shape of the transfer function $T(k)$ in
 * $P(k,z) = A k^{n_s} T(k)^2 B(k) D(z)^2$.
 *
 * %NCM_POWSPEC_ANALYTIC_SHAPE_POWER_LAW has no scale of its own and is the
 * case that isolates everything else.
 *
 * %NCM_POWSPEC_ANALYTIC_SHAPE_BBKS is the default and the one that resembles a
 * real matter power spectrum, with $q = k/k_\mathrm{eq}$,
 * \begin{equation}
 *   T(q) = \frac{\ln(1 + 2.34q)}{2.34q}
 *          \Big[1 + 3.89q + (16.1q)^2 + (5.46q)^3 + (6.71q)^4\Big]^{-1/4} .
 * \end{equation}
 * It carries both features a $k$-quadrature benchmark cares about: a turnover,
 * and a $k^{n_s-4}\ln^2 k$ tail rather than a pure power law.
 *
 * %NCM_POWSPEC_ANALYTIC_SHAPE_RATIONAL is a deliberate contrast: it has a
 * turnover but a clean, log-free $k^{n_s-4}$ tail. It is *not* a good
 * imitation of a matter power spectrum -- fitted against one it is off by a
 * factor of five across $k \in [10^{-5}, 10]\,\mathrm{Mpc}^{-1}$ -- and exists
 * so that the effect of the tail's slope can be measured rather than assumed.
 *
 */
typedef enum _NcmPowspecAnalyticShape /*< prefix=NCM_POWSPEC_ANALYTIC_SHAPE >*/
{
  NCM_POWSPEC_ANALYTIC_SHAPE_POWER_LAW = 0,
  NCM_POWSPEC_ANALYTIC_SHAPE_BBKS,
  NCM_POWSPEC_ANALYTIC_SHAPE_RATIONAL,
  /*< private >*/
  NCM_POWSPEC_ANALYTIC_SHAPE_LEN,
} NcmPowspecAnalyticShape;

/**
 * NcmPowspecAnalyticGrowth:
 * @NCM_POWSPEC_ANALYTIC_GROWTH_NONE: $D(z) = 1$
 * @NCM_POWSPEC_ANALYTIC_GROWTH_LCDM: flat matter-plus-$\Lambda$ linear growth
 * @NCM_POWSPEC_ANALYTIC_GROWTH_RATIONAL: $D(a) \propto a(1 + (a/a_t)^3)^{-1/3}$
 *
 * Shape of the growth factor $D(z)$, normalized so that $D(0) = 1$.
 *
 * %NCM_POWSPEC_ANALYTIC_GROWTH_LCDM is the exact linear growth of a flat
 * universe with matter and a cosmological constant,
 * \begin{equation}
 *   D(a) \propto a \; {}_2F_1\!\left(\tfrac13, 1; \tfrac{11}{6};
 *                     -\frac{1-\Omega_m}{\Omega_m} a^3\right) ,
 * \end{equation}
 * which is a closed form and therefore reproducible in arbitrary precision. It
 * omits radiation, so it departs from #NcGrowthFunc by about $4\times10^{-3}$
 * at $z = 20$ and $2.5\times10^{-6}$ at $z = 0.1$.
 *
 * %NCM_POWSPEC_ANALYTIC_GROWTH_RATIONAL has the same two limits -- $D \propto
 * a$ in matter domination, $D \to$ constant at late times -- with no special
 * function, and is the contrast case.
 *
 */
typedef enum _NcmPowspecAnalyticGrowth /*< prefix=NCM_POWSPEC_ANALYTIC_GROWTH >*/
{
  NCM_POWSPEC_ANALYTIC_GROWTH_NONE = 0,
  NCM_POWSPEC_ANALYTIC_GROWTH_LCDM,
  NCM_POWSPEC_ANALYTIC_GROWTH_RATIONAL,
  /*< private >*/
  NCM_POWSPEC_ANALYTIC_GROWTH_LEN,
} NcmPowspecAnalyticGrowth;

NcmPowspecAnalytic *ncm_powspec_analytic_new (NcmPowspecAnalyticShape shape, NcmPowspecAnalyticGrowth growth);
NcmPowspecAnalytic *ncm_powspec_analytic_new_full (NcmPowspecAnalyticShape shape, NcmPowspecAnalyticGrowth growth, gdouble amplitude, gdouble n_s, gdouble k_eq, gdouble Omega_m);

NcmPowspecAnalytic *ncm_powspec_analytic_ref (NcmPowspecAnalytic *psa);
void ncm_powspec_analytic_free (NcmPowspecAnalytic *psa);
void ncm_powspec_analytic_clear (NcmPowspecAnalytic **psa);

NcmPowspecAnalyticShape ncm_powspec_analytic_get_shape (NcmPowspecAnalytic *psa);
NcmPowspecAnalyticGrowth ncm_powspec_analytic_get_growth (NcmPowspecAnalytic *psa);

gdouble ncm_powspec_analytic_get_amplitude (NcmPowspecAnalytic *psa);
gdouble ncm_powspec_analytic_get_n_s (NcmPowspecAnalytic *psa);
gdouble ncm_powspec_analytic_get_k_eq (NcmPowspecAnalytic *psa);
gdouble ncm_powspec_analytic_get_Omega_m (NcmPowspecAnalytic *psa);

gdouble ncm_powspec_analytic_eval_transfer (NcmPowspecAnalytic *psa, const gdouble k);
gdouble ncm_powspec_analytic_eval_bao (NcmPowspecAnalytic *psa, const gdouble k);
gdouble ncm_powspec_analytic_eval_growth (NcmPowspecAnalytic *psa, const gdouble z);

G_END_DECLS

#endif /* _NCM_POWSPEC_ANALYTIC_H_ */

