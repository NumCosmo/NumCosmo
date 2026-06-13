"""Mock catalog generation for halos, clusters and galaxy members.

This module provides :class:`MockGenerator`, which samples synthetic catalogs of
halos and clusters (optionally from a halo mass function) together with galaxy
members following a Halo Occupation Distribution. It also provides small,
well-defined selection models (:class:`CompletenessModel`, :class:`PurityModel`)
used to inject detection incompleteness and impurity into the mocks.

Geometry and cosmology helpers live in :mod:`numcosmo_py.catalog._geometry` and
:mod:`numcosmo_py.catalog._cosmology`, which wrap NumCosmo primitives. The only
heavy/uncommon dependency, ``scipy``, is imported lazily inside the one method
that needs it, so importing this module -- and therefore
:mod:`numcosmo_py.catalog` -- stays cheap.
"""

from __future__ import annotations

from typing import Protocol

import numpy as np
import numpy.typing as npt
from astropy.table import Table, vstack

from numcosmo_py import Nc, Ncm

from . import sky_match
from . import _cosmology, _geometry


class CompletenessModel(Protocol):
    """Detection-completeness model.

    A callable returning, for each object, the probability (0 to 1) that it is
    detected, as a function of ``log_mass`` (natural log of the mass) and ``z``.
    """

    def __call__(
        self, log_mass: npt.NDArray[np.float64], z: npt.NDArray[np.float64]
    ) -> npt.NDArray[np.float64]: ...


class PurityModel(Protocol):
    """Detection-purity model.

    A callable returning, for each object, the probability (0 to 1) that a
    detection is real, as a function of ``log_mass`` (natural log of the mass)
    and ``z``. Lower values inject more spurious (fake) detections.
    """

    def __call__(
        self, log_mass: npt.NDArray[np.float64], z: npt.NDArray[np.float64]
    ) -> npt.NDArray[np.float64]: ...


class ConstantCompleteness:
    """Completeness model with a single, mass- and redshift-independent value.

    :param completeness: detection probability applied to every object
        (default ``1.0``, i.e. a fully complete catalog).
    """

    def __init__(self, completeness: float = 1.0) -> None:
        """Create a constant completeness model."""
        self.completeness = completeness

    def __call__(
        self, log_mass: npt.NDArray[np.float64], z: npt.NDArray[np.float64]
    ) -> npt.NDArray[np.float64]:
        """Return the constant detection probability for each object."""
        del z
        return np.full_like(np.asarray(log_mass, dtype=float), self.completeness)


class ConstantPurity:
    """Purity model with a single, mass- and redshift-independent value.

    :param purity: probability that a detection is real, applied to every object
        (default ``1.0``, i.e. a fully pure catalog).
    """

    def __init__(self, purity: float = 1.0) -> None:
        """Create a constant purity model."""
        self.purity = purity

    def __call__(
        self, log_mass: npt.NDArray[np.float64], z: npt.NDArray[np.float64]
    ) -> npt.NDArray[np.float64]:
        """Return the constant purity for each object."""
        del z
        return np.full_like(np.asarray(log_mass, dtype=float), self.purity)


def identity_scaling_relation(
    log_mass_obs: npt.NDArray[np.float64], z: npt.NDArray[np.float64]
) -> npt.NDArray[np.float64]:
    """Trivial mass-observable scaling relation returning the observed mass.

    :param log_mass_obs: observed natural-log mass.
    :param z: redshift (unused).
    :return: ``log_mass_obs`` unchanged.
    """
    del z
    return log_mass_obs


class MockGenerator:
    """Class to generate clusters, halos and galaxy members."""

    C_LIGHT = Ncm.C.c()

    def __init__(
        self,
        cosmo: Nc.HICosmo,
        halo_set_size: int | None = 200,
        cluster_set_size: int | None = 100,
        ra_interval: tuple[float, float] = (-10, 10),
        dec_interval: tuple[float, float] = (-10, 10),
        z_interval: tuple[float, float] = (0.2, 0.5),
        cluster_mass_interval: tuple[float, float] = (1e14, 1e15),
        halo_mass_interval: tuple[float, float] = (1e13, 1e14),
        hmf: Nc.HaloMassFunction | None = None,
        cluster_m: Nc.ClusterMass | None = None,
        seed: int | None = 42,
        stochastic_central: bool = True,
    ) -> None:
        """Create a new MockGenerator object."""
        self.cosmo = cosmo
        self.halo_set_size = halo_set_size
        self.cluster_set_size = cluster_set_size
        self.ra_interval = ra_interval
        self.dec_interval = dec_interval
        self.z_interval = z_interval
        self.cluster_mass_interval = cluster_mass_interval
        self.halo_mass_interval = halo_mass_interval
        self.hmf = hmf
        self.cluster_m = cluster_m
        self.seed = seed
        self.stochastic_central = stochastic_central
        self._dist: Nc.Distance | None = None
        if self.seed is not None:
            np.random.seed(self.seed)

    @property
    def ra_min(self) -> float:
        """:return float: The minimum RA of the mock catalog."""
        return self.ra_interval[0]

    @property
    def ra_max(self) -> float:
        """:return float: The maximum RA of the mock catalog."""
        return self.ra_interval[1]

    @property
    def dec_min(self) -> float:
        """:return float: The minimum DEC of the mock catalog."""
        return self.dec_interval[0]

    @property
    def dec_max(self) -> float:
        """:return float: The maximum DEC of the mock catalog."""
        return self.dec_interval[1]

    @property
    def z_min(self) -> float:
        """:return float: The minimum redshift of the mock catalog."""
        return self.z_interval[0]

    @property
    def z_max(self) -> float:
        """:return float: The maximum redshift of the mock catalog."""
        return self.z_interval[1]

    @property
    def cluster_mass_min(self) -> float:
        """:return float: The minimum cluster mass of the mock catalog."""
        return self.cluster_mass_interval[0]

    @property
    def cluster_mass_max(self) -> float:
        """:return float: The maximum cluster mass of the mock catalog."""
        return self.cluster_mass_interval[1]

    @property
    def halo_mass_min(self) -> float:
        """:return float: The minimum halo mass of the mock catalog."""
        return self.halo_mass_interval[0]

    @property
    def halo_mass_max(self) -> float:
        """:return float: The maximum halo mass of the mock catalog."""
        return self.halo_mass_interval[1]

    def sky_area(self):
        """Compute the survey footprint area in square degrees."""
        return _geometry.rectangular_sky_area(
            self.ra_min, self.ra_max, self.dec_min, self.dec_max
        )

    def generate_cluster_positions(self):
        """Generate random cluster positions within the specified RA, Dec, and redshift
        intervals.

        :return tuple[ndarray, ndarray, ndarray]: cluster_ra: cluster right ascension
            (RA), cluster_dec: cluster declination (Dec), cluster_z: redshift (z)
        """

        cluster_ra = np.random.uniform(self.ra_min, self.ra_max, self.cluster_set_size)

        cluster_sin_dec = np.random.uniform(
            np.sin(np.radians(self.dec_min)),
            np.sin(np.radians(self.dec_max)),
            self.cluster_set_size,
        )

        cluster_dec = np.degrees(np.arcsin(cluster_sin_dec))

        cluster_z = np.random.uniform(self.z_min, self.z_max, self.cluster_set_size)

        return cluster_ra, cluster_dec, cluster_z

    def generate_halos_from_hmf(self, completeness_model=None):
        """Generate a synthetic halo catalog by sampling a Halo Mass Function (HMF).

        Orchestrates the stages: draw the halo population from the HMF
        (:meth:`_sample_halo_abundance`), assign sky positions
        (:meth:`_assign_positions`), compute geometry/R200c, apply the detection
        completeness (:meth:`_apply_completeness`) and package the result.

        :param completeness_model: optional callable ``f(logM, z)`` returning a
            detection probability in [0, 1]. If ``None`` every halo is detected.
        :return: astropy Table with IDs, coordinates (RA, Dec, x1/x2/x3), redshift,
            log-mass (true and observed), R200c, comoving radius and detection flag.
        """
        # Core NumCosmo step: realize the halo population (true z, true/obs lnM).
        halo_z, halo_logm, halo_mobs = self._sample_halo_abundance()
        halo_id = np.arange(200000, 200000 + self.halo_set_size)

        # Decoration: positions, geometry, derived properties, detection.
        halo_ra, halo_dec, halo_r = self._assign_positions(self.halo_set_size, halo_z)
        halo_x1, halo_x2, halo_x3 = self.get_3D_coordinates(halo_ra, halo_dec, halo_r)
        halos_R200c = self.get_R200c(np.exp(halo_logm), halo_z)
        is_detected = self._apply_completeness(halo_logm, halo_z, completeness_model)

        halos = Table(
            [
                halo_id,
                halo_ra,
                halo_dec,
                halo_z,
                halo_logm,
                halo_mobs,
                halos_R200c,
                halo_x1,
                halo_x2,
                halo_x3,
                halo_r,
                is_detected,
            ],
            names=(
                "halo_id",
                "RA",
                "DEC",
                "z",
                "Mass",
                "Mass_obs",
                "R200c",
                "x1",
                "x2",
                "x3",
                "halo_r",
                "is_detected",
            ),
            dtype=(
                int,
                float,
                float,
                float,
                float,
                float,
                float,
                float,
                float,
                float,
                float,
                int,
            ),
        )
        return halos[
            (halos["Mass"] >= np.log(self.halo_mass_min))
            & (halos["Mass"] <= np.log(self.halo_mass_max))
        ]

    def _sample_halo_abundance(self):
        """Sample the halo population from the HMF via NumCosmo's cluster abundance.

        This is the core step: it realizes the number of objects over the survey
        volume and mass range and returns their true redshift, true log-mass and
        observed log-mass. Also sets :attr:`halo_set_size`.

        :return: tuple ``(z_true, lnM_true, lnM_obs)`` numpy arrays.
        """
        sky_area_rad = self.sky_area() * (np.pi / 180) ** 2  # deg^2 -> steradians
        cluster_z = Nc.ClusterRedshiftNodist(z_min=self.z_min, z_max=self.z_max)

        self.hmf.prepare(self.cosmo)
        self.hmf.set_area(sky_area_rad)

        # Cluster abundance integrates the HMF over the volume and mass range.
        cad = Nc.ClusterAbundance.new(self.hmf, None)
        cad.set_area(sky_area_rad)
        cad.prepare(self.cosmo, cluster_z, self.cluster_m)

        rng = Ncm.RNG.seeded_new(None, self.seed)
        mset = Ncm.MSet.new_array([self.cosmo, self.cluster_m, cluster_z])
        ncount = Nc.DataClusterNCount.new(
            cad,
            "NcClusterRedshiftNodist",
            "Nc" + str(type(self.cluster_m)).rsplit(".", maxsplit=1)[-1].strip("'>"),
        )
        ncount.init_from_sampling(mset, sky_area_rad, rng)

        self.halo_set_size = ncount.get_lnM_true().len()
        z_true = np.array(ncount.get_z_true().dup_array())
        lnm_true = np.array(ncount.get_lnM_true().dup_array())
        lnm_obs = np.array(ncount.get_lnM_obs().dup_array())
        return z_true, lnm_true, lnm_obs

    def _assign_positions(self, n, z):
        """Assign random sky positions within the footprint and comoving radii.

        Samples RA uniformly and DEC uniformly in ``sin(DEC)`` (uniform on the
        sphere patch), and computes the comoving distance for each redshift.

        :param n: number of objects.
        :param z: per-object redshift array.
        :return: tuple ``(ra, dec, r)`` numpy arrays (degrees, degrees, Mpc).
        """
        ra = np.random.uniform(self.ra_min, self.ra_max, n)
        sin_dec = np.random.uniform(
            np.sin(np.radians(self.dec_min)),
            np.sin(np.radians(self.dec_max)),
            n,
        )
        dec = np.degrees(np.arcsin(sin_dec))
        r = np.array(self.dist().comoving_array(self.cosmo, z)) * self.cosmo.RH_Mpc()
        return ra, dec, r

    def _apply_completeness(self, log_mass, z, completeness_model):
        """Realize the detection flag from a completeness model.

        :param log_mass: per-object natural-log mass.
        :param z: per-object redshift.
        :param completeness_model: callable ``f(logM, z) -> probability`` or ``None``.
        :return: integer detection flag array (all ones when no model is given).
        """
        if completeness_model is None:
            return np.ones_like(log_mass)
        probs = completeness_model(log_mass, z)
        return np.random.binomial(n=1, p=probs)

    def generate_clusters_from_halos(
        self, halos, D_DIM=2.0, purity_model=None, scaling_relation=None
    ):
        """Generate clusters from an HMF-sampled halo catalog.

        Real detections are built by jittering each detected halo's position
        (:meth:`_clusters_from_detected_halos`). If a ``purity_model`` is given,
        spurious detections are sampled and appended (:meth:`_inject_fake_clusters`).

        :param Table halos: halo catalog (from :meth:`generate_halos_from_hmf`).
        :param float D_DIM: position-jitter half-width applied to detected halos.
        :param purity_model: optional callable ``f(logM, z) -> purity`` injecting fakes.
        :param scaling_relation: callable ``f(logM_obs, z) -> logM_true`` for fakes.
        :return: astropy Table of clusters (``parent_id == 0`` for fakes).
        """
        parts = self._clusters_from_detected_halos(halos, D_DIM)
        n_real = len(parts["z"])

        if purity_model is None:
            self.cluster_set_size = n_real
        else:
            fakes, n_fake = self._inject_fake_clusters(purity_model, scaling_relation)
            self.cluster_set_size = n_real + n_fake
            parts = {key: np.append(parts[key], fakes[key]) for key in parts}

        cluster_R200c = self.get_R200c(np.exp(parts["logm_true"]), parts["z"])
        cluster_id = np.array([int(i + 100000) for i in range(self.cluster_set_size)])

        return Table(
            [
                cluster_id,
                parts["ra"],
                parts["dec"],
                parts["z"],
                parts["logm_obs"],
                parts["logm_true"],
                cluster_R200c,
                parts["x1"],
                parts["x2"],
                parts["x3"],
                parts["r"],
                parts["parent_id"],
            ],
            names=(
                "cluster_id",
                "RA",
                "DEC",
                "z",
                "Mass_obs",
                "Mass",
                "R200c",
                "x1",
                "x2",
                "x3",
                "cluster_r",
                "parent_id",
            ),
            dtype=(
                int,
                float,
                float,
                float,
                float,
                float,
                float,
                float,
                float,
                float,
                float,
                int,
            ),
        )

    def _clusters_from_detected_halos(self, halos, D_DIM):
        """Build real cluster detections by jittering detected halo positions.

        Each detected halo (``is_detected == 1``) gets a cluster whose Cartesian
        position is offset by a uniform jitter of half-width ``D_DIM``; RA/DEC/z are
        recomputed from it and ``parent_id`` is set to the halo id.

        :return: dict of per-cluster arrays (keys: ra, dec, z, r, x1, x2, x3,
            logm_obs, logm_true, parent_id).
        """
        detected = halos[halos["is_detected"] == 1]
        n = len(detected["z"])

        x1 = detected["x1"] + np.random.uniform(-D_DIM, D_DIM, n)
        x2 = detected["x2"] + np.random.uniform(-D_DIM, D_DIM, n)
        x3 = detected["x3"] + np.random.uniform(-D_DIM, D_DIM, n)

        ra = np.degrees(np.arctan2(x2, x1))
        dec = np.degrees(np.arcsin(x3 / detected["halo_r"]))
        r = np.sqrt(x1**2 + x2**2 + x3**2)
        z = np.array(
            [self.dist().inv_comoving(self.cosmo, ri / self.cosmo.RH_Mpc()) for ri in r]
        )
        return {
            "ra": ra,
            "dec": dec,
            "z": z,
            "r": r,
            "x1": x1,
            "x2": x2,
            "x3": x3,
            "logm_obs": np.asarray(detected["Mass_obs"]),
            "logm_true": np.asarray(detected["Mass"]),
            "parent_id": np.asarray(detected["halo_id"]),
        }

    def _inject_fake_clusters(self, purity_model, scaling_relation):
        """Sample spurious (fake) cluster detections implied by an impurity model.

        The fake-detection rate is ``hmf * (1/purity - 1)``; its integral over the
        ``(z, lnM)`` range gives the expected number of fakes (Poisson realized),
        which are rejection-sampled and placed uniformly in the footprint.

        :return: tuple ``(fakes, n_fake)`` where ``fakes`` is a dict of per-cluster
            arrays (same keys as :meth:`_clusters_from_detected_halos`, ``parent_id``
            all zero).
        """
        from scipy.integrate import dblquad

        sky_area_rad = self.sky_area() * (np.pi / 180) ** 2  # deg^2 -> steradians

        if self.cluster_m is None or isinstance(self.cluster_m, Nc.ClusterMassNodist):
            self.hmf.set_area(sky_area_rad)

            def hmf_fake(logm, z):
                return self.hmf.d2n_dzdlnM(self.cosmo, logm, z) * (
                    1 / purity_model(logm, z) - 1
                )

            def hmf_fake_pdf(logm, z):
                return np.log(
                    self.hmf.d2n_dzdlnM(self.cosmo, logm, z)
                    * (1 / purity_model(logm, z) - 1)
                    / mean_fake_clusters_size
                )

        else:
            clusterz = Nc.ClusterRedshiftNodist(z_min=self.z_min, z_max=self.z_max)
            cad = Nc.ClusterAbundance.new(self.hmf, None)
            cad.set_area(sky_area_rad)
            cad.prepare(self.cosmo, clusterz, self.cluster_m)

            def hmf_fake(logm, z):
                return cad.lnM_p_d2n(
                    self.cosmo, clusterz, self.cluster_m, [logm], None, z
                ) * (1 / purity_model(logm, z) - 1)

            def hmf_fake_pdf(logm, z):
                return np.log(
                    cad.lnM_p_d2n(self.cosmo, clusterz, self.cluster_m, [logm], None, z)
                    * (1 / purity_model(logm, z) - 1)
                    / mean_fake_clusters_size
                )

        mean_fake_clusters_size = dblquad(
            hmf_fake,
            self.z_min,
            self.z_max,
            np.log(self.cluster_mass_min),
            np.log(self.cluster_mass_max),
        )[0]
        n_fake = np.random.poisson(lam=mean_fake_clusters_size)

        if n_fake == 0:
            empty = np.array([])
            keys = ("ra", "dec", "z", "r", "x1", "x2", "x3", "logm_obs", "logm_true")
            fakes = {key: empty for key in keys}
            fakes["parent_id"] = empty
            return fakes, 0

        # ARRUMAR ESSE SAMPLER -- faster unbinned alternative to MCMC.
        def calculate_approx_max(target_func):
            z_grid = np.linspace(self.z_min, self.z_max, 50)
            logm_grid = np.linspace(
                np.log(self.cluster_mass_min), np.log(self.cluster_mass_max), 50
            )
            zz, mm = np.meshgrid(z_grid, logm_grid)
            # List comprehension: target_func may not be fully vectorized.
            vals = np.array(
                [target_func(m, z) for m, z in zip(mm.flatten(), zz.flatten())]
            )
            return np.max(vals)

        def sample_fakes_rejection(target_func, n_to_generate):
            sampled = []
            max_val = calculate_approx_max(target_func)
            while len(sampled) < n_to_generate:
                m_test = np.random.uniform(
                    np.log(self.cluster_mass_min), np.log(self.cluster_mass_max)
                )
                z_test = np.random.uniform(self.z_min, self.z_max)
                if np.random.rand() < target_func(m_test, z_test) / max_val:
                    sampled.append([m_test, z_test])
            return np.array(sampled)

        sampled = sample_fakes_rejection(hmf_fake_pdf, n_fake)
        logm_fake = sampled[:, 0]
        z_fake = sampled[:, 1]
        logm_true_fake = scaling_relation(logm_fake, z_fake)

        ra_fake = np.random.uniform(self.ra_min, self.ra_max, n_fake)
        sin_dec_fake = np.random.uniform(
            np.sin(np.radians(self.dec_min)),
            np.sin(np.radians(self.dec_max)),
            n_fake,
        )
        dec_fake = np.degrees(np.arcsin(sin_dec_fake))
        r_fake = (
            np.array(self.dist().comoving_array(self.cosmo, z_fake))
            * self.cosmo.RH_Mpc()
        )
        x1_fake, x2_fake, x3_fake = self.get_3D_coordinates(ra_fake, dec_fake, r_fake)

        fakes = {
            "ra": ra_fake,
            "dec": dec_fake,
            "z": z_fake,
            "r": r_fake,
            "x1": x1_fake,
            "x2": x2_fake,
            "x3": x3_fake,
            "logm_obs": logm_fake,
            "logm_true": logm_true_fake,
            "parent_id": np.zeros(n_fake),
        }
        return fakes, n_fake

    def generate_cluster_logm(self):
        """Generate random log10(masses) within the specified mass interval."""

        cluster_logm = np.random.uniform(
            np.log10(self.cluster_mass_min),
            np.log10(self.cluster_mass_max),
            self.cluster_set_size,
        )

        return cluster_logm

    def get_3D_coordinates(self, RA, DEC, R):
        """Get positions in 3D coordinates (NumCosmo astro convention)."""
        return _geometry.spherical_to_cartesian(RA, DEC, R)

    def rho_crit(self, z):
        """Critical density of the Universe at redshift ``z`` (solar masses / Mpc^3)."""
        return _cosmology.critical_density(self.cosmo, z)

    def get_R200c(self, mass, z):
        """Calculate R_{200c} for a given mass and redshift."""
        return _cosmology.r200c(self.cosmo, mass, z)

    def dist(self, zf=100.0):
        """Return a prepared NumCosmo distance object, built once and reused."""
        if self._dist is None:
            self._dist = Nc.Distance.new(zf)
            self._dist.compute_inv_comoving(True)
            self._dist.prepare(self.cosmo)
        return self._dist

    def generate_clusters(self):
        """Generate clusters.

        :return Table: clusters: Astropy Table with generated clusters"""

        # Generate cluster positions, redshifts, and masses
        cluster_ra, cluster_dec, cluster_z = self.generate_cluster_positions()
        cluster_logm = self.generate_cluster_logm()

        clusters_R200c = self.get_R200c(10**cluster_logm, cluster_z)

        dist = self.dist()
        cluster_r = (
            np.array(dist.comoving_array(self.cosmo, cluster_z)) * self.cosmo.RH_Mpc()
        )

        # Compute the 3D coordinates of the clusters
        cluster_x1, cluster_x2, cluster_x3 = self.get_3D_coordinates(
            cluster_ra, cluster_dec, cluster_r
        )

        # Create the cluster ID
        cluster_id = np.array([int(i + 100000) for i in range(len(cluster_z))])

        # Table with cluster properties
        clusters = Table(
            [
                cluster_id,
                cluster_ra,
                cluster_dec,
                cluster_z,
                cluster_logm,
                clusters_R200c,
                cluster_x1,
                cluster_x2,
                cluster_x3,
                cluster_r,
            ],
            names=(
                "cluster_id",
                "RA",
                "DEC",
                "z",
                "Mass",
                "R200c",
                "x1",
                "x2",
                "x3",
                "cluster_r",
            ),
            dtype=(int, float, float, float, float, float, float, float, float, float),
        )

        return clusters

    def generate_halos(self, D_DIM=2.0, clusters=None):
        """Generate halos.

        :param float D_DIM: size of the region around each cluster where the main halos
            are generated.

        :return Table: halos: Astropy Table with generated halos."""

        if clusters is None:
            clusters = self.generate_clusters()

        dist = self.dist()

        cluster_x1 = clusters["x1"]
        cluster_x2 = clusters["x2"]
        cluster_x3 = clusters["x3"]
        cluster_r = clusters["cluster_r"]
        cluster_logm = clusters["Mass"]

        # Generate halo positions, redshifts, and masses around the clusters
        halo_x1 = cluster_x1 + np.random.uniform(-D_DIM, D_DIM, self.cluster_set_size)
        halo_x2 = cluster_x2 + np.random.uniform(-D_DIM, D_DIM, self.cluster_set_size)
        halo_x3 = cluster_x3 + np.random.uniform(-D_DIM, D_DIM, self.cluster_set_size)
        halo_ra = np.degrees(np.arctan2(halo_x2, halo_x1))
        halo_dec = np.degrees(np.arcsin(halo_x3 / cluster_r))
        halo_r = np.sqrt(halo_x1**2 + halo_x2**2 + halo_x3**2)
        halo_z = [
            dist.inv_comoving(self.cosmo, r / self.cosmo.RH_Mpc()) for r in halo_r
        ]

        # For the halo masses we use the cluster's masses added a Gaussian noise
        halo_logm = cluster_logm + np.random.normal(0, 0.1, self.cluster_set_size)

        # Add more halos randomly
        DELTA_OBJECTS = self.halo_set_size - self.cluster_set_size

        halo_ra = np.append(
            halo_ra, np.random.uniform(self.ra_min, self.ra_max, DELTA_OBJECTS)
        )
        halo_dec = np.append(
            halo_dec, np.random.uniform(self.dec_min, self.dec_max, DELTA_OBJECTS)
        )
        halo_z = np.append(
            halo_z, np.random.uniform(self.z_min, self.z_max, DELTA_OBJECTS)
        )
        halo_logm = np.append(
            halo_logm,
            np.random.uniform(
                np.log10(self.halo_mass_min),
                np.log10(self.halo_mass_max),
                DELTA_OBJECTS,
            ),
        )
        # Compute R200 for halos and clusters
        halos_R200c = self.get_R200c(10**halo_logm, halo_z)

        # Create the halo ID
        halo_id = np.array([int(i + 200000) for i in range(len(halo_z))])

        # Table with halos properties
        halos = Table(
            [
                halo_id,
                halo_ra,
                halo_dec,
                halo_z,
                halo_logm,
                halos_R200c,
            ],  # halo_x1, halo_x2, halo_x3, halo_r
            names=(
                "halo_id",
                "RA",
                "DEC",
                "z",
                "Mass",
                "R200c",
            ),  # , "x1", "x2", "x3", "halo_r"
            dtype=(int, float, float, float, float, float),
        )  # float, float, float, float

        return halos

    def get_galaxy_coords(
        self, ra_c: float, dec_c: float, sep_angular_rad: float, phi_rad: float
    ) -> tuple[float, float]:
        """
        Compute the RA and DEC of a galaxy given the angular separation, position angle
        and center coordinates of the cluster/halo.

        :param float ra_c: cluster/halo central right ascension angle in degrees.
        :param float dec_c: cluster/halo central declination angle in degrees.
        :param float sep_angular_rad: galaxy angular separation from the center in
            radians.
        :param float phi_rad: galaxy position angle (angle from north to east) in
            radians.

        :return tuple[float, float]: ra: galaxy right ascension angle in degrees, dec:
            galaxy declination angle in degrees
        """

        # Convert center coordinates to radians
        ra_c_rad = np.radians(ra_c)
        dec_c_rad = np.radians(dec_c)

        # Compute the declination (delta)
        # sin(delta) = sin(dec_c)cos(sep) + cos(dec_c)sin(sep)cos(phi)
        sin_dec = np.sin(dec_c_rad) * np.cos(sep_angular_rad) + np.cos(
            dec_c_rad
        ) * np.sin(sep_angular_rad) * np.cos(phi_rad)

        dec_gal_rad = np.arcsin(sin_dec)

        # Compute the right ascension (alpha)
        # We use the y and x terms for the arctan2(y, x)
        y = np.sin(phi_rad) * np.sin(sep_angular_rad) * np.cos(dec_c_rad)
        x = np.cos(sep_angular_rad) - np.sin(dec_c_rad) * np.sin(dec_gal_rad)

        ra_gal_rad = ra_c_rad + np.arctan2(y, x)

        # Convert back to degrees and normalize RA to [0, 360]
        ra_gal = np.degrees(ra_gal_rad) % 360
        dec_gal = np.degrees(dec_gal_rad)

        return ra_gal, dec_gal

    def HOD_model(
        self,
        m_halo,
        logMmin=12.72,
        sigma_logM=0.26,
        alpha=1.15,
        logM1=13.93,
        logM0=12.7,
        stochastic_central: bool | None = None,
    ) -> tuple[int, int]:
        """
        Standard HOD model (Zheng et al. 2007)

        :param float m_halo: halo mass in solar masses.
        :param float logMmin: minimum log10(Mass) for central occupation.
        :param float sigma_logM: width of the transition region.
        :param float alpha: power law index for satellite occupation.
        :param float logM1: log10(Mass) at which satellite occupation starts.
        :param float logM0: log10(Mass) at which central occupation starts.
        :param stochastic_central: if ``True`` the central is a Bernoulli draw with
            probability ``<N_cen>``; if ``False`` it is present deterministically
            when ``<N_cen> >= 0.5``. Defaults to the instance's ``stochastic_central``.

        :return tuple[int, int]: n_cen: number of central galaxies (0 or 1), n_sat:
            number of satellite galaxies "
        """
        if stochastic_central is None:
            stochastic_central = self.stochastic_central

        # 1. Mean Central Occupancy (Error function)
        # 0.5 * (1 + erf(arg)) is the unit-normal CDF at arg * sqrt(2), i.e.
        # 0.5 + normal_gaussian_integral(0, arg * sqrt(2)).
        arg = (np.log10(m_halo) - logMmin) / sigma_logM
        mean_n_cen = 0.5 + Ncm.util_normal_gaussian_integral(0.0, arg * np.sqrt(2.0))
        if stochastic_central:
            n_cen = int(np.random.binomial(1, np.clip(mean_n_cen, 0.0, 1.0)))
        else:
            n_cen = int(mean_n_cen >= 0.5)

        # 2. Mean Satellite Occupancy (Power Law)
        n_sat = 0
        if n_cen == 1:
            # Satellites only exist if a central exists
            diff = m_halo - 10**logM0
            mean_n_sat = (np.maximum(0, diff) / 10**logM1) ** alpha
            if mean_n_sat > 0:
                n_sat = int(np.random.poisson(mean_n_sat))

        return n_cen, n_sat

    def generate_galaxies(self, catalog, object_type):
        """Generate galaxy members for each object in ``catalog`` via the HOD model."""

        if object_type == "halo":
            obj_id_base = 200000
            gal_id_base = 2000000
        elif object_type == "cluster":
            obj_id_base = 100000
            gal_id_base = 1000000
        else:
            raise ValueError("object_type deve ser 'halo' ou 'cluster'")

        dist = self.dist()
        all_galaxies = []

        for i in range(len(catalog)):
            n_cen, n_sat = self.HOD_model(10 ** catalog["Mass"][i])
            total_gals = n_cen + n_sat
            if total_gals == 0:
                continue

            object_ra = catalog["RA"][i]
            object_dec = catalog["DEC"][i]
            object_z = catalog["z"][i]
            object_r200 = catalog["R200c"][i]

            # Generate local coordinates
            # Central is at (0,0,0), Satellites follow NFW/Uniform profile
            galaxy_R = np.zeros(total_gals)

            if n_sat > 0:
                # Satellites distributed within R200
                galaxy_R[n_cen:] = object_r200 * np.random.uniform(0, 1, n_sat) ** (
                    1 / 3
                )

            # galaxy_R = object_r200 * np.random.uniform(0, 1, total_gals)**(1/3)
            cos_theta = np.random.uniform(-1, 1, total_gals)
            galaxy_theta = np.arccos(cos_theta)
            galaxy_phi = np.random.uniform(0, 2 * np.pi, total_gals)

            Hz = self.cosmo.H(object_z)

            distancia_radial = galaxy_R * np.cos(galaxy_theta)
            delta_z = (Hz / self.C_LIGHT) * distancia_radial
            z_galaxy = object_z + delta_z

            dA = dist.angular_diameter(self.cosmo, object_z) * self.cosmo.RH_Mpc()
            sep_angular = galaxy_R / dA

            ra_g, dec_g = self.get_galaxy_coords(
                object_ra, object_dec, sep_angular, galaxy_phi
            )

            ra_corrected = ra_g % 360.0
            ra_corrected = np.where(
                ra_corrected > 180, ra_corrected - 360, ra_corrected
            )
            ra_corrected = np.clip(ra_corrected, self.ra_min, self.ra_max)

            temp_table = Table()

            temp_table[f"{object_type}_id"] = np.full(total_gals, i + obj_id_base)
            temp_table["RA"] = ra_corrected
            temp_table["DEC"] = dec_g
            temp_table["z"] = z_galaxy
            temp_table["is_central"] = [
                True if j < n_cen else False for j in range(total_gals)
            ]
            temp_table[f"{object_type}_RA"] = np.full(total_gals, object_ra)
            temp_table[f"{object_type}_DEC"] = np.full(total_gals, object_dec)
            temp_table[f"{object_type}_z"] = np.full(total_gals, object_z)
            temp_table[f"{object_type}_mass"] = np.full(total_gals, catalog["Mass"][i])
            all_galaxies.append(temp_table)

        if all_galaxies:
            galaxy_catalog = vstack(all_galaxies)
            galaxy_catalog["galaxy_id"] = np.arange(len(galaxy_catalog)) + gal_id_base

            ordered_cols = [
                "galaxy_id",
                f"{object_type}_id",
                "is_central",
                "RA",
                "DEC",
                "z",
                f"{object_type}_z",
                f"{object_type}_mass",
            ]

            galaxy_catalog = galaxy_catalog[ordered_cols]

            print(
                f"Generated {len(galaxy_catalog)} galaxies "
                f"from {len(catalog)} {object_type}(s)."
            )

        else:
            galaxy_catalog = Table()
            print("No galaxies were generated. Check your HOD mass thresholds.")

        return galaxy_catalog

    def match_galaxies_to_objects(self, object_catalog, galaxy_catalog, object_type):
        """Match galaxies to objects using SkyMatch.

        :param Table object_catalog: Astropy Table with object properties.
        :param Table galaxy_catalog: Astropy Table with galaxy properties.

        :return Table: halo_galaxies: Updated Astropy Table with matched galaxies to
            objects.
        """

    def get_common_galaxies_by_distance(
        self, object_catalog, galaxy_catalog, object_type
    ):
        """Return galaxies lying within R200c of each object, matched by distance."""

        object_coordinates = {"RA": "RA", "DEC": "DEC", "z": "z"}
        galaxy_coordinates = {"RA": "RA", "DEC": "DEC", "z": "z"}

        obj_id = f"{object_type}_id"
        obj_mass = f"{object_type}_mass"
        obj_properties = {obj_id: obj_id, "Mass": obj_mass, "R200c": "R200c"}

        obj_m = sky_match.SkyMatch(
            query_data=object_catalog,
            query_coordinates=object_coordinates,
            match_data=galaxy_catalog,
            match_coordinates=galaxy_coordinates,
        )
        matched_gal = obj_m.match_3d(self.cosmo, 20)

        matched_gal_table = matched_gal.to_table_complete(
            query_properties=obj_properties
        )

        obj_z = f"{object_type}_z"

        matched_index = []
        matched_obj_id = []
        matched_obj_z = []
        matched_obj_mass = []
        for h7 in matched_gal_table:

            mask = h7["distances"] <= h7["R200c"]

            indices = h7["Index_matched"][mask]

            for idx in indices:
                matched_index.append(idx)
                matched_obj_id.append(h7[obj_id])
                matched_obj_z.append(h7["z"])
                matched_obj_mass.append(h7[obj_mass])

        # ``matched_index`` holds positional indices into ``galaxy_catalog``, so
        # gathering those rows reproduces the previous index-based pandas merge.
        gathered = galaxy_catalog[np.array(matched_index, dtype=int)]

        common_galaxies = Table()
        common_galaxies["galaxy_id"] = gathered["galaxy_id"]
        common_galaxies[obj_id] = matched_obj_id
        common_galaxies["RA"] = gathered["RA"]
        common_galaxies["DEC"] = gathered["DEC"]
        common_galaxies["z"] = gathered["z"]
        common_galaxies[obj_z] = matched_obj_z
        common_galaxies[obj_mass] = matched_obj_mass
        common_galaxies["is_central"] = np.full(len(gathered), False)

        return common_galaxies
