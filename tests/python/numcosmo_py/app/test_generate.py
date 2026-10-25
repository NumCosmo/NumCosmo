#!/usr/bin/env python
#
# test_py_generate.py
#
# Wed Jan 31 10:04:00 2024
# Copyright  2023  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# test_py_generate.py
# Copyright (C) 2024 Sandro Dias Pinto Vitenti <vitenti@uel.br>
#
# numcosmo is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# numcosmo is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program.  If not, see <http://www.gnu.org/licenses/>.

"""Test Generate classes."""

from typing import cast
from pathlib import Path
import pytest
import numpy as np
import typer

pytest.importorskip("astropy")
pytest.importorskip("getdist")
# flake8: noqa: E402


from numcosmo_py import Nc, Ncm
import numcosmo_py.app.generate as gen
import numcosmo_py.datasets.dist_mu_obs as dist_mu_obs
import numcosmo_py.datasets.hicosmo as hicosmo

Ncm.cfg_init()

pytestmark = pytest.mark.app


# JPAS Forecast
def test_generate_jpas(tmp_path: Path):
    """Test JPAS forecast generation with YAML file."""
    exp_file = tmp_path / "jpas_forecast.yaml"
    _ = gen.GenerateJpasForecast(experiment=exp_file.absolute())

    assert exp_file.exists()


def test_generate_jpas_invalid_suffix(tmp_path: Path):
    """Test JPAS forecast rejects invalid file suffix."""
    exp_file = tmp_path / "jpas_forecast.txt"
    with pytest.raises(ValueError, match="Invalid experiment file suffix: .txt"):
        _ = gen.GenerateJpasForecast(experiment=exp_file.absolute())


# Native Planck
@pytest.mark.planck_data
def test_generate_planck_native(tmp_path: Path):
    """Native Planck experiment generates all-native, clik-free serialized blocks."""
    from numcosmo_py.experiments.planck_lite import find_baseline_file
    from numcosmo_py.experiments.planck_commander import COMMANDER_RELPATH
    from numcosmo_py.experiments.planck18 import Planck18Types

    if find_baseline_file(COMMANDER_RELPATH) is None:
        pytest.skip("Planck baseline data not found")

    exp_file = tmp_path / "planck_native.yaml"
    _ = gen.GeneratePlanck(
        experiment=exp_file.absolute(),
        native=True,
        data_type=Planck18Types.TT,
        include_lens_lkl=True,
    )
    assert exp_file.exists()

    dset_file = exp_file.with_suffix(".dataset.gvar")
    assert dset_file.exists()

    # The serialized dataset must hold only native NumCosmo objects (no clik
    # DataPlanckLKL), so it reloads without the clik data or the PLC library.
    ser = Ncm.Serialize.new(Ncm.SerializeOpt.KEEP_SHARED)
    dset = ser.from_binfile(str(dset_file.absolute()))
    names = {dset.get_data(i).__class__.__name__ for i in range(dset.get_ndata())}
    assert names == {
        "DataPlanckSimall",
        "DataPlanckCommander",
        "DataPlanckSmica",
        "DataPlanckLensing",
    }

    # The reloaded blocks self-configure the shared (serialized) Boltzmann and
    # evaluate, with no clik data or PLC library involved.
    import numpy as np

    exp = ser.dict_str_from_yaml_file(str(exp_file.absolute()))
    mset = exp.peek("model-set")
    assert isinstance(mset, Ncm.MSet)
    for i in range(dset.get_ndata()):
        dset.get_data(i).prepare(mset)
    assert np.isfinite(dset.m2lnL_val(mset))


@pytest.mark.planck_data
def test_build_planck_release(tmp_path: Path):
    """The publish command writes serialized native Planck release objects."""
    from numcosmo_py.experiments.planck_lite import find_baseline_file
    from numcosmo_py.experiments.planck_commander import COMMANDER_RELPATH

    if find_baseline_file(COMMANDER_RELPATH) is None:
        pytest.skip("Planck baseline data not found")

    out_dir = tmp_path / "planck_release"
    _ = gen.BuildPlanckRelease(output_dir=out_dir.absolute())

    written = list(out_dir.glob("planck_native_*.gvar"))
    assert written


# QSpline
def test_generate_qspline(tmp_path: Path):
    """Test QSpline generation with all data types."""
    exp_file = tmp_path / "qspline.yaml"
    bao_id = hicosmo.BAOID.ALL_COMBINED_JUN_2025
    sne_id = hicosmo.SNIaID.COV_DES_Y5_STAT_SYS
    h_id = hicosmo.HID.ALL_COMBINED_APR_2025
    _ = gen.GenerateQSpline(
        experiment=exp_file.absolute(),
        include_snia=sne_id,
        include_hubble=h_id,
        include_bao=bao_id,
    )

    assert exp_file.exists()


def test_generate_qspline_no_data(tmp_path: Path):
    """Test QSpline rejects experiment with no data."""
    exp_file = tmp_path / "qspline.yaml"
    with pytest.raises(ValueError, match="No data included in the experiment."):
        _ = gen.GenerateQSpline(experiment=exp_file.absolute())


def test_generate_qspline_invalid_suffix(tmp_path: Path):
    """Test QSpline rejects invalid file suffix."""
    exp_file = tmp_path / "qspline.txt"
    bao_id = hicosmo.BAOID.ALL_COMBINED_JUN_2025
    with pytest.raises(ValueError, match="Invalid experiment file suffix: .txt"):
        _ = gen.GenerateQSpline(experiment=exp_file.absolute(), include_bao=bao_id)


# XCDM
def test_generate_xcdm(tmp_path: Path):
    """Test XCDM generation with multiple data sets."""
    exp_file = tmp_path / "xcdm.yaml"
    bao_ids = [hicosmo.BAOID.ALL_COMBINED_JUN_2025, hicosmo.BAOID.ALL_COMBINED_JAN_2023]
    sne_ids = [hicosmo.SNIaID.COV_DES_Y5_STAT_SYS, hicosmo.SNIaID.COV_DES_Y5_STATONLY]
    h_ids = [hicosmo.HID.ALL_COMBINED_APR_2025, hicosmo.HID.ALL_COMBINED_JAN_2023]

    for bao_id, sne_id, h_id in zip(bao_ids, sne_ids, h_ids):
        _ = gen.GenerateXCDM(
            experiment=exp_file.absolute(),
            include_snia=sne_id,
            include_hubble=h_id,
            include_bao=bao_id,
        )

    assert exp_file.exists()


def test_generate_xcdm_no_data(tmp_path: Path):
    """Test XCDM rejects experiment with no data."""
    exp_file = tmp_path / "xcdm.yaml"
    with pytest.raises(ValueError, match="No data included in the experiment."):
        _ = gen.GenerateXCDM(experiment=exp_file.absolute())


def test_generate_xcdm_no_valid_bao_id(tmp_path: Path):
    """Test XCDM rejects invalid BAO ID."""
    exp_file = tmp_path / "xcdm.yaml"
    bao_id = cast(hicosmo.BAOID, str(-1))
    with pytest.raises(ValueError, match="Unknown BAO data set id: -1"):
        _ = gen.GenerateXCDM(experiment=exp_file.absolute(), include_bao=bao_id)


def test_generate_xcdm_no_valid_hubble_id(tmp_path: Path):
    """Test XCDM rejects invalid Hubble ID."""
    exp_file = tmp_path / "xcdm.yaml"
    h_id = cast(hicosmo.HID, str(-1))
    with pytest.raises(ValueError, match="Unknown Hubble data set id: -1"):
        _ = gen.GenerateXCDM(experiment=exp_file.absolute(), include_hubble=h_id)


def test_generate_xcdm_invalid_suffix(tmp_path: Path):
    """Test XCDM rejects invalid file suffix."""
    exp_file = tmp_path / "xcdm.txt"
    bao_id = hicosmo.BAOID.ALL_COMBINED_JUN_2025
    with pytest.raises(ValueError, match="Invalid experiment file suffix: .txt"):
        _ = gen.GenerateXCDM(experiment=exp_file.absolute(), include_bao=bao_id)


# DEWSpline
# QSpline
def test_generate_dewspline(tmp_path: Path):
    """Test DEWSpline generation with all data types."""
    exp_file = tmp_path / "dewspline.yaml"
    bao_id = hicosmo.BAOID.ALL_COMBINED_JUN_2025
    sne_id = hicosmo.SNIaID.COV_DES_Y5_STAT_SYS
    h_id = hicosmo.HID.ALL_COMBINED_APR_2025
    _ = gen.GenerateDEWSpline(
        experiment=exp_file.absolute(),
        include_snia=sne_id,
        include_hubble=h_id,
        include_bao=bao_id,
    )

    assert exp_file.exists()


@pytest.mark.parametrize("use_as_drag", [False, True])
def test_generate_dewspline_use_as_drag(tmp_path: Path, use_as_drag: bool):
    """With use_as_drag the BAO sound horizon is the free asdrag parameter."""
    exp_file = tmp_path / "dewspline.yaml"
    _ = gen.GenerateDEWSpline(
        experiment=exp_file.absolute(),
        include_bao=hicosmo.BAOID.ALL_COMBINED_JUN_2025,
        use_as_drag=use_as_drag,
    )

    # The experiment refers to the dataset saved beside it, which must be loaded first.
    ser = Ncm.Serialize.new(Ncm.SerializeOpt.KEEP_SHARED)
    _ = ser.from_binfile(exp_file.with_suffix(".dataset.gvar").absolute().as_posix())
    exp = ser.dict_str_from_yaml_file(exp_file.absolute().as_posix())
    mset = exp.peek("model-set")
    assert isinstance(mset, Ncm.MSet)
    cosmo = mset.peek(Nc.HICosmo.id())
    assert isinstance(cosmo, Nc.HICosmoDEWSpline)
    assert cosmo.get_use_as_drag() == use_as_drag
    assert (mset.fparam_get_pi_by_name("NcHICosmo:asdrag") is not None) == use_as_drag


CURVATURE_PRIORS = [
    gen.CurvaturePriorType.NONE,
    gen.CurvaturePriorType.MEAN_KAPPA,
    gen.CurvaturePriorType.LP_KAPPA,
    gen.CurvaturePriorType.LP_D2,
]


@pytest.mark.parametrize("curvature_prior", CURVATURE_PRIORS)
def test_generate_dewspline_curvature_prior(
    tmp_path: Path, curvature_prior: "gen.CurvaturePriorType"
):
    """Test DEWSpline generation across curvature prior options."""
    exp_file = tmp_path / "dewspline.yaml"
    _ = gen.GenerateDEWSpline(
        experiment=exp_file.absolute(),
        include_bao=hicosmo.BAOID.ALL_COMBINED_JUN_2025,
        curvature_prior=curvature_prior,
        curvature_sigma=1.5,
        curvature_p=12.0,
    )
    assert exp_file.exists()
    assert exp_file.with_suffix(".functions.yaml").exists()


@pytest.mark.parametrize("curvature_prior", CURVATURE_PRIORS)
def test_generate_qspline_curvature_prior(
    tmp_path: Path, curvature_prior: "gen.CurvaturePriorType"
):
    """Test QSpline generation across curvature prior options."""
    exp_file = tmp_path / "qspline.yaml"
    _ = gen.GenerateQSpline(
        experiment=exp_file.absolute(),
        include_bao=hicosmo.BAOID.ALL_COMBINED_JUN_2025,
        curvature_prior=curvature_prior,
        curvature_sigma=1.5,
        curvature_p=12.0,
    )
    assert exp_file.exists()
    assert exp_file.with_suffix(".functions.yaml").exists()


@pytest.mark.parametrize(
    "generator, func_name",
    [
        (gen.GenerateDEWSpline, "NcHICosmoDEWSpline:lp_w2_ls"),
        (gen.GenerateQSpline, "NcHICosmoQSpline:lp_q2_ls"),
    ],
)
def test_generate_spline_default_curvature_prior(
    tmp_path: Path, generator: type, func_name: str
):
    """The default curvature prior is the L10 norm of the second derivative.

    It is a Gaussian prior on the residual functions, with one residual per
    quadrature node, and the same value as the prior on the norm.
    """
    exp_file = tmp_path / "spline.yaml"
    _ = generator(
        experiment=exp_file.absolute(),
        include_bao=hicosmo.BAOID.ALL_COMBINED_JUN_2025,
    )

    ser = Ncm.Serialize.new(Ncm.SerializeOpt.KEEP_SHARED)
    _ = ser.from_binfile(exp_file.with_suffix(".dataset.gvar").absolute().as_posix())
    exp = ser.dict_str_from_yaml_file(exp_file.absolute().as_posix())
    likelihood = exp.peek("likelihood")
    assert isinstance(likelihood, Ncm.Likelihood)
    priors = [likelihood.priors_peek_f(i) for i in range(likelihood.priors_length_f())]
    curvature = [
        prior
        for prior in priors
        if isinstance(prior, Ncm.PriorGaussFunc)
        and prior.props.mean_func.props.full_name == func_name
    ]
    assert len(curvature) == 1
    assert curvature[0].props.variable == 10.0

    mset = exp.peek("model-set")
    assert isinstance(mset, Ncm.MSet)
    cosmo = mset.peek(Nc.HICosmo.id())
    assert curvature[0].props.dimension == cosmo.d2_lp_residuals_len(None, 10)

    norm = Ncm.MSetFuncList.new(func_name.removesuffix("_ls"), None)
    f = np.array(curvature[0].eval_array(mset, None))
    expected = norm.eval1(mset, 10.0) / curvature[0].props.sigma
    assert np.sum(f**2) == pytest.approx(expected**2, rel=1.0e-10)


# LOCAL_* priors build a data-driven weight from the expected Fisher, so they
# need Gaussian-family data (Hubble) rather than the empirical-fit BAO alone.
LOCAL_CURVATURE_PRIORS = [
    gen.CurvaturePriorType.LOCAL_KAPPA,
    gen.CurvaturePriorType.LOCAL_D2,
]


@pytest.mark.parametrize("curvature_prior", LOCAL_CURVATURE_PRIORS)
@pytest.mark.parametrize("generator", [gen.GenerateDEWSpline, gen.GenerateQSpline])
def test_generate_spline_local_curvature_prior(
    tmp_path: Path, generator, curvature_prior: "gen.CurvaturePriorType"
):
    """Local (data-driven) curvature priors build the weight and the wlp_* func."""
    exp_file = tmp_path / "spline.yaml"
    _ = generator(
        experiment=exp_file.absolute(),
        include_hubble=hicosmo.HID.ALL_COMBINED_APR_2025,
        curvature_prior=curvature_prior,
        curvature_sigma=1.5,
        curvature_p=2.0,
        curvature_ref_factor=0.5,
    )
    assert exp_file.exists()
    assert exp_file.with_suffix(".functions.yaml").exists()


@pytest.mark.parametrize(
    "knots", [gen.KnotPlacement.UNIFORM, gen.KnotPlacement.CHEBYSHEV]
)
@pytest.mark.parametrize("generator", [gen.GenerateDEWSpline, gen.GenerateQSpline])
def test_generate_spline_knot_placement(
    tmp_path: Path, generator, knots: "gen.KnotPlacement"
):
    """Explicit knot placement constructs the model through the knots= path."""
    exp_file = tmp_path / "spline.yaml"
    _ = generator(
        experiment=exp_file.absolute(),
        include_bao=hicosmo.BAOID.ALL_COMBINED_JUN_2025,
        knots=knots,
    )
    assert exp_file.exists()


def test_add_curvature_prior_local_requires_weight():
    """A LOCAL_* prior without a precomputed weight spline is rejected."""
    likelihood = Ncm.Likelihood.new(Ncm.Dataset.new())
    with pytest.raises(ValueError, match="requires a precomputed weight spline"):
        gen._add_curvature_prior(
            likelihood,
            cosmo=Nc.HICosmoDEWSpline.new(6, 2.33),
            namespace="NcHICosmoDEWSpline",
            d2_name="lp_w2",
            prior_type=gen.CurvaturePriorType.LOCAL_D2,
            sigma=1.0,
            p=2.0,
            weight=None,
        )


@pytest.mark.parametrize(
    "generator, n_base",
    [(gen.GenerateQSpline, 2), (gen.GenerateDEWSpline, 2)],
)
def test_generate_band_nodes(tmp_path: Path, generator, n_base: int):
    """band_nodes adds that many reconstruction-band functions to the catalog."""
    exp_file = tmp_path / "band.yaml"
    band_nodes = 7
    _ = generator(
        experiment=exp_file.absolute(),
        include_bao=hicosmo.BAOID.ALL_COMBINED_JUN_2025,
        band_nodes=band_nodes,
    )
    ser = Ncm.Serialize.new(Ncm.SerializeOpt.KEEP_SHARED)
    funcs = ser.array_from_yaml_file(
        exp_file.with_suffix(".functions.yaml").absolute().as_posix()
    )
    # n_base nvar=0 diagnostics + band_nodes bound-z functions.
    assert funcs.len() == n_base + band_nodes


def test_generate_dewspline_no_data(tmp_path: Path):
    """Test DEWSpline rejects experiment with no data."""
    exp_file = tmp_path / "qspline.yaml"
    with pytest.raises(ValueError, match="No data included in the experiment."):
        _ = gen.GenerateDEWSpline(experiment=exp_file.absolute())


def test_generate_dewspline_invalid_suffix(tmp_path: Path):
    """Test DEWSpline rejects invalid file suffix."""
    exp_file = tmp_path / "dewspline.txt"
    bao_id = hicosmo.BAOID.ALL_COMBINED_JUN_2025
    with pytest.raises(ValueError, match="Invalid experiment file suffix: .txt"):
        _ = gen.GenerateDEWSpline(experiment=exp_file.absolute(), include_bao=bao_id)


# Native Planck: paths that need no local plc_3.0 tree (see fixtures_planck.py).
def test_generate_planck_from_release_requires_native(tmp_path: Path):
    """--from-release only makes sense together with --native."""
    from numcosmo_py.experiments.planck18 import Planck18Types

    exp_file = tmp_path / "planck.yaml"
    with pytest.raises(ValueError, match="--from-release requires --native"):
        _ = gen.GeneratePlanck(
            experiment=exp_file.absolute(),
            data_type=Planck18Types.TT,
            from_release=True,
        )


def test_build_planck_release_without_data(tmp_path: Path, monkeypatch):
    """Rebuilding the release with no clik tree present is an explicit error."""
    import numcosmo_py.experiments.planck_native_release as pnr

    monkeypatch.setattr(pnr, "find_baseline_file", lambda relpath: None)

    with pytest.raises(ValueError, match="No Planck clik data found"):
        _ = gen.BuildPlanckRelease(output_dir=(tmp_path / "release").absolute())


def test_build_planck_release_writes_available_ids(tmp_path: Path, monkeypatch):
    """The publish command writes the ids whose source data it can find."""
    import numcosmo_py.experiments.planck_native_release as pnr
    from numcosmo_py.experiments.planck_commander import COMMANDER_RELPATH
    from python.fixtures_planck import make_commander_cldf

    clik = make_commander_cldf(tmp_path)
    monkeypatch.setattr(
        pnr,
        "find_baseline_file",
        lambda relpath: clik if relpath == COMMANDER_RELPATH else None,
    )

    out_dir = tmp_path / "release"
    _ = gen.BuildPlanckRelease(output_dir=out_dir.absolute())

    assert [p.name for p in out_dir.glob("planck_native_*.gvar")] == [
        "planck_native_pr3_commander.gvar"
    ]


def test_generate_planck_native_from_release(tmp_path: Path, monkeypatch):
    """--native --from-release writes an experiment holding only native blocks."""
    import numcosmo_py.experiments.planck_native_release as pnr
    from numcosmo_py.experiments.planck18 import Planck18Types
    from numcosmo_py.experiments.planck_commander import build_commander
    from numcosmo_py.experiments.planck_simall import build_simall
    from numcosmo_py.experiments.planck_smica import build_smica_tt
    from python.fixtures_planck import (
        make_commander_cldf,
        make_simall_cldf,
        make_smica_cldf,
    )

    trees = {
        pnr.PlanckReleaseId.PR3_SIMALL_EE: (make_simall_cldf(tmp_path), build_simall),
        pnr.PlanckReleaseId.PR3_COMMANDER: (
            make_commander_cldf(tmp_path),
            build_commander,
        ),
        pnr.PlanckReleaseId.PR3_PLIK_TT: (make_smica_cldf(tmp_path), build_smica_tt),
    }

    def _load(rid, pb=None, cache_dir=None):
        del cache_dir
        clik, builder = trees[rid]
        return builder(clik, pb)

    monkeypatch.setattr(pnr, "load_planck_release", _load)

    exp_file = tmp_path / "planck_native_release.yaml"
    _ = gen.GeneratePlanck(
        experiment=exp_file.absolute(),
        native=True,
        from_release=True,
        data_type=Planck18Types.TT,
    )

    assert exp_file.exists()
    dset_file = exp_file.with_suffix(".dataset.gvar")
    assert dset_file.exists()

    ser = Ncm.Serialize.new(Ncm.SerializeOpt.KEEP_SHARED)
    dset = ser.from_binfile(str(dset_file.absolute()))
    names = {dset.get_data(i).__class__.__name__ for i in range(dset.get_ndata())}
    assert names == {"DataPlanckSimall", "DataPlanckCommander", "DataPlanckSmica"}


# Bianchi
def test_generate_bianchi(tmp_path: Path):
    """Test Bianchi I generation from the CF4 SNIa catalog with a redshift cut."""
    exp_file = tmp_path / "bianchi.yaml"
    cat_file = tmp_path / "cf4.obj"
    _ = gen.GenerateBianchi(
        experiment=exp_file.absolute(),
        catalog=dist_mu_obs.DistMuObsCatalogID.CF4_SNIA,
        z_max=0.02,
        catalog_out=cat_file.absolute(),
    )
    assert exp_file.exists()
    assert exp_file.with_suffix(".functions.yaml").exists()
    assert cat_file.exists()

    ser = Ncm.Serialize.new(Ncm.SerializeOpt.KEEP_SHARED)
    experiment = ser.dict_str_from_yaml_file(exp_file.absolute().as_posix())
    likelihood = cast(Ncm.Likelihood, experiment.get("likelihood"))
    data = cast(Nc.DataHCosmoBIDistMu, likelihood.peek_dataset().get_data(0))
    z = data.peek_z()
    assert 0 < z.len() < 1008
    assert max(z.get(i) for i in range(z.len())) <= 0.02

    # The orientation is bounded to one fundamental domain by default
    mset = cast(Ncm.MSet, experiment.get("model-set"))
    cosmo = cast(Nc.HCosmoBI, mset.peek(Nc.HCosmoBI.id()))
    assert cosmo.param_get_upper_bound(cosmo.param_index_from_name("u1")[1]) == 0.5
    assert cosmo.param_get_upper_bound(cosmo.param_index_from_name("u2")[1]) == 0.25
    assert cosmo.param_get_upper_bound(cosmo.param_index_from_name("u3")[1]) == 1.0

    # The written catalog loads as a NcDistMuObs and gives the same experiment
    exp_file2 = tmp_path / "bianchi2.yaml"
    _ = gen.GenerateBianchi(
        experiment=exp_file2.absolute(),
        data_file=cat_file.absolute(),
        mu_col="DM_av",
        z_max=0.02,
        fundamental_domain=False,
    )
    assert exp_file2.exists()
    ser2 = Ncm.Serialize.new(Ncm.SerializeOpt.KEEP_SHARED)
    experiment2 = ser2.dict_str_from_yaml_file(exp_file2.absolute().as_posix())
    mset2 = cast(Ncm.MSet, experiment2.get("model-set"))
    cosmo2 = cast(Nc.HCosmoBI, mset2.peek(Nc.HCosmoBI.id()))
    assert cosmo2.param_get_upper_bound(cosmo2.param_index_from_name("u1")[1]) == 1.0
    assert cosmo2.param_get_upper_bound(cosmo2.param_index_from_name("u2")[1]) == 1.0


def test_generate_bianchi_sigma_v_per_sample(
    tmp_path: Path, edd_cf4_excerpt_file: Path
):
    """Test Bianchi I generation with a velocity dispersion per sample."""
    exp_file = tmp_path / "bianchi_svs.yaml"
    _ = gen.GenerateBianchi(
        experiment=exp_file.absolute(),
        edd_table=edd_cf4_excerpt_file.absolute(),
        sigma_v=250.0,
        fit_sigma_v=True,
        sigma_v_per_sample=True,
        sample_rule=dist_mu_obs.SampleRule.SMALLEST_ERROR,
    )
    ser = Ncm.Serialize.new(Ncm.SerializeOpt.KEEP_SHARED)
    experiment = ser.dict_str_from_yaml_file(exp_file.absolute().as_posix())
    mset = cast(Ncm.MSet, experiment.get("model-set"))
    nuisance = cast(Nc.DistMuNuisance, mset.peek(Nc.DistMuNuisance.id()))
    n = nuisance.get_nsamples()
    assert n > 0
    assert nuisance.get_nsigma_v() == n
    for s in range(n):
        assert nuisance.get_sigma_v_sample(s) == 250.0
        assert (
            nuisance.param_get_ftype(nuisance.param_index_from_name(f"sigmav_s_{s}")[1])
            == Ncm.ParamType.FREE
        )
    # Sample -1 keeps the common value, which stays fixed since every object has a sample
    assert nuisance.get_sigma_v_sample(-1) == 250.0
    assert (
        nuisance.param_get_ftype(nuisance.param_index_from_name("sigmav")[1])
        == Ncm.ParamType.FIXED
    )

    with pytest.raises(ValueError, match="needs samples"):
        _ = gen.GenerateBianchi(
            experiment=(tmp_path / "bianchi_svs2.yaml").absolute(),
            edd_table=edd_cf4_excerpt_file.absolute(),
            sigma_v_per_sample=True,
        )


def test_generate_bianchi_velocity(tmp_path: Path):
    """Test Bianchi I generation with the observer velocity free."""
    exp_file = tmp_path / "bianchi_v.yaml"
    _ = gen.GenerateBianchi(
        experiment=exp_file.absolute(),
        catalog=dist_mu_obs.DistMuObsCatalogID.CF4_SNIA,
        z_min=0.01,
        z_max=0.03,
        fit_velocity=True,
    )
    ser = Ncm.Serialize.new(Ncm.SerializeOpt.KEEP_SHARED)
    experiment = ser.dict_str_from_yaml_file(exp_file.absolute().as_posix())
    mset = cast(Ncm.MSet, experiment.get("model-set"))
    velocity = cast(Nc.ObserverVelocity, mset.peek(Nc.ObserverVelocity.id()))
    assert velocity is not None
    assert all(
        velocity.param_get_ftype(velocity.param_index_from_name(n)[1])
        == Ncm.ParamType.FREE
        for n in ("vx", "vy", "vz")
    )
    assert velocity.speed() == 0.0
    # The distance table reaches beyond the largest redshift by the largest beta
    likelihood = cast(Ncm.Likelihood, experiment.get("likelihood"))
    data = cast(Nc.DataHCosmoBIDistMu, likelihood.peek_dataset().get_data(0))
    dist = data.props.dist
    assert dist is not None
    assert dist.props.zmax > 1.05 * 0.03 + Nc.OBSERVER_VELOCITY_MAX_V / 299792.458

    # Too small a lower redshift limit for the velocity bound
    with pytest.raises(ValueError, match="raise z_min"):
        _ = gen.GenerateBianchi(
            experiment=(tmp_path / "bianchi_v2.yaml").absolute(),
            catalog=dist_mu_obs.DistMuObsCatalogID.CF4_SNIA,
            z_max=0.03,
            fit_velocity=True,
        )


def test_generate_bianchi_sources(tmp_path: Path):
    """Test Bianchi I generation rejects zero or two catalog sources."""
    exp_file = tmp_path / "bianchi.yaml"
    with pytest.raises(typer.BadParameter, match="exactly one of"):
        _ = gen.GenerateBianchi(experiment=exp_file.absolute())
    with pytest.raises(typer.BadParameter, match="exactly one of"):
        _ = gen.GenerateBianchi(
            experiment=exp_file.absolute(),
            catalog=dist_mu_obs.DistMuObsCatalogID.CF4_SNIA,
            data_file=exp_file.absolute(),
        )


def test_generate_bianchi_invalid_suffix(tmp_path: Path):
    """Test Bianchi I generation rejects invalid file suffix."""
    exp_file = tmp_path / "bianchi.txt"
    with pytest.raises(ValueError, match="Invalid experiment file suffix: .txt"):
        _ = gen.GenerateBianchi(
            experiment=exp_file.absolute(),
            catalog=dist_mu_obs.DistMuObsCatalogID.CF4_SNIA,
        )


def test_generate_bianchi_edd_table(tmp_path: Path, edd_cf4_excerpt_file: Path):
    """Test Bianchi I generation from an EDD export, with a catalog written out."""
    exp_file = tmp_path / "bianchi_edd.yaml"
    cat_file = tmp_path / "edd.gvar"
    _ = gen.GenerateBianchi(
        experiment=exp_file.absolute(),
        edd_table=edd_cf4_excerpt_file,
        edd_downloaded="2026-10-17",
        omegam_mean=None,
        ordering=False,
        catalog_out=cat_file.absolute(),
    )
    assert exp_file.exists()
    ser = Ncm.Serialize.new(Ncm.SerializeOpt.KEEP_SHARED)
    obs = ser.from_file(cat_file.absolute().as_posix())
    assert isinstance(obs, Nc.DistMuObs)
    assert obs.len() == 6


def test_generate_bianchi_nuisance(tmp_path: Path, edd_cf4_excerpt_file: Path):
    """Test Bianchi I generation with the velocity dispersion and per-sample offsets."""
    exp_file = tmp_path / "bianchi_nuisance.yaml"
    _ = gen.GenerateBianchi(
        experiment=exp_file.absolute(),
        edd_table=edd_cf4_excerpt_file,
        sigma_v=300.0,
        fit_sigma_v=True,
        sample_rule=dist_mu_obs.SampleRule.SMALLEST_ERROR,
        fit_offsets=True,
        reference_sample="0",
        omegam_mean=None,
    )
    ser = Ncm.Serialize.new(Ncm.SerializeOpt.KEEP_SHARED)
    experiment = ser.dict_str_from_yaml_file(exp_file.absolute().as_posix())
    mset = cast(Ncm.MSet, experiment.get("model-set"))
    nuisance = cast(Nc.DistMuNuisance, mset.peek(Nc.DistMuNuisance.id()))
    assert nuisance.get_sigma_v() == 300.0
    assert nuisance.get_nsamples() >= 1
    assert mset.fparams_len() == 7 + 1 + nuisance.get_nsamples() - 1
    with pytest.raises(ValueError, match="not among the samples"):
        gen.GenerateBianchi(
            experiment=(tmp_path / "b3.yaml").absolute(),
            edd_table=edd_cf4_excerpt_file,
            sample_rule=dist_mu_obs.SampleRule.SMALLEST_ERROR,
            fit_offsets=True,
            reference_sample="nosuchsurvey",
        )
    with pytest.raises(ValueError, match="needs a reference_sample"):
        gen.GenerateBianchi(
            experiment=(tmp_path / "b2.yaml").absolute(),
            edd_table=edd_cf4_excerpt_file,
            sample_rule=dist_mu_obs.SampleRule.HEMISPHERE,
            fit_offsets=True,
        )
    likelihood = cast(Ncm.Likelihood, experiment.get("likelihood"))
    data = cast(Nc.DataHCosmoBIDistMu, likelihood.peek_dataset().get_data(0))
    samples = data.peek_samples()
    assert (
        max(samples.get(i) for i in range(samples.len())) == nuisance.get_nsamples() - 1
    )
