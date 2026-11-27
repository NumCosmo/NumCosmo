"""Tests of the tool that updates the p(z) objects of a serialized WL catalog."""

import math
from pathlib import Path

import pytest
from gi.repository import GLib, GObject

from numcosmo_py import Nc, Ncm
from tools import nc_galaxy_pz_update_catalog as tool

Ncm.cfg_init()

# Rows of the test catalog, and those that have a p(z).
N_ROWS = 3
ROWS_WITH_PZ = [0, 2]

# Points of each p(z) table.
N_POINTS = 40


def table(row: int) -> tuple[list[float], list[float]]:
    """Return the table of a row: a Gaussian whose mean depends on the row."""
    z = [0.3 + 0.01 * j for j in range(N_POINTS)]
    p = [math.exp(-0.5 * ((z_j - 0.45 - 0.05 * row) / 0.05) ** 2) for z_j in z]
    return z, p


def spline(row: int) -> Ncm.Spline:
    """Return the table of a row as the NcmSpline of a catalog written before."""
    z, p = table(row)
    return Ncm.SplineCubicNotaknot.new_full(
        Ncm.Vector.new_array(z), Ncm.Vector.new_array(p), True
    )


def catalog_variant(entries: dict[int, GObject.Object]) -> GLib.Variant:
    """Return a serialized catalog whose ``pz`` property holds ``entries``.

    The catalog is serialized without p(z) and the dictionary is written into
    the GVariant, so the result does not depend on the type that the ``pz``
    column of the library accepts.
    """
    obs = Nc.GalaxyWLObs.new(
        Nc.GalaxyWLObsEllipConv.TRACE_DET,
        Nc.WLEllipticityFrame.CELESTIAL,
        N_ROWS,
        ["ra", "dec", "z"],
        [Ncm.CatalogColType.DOUBLE] * len(["ra", "dec", "z"]),
    )
    for i in range(N_ROWS):
        obs.set("ra", i, 0.1 * i)
        obs.set("dec", i, -0.2 * i)
        obs.set("z", i, 0.5 + 0.1 * i)

    ser = Ncm.Serialize.new(Ncm.SerializeOpt.NONE)
    builder = GLib.VariantBuilder.new(GLib.VariantType.new(tool.PZ_DICT_TYPE))
    for row, obj in entries.items():
        builder.add_value(
            GLib.Variant.new_dict_entry(
                GLib.Variant.new_int32(row), ser.to_variant(obj)
            )
        )

    return tool.replace_pz(ser.to_variant(obs), builder.end())


def write(path: Path, variant: GLib.Variant) -> Path:
    """Write a serialized object as NcmSerialize writes a binary file."""
    path.write_bytes(variant.get_data_as_bytes().get_data())
    return path


def without(catalog: GLib.Variant, removed: str) -> GLib.Variant:
    """Return a serialized object without one of its properties."""
    builder = GLib.VariantBuilder.new(GLib.VariantType.new("a{sv}"))
    for name, value in tool.properties(catalog):
        if name != removed:
            builder.add_value(
                GLib.Variant.new_dict_entry(
                    GLib.Variant.new_string(name), GLib.Variant.new_variant(value)
                )
            )
    return GLib.Variant.new_tuple(catalog.get_child_value(0), builder.end())


@pytest.fixture(name="old_catalog")
def fixture_old_catalog(tmp_path: Path) -> Path:
    """Write a catalog whose p(z) are NcmSpline objects."""
    entries = {row: spline(row) for row in ROWS_WITH_PZ}
    return write(tmp_path / "old.gvar", catalog_variant(entries))


def test_update_replaces_the_splines(
    old_catalog: Path, tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    """Each NcmSpline becomes the NcGalaxyPzBSplineSqrt of its table."""
    output = tmp_path / "new.gvar"

    assert tool.main([str(old_catalog), str(output)]) == 0
    assert "2 entries replaced, 0 copied" in capsys.readouterr().out

    old = tool.load_catalog(old_catalog)
    new = tool.load_catalog(output)
    old_props = dict(tool.properties(old))
    new_props = dict(tool.properties(new))

    assert list(old_props) == list(new_props)
    for name, value in old_props.items():
        if name != "pz":
            assert value.equal(new_props[name])

    new_pz = new_props["pz"]
    assert tool.pz_rows(new_pz) == ROWS_WITH_PZ

    ser = Ncm.Serialize.new(Ncm.SerializeOpt.NONE)
    for i, row in enumerate(ROWS_WITH_PZ):
        obj_variant = new_pz.get_child_value(i).get_child_value(1)
        gpz = ser.from_variant(obj_variant)
        z, p = table(row)
        ref = Nc.GalaxyPzBSplineSqrt.new_from_table(
            Ncm.Vector.new_array(z), Ncm.Vector.new_array(p)
        )

        assert tool.type_name(obj_variant) == "NcGalaxyPzBSplineSqrt"
        assert isinstance(gpz, Nc.GalaxyPzBSplineSqrt)
        assert gpz.get_range() == (z[0], z[-1])
        assert gpz.get_norm() == ref.get_norm()
        for z_j in z:
            assert gpz.eval_p(z_j + 0.004) == ref.eval_p(z_j + 0.004)


def test_update_of_an_updated_catalog_copies_it(
    old_catalog: Path, tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    """A second run finds only NcGalaxyPz entries and writes the same bytes."""
    first = tmp_path / "first.gvar"
    second = tmp_path / "second.gvar"

    assert tool.main([str(old_catalog), str(first)]) == 0
    capsys.readouterr()
    assert tool.main([str(first), str(second)]) == 0
    assert "0 entries replaced, 2 copied" in capsys.readouterr().out
    assert first.read_bytes() == second.read_bytes()


def test_existing_output_needs_force(old_catalog: Path, tmp_path: Path) -> None:
    """An existing output file is overwritten only with --force."""
    output = tmp_path / "new.gvar"
    output.write_bytes(b"")

    with pytest.raises(SystemExit):
        tool.main([str(old_catalog), str(output)])
    assert output.read_bytes() == b""

    assert tool.main([str(old_catalog), str(output), "--force"]) == 0
    assert tool.is_update_of(tool.load_catalog(output), tool.load_catalog(old_catalog))


def test_output_must_not_be_the_input(old_catalog: Path) -> None:
    """The input file is never overwritten."""
    before = old_catalog.read_bytes()

    with pytest.raises(SystemExit):
        tool.main([str(old_catalog), str(old_catalog), "--force"])
    assert old_catalog.read_bytes() == before


def test_file_of_another_object(tmp_path: Path) -> None:
    """A file that does not hold a NcGalaxyWLObs is rejected."""
    ser = Ncm.Serialize.new(Ncm.SerializeOpt.NONE)
    path = write(tmp_path / "vector.gvar", ser.to_variant(Ncm.Vector.new_array([1.0])))

    with pytest.raises(ValueError, match="holds a NcmVector, not a NcGalaxyWLObs"):
        tool.main([str(path), str(tmp_path / "new.gvar")])


def test_catalog_without_pz(old_catalog: Path, tmp_path: Path) -> None:
    """A catalog without the pz property is rejected."""
    path = write(tmp_path / "no_pz.gvar", without(tool.load_catalog(old_catalog), "pz"))

    with pytest.raises(ValueError, match="has no pz property"):
        tool.main([str(path), str(tmp_path / "new.gvar")])


def test_entry_of_another_type(tmp_path: Path) -> None:
    """An entry that is neither a NcmSpline nor a NcGalaxyPz is rejected."""
    entries: dict[int, GObject.Object] = {1: Ncm.Vector.new_array([1.0])}
    path = write(tmp_path / "bad.gvar", catalog_variant(entries))

    with pytest.raises(ValueError, match="pz entry 1 is a NcmVector"):
        tool.main([str(path), str(tmp_path / "new.gvar")])


def test_pz_of_another_variant_type() -> None:
    """A pz value that is not a dictionary of objects by row is rejected."""
    with pytest.raises(ValueError, match="pz has type as"):
        tool.update_pz(GLib.Variant("as", ["a"]))


def test_is_update_of(old_catalog: Path, tmp_path: Path) -> None:
    """Only the input with every pz entry updated counts as its update."""
    old = tool.load_catalog(old_catalog)
    output = tmp_path / "new.gvar"
    assert tool.main([str(old_catalog), str(output)]) == 0
    new = tool.load_catalog(output)
    assert tool.is_update_of(new, old)
    # The entries were not updated.
    assert not tool.is_update_of(old, old)
    # Another object name.
    renamed = GLib.Variant.new_tuple(
        GLib.Variant.new_string("NcGalaxyWLObs[S7]"), new.get_child_value(1)
    )
    assert not tool.is_update_of(renamed, old)
    # A property is missing.
    assert not tool.is_update_of(without(new, "len"), old)
    # A property other than pz has another value.
    changed = GLib.VariantBuilder.new(GLib.VariantType.new("a{sv}"))
    for name, value in tool.properties(new):
        changed.add_value(
            GLib.Variant.new_dict_entry(
                GLib.Variant.new_string(name),
                GLib.Variant.new_variant(
                    GLib.Variant.new_uint32(N_ROWS + 1) if name == "len" else value
                ),
            )
        )
    assert not tool.is_update_of(
        GLib.Variant.new_tuple(new.get_child_value(0), changed.end()), old
    )
    # The rows of pz differ.
    one_row = write(tmp_path / "one.gvar", catalog_variant({0: spline(0)}))
    fewer = tmp_path / "fewer.gvar"
    assert tool.main([str(one_row), str(fewer)]) == 0
    assert not tool.is_update_of(tool.load_catalog(fewer), old)


def test_main_checks_the_written_file(
    old_catalog: Path, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """The run fails when the file written is not the input with pz updated."""
    monkeypatch.setattr(tool, "is_update_of", lambda new, old: False)

    with pytest.raises(RuntimeError, match="with pz updated"):
        tool.main([str(old_catalog), str(tmp_path / "new.gvar")])
