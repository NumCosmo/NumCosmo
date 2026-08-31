#!/usr/bin/env python
"""Freeze pre-NcBBN *catalogs* as fixtures.

The BBN migration took NcHICosmoDE from eight scalar parameters to seven, so
every catalog written before it carries an NcHICosmoDEReparamOk whose length is
8 -- and reading one aborted inside ncm_vector_memcpy() until models learned to
rebuild a reparametrization at their own length.

Unlike tests/tools/make_bbn_compat_fixtures.py, this script does not have to run
against the old library: a catalog's model-set is serialized text, so the old
shape is written out here literally. The three containers below are the three
this repository has ever produced, transcribed from real files under ~/data:

  old_reparam_vardict.mc.fits   HDU0 holds a {'model-set', 'functions'} vardict
                                (MSETFMT = 'vardict'), the current container
  old_reparam_object.mc.fits    HDU0 holds the bare model-set object
                                (MSETFMT = 'gvariant'), the container before it
  old_reparam_sidecar.mc.fits   empty HDU0, model-set in a `.mset' GKeyFile
    + old_reparam_sidecar.mc.mset

The table extension is written by the current library -- only HDU0 and the
sidecar are old -- so the fixtures stay small and the columns stay consistent
with the model-set's free parameters.
"""

import os
import tempfile

import numpy as np
from astropy.io import fits
from gi.repository import GLib

from numcosmo_py import Nc, Ncm

Ncm.cfg_init()

OUT = os.path.join(
    os.path.dirname(os.path.abspath(__file__)),
    "..",
    "..",
    "data",
    "truth_tables",
    "bbn",
)

NROWS = 40
NCHAINS = 4

# What the old library wrote and the current one does not. Yp was a parameter of
# the cosmology; the reparametrization counted it.
OLD_EDITS = (
    ("'length': <uint32 7>", "'length': <uint32 8>"),
    ("'ENnu': <", "'Yp': <0.23999999999999999>, 'ENnu': <"),
    ("'ENnu-fit': <", "'Yp-fit': <false>, 'ENnu-fit': <"),
    # Old cosmologies had no NcBBN submodel at all.
    (
        "'submodel-array': <[('NcBBNParthenope', {'sparam-array': <@a{i(sa{sv})} {}>, "
        "'submodel-array': <@a(sa{sv}) []>, 'table': <0>})]>",
        "'submodel-array': <@a(sa{sv}) []>",
    ),
)

OLD_SIDECAR_EDITS = (
    ("'length': <uint32 7>", "'length': <uint32 8>"),
    ("\nENnu=", "\nYp=0.23999999999999999\nENnu="),
    ("\nENnu-fit=", "\nYp-fit=false\nENnu-fit="),
)


def make_old(text, edits):
    """Apply @edits to @text, each exactly once."""
    for old, new in edits:
        if text.count(old) != 1:
            raise RuntimeError(f"anchor {old!r} appears {text.count(old)} times")

        text = text.replace(old, new)

    return text


def build_mset():
    """The current-library equivalent of the cosmology the fixtures carry."""
    cosmo = Nc.HICosmoDEXcdm.new()
    cosmo.set_reparam(Nc.HICosmoDEReparamOk.new(cosmo.len()))
    cosmo.param_set_by_name("H0", 67.81)
    cosmo.param_set_by_name("Omegac", 0.2312)
    cosmo.param_set_by_name("Omegab", 0.0486)

    mset = Ncm.MSet.new_array([cosmo])
    mset.param_set_ftype(
        Nc.HICosmo.id(), Nc.HICosmoDESParams.OMEGA_X, Ncm.ParamType.FREE
    )
    mset.param_set_ftype(Nc.HICosmo.id(), Nc.HICosmoDEXCDMSParams.W, Ncm.ParamType.FREE)
    mset.prepare_fparam_map()

    return mset


def build_catalog(mset, path):
    """Write a small two-parameter chain, the way a real run would."""
    mcat = Ncm.MSetCatalog(
        mset=mset,
        nadd_vals=1,
        nadd_val_names=["m2lnL"],
        nadd_val_symbols=["-2\\ln(L)"],
        nchains=NCHAINS,
        m2lnp_var=0,
        weighted=False,
    )
    mcat.set_run_type("Fixture")
    mcat.set_file(path)

    rng = np.random.default_rng(1234)

    for _ in range(NROWS):
        omegak, w = rng.normal(scale=[0.02, 0.05]) + [0.0, -1.0]
        mset.fparams_set_vector(Ncm.Vector.new_array([omegak, w]))
        mcat.add_from_mset_array(mset, [omegak**2 + (w + 1.0) ** 2])

    mcat.sync(True)


def old_mset_text(mset):
    """The mset as a pre-migration library would have serialized it."""
    ser = Ncm.Serialize.new(0)

    return make_old(ser.to_string(mset, True), OLD_EDITS)


def old_sidecar_text(mset):
    """The mset as a pre-migration `.mset' GKeyFile."""
    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, "sidecar.mset")
        ser = Ncm.Serialize.new(0)
        mset.save(ser, path, True)

        # Old cosmologies had no NcBBN submodel, so no group for it either.
        keyfile = GLib.KeyFile()
        keyfile.load_from_file(path, GLib.KeyFileFlags.KEEP_COMMENTS)
        keyfile.remove_comment("NcBBN", None)
        keyfile.remove_group("NcBBN")
        text = keyfile.to_data()[0]

    return make_old(text, OLD_SIDECAR_EDITS)


def blob(variant):
    """A GVariant's serialized bytes, as the FITS image HDU0 stores them."""
    return np.frombuffer(variant.get_data_as_bytes().get_data(), dtype=np.uint8)


def write_fits(template, path, data, fmt):
    """Copy @template, replacing HDU0 with @data (or emptying it if %None)."""
    with fits.open(template) as hdul:
        primary = fits.PrimaryHDU(data=data)

        if fmt is not None:
            primary.header["MSETFMT"] = (fmt, "Format of the data stored in this HDU.")

        fits.HDUList([primary] + [h.copy() for h in hdul[1:]]).writeto(
            path, overwrite=True
        )


def main():
    """Write the three fixtures."""
    mset = build_mset()

    # Before build_catalog(), which walks the free parameters over a chain.
    mset_variant = GLib.Variant.parse(None, old_mset_text(mset))
    sidecar = old_sidecar_text(mset)
    vardict = GLib.Variant("a{sv}", {"model-set": mset_variant})

    with tempfile.TemporaryDirectory() as tmp:
        template = os.path.join(tmp, "template.mc.fits")
        build_catalog(mset, template)

        write_fits(
            template,
            os.path.join(OUT, "old_reparam_vardict.mc.fits"),
            blob(vardict),
            "vardict",
        )
        write_fits(
            template,
            os.path.join(OUT, "old_reparam_object.mc.fits"),
            blob(mset_variant),
            "gvariant",
        )
        write_fits(
            template, os.path.join(OUT, "old_reparam_sidecar.mc.fits"), None, None
        )

    with open(
        os.path.join(OUT, "old_reparam_sidecar.mc.mset"), "w", encoding="utf-8"
    ) as f:
        f.write(sidecar)

    for name in sorted(os.listdir(OUT)):
        if name.startswith("old_reparam"):
            print(f"{name:34s} {os.path.getsize(os.path.join(OUT, name)):6d} bytes")


if __name__ == "__main__":
    main()
