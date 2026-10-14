#!/usr/bin/env python
#
# test_py_xcor_kernel_serialization.py
#
# Thu Apr 17 2026
# Copyright  2026  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# test_py_xcor_kernel_serialization.py
# Copyright (C) 2026 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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

"""XcorKernel serialization through the paths only Python reaches.

That a serialized copy is the same kernel (type, knobs, range, Limber window, both
serializer options, re-preparation) is checked in C, test_nc_xcor_kernel_common.c
(serialize). What stays here: every kernel type resolves by name in a fresh process, and
every kernel of the LSST-Y1 fixtures rebuilds from its serialized string.
"""

import pytest

from numcosmo_py import Ncm, Nc
from numcosmo_py.cosmology import Cosmology

pytestmark = pytest.mark.xcor

Ncm.cfg_init()
pytest_plugins = ["python.fixtures_xcor"]


_REGISTRY_CHILD = """
import sys
from numcosmo_py import Ncm
from gi.repository import GObject

Ncm.cfg_init()

missing = []
for name in sys.argv[1:]:
    try:
        GObject.type_from_name(name)
    except RuntimeError:
        missing.append(name)
print(",".join(missing))
"""


def test_every_kernel_type_is_registered_in_a_fresh_process() -> None:
    """Concrete kernel types must resolve by name after cfg_init() alone.

    ncm_serialize_from_name_params() resolves a serialized name with
    g_type_from_name(), which only succeeds once the GType has been realized in
    the process. ncm_cfg_register_objects() realizes it at startup; otherwise
    merely touching the class from Python realizes it too, which is why
    dup_obj-based tests -- and even a to_string/from_string round-trip on a live
    instance -- pass for a kernel that was never registered. Reading a YAML
    experiment in a fresh process does neither, so it aborts with "object `X' is
    not registered". This ran green while NcXcorKernelCMBISW was missing from
    ncm_cfg_register_objects(), so it must not construct anything.
    """
    import subprocess
    import sys

    # Names are collected here, where realizing the classes is harmless; the
    # child only ever looks them up.
    names = sorted(
        "Nc" + n
        for n in dir(Nc)
        if n.startswith("XcorKernel")
        and not n.endswith(("SParams", "VParams", "Class", "Impl", "Integrand"))
        and n not in ("XcorKernel", "XcorKernelComponent")
    )
    assert "NcXcorKernelCMBISW" in names, "ISW kernel dropped from the sweep"

    proc = subprocess.run(
        [sys.executable, "-c", _REGISTRY_CHILD, *names],
        capture_output=True,
        text=True,
        check=False,
        timeout=600,
    )

    assert proc.returncode == 0, proc.stderr
    missing = [n for n in proc.stdout.strip().split(",") if n]
    assert not missing, (
        f"not passed to ncm_cfg_register_obj(): {missing} -- these cannot be "
        f"built from a serialized name in a fresh process"
    )


def test_kernel_reconstructible_by_type_name(
    kernel: Nc.XcorKernel, cosmology: Cosmology
) -> None:
    """Every kernel type must be reconstructible from its serialized name.

    Ncm.Serialize.dup_obj() resolves the GType from the live instance, so it
    works whether or not the type was passed to ncm_cfg_register_obj(). Going
    through a *string* does not: ncm_serialize_from_name_params() looks the name
    up in that registry and aborts if it is absent. That is the path YAML
    experiment files, the CLI and ncm_obj_array take, so a kernel missing from
    ncm_cfg_register_objects() is unusable there while every dup_obj-based test
    still passes -- which is exactly how NcXcorKernelCMBISW stayed broken.
    """
    ser = Ncm.Serialize.new(Ncm.SerializeOpt.KEEP_SHARED)

    # The GType name keeps the Nc namespace that introspection strips, and it
    # is the name the registry is keyed on.
    gtype_name = kernel.__gtype__.name
    serialized = ser.to_string(kernel, False)
    assert serialized.startswith(gtype_name)

    rebuilt = ser.from_string(serialized)

    assert rebuilt is not kernel
    assert isinstance(rebuilt, type(kernel))
    rebuilt.prepare(cosmology.cosmo)
