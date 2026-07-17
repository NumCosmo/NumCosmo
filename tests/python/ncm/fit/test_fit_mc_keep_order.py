#!/usr/bin/env python
#
# test_fit_mc_keep_order.py
#
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
# with this program. If not, see <http://www.gnu.org/licenses/>.

"""Regression tests for ``NcmFitMC``'s ``keep-order`` property.

A full implementation of ``keep_order`` (field, GObject property, and a
dedicated ordered-evaluation code path) existed years ago but was silently
dropped by the OpenMP-parallelism rewrite in commit 83052fb5 -- only the
public header declaration ``ncm_fit_mc_keep_order()`` survived, with no
implementation anywhere in the codebase (a call to it would fail at
runtime). Found while trying to build a row-paired Monte Carlo comparison
between two ``NcmData`` engines: without ``keep_order``, ``sample_index``
(hence catalog row placement) is assigned in whatever order threads happen
to *race* into the resample critical section -- not in loop-iteration
order -- so which physical realization lands in which row is scheduling-
dependent and not reproducible run to run, even for the same seed.

These tests use the same trivial ``Ncm.ModelMVND``/toy-``Ncm.Data`` style
already established in ``test_fit_mc.py``, plus one small custom
``Ncm.Data`` subclass (``_OrderProbeData``) that encodes, directly in the
fitted parameter value, which thread produced a realization and in what
order that thread called ``resample()`` -- letting the tests decode row
order from the catalog alone. Its very first ``resample()`` call sleeps for
much longer than the (near-instant) fit itself, so a second thread reliably
races ahead through its *entire* remaining chunk before the first thread
wakes -- a deliberately engineered, large-margin divergence, not a
timing-window race that could occasionally pass by luck.
"""

import threading
import time
import numpy as np
import pytest

from numcosmo_py import Ncm

Ncm.cfg_init()

# FitMC parallelism is OpenMP; these tests set nthreads>1, so they exercise
# the OpenMP-parallel path in the dedicated omp lane. See TESTING.md.
pytestmark = pytest.mark.omp

MEAN = 1.2543
SLOW_SLEEP_S = 1.5  # far larger than the trivial fit's own (sub-ms) cost
N_PER_CHUNK = 4
NTHREADS = 2
N_TOTAL = N_PER_CHUNK * NTHREADS


def _build_mvnd_mc(seed):
    """The same trivial MVND fixture style as test_fit_mc.py."""
    rng = Ncm.RNG.seeded_new(None, 1234)
    model_mvnd = Ncm.ModelMVND.new(2)
    data_mvnd = [
        Ncm.DataGaussCovMVND.new_full(2, 1.0e-1, 2.0e-1, 30.0, MEAN, MEAN, rng)
        for _ in range(4)
    ]
    mset = Ncm.MSet.new_array([model_mvnd])
    mset.param_set_all_ftype(Ncm.ParamType.FREE)
    mset.prepare_fparam_map()
    mset.fparams_set_array(np.ones(2) * MEAN)
    likelihood = Ncm.Likelihood.new(Ncm.Dataset.new_array(data_mvnd))
    fit = Ncm.Fit.factory(
        Ncm.FitType.GSL_MMS, None, likelihood, mset, Ncm.FitGradType.NUMDIFF_CENTRAL
    )
    mc = Ncm.FitMC.new(fit, Ncm.FitMCResampleType.FROM_MODEL, Ncm.FitRunMsgs.NONE)
    mc.set_rng(Ncm.RNG.seeded_new("mt19937", seed))
    return mc


def test_keep_order_property_round_trips():
    """The property (not just the setter) is restored and defaults to False."""
    mc = _build_mvnd_mc(1)
    assert mc.props.keep_order is False
    mc.keep_order(True)
    assert mc.props.keep_order is True
    mc.props.keep_order = False
    assert mc.props.keep_order is False


def test_keep_order_matches_serial_baseline_regardless_of_nthreads():
    """With keep_order=True, the catalog is bit-identical between a serial
    (nthreads=1, trivially ordered) run and a multi-threaded run at the same
    seed -- this direction is deterministic by construction (the "omp
    ordered" region forces resample() to execute in loop-iteration order
    regardless of which thread executes which iteration), so it should
    never be flaky."""
    n = 40
    seed = 777

    mc1 = _build_mvnd_mc(seed)
    mc1.keep_order(True)
    mc1.set_nthreads(1)
    mc1.start_run()
    mc1.run(n)
    mc1.end_run()
    cat1 = mc1.peek_catalog()
    rows1 = np.array([cat1.peek_row(i).dup_array() for i in range(cat1.len())])

    mc4 = _build_mvnd_mc(seed)
    mc4.keep_order(True)
    mc4.set_nthreads(4)
    mc4.start_run()
    mc4.run(n)
    mc4.end_run()
    cat4 = mc4.peek_catalog()
    rows4 = np.array([cat4.peek_row(i).dup_array() for i in range(cat4.len())])

    np.testing.assert_array_equal(rows1, rows4)


class _OrderProbeData(Ncm.Data):
    """Toy Ncm.Data whose fitted parameter, after resample(), encodes
    (thread_slot, call_count_within_thread) -- see module docstring."""

    __gtype_name__ = "TestFitMCOrderProbeData"

    # Class-level (not instance-level): each OpenMP thread gets its OWN
    # duplicated NcmData instance (created once, reused for that thread's
    # whole chunk of iterations), so the registry mapping thread identity to
    # slot numbers must be shared across those instances, hence a class
    # attribute. Safe without extra locking: do_resample() is only ever
    # entered by one thread at a time (the whole point of this test --
    # either the "critical" or the "ordered" region around ncm_fit_mc's
    # resample step serializes it).
    #
    # Keyed by threading.get_ident(), NOT id(self): PyGObject's per-thread
    # wrapper resolution for a Python-subclassed GObject called back into
    # from multiple native (OpenMP-spawned) threads was empirically found
    # to sometimes hand back the SAME Python wrapper identity (id(self)) to
    # what are genuinely two different underlying C NcmData instances (one
    # per thread) -- verified independently via threading.get_ident(),
    # which reliably reports two distinct OS thread ids. Whatever the exact
    # PyGObject mechanism, id(self) is not a trustworthy per-thread key
    # here; threading.get_ident() is.
    _thread_slots: dict = {}
    _next_slot = [0]

    @staticmethod
    def reset_registry():
        """Must be called before each test using this class: OS thread ids
        can be recycled once a previous test's threads have exited, so a
        registry surviving across tests would risk a stale key colliding
        with an unrelated thread in a later test."""
        _OrderProbeData._thread_slots = {}
        _OrderProbeData._next_slot = [0]

    def __init__(self):
        super().__init__()
        self.set_init(True)
        self.call_count = 0
        self.target = 0.0

    def do_get_length(self):  # pylint: disable-msg=arguments-differ
        return 1

    def do_begin(self):  # pylint: disable-msg=arguments-differ
        pass

    def do_prepare(self, mset):  # pylint: disable-msg=arguments-differ
        pass

    def do_resample(self, mset, rng):  # pylint: disable-msg=arguments-differ
        key = threading.get_ident()

        if key not in _OrderProbeData._thread_slots:
            _OrderProbeData._thread_slots[key] = _OrderProbeData._next_slot[0]
            _OrderProbeData._next_slot[0] += 1

        thread_slot = _OrderProbeData._thread_slots[key]

        self.target = float(thread_slot * 100 + self.call_count)
        self.call_count += 1

        # Punish whichever thread happens to win the very first resample
        # race overall: sleeping here (still holding the section that
        # serializes resample calls) delays only *this* thread's return to
        # compete for its own next iteration -- any other thread already
        # mid-chunk races ahead through its remaining iterations in the
        # meantime, uncontested.
        if thread_slot == 0 and self.call_count == 1:
            time.sleep(SLOW_SLEEP_S)

    def do_m2lnL_val(self, mset):  # pylint: disable-msg=arguments-differ
        model = mset.peek(Ncm.ModelMVND.id())
        param = model.orig_vparam_get(0, 0)

        return (param - self.target) ** 2


@pytest.fixture(autouse=True)
def _reset_order_probe_registry():
    _OrderProbeData.reset_registry()
    yield


def _build_probe_mc(seed):
    # ncm_fit_mc_set_nthreads()/mc->nthreads does NOT itself call
    # omp_set_num_threads() -- the real OpenMP thread count is controlled
    # only by the OpenMP runtime's own state (env var OMP_NUM_THREADS or an
    # explicit omp_set_num_threads() call), independent of this property
    # (see ncm-fit-mc-nthreads-noop). This test needs the ACTUAL thread
    # count pinned to exactly NTHREADS regardless of the ambient
    # OMP_NUM_THREADS the test suite happens to run under, so it calls the
    # real OpenMP setter directly.
    Ncm.cfg_set_openmp_nthreads(NTHREADS)

    model_mvnd = Ncm.ModelMVND.new(1)
    mset = Ncm.MSet.new_array([model_mvnd])
    mset.param_set_all_ftype(Ncm.ParamType.FREE)
    mset.prepare_fparam_map()
    mset.fparams_set_array(np.zeros(1))

    dset = Ncm.Dataset.new()
    dset.append_data(_OrderProbeData())
    likelihood = Ncm.Likelihood.new(dset)
    fit = Ncm.Fit.factory(
        Ncm.FitType.GSL_MMS, None, likelihood, mset, Ncm.FitGradType.NUMDIFF_CENTRAL
    )
    mc = Ncm.FitMC.new(fit, Ncm.FitMCResampleType.FROM_MODEL, Ncm.FitRunMsgs.NONE)
    mc.set_rng(Ncm.RNG.seeded_new("mt19937", seed))
    mc.set_nthreads(NTHREADS)

    return mc


def _decode_rows(mc):
    """GSL_MMS (derivative-free Nelder-Mead) does not converge to the exact
    target value even for a perfect quadratic within a few evaluations --
    round (not truncate) to the nearest integer before decoding, since a
    slightly-under-converged value (e.g. 0.997) would otherwise floor to
    the wrong bucket."""
    cat = mc.peek_catalog()
    values = [round(cat.peek_row(i).get(1)) for i in range(cat.len())]

    return [(v // 100, v % 100) for v in values]


def _skip_unless_truly_parallel(slots):
    """Ncm.cfg_set_openmp_nthreads() cannot exceed OMP_THREAD_LIMIT (a hard
    ceiling per the OpenMP spec) -- under the fast test lane's thread-pinned
    environment (OMP_THREAD_LIMIT=1), this test's request for NTHREADS real
    threads is silently capped to 1, collapsing every row to the same
    thread_slot. That is an environment constraint, not a finding about
    ordering, so skip rather than fail when it happens (the dedicated OMP
    lane, which runs with a high thread limit, is where this test is
    actually meaningful)."""
    if len(set(slots)) < 2:
        pytest.skip(
            "environment does not provide >=2 real OpenMP threads "
            "(OMP_THREAD_LIMIT capped below NTHREADS) -- run in the OMP lane"
        )


def test_without_keep_order_resample_order_is_scrambled():
    """Deliberately forces a large, reliable divergence from loop-iteration
    order (see _OrderProbeData): without keep_order, the thread that is NOT
    slowed down races through its entire chunk before the slowed-down
    thread's first realization even finishes, so the two threads'
    call-count sequences (which, in correct loop-iteration order, would
    each form one contiguous block of N_PER_CHUNK rows) end up interleaved
    instead."""
    mc = _build_probe_mc(2024)
    mc.keep_order(False)
    mc.start_run()
    mc.run(N_TOTAL)
    mc.end_run()

    decoded = _decode_rows(mc)
    slots = [slot for slot, _count in decoded]
    _skip_unless_truly_parallel(slots)

    # Correct (loop-ordered) placement would be a single contiguous block
    # per thread_slot (e.g. [0,0,0,0,1,1,1,1] or [1,1,1,1,0,0,0,0]) -- i.e.
    # at most one "switch" between slots as row index increases. The forced
    # delay on thread_slot 0's first call guarantees more switches than
    # that: thread_slot 1 finishes its whole chunk while slot 0 is still
    # asleep on row 0, so slot 0 cannot occupy a single contiguous block.
    n_switches = sum(1 for a, b in zip(slots, slots[1:]) if a != b)
    assert n_switches > 1, f"expected scrambled row order, got slots={slots}"


def test_with_keep_order_resample_order_matches_loop_order():
    """The identical setup (same forced delay) as the test above, but with
    keep_order=True: the delay can no longer scramble anything, since the
    "omp ordered" region forces every thread's resample() call to happen in
    strict loop-iteration order regardless of how long any one call takes.
    Each thread's call-count sequence must land in one contiguous block."""
    mc = _build_probe_mc(2024)
    mc.keep_order(True)
    mc.start_run()
    mc.run(N_TOTAL)
    mc.end_run()

    decoded = _decode_rows(mc)
    slots = [slot for slot, _count in decoded]
    counts = [count for _slot, count in decoded]
    _skip_unless_truly_parallel(slots)

    n_switches = sum(1 for a, b in zip(slots, slots[1:]) if a != b)
    assert n_switches == 1, f"expected exactly one contiguous block boundary, got slots={slots}"

    # Within each contiguous block, call_count must be strictly increasing
    # from 0 -- i.e. rows appear in the exact order that thread issued them.
    for block_start in (0, N_PER_CHUNK):
        block_counts = counts[block_start:block_start + N_PER_CHUNK]
        assert block_counts == list(range(N_PER_CHUNK)), (
            f"expected in-order call_count 0..{N_PER_CHUNK - 1} within the "
            f"block, got {block_counts}"
        )


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
