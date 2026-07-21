#!/bin/sh
#
# detect_omp_threads.sh
#
# Sets OMP_NUM_THREADS/OMP_THREAD_LIMIT to the number of PHYSICAL CPU cores
# available right now, then execs the real command. Deliberately physical,
# not logical/hyperthreaded: this is used for CPU-bound floating-point OpenMP
# work, where hyperthread siblings share the same execution units and add
# scheduling overhead/contention rather than real throughput -- using the
# logical (hyperthreaded) count there oversubscribes the real hardware by
# ~2x for no gain.
#
# This is deliberately run fresh on every `meson test` invocation (as a test
# wrapper) rather than computed once at `meson setup` time: a builddir can be
# configured on one machine (or container/cgroup) and later run in a
# differently-sized one, and a value baked in at configure time would then be
# wrong -- either idling most of the machine or (worse) oversubscribing a
# smaller one.
#
# Prefers, in order:
#   - /proc/cpuinfo (Linux): counts unique (physical id, core id) pairs, i.e.
#     physical cores, not hyperthreads. Falls through to the logical count
#     below if the topology fields are absent/degenerate (common in VMs that
#     don't expose real hyperthread topology to the guest).
#   - sysctl -n hw.physicalcpu (macOS): the physical-core equivalent of
#     hw.logicalcpu.
#   - nproc / getconf _NPROCESSORS_ONLN: logical CPU count, used only if
#     neither of the above worked -- better than nothing, but does not
#     distinguish hyperthreads from real cores.
#   - 4, if nothing above works.

set -eu

n=""

if [ -r /proc/cpuinfo ]; then
    n=$(awk '
        /^physical id/ { p = $NF }
        /^core id/     { print p "," $NF }
    ' /proc/cpuinfo | sort -u | wc -l | tr -d ' ')
fi

if { [ -z "$n" ] || [ "$n" = "0" ]; } && command -v sysctl >/dev/null 2>&1; then
    n=$(sysctl -n hw.physicalcpu 2>/dev/null || true)
fi

if [ -z "$n" ] || [ "$n" = "0" ]; then
    if command -v nproc >/dev/null 2>&1; then
        n=$(nproc 2>/dev/null || true)
    elif command -v getconf >/dev/null 2>&1; then
        n=$(getconf _NPROCESSORS_ONLN 2>/dev/null || true)
    fi
fi

case "$n" in
    '' | *[!0-9]* | 0) n=4 ;;
esac

export OMP_NUM_THREADS="$n"
export OMP_MAX_ACTIVE_LEVELS=1
export OMP_THREAD_LIMIT="$n"

exec "$@"
