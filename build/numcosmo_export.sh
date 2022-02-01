# NumCosmo buildir exports necessary to use 
# the library without installation

export NUMCOSMO_DIR=$(cd "/home/cinthia/NumCosmo/build/.."; pwd)
export NUMCOSMO_BUILD_DIR=$(cd "/home/cinthia/NumCosmo/build"; pwd)
export NUMCOSMO_DATA_DIR="${NUMCOSMO_DIR}"

export CPATH="${CPATH:+$CPATH:}${NUMCOSMO_DIR}:${NUMCOSMO_BUILD_DIR}"
export LIBRARY_PATH="${LIBRARY_PATH:+$LIBRARY_PATH:}${NUMCOSMO_BUILD_DIR}/numcosmo/.libs"
export LD_LIBRARY_PATH="${LD_LIBRARY_PATH:+$LD_LIBRARY_PATH:}${NUMCOSMO_BUILD_DIR}/numcosmo/.libs"
export GI_TYPELIB_PATH="${GI_TYPELIB_PATH:+$GI_TYPELIB_PATH:}${NUMCOSMO_BUILD_DIR}/numcosmo"
