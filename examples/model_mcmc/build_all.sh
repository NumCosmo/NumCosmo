#!/usr/bin/sh


# These commands are used to build the library libnumcosmo_curve.la which
# already include numcosmo as a dependency.

libtool --mode=compile gcc -Wall -c nc_curve.c -o nc_curve.lo \
  -I${NUMCOSMO_DIR} -I${NUMCOSMO_BUILD_DIR} -I${NUMCOSMO_BUILD_DIR}/numcosmo \
  `pkg-config --cflags glib-2.0`

libtool --mode=compile gcc -Wall -c nc_curve_linear.c -o nc_curve_linear.lo  \
  -I${NUMCOSMO_DIR} -I${NUMCOSMO_BUILD_DIR} -I${NUMCOSMO_BUILD_DIR}/numcosmo \
  `pkg-config --cflags glib-2.0`

libtool --mode=compile gcc -Wall -c nc_data_curve.c -o nc_data_curve.lo \
  -I${NUMCOSMO_DIR} -I${NUMCOSMO_BUILD_DIR} -I${NUMCOSMO_BUILD_DIR}/numcosmo \
  `pkg-config --cflags glib-2.0`

libtool --mode=link gcc -Wall -o libnumcosmo_curve.la -rpath `pwd` \
  nc_curve.lo nc_curve_linear.lo nc_data_curve.lo                  \
  ${NUMCOSMO_BUILD_DIR}/numcosmo/libnumcosmo.la

libtool --mode=link gcc -Wall curve_generate.c -o curve_generate   \
  libnumcosmo_curve.la                                             \
  -I${NUMCOSMO_DIR} -I${NUMCOSMO_BUILD_DIR} -I${NUMCOSMO_BUILD_DIR}/numcosmo \
  `pkg-config --cflags glib-2.0`

