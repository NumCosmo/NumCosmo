#!/usr/bin/sh


libtool --mode=compile gcc -Wall -c nc_curve.c -o nc_curve.lo \
  -I${NUMCOSMO_DIR} -I${NUMCOSMO_BUILD_DIR} `pkg-config --cflags glib-2.0`

libtool --mode=link gcc -Wall nc_curve.lo -o libnumcosmo_curve.la -shared \
  -I${NUMCOSMO_DIR} -I${NUMCOSMO_BUILD_DIR} `pkg-config --cflags glib-2.0`
