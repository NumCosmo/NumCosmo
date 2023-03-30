## Process this file with automake to produce Makefile.in

# gcov support
include $(top_srcdir)/aminclude_static.am

SUBDIRS = numcosmo tools data docs tests examples notebooks

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = numcosmo.pc

ACLOCAL_AMFLAGS = -I m4
DISTCHECK_CONFIGURE_FLAGS=--enable-gtk-doc --enable-introspection

numcosmodocdir = ${docdir}-$(VERSION)
numcosmodoc_DATA = \
	README.md  \
	README     \
	COPYING    \
	AUTHORS    \
	ChangeLog  \
	NEWS

if HAVE_PYTHON
include numcosmo_py.mk
endif # HAVE_PYTHON

EXTRA_DIST = \
	INSTALL                \
	$(numcosmodoc_DATA)    \
	numcosmo-valgrind.supp \
        tap-driver.sh          \
        tap-test

AM_CFLAGS = 

# Copy all the spec files. Of cource, only one is actually used.
dist-hook:
	for specfile in *.spec; do \
		if test -f $$specfile; then \
			cp -p $$specfile $(distdir); \
		fi \
	done

