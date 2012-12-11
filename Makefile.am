## Process this file with automake to produce Makefile.in
## Created by Anjuta

SUBDIRS = numcosmo tools darkenergy data docs tests examples

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = numcosmo.pc

ACLOCAL_AMFLAGS = -I m4
DISTCHECK_CONFIGURE_FLAGS=--enable-gtk-doc --enable-introspection

numcosmodocdir = ${docdir}-$(VERSION)
numcosmodoc_DATA = \
	README     \
	COPYING    \
	AUTHORS    \
	ChangeLog  \
	NEWS

#DISTCLEANFILES = intltool-extract intltool-merge intltool-update

EXTRA_DIST = \
	INSTALL                   \
	$(numcosmodoc_DATA)       \
	numcosmo.pc.in

AM_CFLAGS = 

# Copy all the spec files. Of cource, only one is actually used.
dist-hook:
	for specfile in *.spec; do \
		if test -f $$specfile; then \
			cp -p $$specfile $(distdir); \
		fi \
	done
