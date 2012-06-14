## Process this file with automake to produce Makefile.in

## Created by Anjuta

INCLUDES = 

AM_CFLAGS = \
	 -Wall \
	 -DG_LOG_DOMAIN=\"COSMOCLUSTERLIB\" \
	 $(COSMOLIB_CFLAGS) \
	 $(GLIB_CFLAGS) \
	 $(GSL_CFLAGS) \
	 $(SUNDIALS_CFLAGS) \
	 -I$(top_srcdir)

SUBDIRS = 

lib_LTLIBRARIES = \
	libcosmocluster.la

libcosmocluster_la_SOURCES = \
	catalog_parser.c \
	catalog_parser.h \
	growth_function.c \
	growth_function.h \
	mass_function.c \
	mass_function.h \
	multiplicity_function.c \
	multiplicity_function.h \
	read_matrix.c \
	read_matrix.h \
	transfer_function.c \
	transfer_function.h \
	variance.c \
	variance.h \
	window.c \
	window.h \
	cosmoclusterlib.h \
	data_cluster_abundance.c \
	data_cluster_abundance.h \
	cluster_abundance_distribution.c \
	cluster_abundance_distribution.h \
        bias_function_type.c \
        bias_function_type.h \
        mean_bias_function.c \
        mean_bias_function.h \
	mass_distribution.c \
	mass_distribution.h \
	print_data.c \
	print_data.h

libcosmocluster_la_LIBADD = \
	$(SUNDIALS_LIBS) \
	$(GLIB_LIBS) \
	$(GSL_LIBS) \
	$(COSMOLIB_LIBS)

libcosmocluster_la_LDFLAGS = -version-info $(COSMOCLUSTERLIB_LT_VERSION)

AM_CPPFLAGS = \
	$(SUNDIALS_CFLAGS) \
	$(COSMOLIB_CFLAGS) \
	$(GLIB_CFLAGS) \
	$(GSL_CFLAGS)

cosmocluster_include_HEADERS = \
	catalog_parser.h \
	growth_function.h \
	mass_function.h \
	multiplicity_function.h \
	read_matrix.h \
	transfer_function.h \
	variance.h \
	window.h \
	cluster_abundance_distribution.h \
	cosmoclusterlib.h \
	macros.h \
	data_cluster_abundance.h \
	mass_distribution.h \
	bias_function_type.h \
	mean_bias_function.h \
	config.h \
	print_data.h

cosmocluster_includedir = \
	$(pkgincludedir)

EXTRA_DIST = \
	$(cosmocluster_include_HEADERS) \
	config.h.in

-include $(INTROSPECTION_MAKEFILE)
INTROSPECTION_GIRS =
INTROSPECTION_SCANNER_ARGS = --add-include-path=$(srcdir) -I$(srcdir)/../ -I$(top_builddir) \
  --identifier-prefix=Ccl --warn-all --verbose

INTROSPECTION_COMPILER_ARGS = --includedir=$(srcdir)

if HAVE_INTROSPECTION

introspection_sources = cosmoclusterlib.h $(libcosmocluster_la_SOURCES)

Cosmoclusterlib-1.0.gir: libcosmocluster.la
Cosmoclusterlib_1_0_gir_INCLUDES = GObject-2.0 Cosmolib-1.0
Cosmoclusterlib_1_0_gir_CFLAGS = $(INCLUDES)
Cosmoclusterlib_1_0_gir_LIBS = libcosmocluster.la
Cosmoclusterlib_1_0_gir_FILES = $(introspection_sources)
INTROSPECTION_GIRS += Cosmoclusterlib-1.0.gir

girdir = $(datadir)/gir-1.0
gir_DATA = $(INTROSPECTION_GIRS)   
  
typelibdir = $(libdir)/girepository-1.0
typelib_DATA = $(INTROSPECTION_GIRS:.gir=.typelib)
      
CLEANFILES = $(gir_DATA) $(typelib_DATA)
  
endif
