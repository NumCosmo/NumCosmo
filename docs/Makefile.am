## Process this file with automake to produce Makefile.in

# Local xml files
XML_FILES= \
	introduction.xml \
	darkenergy.xml   \
	refs.xml

# Generated xml files
XML_GEN_FILES= \
	references.xml

# We require automake 1.6 at least.
AUTOMAKE_OPTIONS = 1.6

# This is a blank Makefile.am for using gtk-doc.
# Copy this to your project's API docs directory and modify the variables to
# suit your project. See the GTK+ Makefiles in gtk+/docs/reference for examples
# of using the various options.

# The name of the module, e.g. 'glib'.
DOC_MODULE=numcosmo

# Uncomment for versioned docs and specify the version of the module, e.g. '2'.
#DOC_MODULE_VERSION=2

# The top-level SGML file. You can change this if you want to.
DOC_MAIN_SGML_FILE=$(DOC_MODULE)-docs.sgml

$(DOC_MAIN_SGML_FILE): manual_head.xml

# Directories containing the source code, relative to $(srcdir).
# gtk-doc will search all .c and .h files beneath these paths
# for inline comments documenting functions and macros.
# e.g. DOC_SOURCE_DIR=../../../gtk ../../../gdk
DOC_SOURCE_DIR=$(top_srcdir)/numcosmo \
	       $(top_srcdir)/numcosmo/lss \
	       $(top_srcdir)/numcosmo/abc \
	       $(top_srcdir)/numcosmo/math \
	       $(top_srcdir)/numcosmo/model \
	       $(top_srcdir)/numcosmo/data \
	       $(top_srcdir)/numcosmo/perturbations \
	       $(top_srcdir)/numcosmo/xcor \
	       $(top_builddir)/numcosmo

DOC_BUILD_DIR=$(top_builddir)/numcosmo
DOC_TARGET_DIR=$(top_builddir)/docs

# Extra options to pass to gtkdoc-scangobj. Not normally needed.
SCANGOBJ_OPTIONS=

# Extra options to supply to gtkdoc-scan.
# e.g. SCAN_OPTIONS=--deprecated-guards="GTK_DISABLE_DEPRECATED"
SCAN_OPTIONS=--rebuild-types --rebuild-sections --ignore-decorators="NCM_INLINE"

# Extra options to supply to gtkdoc-mkdb.
# e.g. MKDB_OPTIONS=--xml-mode --output-format=xml
MKDB_OPTIONS=--sgml-mode --output-format=xml --name-space=nc --name-space=ncm \
     --ignore-files "${IGNORE_HFILES}"

# Extra options to supply to gtkdoc-mktmpl
# e.g. MKTMPL_OPTIONS=--only-section-tmpl
MKTMPL_OPTIONS=

# Extra options to supply to gtkdoc-mkhtml
MKHTML_OPTIONS=

# Extra options to supply to gtkdoc-fixref. Not normally needed.
# e.g. FIXXREF_OPTIONS=--extra-dir=../gdk-pixbuf/html --extra-dir=../gdk/html
FIXXREF_OPTIONS=--extra-dir=$(GLIB_PREFIX)/share/gtk-doc/html

# Used for dependencies. The docs will be rebuilt if any of these change.
# e.g. HFILE_GLOB=$(top_srcdir)/gtk/*.h
# e.g. CFILE_GLOB=$(top_srcdir)/gtk/*.c
HFILE_GLOB=$(top_srcdir)/numcosmo/*.h \
	   $(top_srcdir)/numcosmo/lss/*.h \
	   $(top_srcdir)/numcosmo/abc/*.h \
	   $(top_srcdir)/numcosmo/math/*.h \
	   $(top_srcdir)/numcosmo/model/*.h \
	   $(top_srcdir)/numcosmo/data/*.h \
	   $(top_srcdir)/numcosmo/perturbations/*.h

CFILE_GLOB=$(top_srcdir)/numcosmo/*.c \
	   $(top_srcdir)/numcosmo/lss/*.c \
	   $(top_srcdir)/numcosmo/abc/*.c \
	   $(top_srcdir)/numcosmo/math/*.c \
	   $(top_srcdir)/numcosmo/model/*.c \
	   $(top_srcdir)/numcosmo/data/*.c \
	   $(top_srcdir)/numcosmo/perturbations/*.c

# Extra header to include when scanning, which are not under DOC_SOURCE_DIR
# e.g. EXTRA_HFILES=$(top_srcdir}/contrib/extra.h
EXTRA_HFILES=

if ! HAVE_LIBNLOPT
ignore_nlopt = \
	ncm_fit_nlopt_enum.h \
	ncm_fit_nlopt.h
endif

# Header files or dirs to ignore when scanning. Use base file/dir names
# e.g. IGNORE_HFILES=gtkdebug.h gtkintl.h private_code
IGNORE_HFILES= \
              config.h                  \
              config_extra.h            \
              ncm_gsl_blas_types.h      \
              ncm_fit_nlopt_enum_meta.h \
              ncm_sphere_map_block.c    \
              ncm_flapack.h             \
              nc_hipert_private.h       \
              gsl_rstat.h               \
              gsl_cqp.h                 \
              initial_point.h           \
              rquantile.c               \
              $(ignore_nlopt)           \
              Faddeeva.h                \
              levmar                    \
              libcuba                   \
              plc			\
              class			\
              class/source		\
              class/hyrec		\
              class/include             \
              toeplitz                  \
              sundials                  \
              rcm.h                     \
              rcm.c

# Images to copy into HTML directory.
# e.g. HTML_IMAGES=$(top_srcdir)/gtk/stock-icons/stock_about_24.png
HTML_IMAGES=

# Extra SGML files that are included by $(DOC_MAIN_SGML_FILE).
# e.g. content_files=running.sgml building.sgml changes-2.0.sgml
content_files=$(XML_FILES) $(XML_GEN_FILES) references.bib container.css

# SGML files where gtk-doc abbrevations (#GtkWidget) are expanded
# These files must be listed here *and* in content_files
# e.g. expand_content_files=running.sgml
expand_content_files=

# CFLAGS and LDFLAGS for compiling gtkdoc-scangobj with your library.
# Only needed if you are using gtkdoc-scangobj to dynamically query widget
# signals and properties.
# e.g. GTKDOC_CFLAGS=-I$(top_srcdir) -I$(top_builddir) $(GTK_DEBUG_FLAGS)
# e.g. GTKDOC_LIBS=$(top_builddir)/gtk/$(gtktargetlib)
GTKDOC_CFLAGS=$(GLIB_CFLAGS)
GTKDOC_LIBS=$(top_builddir)/numcosmo/libnumcosmo.la

# This includes the standard gtk-doc make rules, copied by gtkdocize.
include $(top_srcdir)/gtk-doc.make

setup-build.stamp: $(addprefix $(top_srcdir)/docs/, $(content_files)) $(expand_content_files)

# Other files to distribute
# e.g. EXTRA_DIST += version.xml.in
EXTRA_DIST += $(XML_FILES) $(XML_GEN_FILES) references.bib container.css addmathjax.pl xmltogtk.pl

# Files not to distribute
# for --rebuild-types in $(SCAN_OPTIONS), e.g. $(DOC_MODULE).types
# for --rebuild-sections in $(SCAN_OPTIONS) e.g. $(DOC_MODULE)-sections.txt
DISTCLEANFILES =

# Comment this out if you want 'make check' to test you doc status
# and run some sanity checks
if ENABLE_GTK_DOC
#TESTS_ENVIRONMENT = cd $(top_builddir)/dsdfocs && \
TESTS_ENVIRONMENT = \
  DOC_MODULE=$(DOC_MODULE) DOC_MAIN_SGML_FILE=$(DOC_MAIN_SGML_FILE) \
  SRCDIR=$(abs_srcdir) BUILDDIR=$(abs_builddir)
#TESTS = $(GTKDOC_CHECK)
endif

########################################################################

man_MANS := darkenergy.1

EXTRA_DIST += $(man_MANS)

if ENABLE_MAN

%.1: %.xml
	@@XSLTPROC@ -nonet http://docbook.sourceforge.net/release/xsl/current/manpages/docbook.xsl $<

dist-local-check-mans-enabled:
	if grep "Man generation disabled" $(man_MANS) >/dev/null; then $(RM) $(man_MANS); fi
else

$(man_MANS):
	@echo Man generation disabled.  Creating dummy $@.  Configure with --enable-man to enable it.
	@echo Man generation disabled.  Remove this file, configure with --enable-man, and rebuild > $@

dist-local-check-mans-enabled:
	@echo "*** --enable-man must be used in order to make dist"
	false

endif

MAINTAINERCLEANFILES = $(man_MANS) $(BUILT_SOURCES)

BUILT_EXTRA_DIST = $(man_MANS)

CLEANFILES += bibfile.stamp addmathjax.stamp

dist-hook-local: dist-local-check-mans-enabled all-local

numcosmo-docs-clean: clean
	cd $(srcdir) && rm -rf xml html

bibfile.stamp: $(top_srcdir)/docs/references.bib
	-$(GTK_DOC_V_SETUP)if test "$(abs_srcdir)" != "$(abs_builddir)" ; then \
	    cp -p $(top_srcdir)/docs/references.bib $(DOC_TARGET_DIR); \
	fi
	$(AM_V_at)touch bibfile.stamp

%.xml: %.tex bibfile.stamp
	$(AM_V_at)echo "  DOC   Building XML $@ from $<"
	$(AM_V_at)rm -f $*.aux $*.bbl 
	$(AM_V_at)latex --interaction=nonstopmode $< >  $(@F).htlatex.log 2>&1
	$(AM_V_at)bibtex8 $(*F)                      >> $(@F).htlatex.log 2>&1
	$(AM_V_at)latex --interaction=nonstopmode $< >> $(@F).htlatex.log 2>&1
	$(AM_V_at)latex --interaction=nonstopmode $< >> $(@F).htlatex.log 2>&1
	$(AM_V_at)htlatex $< "xhtml,docbook" " -cunihtf" "-cdocbk" "--interaction=nonstopmode" >> $(@F).htlatex.log 2>&1
	$(AM_V_at)$(top_srcdir)/docs/xmltogtk.pl $(DOC_TARGET_DIR)/$(@F) > $(@F).xmltogtk.log 2>&1
	$(AM_V_at)if test ! "$(DOC_TARGET_DIR)/$(@F)" -ef "$(<D)/$(@F)"; then cp $(DOC_TARGET_DIR)/$(@F) $(<D)/$(@F); fi

if ENABLE_GTK_DOC
HTML_BUILD_STAMP += addmathjax.stamp
endif

addmathjax.stamp: html-build.stamp
	$(AM_V_at)echo "  DOC   Adding support to MathJax"
	$(AM_V_at)$(top_srcdir)/docs/addmathjax.pl $(DOC_TARGET_DIR)/html/*.html > $(@F).addmathjax.log 2>&1
	$(AM_V_at)cp -p $(top_srcdir)/docs/container.css $(DOC_TARGET_DIR)/html/
	$(AM_V_at)touch addmathjax.stamp

-include $(top_srcdir)/git.mk
