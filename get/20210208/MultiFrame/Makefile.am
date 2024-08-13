include $(top_srcdir)/am/doxygen.am
include $(top_srcdir)/am/cppcheck.am

ACLOCAL_AMFLAGS = -I m4

SUBDIRS = src test

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = @MFM_LIB_NAME@.pc

EXTRA_DIST = autogen.sh doc/doxyfile

dist-hook:
	mkdir -p $(distdir)/doc
	cp -p $(srcdir)/doc/*.dox $(distdir)/doc/

dist_noinst_SCRIPTS = autogen.sh cleanup.sh
dist_noinst_DATA = COPYING AUTHORS doc/doxyfile