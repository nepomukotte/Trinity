include $(top_srcdir)/am/doxygen.am
include $(top_srcdir)/am/cppcheck.am

ACLOCAL_AMFLAGS = -I m4

SUBDIRS = src test

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = @CCFG_LIB_NAME@.pc

EXTRA_DIST = autogen.sh doc/doxyfile

dist-hook:
	mkdir -p $(distdir)/doc
	cp -p $(srcdir)/doc/*.dox $(distdir)/doc/

include $(top_srcdir)/am/pkg-linux.am
include $(top_srcdir)/am/pkg-win.am
