include $(top_srcdir)/am/doxygen.am
include $(top_srcdir)/am/cppcheck.am

ACLOCAL_AMFLAGS = -I m4

SUBDIRS = src

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = @CONFIG_WIDGETS_LIB_NAME@.pc

EXTRA_DIST = autogen.sh doc/doxyfile

dist-hook:
	mkdir -p $(distdir)/res/images
	cp -rp $(srcdir)/res/images/* $(distdir)/res/images/

include $(top_srcdir)/am/pkg-linux.am
include $(top_srcdir)/am/pkg-win.am
