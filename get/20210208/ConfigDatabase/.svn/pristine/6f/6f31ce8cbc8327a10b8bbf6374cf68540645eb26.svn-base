include $(top_srcdir)/am/doxygen.am
include $(top_srcdir)/am/cppcheck.am

ACLOCAL_AMFLAGS = -I m4

SUBDIRS = src test

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = @CONFIG_DATABASE_LIB_NAME@.pc

EXTRA_DIST = autogen.sh doc/doxyfile
