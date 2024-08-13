include $(top_srcdir)/am/doxygen.am
include $(top_srcdir)/am/cppcheck.am

ACLOCAL_AMFLAGS = -I m4

SUBDIRS = src test

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = dhsm@DEBUG_SUFFIX@.pc

if SERVANT
pkgconfig_DATA += dhsm-servant@DEBUG_SUFFIX@.pc
endif

if MASTER
pkgconfig_DATA += dhsm-master@DEBUG_SUFFIX@.pc
endif

EXTRA_DIST = autogen.sh doc/doxyfile
