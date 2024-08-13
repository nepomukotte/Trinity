include $(top_srcdir)/am/cppcheck.am

ACLOCAL_AMFLAGS = -I m4

resources =                            \
	res/images/connect_established.png \
	res/images/connect_no.png          \
	res/images/ccfg-db-logo.svg

SUBDIRS = src

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = @DATABASE_WIDGETS_LIB_NAME@.pc

EXTRA_DIST = autogen.sh $(resources)
