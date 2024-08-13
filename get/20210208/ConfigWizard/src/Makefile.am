include $(top_srcdir)/am/autotroll.am
AUTOMAKE_OPTIONS = subdir-objects

AM_CPPFLAGS = -I@top_srcdir@/src -Isrc -I@top_builddir@/src \
 $(UTL_CFLAGS) $(CONFIG_WIDGETS_CFLAGS) \
 $(QT_CPPFLAGS) $(QWT_CPPFLAGS) $(BOOST_CPPFLAGS) $(QWT_PLOT_3D_CPPFLAGS) \
 $(ROOTCFLAGS)
AM_CXXFLAGS = $(QT_CXXFLAGS) $(QWT_PLOT_3D_CXXFLAGS)

bin_PROGRAMS=get-config-wizard@DEBUG_SUFFIX@

moc_sources = \
	get/ConfigWizardMainWindow.moc.cpp \
	get/cfg/ChannelConfigModel.moc.cpp \
	get/geo/GeometryModel.moc.cpp \
	get/geo/EditMappingDialog.moc.cpp \
	get/root/RangeDialog.moc.cpp \
	get/root/HistogramWidget.moc.cpp \
	get/root/HistogramDialog.moc.cpp
ui_files = \
	get/ConfigWizardMainWindow.ui \
	get/geo/EditMappingDialog.ui \
	get/root/RangeDialog.ui
ui_headers = \
	get/ConfigWizardMainWindow.ui.h \
	get/geo/EditMappingDialog.ui.h \
	get/root/RangeDialog.ui.h
gen_sources = $(moc_sources) $(qrc_sources) $(ui_headers)

headers = \
	get/ConfigWizardMainWindow.h \
	get/utl/uid.h \
	get/cfg/ChannelConfigModel.h \
	get/geo/GeometryModel.h \
	get/geo/EditMappingDialog.h \
	get/root/RangeDialog.h \
	get/root/HistogramWidget.h \
	get/root/HistogramDialog.h

get_config_wizard@DEBUG_SUFFIX@_SOURCES= \
	get-config-wizard.cpp \
	get/utl/uid.cpp \
	get/ConfigWizardMainWindow.cpp \
	get/cfg/ChannelConfigModel.cpp \
	get/geo/GeometryModel.cpp \
	get/geo/EditMappingDialog.cpp \
	get/root/RangeDialog.cpp \
	get/root/HistogramWidget.cpp \
	get/root/HistogramDialog.cpp \
	$(headers)
nodist_get_config_wizard@DEBUG_SUFFIX@_SOURCES = $(gen_sources)
get_config_wizard@DEBUG_SUFFIX@_LDFLAGS = $(QT_LDFLAGS) -L$(ROOTLIBDIR)
get_config_wizard@DEBUG_SUFFIX@_LDADD = $(UTL_LIBS) $(CONFIG_WIDGETS_LIBS) $(ROOTGLIBS)
if HAVE_QTROOT
get_config_wizard@DEBUG_SUFFIX@_LDADD += -lQtRoot -lGQt
endif
get_config_wizard@DEBUG_SUFFIX@_LDADD += $(QT_LIBS)

EXTRA_DIST = $(ui_files)
BUILT_SOURCES = $(gen_sources)