include $(top_srcdir)/am/autotroll.am

AUTOMAKE_OPTIONS = subdir-objects

headers = get/EccClient.h \
	get/CatchAllExceptions.h \
	get/daq/DiscardStorage.h \
	get/daq/CustomStorage.h \
	get/daq/DataRouterWizard.h \
	get/daq/FramePublisher.h \
	get/utl/SlotsLoggingBackend.h \
	get/Console.h \
	get/DataBloc.h \
	get/GetTest.h \
	get/TestController.h \
	get/TestInfo.h \
	get/complexnumber.h \
	get/pixmaps.h \
	get/plot.h \
	get/viewgraph.h \
	get/EditChannelMaskDialog.h \
	get/ChannelMask.h \
	get/TestInfo.h \
	get/plot.h \
	get/viewgraph.h \
	get/PlotCurve.h \
	get/EditTestDialog.h \
	get/TestSummary.h \
	get/WorkspaceManager.h \
	get/WorkspaceWizard.h \
	get/SelectTestDialog.h \
	get/GetControllerException.h \
	get/TestSummaryWidget.h \
	get/WorkspaceModel.h \
	get/Channel.h \
	get/CoBoEvent.h \
	get/ControllerOptions.h \
	get/QwtSymbolDialog.h

sources = get/EccClient.cpp \
	get/daq/DiscardStorage.cpp \
	get/daq/CustomStorage.cpp \
	get/daq/DataRouterWizard.cpp \
	get/daq/FramePublisher.cpp \
	get/utl/SlotsLoggingBackend.cpp \
	get/EditChannelMaskDialog.cpp \
	get/ChannelMask.cpp \
	get/Console.cpp \
	get/DataBloc.cpp \
	get/GetTest.cpp \
	get/TestController.cpp \
	get/TestInfo.cpp \
	get/main.cpp \
	get/plot.cpp \
	get/EditTestDialog.cpp \
	get/viewgraph.cpp \
	get/PlotCurve.cpp \
	get/TestSummary.cpp \
	get/WorkspaceManager.cpp \
	get/WorkspaceWizard.cpp \
	get/SelectTestDialog.cpp \
	get/GetControllerException.cpp \
	get/TestSummaryWidget.cpp \
	get/WorkspaceModel.cpp \
	get/Channel.cpp \
	get/CoBoEvent.cpp \
	get/ControllerOptions.cpp \
	get/QwtSymbolDialog.cpp

moc_sources = \
	get/daq/DataRouterWizard.moc.cpp \
	get/daq/FramePublisher.moc.cpp \
	get/utl/SlotsLoggingBackend.moc.cpp \
	get/viewgraph.moc.cpp \
	get/GetTest.moc.cpp \
	get/Console.moc.cpp \
	get/DataBloc.moc.cpp \
	get/plot.moc.cpp \
	get/TestController.moc.cpp \
	get/EditChannelMaskDialog.moc.cpp \
	get/EditTestDialog.moc.cpp \
	get/WorkspaceWizard.moc.cpp \
	get/SelectTestDialog.moc.cpp \
	get/TestSummaryWidget.moc.cpp \
	get/QwtSymbolDialog.moc.cpp

qrc_files = GetController.qrc
qrc_sources = GetController.qrc.cpp

ui_files = \
	get/ControlWidget.ui \
	get/GetTest.ui \
	get/GetTest.ui \
	get/EditTestDialog.ui \
	get/TestSummaryWidget.ui

ui_headers = \
	get/ControlWidget.ui.h \
	get/GetTest.ui.h \
	get/GetTest.ui.h \
	get/EditTestDialog.ui.h \
	get/TestSummaryWidget.ui.h

gen_sources = $(moc_sources) $(qrc_sources) $(ui_headers)

AM_CPPFLAGS  = -I@top_srcdir@/src -I@top_builddir@/src -I@top_builddir@/get -Iget
AM_CPPFLAGS += $(QT_CPPFLAGS)
AM_CPPFLAGS += $(QWT_CPPFLAGS)
AM_CPPFLAGS += $(ICE_CPPFLAGS)
AM_CPPFLAGS += $(ICESTORM_CPPFLAGS)
AM_CPPFLAGS += $(BOOST_CPPFLAGS)
AM_CPPFLAGS += $(COMPOUND_CONFIG_CFLAGS)
AM_CPPFLAGS += $(CONFIG_WIDGETS_CFLAGS)

AM_CXXFLAGS = $(QT_CXXFLAGS)

bin_PROGRAMS = GetController@DEBUG_SUFFIX@
GetController@DEBUG_SUFFIX@_SOURCES        = $(sources)
nodist_GetController@DEBUG_SUFFIX@_SOURCES = $(gen_sources)
GetController@DEBUG_SUFFIX@_LDFLAGS        = $(QT_LDFLAGS) $(QWT_LDFLAGS) $(ICE_LDFLAGS) $(BOOST_LDFLAGS)
GetController@DEBUG_SUFFIX@_LDADD  = $(QT_LIBS) $(QWT_LIBS)
GetController@DEBUG_SUFFIX@_LDADD += $(GET_BENCH_LIBS)
GetController@DEBUG_SUFFIX@_LDADD += $(ICE_LIBS)
GetController@DEBUG_SUFFIX@_LDADD += $(ICESTORM_LIBS)
GetController@DEBUG_SUFFIX@_LDADD += $(CONFIG_WIDGETS_LIBS)
GetController@DEBUG_SUFFIX@_LDADD += $(COMPOUND_CONFIG_LIBS)
GetController@DEBUG_SUFFIX@_LDADD += $(CONFIG_WIDGETS_LIBS)
GetController@DEBUG_SUFFIX@_LDADD += $(BOOST_SYSTEM_LIB) 
GetController@DEBUG_SUFFIX@_LDADD += $(BOOST_THREAD_LIB)
GetController@DEBUG_SUFFIX@_LDADD += $(BOOST_PROGRAM_OPTIONS_LIB)

BUILT_SOURCES = $(gen_sources)
EXTRA_DIST = $(qrc_files) $(ui_files) $(resources) $(headers)
