include $(top_srcdir)/am/autotroll.am
AUTOMAKE_OPTIONS = subdir-objects

common_sources = \
	Channel.cpp \
	CoBoEvent.cpp

common_headers = \
	Channel.h \
	CoBoEvent.h

# Library
noinst_LTLIBRARIES = libcobo-event@DEBUG_SUFFIX@.la
libcobo_event@DEBUG_SUFFIX@_la_SOURCES = $(common_sources)
libcobo_event@DEBUG_SUFFIX@_la_LDFLAGS = $(QT_LDFLAGS)
libcobo_event@DEBUG_SUFFIX@_la_LIBADD  = $(QT_LIBS)
	
viewer_sources = \
	CoBoStandardDisplay.cpp \
	CoBoFrameViewer.cpp \
	ChannelMask.cpp \
	EditChannelMaskDialog.cpp \
	PlotCurve.cpp \
	plot.cpp \
	viewgraph.cpp \
	QwtSymbolDialog.cpp
if BUILD_3DVIEWER
viewer_sources += CoBoSurfDisplay.cpp
endif

viewer_headers = \
	CoBoDisplay.h \
	CoBoStandardDisplay.h \
	CoBoFrameViewer.h \
	ChannelMask.h \
	EditChannelMaskDialog.h \
	PlotCurve.h \
	plot.h \
	viewgraph.h \
	complexnumber.h \
	QwtSymbolDialog.h \
	ana/algorithm.h
if BUILD_3DVIEWER
viewer_headers += CoBoSurfDisplay.h
endif

moc_sources = \
	CoBoStandardDisplay.moc.cpp \
	EditChannelMaskDialog.moc.cpp \
	CoBoFrameViewer.moc.cpp \
	plot.moc.cpp \
	viewgraph.moc.cpp \
	QwtSymbolDialog.moc.cpp
if BUILD_3DVIEWER
moc_sources += CoBoSurfDisplay.moc.cpp
endif

qrc_files = ../CoBoFrameViewer.qrc
qrc_sources = ../CoBoFrameViewer.qrc.cpp

gen_sources = $(moc_sources) $(qrc_sources)

AM_CPPFLAGS = -I@top_srcdir@/src -I@top_builddir@/src \
 $(QT_CPPFLAGS) $(QWT_CPPFLAGS) $(BOOST_CPPFLAGS) $(MFM_CFLAGS) $(QWT_PLOT_3D_CPPFLAGS)
AM_CXXFLAGS = $(QT_CXXFLAGS) $(QWT_PLOT_3D_CXXFLAGS)

bin_PROGRAMS=cobo-frame-viewer@DEBUG_SUFFIX@
cobo_frame_viewer@DEBUG_SUFFIX@_SOURCES        = $(viewer_sources) cobo-frame-viewer.cpp
nodist_cobo_frame_viewer@DEBUG_SUFFIX@_SOURCES = $(gen_sources)
cobo_frame_viewer@DEBUG_SUFFIX@_LDFLAGS        = $(QT_LDFLAGS) $(QWT_LDFLAGS) $(BOOST_LDFLAGS) $(MFM_LDFLAGS) $(QWT_PLOT_3D_LDFLAGS)
cobo_frame_viewer@DEBUG_SUFFIX@_LDADD		   = libcobo-event@DEBUG_SUFFIX@.la
cobo_frame_viewer@DEBUG_SUFFIX@_LDADD          += $(QT_LIBS) $(QWT_LIBS) $(MFM_LIBS)
cobo_frame_viewer@DEBUG_SUFFIX@_LDADD          += $(BOOST_PROGRAM_OPTIONS_LIB) $(BOOST_FILESYSTEM_LIB) $(BOOST_SYSTEM_LIB)
cobo_frame_viewer@DEBUG_SUFFIX@_LDADD          += $(QWT_PLOT_3D_LIBS)

install-exec-hook:
	cd $(DESTDIR)$(bindir) && \
	$(LN_S) -f cobo-frame-viewer@DEBUG_SUFFIX@$(EXEEXT) cobo-frame-graw2disp@DEBUG_SUFFIX@$(EXEEXT)

BUILT_SOURCES = $(gen_sources)
EXTRA_DIST = $(qrc_files) $(common_headers) $(viewer_headers)
