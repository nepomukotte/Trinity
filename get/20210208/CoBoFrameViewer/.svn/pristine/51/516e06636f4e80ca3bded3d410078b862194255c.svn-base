
AM_CPPFLAGS = -I@top_srcdir@/src -I@top_builddir@/src \
 $(QT_CPPFLAGS) $(BOOST_CPPFLAGS) $(MFM_CFLAGS)
AM_CXXFLAGS = $(QT_CXXFLAGS)

bin_PROGRAMS = cobo-frame-graw2gmat@DEBUG_SUFFIX@
cobo_frame_graw2gmat@DEBUG_SUFFIX@_SOURCES        = cobo-frame-graw2gmat.cpp
cobo_frame_graw2gmat@DEBUG_SUFFIX@_LDFLAGS        = $(QT_LDFLAGS) $(QWT_LDFLAGS) $(BOOST_LDFLAGS) $(MFM_LDFLAGS)
cobo_frame_graw2gmat@DEBUG_SUFFIX@_LDADD		  = ../get/libcobo-event@DEBUG_SUFFIX@.la
cobo_frame_graw2gmat@DEBUG_SUFFIX@_LDADD          += $(QT_LIBS) $(QWT_LIBS)
cobo_frame_graw2gmat@DEBUG_SUFFIX@_LDADD          += $(BOOST_PROGRAM_OPTIONS_LIB) $(BOOST_FILESYSTEM_LIB) $(BOOST_SYSTEM_LIB) $(MFM_LIBS)