AM_CPPFLAGS = -I@top_srcdir@/src -I@top_srcdir@/src/root $(QT_CPPFLAGS) $(BOOST_CPPFLAGS) $(MFM_CFLAGS) $(CCFITS_CFLAGS)

bin_PROGRAMS = 

if HAVE_CCFITS
bin_PROGRAMS += cobo-frame-graw2fits@DEBUG_SUFFIX@
cobo_frame_graw2fits@DEBUG_SUFFIX@_SOURCES = cobo-frame-graw2fits.cpp graw2fits.h graw2fits.cpp
cobo_frame_graw2fits@DEBUG_SUFFIX@_LDFLAGS = $(CCFITS_LDFLAGS) $(BOOST_LDFLAGS) $(QT_LDFLAGS)
cobo_frame_graw2fits@DEBUG_SUFFIX@_LDADD = ../get/libcobo-event@DEBUG_SUFFIX@.la
cobo_frame_graw2fits@DEBUG_SUFFIX@_LDADD += $(MFM_LIBS) $(CCFITS_LIBS)
cobo_frame_graw2fits@DEBUG_SUFFIX@_LDADD += $(BOOST_PROGRAM_OPTIONS_LIB) $(BOOST_FILESYSTEM_LIB) $(BOOST_SYSTEM_LIB) $(QT_LIBS)
endif