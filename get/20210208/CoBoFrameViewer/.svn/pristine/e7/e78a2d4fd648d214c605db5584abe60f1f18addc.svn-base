
AM_CPPFLAGS = -I@top_srcdir@/src -I@top_srcdir@/src/root $(QT_CPPFLAGS) $(ROOTCFLAGS) $(BOOST_CPPFLAGS) $(MFM_CFLAGS)

dict_sources = 
dict_headers = 

SUFFIXES = Dict.cpp Dict.h
%Dict.cpp %Dict.h : %.h
	$(ROOTCINT) -f $*Dict.cpp -c $<

if BUILD_GRAW2ROOT
dict_sources += T2KEventDict.cpp
if !HAVE_ROOT6
dict_headers += T2KEventDict.h
endif
if HAVE_ROOT6
dict_headers += T2KEventDict_rdict.pcm
endif

lib_LTLIBRARIES = libcobo-frame-graw2root@DEBUG_SUFFIX@.la
libcobo_frame_graw2root@DEBUG_SUFFIX@_ladir = $(includedir)/get
libcobo_frame_graw2root@DEBUG_SUFFIX@_la_HEADERS = T2KEvent.h
libcobo_frame_graw2root@DEBUG_SUFFIX@_la_SOURCES = T2KEvent.cpp
nodist_libcobo_frame_graw2root@DEBUG_SUFFIX@_la_HEADERS = $(dict_headers)
nodist_libcobo_frame_graw2root@DEBUG_SUFFIX@_la_SOURCES = $(dict_sources)
libcobo_frame_graw2root@DEBUG_SUFFIX@_la_LIBADD = ../root/libget-root@DEBUG_SUFFIX@.la
libcobo_frame_graw2root@DEBUG_SUFFIX@_la_LIBADD += -L$(ROOTLIBDIR) $(ROOTLIBS)

bin_PROGRAMS = cobo-frame-graw2root@DEBUG_SUFFIX@
cobo_frame_graw2root@DEBUG_SUFFIX@_SOURCES = cobo-frame-graw2root.cpp graw2root.h graw2root.cpp
cobo_frame_graw2root@DEBUG_SUFFIX@_LDFLAGS = $(BOOST_LDFLAGS) $(QT_LDFLAGS)
cobo_frame_graw2root@DEBUG_SUFFIX@_LDADD  = libcobo-frame-graw2root@DEBUG_SUFFIX@.la
cobo_frame_graw2root@DEBUG_SUFFIX@_LDADD  += ../root/libget-root@DEBUG_SUFFIX@.la
cobo_frame_graw2root@DEBUG_SUFFIX@_LDADD  += ../get/libcobo-event@DEBUG_SUFFIX@.la
cobo_frame_graw2root@DEBUG_SUFFIX@_LDADD += $(MFM_LIBS)
cobo_frame_graw2root@DEBUG_SUFFIX@_LDADD += -L$(ROOTLIBDIR) $(ROOTLIBS)
cobo_frame_graw2root@DEBUG_SUFFIX@_LDADD += $(BOOST_PROGRAM_OPTIONS_LIB) $(BOOST_FILESYSTEM_LIB) $(BOOST_SYSTEM_LIB)
cobo_frame_graw2root@DEBUG_SUFFIX@_LDADD += $(QT_LIBS) 

endif

BUILT_SOURCES = $(dict_sources) $(dict_headers)
CLEANFILES = $(dict_sources) $(dict_headers)