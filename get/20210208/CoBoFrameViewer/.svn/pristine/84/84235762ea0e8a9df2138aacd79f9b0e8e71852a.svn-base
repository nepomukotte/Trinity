
AM_CPPFLAGS = -I@top_srcdir@/src -I@top_srcdir@/src/root
AM_CPPFLAGS += $(BOOST_CPPFLAGS) $(ROOTCFLAGS) $(UTILS_CFLAGS)

dict_source = 
dict_header = 

bin_PROGRAMS =

if BUILD_GRAW2ROOT

dict_source += DictGuiGET.cxx
dict_header += DictGuiGET.h
if !HAVE_ROOT6
dict_header += DictGuiGET.h
endif
if HAVE_ROOT6
dict_header += DictGuiGET_rdict.pcm
endif
headers = CoBoFrameBrowser.h
sources = CoBoFrameBrowser.cpp

bin_PROGRAMS += cobo-frame-root2disp@DEBUG_SUFFIX@
cobo_frame_root2disp@DEBUG_SUFFIX@_SOURCES = cobo-frame-root2disp.cpp $(headers) $(sources) LinkDefGuiGET.h
nodist_cobo_frame_root2disp@DEBUG_SUFFIX@_SOURCES = $(dict_header) $(dict_source)
cobo_frame_root2disp@DEBUG_SUFFIX@_LDFLAGS = $(BOOST_LDFLAGS)
cobo_frame_root2disp@DEBUG_SUFFIX@_LDADD  = ../root/libget-root@DEBUG_SUFFIX@.la
cobo_frame_root2disp@DEBUG_SUFFIX@_LDADD += $(UTILS_LIBS)
cobo_frame_root2disp@DEBUG_SUFFIX@_LDADD += -L$(ROOTLIBDIR) $(ROOTGLIBS)
cobo_frame_root2disp@DEBUG_SUFFIX@_LDADD += $(BOOST_PROGRAM_OPTIONS_LIB) $(BOOST_SYSTEM_LIB)

# Support for ROOT silent build rules
rcint_verbose = $(rcint_verbose_$(V))
rcint_verbose_ = $(rcint_verbose_$(AM_DEFAULT_VERBOSITY))
rcint_verbose_0 = @echo "RCINT    $@";

# ROOT built sources
DictGuiGET.cxx : $(headers) LinkDefGuiGET.h
	$(rcint_verbose)$(ROOTCINT) -f $@ -c -I@top_srcdir@/src/root -I@srcdir@ -p $^

DictGuiGET.h : DictGuiGET.cxx

endif

BUILT_SOURCES = $(dict_source) $(dict_header)
CLEANFILES =  $(dict_source) $(dict_header)