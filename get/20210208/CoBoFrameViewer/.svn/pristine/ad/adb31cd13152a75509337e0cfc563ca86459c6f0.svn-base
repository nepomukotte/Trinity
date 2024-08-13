
AM_CPPFLAGS = -I$(ROOTINCDIR)

dict_source = 
dict_header = 

if BUILD_GRAW2ROOT

dict_source += DictGET.cxx
if !HAVE_ROOT6
dict_header += DictGET.h
endif
if HAVE_ROOT6
dict_header += DictGET_rdict.pcm
endif

headers = GDataSample.h GDataChannel.h GFrameHeader.h GDataFrame.h LinkDefGET.h
sources = GDataSample.cpp GDataChannel.cpp GFrameHeader.cpp GDataFrame.cpp

# Installation of header files
get_root_includedir  = $(includedir)/get
get_root_include_HEADERS = $(headers)
nodist_get_root_include_HEADERS = $(dict_header)

# Library
lib_LTLIBRARIES = libget-root@DEBUG_SUFFIX@.la
libget_root@DEBUG_SUFFIX@_la_SOURCES = $(sources)
nodist_libget_root@DEBUG_SUFFIX@_la_SOURCES = $(dict_source)

# Support for ROOT silent build rules
rcint_verbose = $(rcint_verbose_$(V))
rcint_verbose_ = $(rcint_verbose_$(AM_DEFAULT_VERBOSITY))
rcint_verbose_0 = @echo "RCINT    $@";

# ROOT built sources
DictGET.cxx : $(headers) LinkDefGET.h
	$(rcint_verbose)$(ROOTCINT) -f $@ -c -I@srcdir@ -p $^

DictGET.h : DictGET.cxx

endif

BUILT_SOURCES = $(dict_source) $(dict_header)
CLEANFILES = $(dict_source) $(dict_header)