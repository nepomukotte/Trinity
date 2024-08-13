include $(top_srcdir)/am/cppcheck.am

AUTOMAKE_OPTIONS = subdir-objects

headers = dhsm/Action.hpp \
	dhsm/CompositeEngine.h \
	dhsm/Engine.h \
	dhsm/Event.h \
	dhsm/Exception.h \
	dhsm/Message.h \
	dhsm/State.h \
	dhsm/StateData.h \
	dhsm/StateMachine.h \
	dhsm/Transition.h

sources = dhsm/Action.cpp \
	dhsm/CompositeEngine.cpp \
	dhsm/Engine.cpp \
	dhsm/Event.cpp \
	dhsm/Exception.cpp \
	dhsm/State.cpp \
	dhsm/StateMachine.cpp \
	dhsm/Transition.cpp

dhsm_servant_headers = dhsm/EngineServant.h \
	dhsm/ServantApplication.h

dhsm_servant_sources = dhsm/EngineServant.cpp \
	dhsm/ServantApplication.cpp

dhsm_master_headers = dhsm/MasterEngine.h \
	dhsm/Semaphore.h \
	dhsm/ServantCallback.h \
	dhsm/ServantDescriptor.h

dhsm_master_sources = dhsm/MasterEngine.cpp \
	dhsm/Semaphore.cpp \
	dhsm/ServantCallback.cpp \
	dhsm/ServantDescriptor.cpp

ice_sources = dhsm/IEngineServant.ice

ice_gen_sources = dhsm/IEngineServant.cpp
ice_gen_headers = dhsm/IEngineServant.h

EXTRA_DIST = $(ice_sources)

dhsm_includedir = $(includedir)
nobase_dhsm_include_HEADERS = $(headers) $(dhsm_servant_headers) $(dhsm_master_headers) $(ice_gen_headers)

AM_CXXFLAGS = $(CXXFLAGS_DHSM) $(UTILS_CFLAGS)

if ICE_ANY
BUILT_SOURCES  = $(ice_gen_sources) $(ice_gen_headers)
CLEANFILES     = $(BUILT_SOURCES)
DISTCLEANFILES = $(CLEANFILES)

DISTRIBUTED_CPPFLAGS = $(ICE_CPPFLAGS) -I@builddir@/dhsm
DISTRIBUTED_CXXFLAGS = $(PTHREAD_CFLAGS)
DISTRIBUTED_LDFLAGS  = $(ICE_LDFLAGS)
DISTRIBUTED_LIBADD   = $(ICE_LIBS) $(PTHREAD_LIBS)

AM_CPPFLAGS = $(DISTRIBUTED_CPPFLAGS)
AM_CXXFLAGS += $(DISTRIBUTED_CXXFLAGS)
endif

lib_LTLIBRARIES	= libdhsm@DEBUG_SUFFIX@.la
libdhsm@DEBUG_SUFFIX@_la_SOURCES = $(sources)
libdhsm@DEBUG_SUFFIX@_la_LDFLAGS = -version-info 0:0:0 -no-undefined
libdhsm@DEBUG_SUFFIX@_la_LIBADD  = $(UTILS_LIBS)

if SERVANT
lib_LTLIBRARIES += libdhsm-servant@DEBUG_SUFFIX@.la
libdhsm_servant@DEBUG_SUFFIX@_la_SOURCES        = $(dhsm_servant_sources)
nodist_libdhsm_servant@DEBUG_SUFFIX@_la_SOURCES = $(ice_gen_sources)
libdhsm_servant@DEBUG_SUFFIX@_la_LDFLAGS        = -version-info 0:0:0 -no-undefined $(DISTRIBUTED_LDFLAGS)
libdhsm_servant@DEBUG_SUFFIX@_la_LIBADD         = libdhsm@DEBUG_SUFFIX@.la $(UTILS_LIBS) $(DISTRIBUTED_LIBADD)
endif

if MASTER
lib_LTLIBRARIES += libdhsm-master@DEBUG_SUFFIX@.la
libdhsm_master@DEBUG_SUFFIX@_la_SOURCES = $(dhsm_master_sources)
libdhsm_master@DEBUG_SUFFIX@_la_LDFLAGS = -version-info 0:0:0 -no-undefined $(DISTRIBUTED_LDFLAGS)
libdhsm_master@DEBUG_SUFFIX@_la_LIBADD  = libdhsm-servant@DEBUG_SUFFIX@.la $(DISTRIBUTED_LIBADD)
endif

include $(top_srcdir)/am/ice.am
