include $(top_srcdir)/am/cppcheck.am

AUTOMAKE_OPTIONS = subdir-objects

headers = dhsm/Camera.h \
	dhsm/Chrono.h \
	dhsm/CruiseControl.h

sources = dhsm/Camera.cpp \
	dhsm/Chrono.cpp \
	dhsm/CruiseControl.cpp \
	SmTest.cpp

master_sources = dhsm/RunControlMaster.cpp \
    dhsm/RunControlStateMachine.cpp

servant_sources = dhsm/RunControlServant.cpp \
    dhsm/RunControlStateMachine.cpp \
    dhsm/RunServantImpl.cpp

dhsm_headers = dhsm/RunControlBehaviour.h \
	dhsm/RunControlMaster.h \
	dhsm/RunControlServant.h \
	dhsm/RunControlStateMachine.h \
	dhsm/RunServantImpl.h


ice_sources = dhsm/RunServant.ice

ice_gen_sources = dhsm/RunServant.cpp
ice_gen_headers = dhsm/RunServant.h

EXTRA_DIST = $(headers) $(dhsm_headers) $(ice_sources)

if ICE_ANY
BUILT_SOURCES  = $(ice_gen_sources) $(ice_gen_headers)
CLEANFILES     = $(BUILT_SOURCES)
DISTCLEANFILES = $(CLEANFILES)

DISTRIBUTED_CPPFLAGS = $(ICE_CPPFLAGS) -I@builddir@/dhsm
DISTRIBUTED_CXXFLAGS = $(PTHREAD_CFLAGS)
else
BUILT_SOURCES  =
CLEANFILES     =
DISTCLEANFILES =
endif

AM_CPPFLAGS = -I@srcdir@/../src $(UTILS_CFLAGS) $(DISTRIBUTED_CPPFLAGS)
AM_CXXFLAGS = $(CXXFLAGS_DHSM) $(DISTRIBUTED_CXXFLAGS)

#TESTS = test_dhsm.sh
check_SCRIPTS = test_dhsm.sh

noinst_PROGRAMS = StateMachineTest
StateMachineTest_SOURCES = $(sources)
StateMachineTest_LDADD   = @builddir@/../src/libdhsm@DEBUG_SUFFIX@.la $(UTILS_LIBS)

if SERVANT
noinst_PROGRAMS += RunControlServant
RunControlServant_SOURCES        = $(servant_sources)
nodist_RunControlServant_SOURCES = $(ice_gen_sources)
RunControlServant_LDFLAGS        = $(ICE_LDFLAGS)
RunControlServant_LDADD          = @builddir@/../src/libdhsm-servant@DEBUG_SUFFIX@.la
RunControlServant_LDADD         += @builddir@/../src/libdhsm@DEBUG_SUFFIX@.la
RunControlServant_LDADD         += $(UTILS_LIBS)
RunControlServant_LDADD         += $(ICE_LIBS) $(PTHREAD_LIBS)
endif

if MASTER
noinst_PROGRAMS += RunControlMaster
RunControlMaster_SOURCES        = $(master_sources)
nodist_RunControlMaster_SOURCES = $(ice_gen_sources)
RunControlMaster_LDFLAGS        = $(ICE_LDFLAGS)
RunControlMaster_LDADD          = @builddir@/../src/libdhsm-master@DEBUG_SUFFIX@.la
RunControlMaster_LDADD         += @builddir@/../src/libdhsm@DEBUG_SUFFIX@.la
RunControlMaster_LDADD         += $(UTILS_LIBS)
RunControlMaster_LDADD         += $(ICE_LIBS) $(PTHREAD_LIBS)
endif

dist-hook:
	mkdir $(distdir)/ref
	cp -p $(srcdir)/ref/*.txt $(distdir)/ref/

distclean-local:
	rm -f $(builddir)/*.txt
	rm -f $(builddir)/*.log

include $(top_srcdir)/am/ice.am
include $(srcdir)/Makefile.vxworks.am
