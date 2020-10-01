CCC = g++       # compiler
AR = ar         # archiver
LD = g++        # linker


ARFLAGS = rcs

LDFLAGS += -Wl,--no-as-needed -fPIC -ggdb3 -O3 -fno-omit-frame-pointer
LDFLAGS += $(shell root-config --libs)

CXXFLAGS += -c -fPIC  -ggdb3 -O3 
CXXFLAGS += $(shell root-config --cflags)

SRCFILES = $(wildcard *.cc)
OBJECTS = $(patsubst %.cc, %.o, ${SRCFILES})


ALLFLAGS = $(CXXFLAGS) -Wall

OutPutOpt     = -o



# rule for any compiling any .cpp file
%.o : %.cpp
	@printf "Compiling $< ... "
	@g++ $(ALLFLAGS) -c $<
	@echo "Done"


EXE = TrinityPerformanceCalculations

all: ${EXE}

TrinityPerformanceCalculations: TrinityPerformanceCalculations.o
	$(LD)  $(LDFLAGS) $^ $(OutPutOpt) $@

clean:
	@rm -f *.o *.a *.so *~ core $(EXE)

