# Makefile for EventBuilder
#
#
# 
#By: Oscar F Romero Matamala

#SourceDir
SRC_DIR := src
OBJ_DIR := obj
DICT_DIR := ${EXACT_DIR}/dict
BIN_DIR := .
EXE := $(BIN_DIR)/EventBuilder
GET_DIR := /usr/share
ARCH := $(shell uname)

# linux flags
ifeq ($(ARCH),Linux)
CXX           = g++ 
CXXFLAGS      =  -g -O3 -w -ggdb3 -fPIC -fno-strict-aliasing  -D_FILE_OFFSET_BITS=64 -D_LARGE_FILE_SOURCE -D_LARGEFILE64_SOURCE
LD            = g++
LDFLAGS       = -g
SOFLAGS       = -shared
endif

# Apple OS X flags
ifeq ($(ARCH),Darwin)
CXX           = g++ 
CXXFLAGS      = -g -w -O3 -fPIC  -fno-strict-aliasing
LD            = g++
LDFLAGS       = -O
SOFLAGS       = -shared
endif

OutPutOpt     = -o

# Root

ROOTCFLAGS   := $(shell root-config --cflags)
ROOTLIBS     := $(shell root-config --libs)
ROOTGLIBS    := $(shell root-config --glibs)
ROOTLDFLAGS  := $(shell root-config --ldflags)
ROOTLIBDIR   := $(shell root-config --libdir)

#GETLIBS = /usr/share/lib	

CLLFLAGS += $(ROOTLDFLAGS)

CXXFLAGS     += $(ROOTCFLAGS)
LIBS          = $(ROOTGLIBS)
LIBS         += -L$(GET_DIR)/lib
LIBS 		+= -L$(DICT_DIR)
LIBS         += -lMinuit -lget-root -lMultiFrame -lExACT
#LIBEXACT := $(DICT_DIR)/libExACT.so


INCLUDEFLAGS  = -Iinclude -I./inc/ -I/usr/local/include/get/ -I$(GET_DIR)/include/ -Isrc -I$(DICT_DIR)/include
CXXFLAGS     += $(INCLUDEFLAGS)
#CXXFLAGS += $(GETLIBS)

ALLFLAGS = $(CXXFLAGS) $(CPPFLAGS) -Wall

SRC := $(wildcard $(SRC_DIR)/*.cpp)
OBJ := $(SRC:$(SRC_DIR)/%.cpp=$(OBJ_DIR)/%.o)
#DICT := $(DICT_DIR)/exactDict.C
RMAP := libExACT.rootmap
#LIB := libExACT.so

CP  := cp

#Recipe

.PHONY: all
all: $(EXE)

$(EXE): $(OBJ) | $(BIN_DIR) 
	$(LD)  $(CXXFLAGS) $^ $(LIBS) $(LDFLAGS) $(OutPutOpt) $@ -D__STDC_LIMIT_MACROS -D__STDC_CONSTANT_MACROS 
	@echo "$@ done"

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp | $(OBJ_DIR)
	$(LD) -c $(CXXFLAGS) ${LIBS} $< -o $@ 

$(BIN_DIR)/$(OBJ_DIR):
		mkdir -p $@

clean:
	rm -rv $(BIN_DIR)/$(OBJ_DIR)
	rm -rv $(EXE)

-include $(OBJ:.o=.d)

.SUFFIXES: .o
