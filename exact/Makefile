
# Makefile for ExACT
#
#
# 
#By: Oscar F Romero Matamala

#SourceDir
SRC_DIR := src
OBJ_DIR := obj
PWD := $(shell pwd)
DICT_DIR := dict
DICT_SRC_DIR := $(DICT_DIR)/src
DICT_OBJ_DIR := $(DICT_DIR)/obj

BIN_DIR := .
EXE := $(BIN_DIR)/ExACT

ARCH := $(shell uname)

# linux flags
ifeq ($(ARCH),Linux)
CXX           = g++ 
CXXFLAGS      =  -g -O3 -Wall -ggdb3 -fPIC -fno-strict-aliasing  -D_FILE_OFFSET_BITS=64 -D_LARGE_FILE_SOURCE -D_LARGEFILE64_SOURCE
LD            = g++
LDFLAGS       = -g
SOFLAGS       = -shared
endif

# Apple OS X flags
ifeq ($(ARCH),Darwin)
CXX           = g++ 
CXXFLAGS      = -g -O3 -Wall -fPIC  -fno-strict-aliasing
LD            = g++
LDFLAGS       = -O
SOFLAGS       = -shared
endif

OutPutOpt     = -o

AR			= ar
PREFIX		= /usr/local
# Root

ROOTCFLAGS   := $(shell root-config --cflags)
ROOTLIBS     := $(shell root-config --libs)
ROOTGLIBS    := $(shell root-config --glibs)
ROOTLDFLAGS  := $(shell root-config --ldflags)
ROOTLIBDIR   := $(shell root-config --libdir)

CLLFLAGS += $(ROOTLDFLAGS)

CXXFLAGS     += $(ROOTCFLAGS)
LIBS          = $(ROOTGLIBS)
LIBS         += -lMinuit -L./ -L./include -L./dict


INCLUDEFLAGS  = -Idict -Iinclude -I./inc/ -I$(PWD)/dict/include -I./dict/cling -I./
CXXFLAGS     += $(INCLUDEFLAGS)

ALLFLAGS = $(CXXFLAGS) $(CPPFLAGS) -Wall

SRC := $(wildcard $(SRC_DIR)/*.cpp)
OBJ := $(SRC:$(SRC_DIR)/%.cpp=$(OBJ_DIR)/%.o)
INC := $(SRC:$(SRC_DIR)/%.cpp=$(BIN_DIR)/include/%.h)

LIB_INC := $(wildcard $(BIN_DIR)/include/*.h)
LIB_OBJ := $(LIB_INC:$(BIN_DIR)/include/%.h=$(OBJ_DIR)/%.o)

LIB_EX	:= $(BIN_DIR)/include/libIExACT.a

DICT := $(DICT_DIR)/exactDict.cpp
DICT_SRC := $(wildcard $(DICT_SRC_DIR)/*.cpp)
DICT_OBJ := $(DICT_SRC:$(DICT_SRC_DIR)/%.cpp=$(DICT_OBJ_DIR)/%.o)
DICT_INC := $(DICT_SRC:$(DICT_SRC_DIR)/%.cpp=$(DICT_DIR)/include/%.h)
DICT_LINK := $(DICT_DIR)/cling/Linkdef.h
PCM := $(DICT_DIR)/exactDict_rdict.pcm
RMAP := $(DICT_DIR)/libExACT.rootmap
LIB := $(DICT_DIR)/libExACT.so

CP  := cp

#Recipe

.PHONY: all
all: $(EXE)

$(EXE): $(DICT) $(LIB) $(OBJ) $(LIB_EX)| $(BIN_DIR) 
	$(LD) $(OBJ) $(CXXFLAGS) $(LIBS) $(LDFLAGS) $(OutPutOpt) $@ -D__STDC_LIMIT_MACROS -D__STDC_CONSTANT_MACROS -lExACT
	@echo "$@ Compilation done"

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp | $(OBJ_DIR)
	$(LD) $(CXXFLAGS) -c $< -o $@

$(DICT_OBJ_DIR)/%.o: $(DICT_SRC_DIR)/%.cpp | $(DICT_OBJ_DIR)
	$(LD) $(CXXFLAGS) -c $< -o $@

$(DICT): $(DICT_INC) $(DICT_LINK) | $(DICT_DIR)
	rootcling -f $@ -rmf $(RMAP) -c $(INCLUDEFLAGS) -p $(notdir $^)

$(LIB): $(DICT) $(DICT_OBJ)
	$(LD) $(OutPutOpt) $(LIB) $(SOFLAGS) $(CXXFLAGS) $(INCLUDEFLAGS) $(LIBS) $^ $(LDFLAGS) $(CLLFLAGS) 
#	cp $(LIB) $(ROOTLIBDIR)
#	cp $(PCM) $(ROOTLIBDIR)
#	cp $(RMAP) $(ROOTLIBDIR)

$(LIB_EX): $(LIB_OBJ)
	$(AR) -rcs $(LIB_EX) $(LIB_OBJ)


$(BIN_DIR)/$(OBJ_DIR):
		mkdir -p $@

$(BIN_DIR)/$(DICT_OBJ_DIR):
		mkdir -p $@

install: $(LIB_EX)
	install -d $(PREFIX)/lib
	install -m 777 $(LIB_EX) $(PREFIX)/lib
	install -d ${PREFIX}/include
	install -m 644 $(LIB_INC) $(PREFIX)/include

clean:
	-rm -rv $(BIN_DIR)/$(OBJ_DIR)
	-rm -rv $(BIN_DIR)/$(DICT_OBJ_DIR)
	-rm -rv $(BIN_DIR)/$(LIB_EX)
	-rm -rv $(EXE)
	-rm -rv $(LIB)
	-rm -rv $(DICT)
	-rm -rv $(PCM)
	-rm -rv $(RMAP)
-include $(OBJ:.o=.d)

.SUFFIXES: .o
