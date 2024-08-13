# Makefile for Commanding
#
#
# 
#By: Oscar F Romero Matamala

#ifeq($(PREFIX),)
PREFIX := /usr/local
#endif
#SourceDir
SRC_DIR := src
OBJ_DIR := obj
PWD := $(shell pwd)

BIN_DIR := .
EXE := $(BIN_DIR)/CTMessaging

ARCH := $(shell uname)

# linux flags
ifeq ($(ARCH),Linux)
CXX           = gcc
CXXFLAGS      = -Wall -std=c99
LD            = gcc
LDFLAGS       = -g
SOFLAGS       = -shared
endif

# Apple OS X flags
ifeq ($(ARCH),Darwin)
CXX           = gcc 
CXXFLAGS      = -g -O3 -Wall -fPIC  -fno-strict-aliasing
LD            = gcc
LDFLAGS       = -O
SOFLAGS       = -shared
endif

OutPutOpt     = -o




INCLUDEFLAGS	= -Idict -Iinclude -I./inc/ -I./
CXXFLAGS	+= $(INCLUDEFLAGS)

LFLAGS 			= -lcurl -lrt -lboost_system

ALLFLAGS = $(CXXFLAGS)

SRC := $(wildcard $(SRC_DIR)/*.c)
OBJ := $(SRC:$(SRC_DIR)/%.c=$(BIN_DIR)/$(OBJ_DIR)/%.o)
INC := $(SRC:$(SRC_DIR)/%.c=$(BIN_DIR)/include/%.h)
INCH := $(wildcard $(BIN_DIR)/include/*.h)
INCHI := $(INCH:$(BIN_DIR)/include/%.h=$(PREFIX)/include/%.h)

CP  := cp

#Recipe

.PHONY: all

all: $(EXE)

$(EXE): $(OBJ) | $(BIN_DIR) 
	$(LD) $(OutPutOpt) $@ $(OBJ) $(LFLAGS)
	@echo "$@ Compilation done"

$(BIN_DIR)/$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c | $(OBJ_DIR)
	$(LD) -c $<  $(CXXFLAGS) -o $(BIN_DIR)/$@  $(LFLAGS)
$(BIN_DIR)/$(OBJ_DIR):
	mkdir -p $@

install: $(EXE)
	install -d $(PREFIX)/bin
	install -m 777 $(EXE) $(PREFIX)/bin
	install -d ${PREFIX}/.login
	install -m 644 ./.login/.users ${PREFIX}/.login
	install -d $(PREFIX)/include
	install -m 644 $(INCH) $(PREFIX)/include

uninstall: 
	-rm -rv $(PREFIX)/bin/$(EXE)
	-rm -rv $(INCHI)

clean:
	-rm -rv $(BIN_DIR)/$(OBJ_DIR)
	-rm -rv $(EXE)

-include $(OBJ:.o=.d)

.SUFFIXES: .o
