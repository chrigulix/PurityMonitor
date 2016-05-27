# $Id: GNUmakefile,v 1.3 2009-05-11 13:03:17 gcosmo Exp $
# --------------------------------------------------------------
# GNUmakefile for examples module.
# --------------------------------------------------------------
CC = g++
PROJECT = PurityMonitor

SOURCES = $(PROJECT).cpp
OBJECTS = $(SOURCES:.cpp=.o)

CPPFLAGS = -std=c++0x -Wall -g -O4 
# Root use
# CPPFLAGS += $(shell $(ROOTSYS)/bin/root-config --cflags)
ROOTCFLAGS   = $(shell $(ROOTSYS)/bin/root-config --cflags)
ROOTLIBS     = $(shell $(ROOTSYS)/bin/root-config --libs)
ROOTGLIBS    = $(shell $(ROOTSYS)/bin/root-config --glibs)
SOFLAGS      = -shared
CPPFLAGS     += $(ROOTCFLAGS)
# CPPFLAGS     += -I/usr/include/eigen3
LIBS         += $(ROOTLIBS) -L. -lThread -L/usr/lib64
GLIBS        += $(ROOTGLIBS)
EXTRALIBS    += $(ROOTLIBS)

# LIBS         += -pthread

all: 
	$(CC) $(SOURCES) $(CPPFLAGS) $(LIBS) -o $(PROJECT)
# $(PROJECT): $(CC) -0 $(PROJECT)

# set up extra flags for explicitly setting mode
# debug:      CXXFLAGS    += -g
# release:    CXXFLAGS    += -O3

.PHONY: clean
clean: 
	rm -rf *.o $(PROJECT)