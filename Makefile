# Makefile for TEUFEL
# gcc 4.9.2
# GNU Make version 4.0

SHELL = /bin/sh

SRCDIR = ./src
TESTDIR = ./tests
EXPLDIR = ./examples
OBJDIR = ./obj
LIBDIR = ./lib

SDDS = ./lib/SDDSToolKit-devel-3.4

HDF5 = $(HOME)/lib/hdf5-1.8.19
HDF5INC= $(HDF5)/include
HDF5LIB= $(HDF5)/lib

INCPATH      = -I $(SDDS) -I $(SRCDIR) -I $(HDF5INC)

CC           = gcc
CXX          = g++
LINK         = g++

CFLAGS       = -O2 -g -Wall
CXXFLAGS     = -O2 -g -Wall -std=c++11
LFLAGS       = -Wl,-rpath,$(HDF5LIB)
LIBS         = -L$(SDDS) -lSDDS1 -L$(HDF5LIB) -lhdf5 -lpugixml -lmdblib -lmdbcommon -llzma -lz -lm

####### Output directory

OBJECTS_DIR   = ./obj

####### Files

#	$(SRCDIR)/cavity.cpp
#	$(SRCDIR)/wave.cpp
SRC = 	$(SRCDIR)/bunch.cpp \
	$(SRCDIR)/fields.cpp \
	$(SRCDIR)/observer.cpp\
	$(SRCDIR)/particle.cpp\
	$(SRCDIR)/simulation.cpp\
	$(SRCDIR)/teufel.cpp\
	$(SRCDIR)/undulator.cpp \
	$(SRCDIR)/vector.cpp

#	$(OBJDIR)/cavity.o \
#	$(OBJDIR)/wave.o
OBJ = 	$(OBJDIR)/bunch.o \
	$(OBJDIR)/fields.o \
	$(OBJDIR)/observer.o \
	$(OBJDIR)/particle.o \
	$(OBJDIR)/simulation.o \
	$(OBJDIR)/undulator.o \
	$(OBJDIR)/vector.o

TARGETOBJ = $(OBJDIR)/teufel.o

TARGET = teufel

TESTS = $(TESTDIR)/teufel.bunch \
	$(TESTDIR)/teufel.magnet \
	$(TESTDIR)/teufel.undulator \
	$(TESTDIR)/teufel.loop \
	$(TESTDIR)/teufel.EcrossB
	
EXAMPLES = $(EXPLDIR)/elbe-u300

####### Implicit rules

.SUFFIXES: .o .c .cpp

.cpp.o:
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o "$@" "$<"

.c.o:
	$(CC) -c $(CFLAGS) $(INCPATH) -o "$@" "$<"

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	@mkdir -p $(@D)
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o $@ $<

####### Build rules

.PHONY : first all tests examples clean

first: $(TARGET)

all: $(TARGET) tests docs

$(TARGET):  $(OBJ)  $(TARGETOBJ)  
	$(LINK) $(LFLAGS) -o $(TARGET) $(OBJ) $(TARGETOBJ) $(LIBS)

#	$(SRCDIR)/cavity.h \
#	$(SRCDIR)/wave.h
$(OBJ): $(SRCDIR)/bunch.h \
	$(SRCDIR)/global.h \
	$(SRCDIR)/fields.h \
	$(SRCDIR)/observer.h \
	$(SRCDIR)/particle.h \
	$(SRCDIR)/simulation.h \
	$(SRCDIR)/undulator.h \
	$(SRCDIR)/vector.h

tests: $(OBJ) 
	$(CXX) $(CXXFLAGS) $(INCPATH) -o $(TESTDIR)/teufel.bunch $(TESTDIR)/teufel.bunch.cpp $(LFLAGS) $(OBJ) $(LIBS)
	$(CXX) $(CXXFLAGS) $(INCPATH) -o $(TESTDIR)/teufel.magnet $(TESTDIR)/teufel.magnet.cpp $(LFLAGS) $(OBJ) $(LIBS)
	$(CXX) $(CXXFLAGS) $(INCPATH) -o $(TESTDIR)/teufel.undulator $(TESTDIR)/teufel.undulator.cpp $(LFLAGS) $(OBJ) $(LIBS)
	$(CXX) $(CXXFLAGS) $(INCPATH) -o $(TESTDIR)/teufel.loop $(TESTDIR)/teufel.loop.cpp $(LFLAGS) $(OBJ) $(LIBS)
	$(CXX) $(CXXFLAGS) $(INCPATH) -o $(TESTDIR)/teufel.EcrossB $(TESTDIR)/teufel.EcrossB.cpp $(LFLAGS) $(OBJ) $(LIBS)

examples: $(OBJ) 
	$(CXX) $(CXXFLAGS) $(INCPATH) -o $(EXPLDIR)/elbe-u300 $(EXPLDIR)/elbe-u300.cpp $(LFLAGS) $(OBJ) $(LIBS)

docs:
	doxygen setup.dox

clean:
	-rm $(OBJ)
	-rm $(TARGETOBJ)
	-rm $(TARGET)
	-rm $(TESTS)
	-rm -rf doc

