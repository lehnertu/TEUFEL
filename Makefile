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

INCPATH      = -I $(SDDS) -I $(SRCDIR)

CC           = gcc
CXX          = g++
LINK         = g++

CFLAGS       = -O2 -g -Wall
CXXFLAGS     = -O2 -g -Wall -std=c++11
LFLAGS       =
LIBS         =  -L$(SDDS) -lSDDS1 -lmdblib -lmdbcommon -llzma -lz -lm

####### Output directory

OBJECTS_DIR   = ./obj

####### Files

#	$(SRCDIR)/cavity.cpp
#	$(SRCDIR)/wave.cpp
SRC = 	$(SRCDIR)/fields.cpp \
	$(SRCDIR)/homogeneousmagnet.cpp \
	$(SRCDIR)/particle.cpp\
	$(SRCDIR)/teufel.cpp\
	$(SRCDIR)/undulator.cpp \
	$(SRCDIR)/vector.cpp

#	$(OBJDIR)/cavity.o \
#	$(OBJDIR)/wave.o
OBJ = 	$(OBJDIR)/fields.o \
	$(OBJDIR)/homogeneousmagnet.o \
	$(OBJDIR)/particle.o \
	$(OBJDIR)/undulator.o \
	$(OBJDIR)/vector.o

TARGETOBJ = $(OBJDIR)/teufel.o

TARGET = teufel

TESTS = $(TESTDIR)/teufel.magnet \
	$(TESTDIR)/teufel.undulator
	
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
$(OBJ): $(SRCDIR)/global.h \
	$(SRCDIR)/fields.h \
	$(SRCDIR)/homogeneousmagnet.h \
	$(SRCDIR)/particle.h \
	$(SRCDIR)/undulator.h \
	$(SRCDIR)/vector.h

tests: $(OBJ) 
	$(CXX) $(CXXFLAGS) $(INCPATH) -o $(TESTDIR)/teufel.magnet $(TESTDIR)/teufel.magnet.cpp $(LFLAGS) $(OBJ) $(LIBS)
	$(CXX) $(CXXFLAGS) $(INCPATH) -o $(TESTDIR)/teufel.undulator $(TESTDIR)/teufel.undulator.cpp $(LFLAGS) $(OBJ) $(LIBS)

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

