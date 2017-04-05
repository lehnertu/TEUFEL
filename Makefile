# Makefile for TEUFEL
# gcc 4.9.2
# GNU Make version 4.0

SHELL = /bin/sh

SRCDIR = ./src
TESTDIR = ./tests
OBJDIR = ./obj
LIBDIR = ./lib

SDDS = ./lib/SDDSToolKit-devel-3.4

INCPATH      = -I $(SDDS) -I $(SRCDIR)

CC           = gcc -fopenmp
CXX          = g++ -fopenmp
LINK         = g++ -fopenmp

CFLAGS       = -O2 -g -Wall
CXXFLAGS     = -std=c++11
LFLAGS       =
LIBS         =  -L$(SDDS) -lSDDS1 -lmdblib -lmdbcommon -llzma -lz -lm

####### Output directory

OBJECTS_DIR   = ./obj

####### Files

SRC = 	$(SRCDIR)/bunch.cpp \
	$(SRCDIR)/cavity.cpp \
	$(SRCDIR)/externalfield.cpp \
	$(SRCDIR)/gen_grid.cpp \
	$(SRCDIR)/homogeneousmagnet.cpp \
	$(SRCDIR)/homogeneouselectricfield.cpp \
	$(SRCDIR)/particle.cpp\
	$(SRCDIR)/teufel.cpp\
	$(SRCDIR)/undulator.cpp \
	$(SRCDIR)/vector.cpp \
	$(SRCDIR)/wave.cpp

OBJ =	$(OBJDIR)/bunch.o \
	$(OBJDIR)/cavity.o \
	$(OBJDIR)/externalfield.o \
	$(OBJDIR)/gen_grid.o \
	$(OBJDIR)/homogeneousmagnet.o \
	$(OBJDIR)/homogeneouselectricfield.o \
	$(OBJDIR)/particle.o \
	$(OBJDIR)/undulator.o \
	$(OBJDIR)/vector.o \
	$(OBJDIR)/wave.o

TARGETOBJ = $(OBJDIR)/teufel.o

TARGET = teufel

TESTS = $(TESTDIR)/teufel.EcrossB \
	$(TESTDIR)/teufel.magnet \
	$(TESTDIR)/teufel.undulator \
	$(TESTDIR)/teufel.radiation \
	$(TESTDIR)/trial \
	$(TESTDIR)/teufel.wave-particle
	

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

.PHONY : first all tests clean

first: $(TARGET)

all: $(TARGET) tests docs

$(TARGET):  $(OBJ)  $(TARGETOBJ)  
	$(LINK) $(LFLAGS) -o $(TARGET) $(OBJ) $(TARGETOBJ) $(LIBS)

$(OBJ): $(SRCDIR)/global.h \
	$(SRCDIR)/bunch.h \
	$(SRCDIR)/cavity.h \
	$(SRCDIR)/externalfield.h \
	$(SRCDIR)/gen_grid.h \
	$(SRCDIR)/homogeneouselectricfield.h \
	$(SRCDIR)/homogeneousmagnet.h \
	$(SRCDIR)/particle.h \
	$(SRCDIR)/undulator.h \
	$(SRCDIR)/vector.h \
	$(SRCDIR)/wave.h

tests: $(OBJ)
	$(CXX) $(CXXFLAGS) $(INCPATH) -o $(TESTDIR)/teufel.EcrossB $(TESTDIR)/teufel.EcrossB.cpp $(LFLAGS) $(OBJ) $(LIBS) 
	$(CXX) $(CXXFLAGS) $(INCPATH) -o $(TESTDIR)/teufel.magnet $(TESTDIR)/teufel.magnet.cpp $(LFLAGS) $(OBJ) $(LIBS)
	$(CXX) $(CXXFLAGS) $(INCPATH) -o $(TESTDIR)/teufel.undulator $(TESTDIR)/teufel.undulator.cpp $(LFLAGS) $(OBJ) $(LIBS)
	$(CXX) $(CXXFLAGS) $(INCPATH) -o $(TESTDIR)/teufel.radiation $(TESTDIR)/teufel.radiation.cpp $(LFLAGS) $(OBJ) $(LIBS)
	$(CXX) $(CXXFLAGS) $(INCPATH) -o $(TESTDIR)/trial $(TESTDIR)/trial.cpp $(LFLAGS) $(OBJ) $(LIBS)
	$(CXX) $(CXXFLAGS) $(INCPATH) -o $(TESTDIR)/teufel.wave-particle $(TESTDIR)/teufel.wave-particle.cpp $(LFLAGS) $(OBJ) $(LIBS)
	

docs:
	doxygen setup.dox

clean:
	-rm $(OBJ)
	-rm $(TARGETOBJ)
	-rm $(TARGET)
	-rm $(TESTS)
	-rm -rf doc

