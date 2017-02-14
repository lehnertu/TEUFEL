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

CC           = gcc
CXX          = g++
LINK         = g++

CFLAGS       = -O2 -g -Wall
CXXFLAGS     = -O2 -g -Wall
LFLAGS       =
LIBS         =  -L$(SDDS) -lSDDS1 -lmdblib -lmdbcommon -llzma -lz -lm

####### Output directory

OBJECTS_DIR   = ./obj

####### Files

SRC = 	$(SRCDIR)/cavity.cpp \
	$(SRCDIR)/externalfield.cpp \
	$(SRCDIR)/homogeneousmagnet.cpp \
	$(SRCDIR)/particle.cpp\
	$(SRCDIR)/teufel.cpp\
	$(SRCDIR)/undulator.cpp \
	$(SRCDIR)/vector.cpp \
	$(SRCDIR)/wave.cpp

OBJ = 	$(OBJDIR)/cavity.o \
	$(OBJDIR)/externalfield.o \
	$(OBJDIR)/homogeneousmagnet.o \
	$(OBJDIR)/particle.o \
	$(OBJDIR)/undulator.o \
	$(OBJDIR)/vector.o \
	$(OBJDIR)/wave.o

TARGETOBJ = $(OBJDIR)/teufel.o

TARGET = teufel

TESTS = $(TESTDIR)/teufel.magnet \
	$(TESTDIR)/teufel.undulator
	

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
	$(SRCDIR)/cavity.h \
	$(SRCDIR)/externalfield.h \
	$(SRCDIR)/homogeneousmagnet.h \
	$(SRCDIR)/particle.h \
	$(SRCDIR)/undulator.h \
	$(SRCDIR)/vector.h \
	$(SRCDIR)/wave.h

tests: $(OBJ) 
	$(CXX) $(CXXFLAGS) $(INCPATH) -o $(TESTDIR)/teufel.magnet $(TESTDIR)/teufel.magnet.cpp $(LFLAGS) $(OBJ) $(LIBS)
	$(CXX) $(CXXFLAGS) $(INCPATH) -o $(TESTDIR)/teufel.undulator $(TESTDIR)/teufel.undulator.cpp $(LFLAGS) $(OBJ) $(LIBS)

docs:
	doxygen setup.dox

clean:
	-rm $(OBJ)
	-rm $(TARGETOBJ)
	-rm $(TARGET)
	-rm $(TESTS)
	-rm -rf doc

