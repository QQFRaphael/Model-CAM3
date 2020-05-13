#CFLAGS = -g -Wall 
#CXXFLAGS = -g -Wall 
CXXFLAGS = -O2 -Wall 
FFLAGS = -g -i4 -r8

LFLAGS = -L/contrib/lib -lnetcdf -lm

INCLUDE_DIRS = -I/contrib/include 
CPPFLAGS = $(INCLUDE_DIRS)

CXX = g++
CC = gcc
FC = f77

SHELL =	/bin/sh

####### Implicit rules

.SUFFIXES: .c

.C.o:
	$(CXX) -c $(CXXFLAGS) $(CPPFLAGS)  $<

.c.o:
	$(CC) -c  $(CFLAGS) $(CPPFLAGS)  $<

.F.o:
	$(FC) -c $(FFLAGS) $(CPPFLAGS)  $<

####### Files

SOURCES = $(wildcard *.C )

OBJECTS = $(SOURCES:.C=.o)

TARGET = nctrans


####### Build rules

all: $(TARGET) 

$(TARGET): $(OBJECTS)
	$(CXX) $(OBJECTS)  -o $(TARGET) $(LFLAGS)

depend:
	@makedepend -v $(SOURCES)  $(CPPFLAGS)  2> /dev/null
clean:
	rm -f *.o  *~  TAGS $(TARGET) 
tags:
	etags *.c *.h
cppout:
	$(CC) -E $(SOURCES) 

# DO NOT DELETE THIS LINE -- make depend depends on it.

