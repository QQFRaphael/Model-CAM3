CFLAGS = -O2 -Wall
CXXFLAGS = -O2 -Wall
FFLAGS = -O2 -i4 -r8

CFLAGS = -g -Wall
CXXFLAGS = -g -Wall
FFLAGS = -g -i4 -r8

LFLAGS = -L /contrib/gcc-2.95.2/lib/gcc-lib/sparc-sun-solaris2.6/2.95.2 -lgcc -L/contrib/lib -lnetcdf -lm 
LFLAGS = -L /contrib-jt/lib -lnetcdf -lm 

INCLUDE_DIRS = -I/contrib-jt/include -I../../scm_init
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

TARGET = ccm2iop

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

############################################################################
# options to pass to $CC to generate dependencies 
CC_DEP_OPTS = -MM

#this rule comes from the Gnu Make info pages
#don't regenerate the dependencies unless the target is "all"
ifeq ($(MAKECMDGOALS),all)
%.d: %.C
	@ echo "generating dependencies for $<"
	@ $(SHELL) -ec '$(CXX) $(CC_DEP_OPTS) $(CPPFLAGS) $< \
               | sed '\''s/\($*\)\.o[ :]*/\1.o $@ : /g'\'' > $@; \
               [ -s $@ ] || rm -f $@'

%.d: %.c
	@ echo "generating dependencies for $<"
	@ $(SHELL) -ec '$(CC) $(CC_DEP_OPTS) $(CPPFLAGS) $< \
               | sed '\''s/\($*\)\.o[ :]*/\1.o $@ : /g'\'' > $@; \
               [ -s $@ ] || rm -f $@'

endif

# include the automatically generated dependencies
ifneq ($(MAKECMDGOALS),clean)
-include $(SOURCES:.C=.d)
endif





