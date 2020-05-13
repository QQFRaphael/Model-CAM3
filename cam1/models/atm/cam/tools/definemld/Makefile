CC = cc
CFLAGS = -O -I$(INC_NETCDF) -DDEBUG
LDFLAGS = -L$(LIB_NETCDF) -lnetcdf

OBJS1 = definemld1x1.o err_exit.o wrap_nc.o
EXE1 = definemld1x1

OBJS2 = definemldbdy.o binf2c.o neighborfill.o wrap_nf.o sm121.o precision.o
EXE2 = definemldbdy

ifeq ($(EXEDIR),$(null))
EXEDIR = .
endif

ifeq ($(LIB_NETCDF),$(null))
LIB_NETCDF := /usr/local/lib
endif

ifeq ($(INC_NETCDF),$(null))
INC_NETCDF := /usr/local/include
endif

UNAMES := $(shell uname -s)

ifeq ($(UNAMES),Linux)
FC = pgf90
FFLAGS = -c -g -Mbounds -DDEBUG
endif

ifeq ($(UNAMES),IRIX64)
FC = f90
CFLAGS += -64 -Xcpluscomm 
FFLAGS = -64 -c -C -I$(INC_NETCDF)
LDFLAGS += -64

ifeq ($(DEBUG),TRUE)
FFLAGS += -g
LDFLAGS += -g
else
FFLAGS += -O2
endif

endif

ifeq ($(UNAMES),AIX)
CFLAGS += -qcpluscmt
FC = xlf90
FFLAGS += -c -qsuffix=f=f90 -I$(INC_NETCDF)
endif

all: $(EXEDIR)/$(EXE1) $(EXEDIR)/$(EXE2)

$(EXEDIR)/$(EXE1): $(OBJS1)
	$(CC) -o $@ $(OBJS1) $(LDFLAGS)

$(EXEDIR)/$(EXE2): $(OBJS2)
	$(FC) -o $@ $(OBJS2) $(LDFLAGS)

clean:
	$(RM) *.o *.mod *.MOD

.SUFFIXES: .f90

.f90.o:
	$(FC) $(FFLAGS) $*.f90

definemld1x1.o: definemld1x1.h
err_exit.o:
wrap_nc.o: definemld1x1.h

definemldbdy.o: definemldbdy.f90 precision.o
binf2c.o: binf2c.f90 precision.o
neighborfill.o: precision.o
wrap_nf.o: wrap_nf.f90
sm121.o: sm121.f90 precision.o
precision.o: precision.f90
