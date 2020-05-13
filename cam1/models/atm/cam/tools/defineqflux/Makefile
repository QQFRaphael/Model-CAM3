CC = cc
CFLAGS = -O2 -I$(INC_NETCDF) -UDEBUG
LDFLAGS = -L$(LIB_NETCDF) -lm -lnetcdf
OBJS1 = check_consistent.o defineqflux.o err_exit.o mksith.o timediddle_mavg.o wrap_nc.o \
	gepp.o backsolve.o printeq.o
EXENAME1 = defineqflux

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
ifeq ($(DEBUG),TRUE)
CFLAGS = -g
LDFLAGS += -g
endif
endif

ifeq ($(UNAMES),IRIX64)
CFLAGS += -64 -Xcpluscomm 
LDFLAGS = -g -L$(LIB_NETCDF) -lm -lnetcdf
LDFLAGS += -64
endif

ifeq ($(UNAMES),AIX)
LDFLAGS = -g -L$(LIB_NETCDF) -lm -lnetcdf
CFLAGS += -qcpluscmt
endif

all: $(EXEDIR)/$(EXENAME1)

$(EXEDIR)/$(EXENAME1): $(OBJS1)
	$(CC) -o $@ $(OBJS1) $(LDFLAGS)

clean:
	$(RM) *.o

check_consistent.o: defineqflux.h
defineqflux.o: defineqflux.h
definesomic.o: definesomic.h
testqflux.o: defineqflux.h
err_exit.o:
mksith.o: defineqflux.h
wrap_nc.o: defineqflux.h
timediddle_mavg.o: defineqflux.h
gepp.o: defineqflux.h
backsolve.o: defineqflux.h
printeq.o: defineqflux.h
