# Makefile to build mkrgridnew on various platforms
# Note: If netcdf library is not built in the standard location, you must set the environment
# variables INC_NETCDF and LIB_NETCDF
#

EXEDIR = .
EXENAME = mkrgridnew
RM = rm

# Check for the NetCDF library and include directories 
ifeq ($(LIB_NETCDF),$(null))
LIB_NETCDF := /usr/local/lib
endif

ifeq ($(INC_NETCDF),$(null))
INC_NETCDF := /usr/local/include
endif

# Determine platform 
UNAMES := $(shell uname -s)
UNAMEM := $(findstring CRAY,$(shell uname -m))

# Architecture-specific flags and rules

#------------------------------------------------------------------------
# Cray 
#------------------------------------------------------------------------

ifeq ($(UNAMEM),CRAY)
FC = f90
FFLAGS = -c -I$(INC_NETCDF)
FFLAGS_FIXED = -c -f fixed
LDFLAGS = -L$(LIB_NETCDF) -lnetcdf
endif

#------------------------------------------------------------------------
# SGI
#------------------------------------------------------------------------

ifeq ($(UNAMES),IRIX64)
FC = f90
FFLAGS = -64 -O2 -c -trapuv -I$(INC_NETCDF)
FFLAGS_FIXED = -64 -O2 -c -trapuv -fixedform
LDFLAGS = -64 -L/usr/local/lib64/r4i4 -lnetcdf
endif

#------------------------------------------------------------------------
# AIX
#------------------------------------------------------------------------

ifeq ($(UNAMES),AIX)
FC = xlf90
FFLAGS = -c -qsuffix=f=f90 -qfree=f90 -I$(INC_NETCDF)
FFLAGS_FIXED = -c -qfixed -qdpc
LDFLAGS = -L/usr/local/lib32/r4i4 -lnetcdf
endif

#------------------------------------------------------------------------
# OSF1
#------------------------------------------------------------------------

ifeq ($(UNAMES),OSF1)
FC = f90
FFLAGS = -c -I$(INC_NETCDF)
FFLAGS_FIXED = -c -fixed
LDFLAGS = -L$(LIB_NETCDF) -lnetcdf
endif

#------------------------------------------------------------------------
# Linux
#------------------------------------------------------------------------

ifeq ($(UNAMES),Linux)

  ifeq ($(USER_FC),$(null))
    FC := pgf90
  else
    FC := $(USER_FC)
  endif

  FFLAGS = -c -I$(INC_NETCDF)
  LDFLAGS = -L$(LIB_NETCDF) -lnetcdf

  ifeq ($(FC),pgf90)
    ifeq ($(DEBUG),TRUE)
      FFLAGS += -g -Ktrap=fp -Mbounds
    else
      FFLAGS += -fast
    endif
    FFLAGS_FIXED = $(FFLAGS) -Mnofreeform
  endif

  ifeq ($(FC),ifort)
    FFLAGS += -traceback
    ifeq ($(DEBUG),TRUE)
      FFLAGS += -g -CB -fpe0 -fpstkchk
    else
      FFLAGS += -O2
    endif
    FFLAGS_FIXED = $(FFLAGS) -fixed
  endif

  ifeq ($(FC),lf95)
    ifeq ($(DEBUG),TRUE)
      FFLAGS += -g --chk e,s,u
      LDFLAGS += -g
    else
      FFLAGS += -O
    endif
    FFLAGS_FIXED = $(FFLAGS) --fix
  endif
endif

#------------------------------------------------------------------------
# Cray X1
#------------------------------------------------------------------------

ifeq ($(UNAMES),UNICOS/mp)
FC      = ftn
FFLAGS  = -c -I$(INC_NETCDF) -O2 -DUNICOSMP
LDFLAGS = -L$(LIB_NETCDF) -lnetcdf
endif

#------------------------------------------------------------------------
# Default rules and macros
#------------------------------------------------------------------------

OBJS := prec.o globals.o fmain.o mkrgrid.o interpolation.o \
	fortfft.o wrap_nf.o special_cases.o

.SUFFIXES:
.SUFFIXES: .F90 .f90 .F .o

.f90.o:
	$(FC) $(FFLAGS) $<

.F90.o:
	$(FC) $(FFLAGS) $<

.F.o:
	$(FC) $(FFLAGS_FIXED) $<

$(EXEDIR)/$(EXENAME): $(OBJS)
	$(FC) -o $@ $(OBJS) $(LDFLAGS)

clean:
	$(RM) -f $(OBJS) *.mod $(EXEDIR)/$(EXENAME)

prec.o:
globals.o: prec.o
fmain.o: mkrgrid.o globals.o
mkrgrid.o: prec.o globals.o interpolation.o
interpolation.o: prec.o globals.o fortfft.o
special_cases.o: prec.o globals.o
wrap_nf.o: prec.o
fortfft.o: prec.o
