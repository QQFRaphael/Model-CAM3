# $Id: makefile,v 1.1.6.1 2002/04/24 03:25:32 erik Exp $ 

ALL: lib 

CFLAGS    =
FFLAGS    =
SOURCEC	  = ESMC_BasicUtil.c
SOURCEF	  = ESMF_BasicUtilMod.F
SOURCEH	  = ESMC_BasicUtil.h
OBJSC     = $(addsuffix .o, $(basename $(SOURCEC)))
OBJSF     = $(addsuffix .o, $(basename $(SOURCEF)))
TEXFILES  = $(addsuffix .tex, $(basename $(SOURCEF)))
MANSEC	  = Timeobj
LIBBASE	  = libesmf
DIRS	  = 
LOCDIR	  = src/Infrastructure/BasicUtil/
CLEANFILES = $(LIBNAME) $(addprefix doc/, $(addsuffix .tex, $(basename $(SOURCEF))))

include ${ESMF_DIR}/build/${ESMF_ARCH}/base






