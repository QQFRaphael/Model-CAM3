# GnuMakefile for cross platform building of SCAM

#include $(CONFIG_FILE)
include config.out

TARGET	= ../scam

##############################################################################
###########################     SOURCE FILES      ############################
##############################################################################
#
# Automagically generate source and object lists by searching for source files
#  in the vpath
#
#------------------------------------------------------------------------
#------------------------------------------------------------------------
# ESMF logic
#------------------------------------------------------------------------
# Point to the ESMF distribution. 
ifeq ($(ESMF_DIR),$(null))
ESMF_DIR = ${ROOTDIR}/models/utils/esmf
endif

ifeq ($(ESMF_BUILD),$(null))
ESMF_BUILD := $(shell (mkdir ./esmf 2>/dev/null;cd ./esmf;pwd))
endif
DEBUG = $(DEBUG_FLAGS)
# Define where the libesmf.a file resides.
ifeq ($(DEBUG), TRUE)
BOPT_ESMF := g
else
BOPT_ESMF := O
endif
LIB_ESMF := $(ESMF_BUILD)/lib/lib$(BOPT_ESMF)
# Define the location of the ESMF .mod files.
MOD_ESMF := $(ESMF_BUILD)/mod/mod$(BOPT_ESMF)

#------------------------------------------------------------------------
# change VPATH to STD_SRC_PATH if you only want to include sources that
# are in the standard directories.
# The inclusion of "$(NETCDF_ROOT)/include" is to allow make to find
#  the header files for netcdf (see the comments regarding dependency 
#  generation below).
#
MODL_DIR  := $(MODL_DIR)
INIT_DIR  := $(INIT_DIR)
STD_SRC_PATH   = $(MODL_DIR):$(INIT_DIR)
# Load dependency search path.
dirs := . $(shell cat Filepath)  $(MOD_ESMF)/$(ESMF_ARCH)

# Set cpp search path, include netcdf
cpp_dirs := $(dirs)
cpp_path := $(foreach dir,$(cpp_dirs),-I$(dir)) # format for command line

VPATH    := $(foreach dir,$(cpp_dirs),$(wildcard $(dir))) 
VPATH    := $(subst $(space),:,$(VPATH))               
VPATH    := $(VPATH):$(ALT_SOURCE_PATH):$(STD_SRC_PATH):$(NETCDF_INC)

# get the filepaths for all of the source files in the VPATH 
C_SRCS_WITH_PATH = $(foreach dir,$(subst :, ,$(VPATH)) $(dirs), $(wildcard $(dir)/*.c))
F77_SRCS_WITH_PATH = $(foreach dir,$(subst :, ,$(VPATH)) $(dirs), $(wildcard $(dir)/*.F))
F90_SRCS_WITH_PATH = $(foreach dir,$(subst :, ,$(VPATH)) $(dirs), $(wildcard $(dir)/*.F90))

# get filenames w/o paths; remove dups by sorting
F90_SRCS = $(sort $(notdir $(F90_SRCS_WITH_PATH))) 
F77_SRCS = $(sort $(notdir $(F77_SRCS_WITH_PATH))) 
C_SRCS = $(sort $(notdir $(C_SRCS_WITH_PATH))) 

SOURCES   := $(shell cat Srcfiles)

Depends: Srcfiles Filepath
	$(ROOTDIR)/models/atm/cam/bld/mkDepends Filepath Srcfiles > $@

Srcfiles: Filepath
	$(ROOTDIR)/models/atm/cam/bld/mkSrcfiles > $@

OBJS      := $(addsuffix .o, $(basename $(SOURCES)))




#OBJS = $(F90_SRCS:.F90=.o) $(F77_SRCS:.F=.o) $(C_SRCS:.c=.o)

# name of object file containing main()
MAIN_OBJ = scam.o

##############################################################################
#####################   COMPILE LINE VARIABLES         #######################
##############################################################################

# uncomment following line to define the preprocessor "DEBUG" 
# DEBUG_FLAG = -DDEBUG

INCLUDES   = $(cpp_path)
INCLUDES   += $(MODULE_INCLUDE_FLAG)$(MOD_ESMF)/$(ESMF_ARCH) -I$(NETCDF_INC)

# CPPFLAGS is also defined in "config.out", thus ":= $(CPPFLAGS)"
# "-DCAM$(CAM_VERS)_" is for backward  compatibility
CPPFLAGS  := $(CPPFLAGS) -DSCM $(DEBUG_FLAG) \
  -DCAM_VERS=$(CAM_VERS) -DCAM$(CAM_VERS)_  \
  -DCAM_MINOR_VERS=$(CAM_MINOR_VERS) \
  -DPLON=1 -DPLEV=$(PLEV) -DPLEVR=$(PLEV) -DPCNST=$(PCNST) -DPCOLS=1\
  -DPTRM=5 -DPTRN=5 -DPTRK=5 -DPLON=1 -DPLAT=1 -DPNATS=$(PNATS) -DNO_SHR_VMATH\
  -DHIDE_SHR_MSG -DHIDE_MPI -DSHELL_MSS 
#  -DHIDE_SHR_MSG -DHIDE_MPI -DSHELL_MSS -DFORTRANUNDERSCORE
#  -DHIDE_SHR_MSG -DSHELL_MSS -DPERGRO=1 $(INCLUDES) -DPGF90

CINCLUDES = $(INCLUDES) $(MISC_INCLUDES)
FINCLUDES = $(INCLUDES)

FPPFLAGS   = $(CPPFLAGS)

LFLAGS     = -L$(NETCDF_LIB) -lnetcdf $(SYS_LIBS)  -L$(LIB_ESMF)/$(ESMF_ARCH) -lesmf


# AIX compiler needs the preprocessor symbols in a special format

ifeq ($(SYSTEM), AIX)
empty     =
comma     = ,
space     = $(empty) $(empty)
FPPFLAGS  = -WF,$(subst $(space),$(comma),$(strip $(CPPFLAGS)))
endif
# --chkglobal for global checking

ifeq ($(FC90),pgf90)
#FFLAGS   = $(cpp_path) -c -CcdRR8 -DHIDE_SHR_MSG --trace --trap -DLINUX -I$(MOD_ESMF)/$(ESMF_ARCH) -xref 
SPEC_FFLAGS := $(F90FLAGS)
#F90FLAGS = $(FFLAGS)
LDFLAGS  = -L$(NETCDF_LIB) -lnetcdf -L$(LIB_ESMF)/$(ESMF_ARCH) -lesmf --staticlink --print-map
ifeq ($(DEBUG),TRUE)
  FFLAGS += -g --chk --pca -O0
endif
endif


##############################################################################
###########################   COMPILATION RULES   ############################
##############################################################################
.SUFFIXES:
.SUFFIXES:.o .h .F .F90 .c

ifeq ($(FC90),pgf90)
iobinary.o: iobinary.F90
	$(FC90) $(F90FLAGS_NODBG) $(FREE_SOURCE_FLAGS)   $(FPPFLAGS)  $(FINCLUDES) $<
initGridCellsMod.o: initGridCellsMod.F90
	$(FC90) -c -r8 -i4  -Mrecursive -Mdalign -g -Mfree $(FREE_SOURCE_FLAGS)   $(FPPFLAGS)  $(FINCLUDES) $<
inicFileMod.o: inicFileMod.F90
	$(FC90) -c -r8 -i4  -Mrecursive -Mdalign -g -Mfree $(FREE_SOURCE_FLAGS)   $(FPPFLAGS)  $(FINCLUDES) $<
endif
ifeq ($(FC90),lf95lkj;lkj)
# lahey fails on binary_io due to writing wrap areas
#binary_io.o: binary_io.F90
#	$(FC90) $(SPEC_FFLAGS) $(FREE_SOURCE_FLAGS)   $(FPPFLAGS)  $(FINCLUDES) $<
wrap_nf.o: wrap_nf.F90
	$(FC90) $(SPEC_FFLAGS)  $(FREE_SOURCE_FLAGS)   $(FPPFLAGS)  $(FINCLUDES) $<
wrap_mpi.o: wrap_mpi.F90
	$(FC90) $(SPEC_FFLAGS)  $(FREE_SOURCE_FLAGS)   $(FPPFLAGS)  $(FINCLUDES) $<
tlc1.o: tlc1.F90
	$(FC90)   -c -CcdRR8 -DHIDE_SHR_MSG -DLINUX  $(FREE_SOURCE_FLAGS)   $(FPPFLAGS)  $(FINCLUDES) $<
endif

ifeq ($(CC),pgcc)
CFLAGS += -fast
else
CFLAGS += -DUSE_GCC
endif


.c.o:
	$(CC)  $(CFLAGS) $(CPPFLAGS)  $(CINCLUDES) $<
.F90.o:
	$(FC90) $(F90FLAGS) $(FREE_SOURCE_FLAGS)   $(FPPFLAGS)  $(FINCLUDES) $<
.F.o:
	$(FC90) $(F77FLAGS)  $(FIXED_SOURCE_FLAGS)  $(FPPFLAGS) $(FINCLUDES) $<

scam.o :  fortran.h scam_fifo.h

# The library is not made to be built in parallel, so we must null passed options and
# specify only one job is to run.
$(LIB_ESMF)/$(ESMF_ARCH)/libesmf.a:
	cd $(ESMF_DIR); \
	$(MAKE) -j 1 BOPT=$(BOPT_ESMF) ESMF_BUILD=$(ESMF_BUILD) ESMF_DIR=$(ESMF_DIR) ESMF_ARCH=$(ESMF_ARCH);

time_manager.o : $(LIB_ESMF)/$(ESMF_ARCH)/libesmf.a

##############################################################################
######################  TARGETS/DEPENDENCIES   ###############################
##############################################################################

all: $(TARGET) 

echocomp:
	echo $(PATH);which pgf90

$(TARGET): $(OBJS)
	$(LINKER) -g $(OBJS) -o $(TARGET) $(LFLAGS)

# if the config file is changed remake everything
#$(OBJS): $(CONFIG_FILE)

#==================================================================

clean: 
	rm -rf ../core *.[ocd] *.mod esmf *.stb 

realclean: clean
	rm -f $(TARGET) $(TOP_DIR)/*/*~ TAGS config.out
	rm -f Rootdir Filepath Srcfiles Depends *.h

# "make dist" will bundle everything up in a tar file, ready for
#   distribution
dist: realclean
	(cd ..;rm -rf ../$(SCAM); mkdir ../$(SCAM); \
	 cp -r cam2 cam_stubs html mymods* scm_init \
	 ../$(SCAM);\
	 cp -r data.dist ../$(SCAM)/data;\
	 mkdir ../$(SCAM)/obj ../$(SCAM)/userdata;\
	 cp obj/makdep* ../$(SCAM)/obj; \
	 cp configure GNUmakefile INSTALL Makefile .scam_defaults scamgui \
	 ../$(SCAM);\
         $(STRIP) ../$(SCAM)/scamgui; \
	mkdir -p ../$(SCAM)/lib/ncarg/database; \
        cp -r $(NCARG_ROOT)/lib/ncarg/database/* ../$(SCAM)/lib/ncarg/database; \
        rm -rf ../$(SCAM)/lib/ncarg/database/rangs; \
        mkdir -p ../$(SCAM)/lib/ncarg/fontcaps; \
        cp -r $(NCARG_ROOT)/lib/ncarg/fontcaps/* ../$(SCAM)/lib/ncarg/fontcaps; \
	rm -rf ../$(SCAM)/cam2/models/atm/cam/src/bfb*;\
	rm -rf ../$(SCAM)/cam2/models/atm/cam/src/obj*;\
        cd ..; \
        tar cvf - $(SCAM) | gzip > $(SCAM).$(SYSTEM).tar.gz)
# "make dist" will bundle everything up in a tar file, ready for
#   distribution
distui: realclean
	(cd ..;rm -rf ../$(SCAM); mkdir ../$(SCAM); \
	 cp -r cam2 cam_stubs html mymods* scm_init ui tools \
	 ../$(SCAM);\
	 cp -r data.dist ../$(SCAM)/data;\
	 mkdir ../$(SCAM)/obj ../$(SCAM)/userdata;\
	 cp obj/makdep* ../$(SCAM)/obj; \
	 cp configure GNUmakefile INSTALL Makefile .scam_defaults scamgui \
	 ../$(SCAM);\
         $(STRIP) ../$(SCAM)/scamgui; \
	mkdir -p ../$(SCAM)/lib/ncarg/database; \
        cp -r $(NCARG_ROOT)/lib/ncarg/database/* ../$(SCAM)/lib/ncarg/database; \
        rm -rf ../$(SCAM)/lib/ncarg/database/rangs; \
        mkdir -p ../$(SCAM)/lib/ncarg/fontcaps; \
        cp -r $(NCARG_ROOT)/lib/ncarg/fontcaps/* ../$(SCAM)/lib/ncarg/fontcaps; \
	rm -rf ../$(SCAM)/cam2/models/atm/cam/src/bfb*;\
	rm -rf ../$(SCAM)/cam2/models/atm/cam/src/obj*;\
        cd ..; \
        tar cvf - $(SCAM) | gzip > $(SCAM).ui.$(SYSTEM).tar.gz)


modsclean:
	rm -f $(ALT_OBJS)

showpath:
	@echo TOP_DIR = $(TOP_DIR)
	@echo VPATH = $(VPATH)

tags:
	etags $(INIT_DIR)/*.[Fch] $(MODL_DIR)/*.[Fh]

depend:
	@echo 
	@echo  "Depend: This target is deprecated - dependencies are now automatically generated and updated."
	@echo

showsrc:
	@echo ***CSRCS***  $(C_SRCS)
	@echo ***F90SRCS***  $(F90_SRCS)
	@echo ***F77SRCS***   $(F77_SRCS) 
	@echo ***SRCS_WITH_PATH***  $(C_SRCS_WITH_PATH) $(F90_SRCS_WITH_PATH)  $(F77_SRCS_WITH_PATH) 

showconfig:
	@echo MODEL_VERSION = $(MODEL_VERSION)
	@echo PLEV = $(PLEV)
	@echo PCONST = $(PCONST)
	@echo VPATH = $(VPATH)
	@echo CPP_PATH = $(cpp_path)
	@echo INCLUDES = $(INCLUDES)
	@echo FINCLUDES = $(FINCLUDES)
	@echo CINCLUDES = $(CINCLUDES)

showobj:
	@echo ***OBJS*** $(OBJS)

.PHONY: all  clean realclean

-include Depends
###########################################################
##########    AUTOMATIC DEPENDENCY GENERATION    ##########
###########################################################
#
# The following rules are modified from the Gnu Make manual which
# contains an explanation of the basics.
# here's an explanation of the modifications:
#  -- uses cpp instead of cc so that it will operate on FORTRAN
#     source files, and put each dependency on a separate line.
#  -- added one more sed command to remove the pathname
#    from the include files
#    that the source depends on; this allows VPATH to properly
#    determine which header file to consider when there is more
#    than one possibility (and incidentally necessitates the presence
#    of the previous command).
############################################################
#ifeq ($(GENDEP), MAKDEP)
#depgen := makdep
#ifeq ($(SYSTEM), Linux)
#depgen := makdep.linux
#endif
#%.d : %.c
#	@echo "Building dependency file $@"
#	@./$(depgen) -f $(cpp_path) $< > $@
#%.d : %.F
#	@echo "Building dependency file $@"
#	@./$(depgen) -f $(cpp_path) $< > $@
#%.d : %.F90
#	@echo "Building dependency file $@"
#	@./$(depgen) -f $(cpp_path) $<  > $@
#endif
#
#ifeq ($(GENDEP), MAKEDEPEND)
#%.d: %.F
#	@ echo "generating dependencies for $<"
#	@ $(SHELL) -ec 'makedepend $(CPPFLAGS) -f- $< \
#               | sed -e '\''s/\/[^ ]*\///g'\''\
#                     -e '\''s/\($*\)\.o :/\1.o $@ : /g'\'' > $@'
#%.d: %.F90
#	@ echo "generating dependencies for $<"
#	@ $(SHELL) -ec 'makedepend $(CPPFLAGS) -f- $< \
#               | sed -e '\''s/\/[^ ]*\///g'\''\
#                     -e '\''s/\($*\)\.o :/\1.o $@ : /g'\'' > $@'
#%.d: %.c
#	@ echo "generating dependencies for $<"
#	@ $(SHELL) -ec 'makedepend $(CPPFLAGS) -f- $< \
#              | sed -e '\''s/\($*\)\.o :/\1.o $@ : /g'\'' > $@'
#endif
#
#ifeq ($(GENDEP), CPP)
#%.d: %.F
#	@ echo "generating dependencies for $<"
#	@ $(SHELL) -ec '$(CPP) -M$(CPPFLAGS) $< \
#               | sed -e '\''s/\/[^ ]*\///g'\''\
#                     -e '\''s/\($*\)\.o :/\1.o $@ : /g'\'' > $@'
#%.d: %.c
#	@ echo "generating dependencies for $<"
#	@ $(SHELL) -ec '$(CPP) -M$(CPPFLAGS) $< \
#              | sed -e '\''s/\($*\)\.o :/\1.o $@ : /g'\'' > $@'
#endif
#
#ifeq ($(GENDEP), CPP_HP)
## on HP-UX, CPP prints the dependencies to stderr, so we need to swap
## stdout and stderr, and then send stderr to /dev/null  with "3>&2 2>&1 1>/dev/null"
#%.d: %.F
#	@ echo "generating dependencies for $<"
#	@ $(SHELL) -ec '$(CPP) -M$(CPPFLAGS) $< \
#               3>&2 2>&1 1>/dev/null \
#               | sed -e '\''s/\/[^ ]*\///g'\''\
#                     -e '\''s/\($*\)\.F\.o :/\1.o $@ : /g'\'' > $@'
#%.d: %.c
#	@ echo "generating dependencies for $<"
#	@ $(SHELL) -ec '$(CPP) -M$(CPPFLAGS) $< \
#              3>&2 2>&1 1>/dev/null \
#              | sed -e '\''s/\($*\)\.o :/\1.o $@ : /g'\'' > $@'
#endif
#
## include the dependencies automatically generated for each source file
## the '-' is to suppress some irrelevant warning messages
#
#ifeq ($(MAKECMDGOALS), all)
#FIND_FILES = $(wildcard $(dir)/*.c $(dir)/*.F $(dir)/*.F90)
#FILES      = $(foreach dir, $(dirs),$(FIND_FILES))
#SOURCES   := $(sort $(notdir $(FILES)))
#DEPS      := $(addsuffix .d, $(basename $(SOURCES)))
#-include $(C_SRCS:.c=.d) $(F90_SRCS:.F90=.d) $(F77_SRCS:.F=.d)
#endif
