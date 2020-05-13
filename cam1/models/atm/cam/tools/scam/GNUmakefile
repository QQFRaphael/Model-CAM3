# SCAM GNUmakefile
# requires Gnu Make
# recursively invokes make in the build directory
#
# $Id: GNUmakefile,v 1.1.6.3 2004/10/18 17:40:38 jmccaa Exp $
#

export TOP_DIR = $(shell pwd)
export SCAM_DIR = $(notdir $(TOP_DIR))
export SCAM_VERS = cam2_0_2_dev18_brnch_sci16
export SCAM = scam-$(SCAM_VERS)
export BUILD_DIR = obj
export GUI_DIR = ui
export CONFIG_FILE = $(TOP_DIR)/$(BUILD_DIR)/config.out
export GUI_CONFIG_FILE = $(TOP_DIR)/$(GUI_DIR)/config.out
export GUI_TARGET = ./scamgui
export SCAM_TARGET = ./scam
SYSTEM=`uname`

all: configure_files
	@cd $(GUI_DIR); $(MAKE) all 
	@cd $(BUILD_DIR); $(MAKE) all

clean: 
	@touch $(CONFIG_FILE)
	@cd $(BUILD_DIR);rm -rf ../core *.[ocd] *.mod esmf *.stb config.out

realclean: 
	@touch $(CONFIG_FILE) $(GUI_CONFIG_FILE)
	@\rm -f scamgui* core.*
	@cd $(GUI_DIR); $(MAKE) "MAKECMDGOALS = $@" $@
	@cd $(BUILD_DIR);touch Rootdir Filepath Srcfiles Depends;$(MAKE) "MAKECMDGOALS = $@" $@
	\rm -f $(CONFIG_FILE) $(GUI_CONFIG_FILE) scm_init/scam_rpc.h
	\rm -f userdata/.scamhist.tmp* userdata/crm_quickstart userdata/*crm.nc

configure_files: $(GUI_CONFIG_FILE) $(CONFIG_FILE)

$(GUI_TARGET): $(GUI_CONFIG_FILE)
	@cd ui; $(MAKE) all

$(SCAM_TARGET): $(GUI_TARGET) $(CONFIG_FILE)
	@cd $(BUILD_DIR); $(MAKE) all

$(CONFIG_FILE):
	@echo '======================================================'
	@echo '      Configuring Single Column Version of CAM        '
	@echo '======================================================'
	./configure

$(GUI_CONFIG_FILE):
	@echo '======================================================'
	@echo '              Configuring qt GUI for SCAM'
	@echo '======================================================'
	@cd ui; ./configure

tgz:
	@cd ..; tar cvf - `find $(SCAM_DIR) -type f ! -name scam ! -name scamgui \
                 ! -name "*.o" ! -name "*.d" ! -name "*.nc" -print` \
	           | gzip > scam-$(SCAM_VERS).dev.tgz
