#! /bin/tcsh -f
#
#=======================================================================
#
#  run-lightning.csh
#
#  Generic batch submission script for lightning using LSF.  
#
#-----------------------------------------------------------------------
# Batch options for machine with LSF batch system. 
# Usage for PathScale compiler (default): 
#   bsub < run-lightning.csh
#-----------------------------------------------------------------------
#
#BSUB -a mpich_gm
#BSUB -x
#BSUB -n 2,4
#BSUB -R "span[ptile=2]"
#BSUB -o cam.o%J
#BSUB -e cam.e%J
#BSUB -q regular
#
#=======================================================================

setenv INC_NETCDF /contrib/2.6/netcdf/3.6.0-p1-pathscale-2.2.1-64/include
setenv LIB_NETCDF /contrib/2.6/netcdf/3.6.0-p1-pathscale-2.2.1-64/lib
set mpich=/contrib/2.6/mpich-gm/1.2.6..14a-pathscale-2.2.1-64
setenv INC_MPI ${mpich}/include
setenv LIB_MPI ${mpich}/lib
set ps=/contrib/2.6/pathscale/2.2.1
setenv PATH ${mpich}/bin:${ps}/bin:${PATH}
setenv LD_LIBRARY_PATH ${ps}/lib/2.2.1:${LD_LIBRARY_PATH}

## ROOT OF CAM DISTRIBUTION - probably needs to be customized.
## Contains the source code for the CAM distribution.
## (the root directory contains the subdirectory "models")
set camroot      = /fis/cgd/...

## ROOT OF CAM DATA DISTRIBUTION - needs to be customized unless running at NCAR.
## Contains the initial and boundary data for the CAM distribution.
## (the root directory contains the subdirectories "atm" and "lnd")
setenv CSMDATA /fis/cgd/cseg/csm/inputdata    

## Default namelist settings:
## $case is the case identifier for this run. It will be placed in the namelist.
## $runtype is the run type: initial, restart, or branch.
## $nelapse is the number of timesteps to integrate, or number of days if negative.
set case         = camrun
set runtype      = initial
set nelapse      = -1

## $wrkdir is a working directory where the model will be built and run.
## $blddir is the directory where model will be compiled.
## $rundir is the directory where the model will be run.
## $cfgdir is the directory containing the CAM configuration scripts.
set wrkdir       = /ptmp/$LOGNAME
set blddir       = $wrkdir/$case/bld
set rundir       = $wrkdir/$case
set cfgdir       = $camroot/models/atm/cam/bld

## Ensure that run and build directories exist
mkdir -p $rundir                || echo "cannot create $rundir" && exit 1
mkdir -p $blddir                || echo "cannot create $blddir" && exit 1

## If an executable doesn't exist, build one.
if ( ! -x $blddir/cam ) then
    cd $blddir                  || echo "cd $blddir failed" && exit 1
    $cfgdir/configure -spmd -nosmp -fc pathf90 -ldflags "-L/opt/gm/lib64 -lgm"  || echo "configure failed" && exit 1
    echo "building CAM in $blddir ..."
    rm -f Depends
    gmake >&! MAKE.out      || echo "CAM build failed: see $blddir/MAKE.out" && exit 1
endif

## Create the namelist
cd $blddir                      || echo "cd $blddir failed" && exit 1
$cfgdir/build-namelist -s -case $case -runtype $runtype -o $rundir/namelist \
 -namelist "&camexp nelapse=$nelapse mss_irt=0 /"  || echo "build-namelist failed" && exit 1

## Run CAM
cd $rundir                      || echo "cd $rundir failed" && exit 1
echo "running CAM in $rundir"
mpirun.lsf -np 2 $blddir/cam < namelist  || echo "CAM run failed" && exit 1

exit 0
