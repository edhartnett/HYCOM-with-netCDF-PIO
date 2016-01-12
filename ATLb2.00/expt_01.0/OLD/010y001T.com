#
set echo
set time = 1
set timestamp
#
# --- check that the C comment command is available.
#
if (`which C | cut -c 1-2` == "no") then
  if (-e ${home}/hycom/ALL/bin/C) then
    set path = ( ${path} ${home}/hycom/ALL/bin )
  else
    echo "Please put the command hycom/ALL/bin/C in your path"
  endif
endif
C
C --- Experiment ATLb2.00 - 01.0
C --- 22 layer HYCOM on 2.00 degree Atlantic region
C
C --- 01.0 - KPP: COADS forcing, Levitus bndry and sur. sal. relax.
C
C --- Preamble, script keys on O/S name.
C
C --- Set parallel configuration, see ../README/README.expt_parallel.
C --- NOMP = number of OpenMP threads, 0 for no OpenMP, 1 for inactive OpenMP
C --- NMPI = number of MPI    tasks,   0 for no MPI
C
setenv OS `uname`
switch ($OS)
case 'SunOS':
case 'Linux':
case 'OSF1':
case 'IRIX64':
case 'AIX':
    setenv NOMP 0
    setenv NMPI 0
    breaksw
case 'unicosmk':
    setenv ACCT `groups | awk '{print $1}'`
#   always TYPE=shmem and NOMP=0 on T3E
    setenv TYPE shmem
    setenv NOMP 0
    setenv NMPI `limit -v | grep "MPP PE limit" | awk '{printf("%d\n",$4)}'`
    breaksw
case 'unicos':
    setenv ACCT `newacct -l | awk '{print $4}'`
    setenv NOMP 0
    setenv NMPI 0
    breaksw
default:
    echo 'Unknown Operating System: ' $OS
    exit (1)
endsw
C
C --- modify NOMP and NMPI based on batch limits
C
if ( $?LSB_MCPU_HOSTS ) then
# LSF batch system
  if      ($NMPI == 0) then
    setenv NOMP `echo $LSB_MCPU_HOSTS | awk '{print $2}'`
  else if ($NOMP == 0) then
    setenv NMPI `echo $LSB_MCPU_HOSTS | awk '{print $2}'`
  else
    setenv NMPI `echo $LSB_MCPU_HOSTS $NOMP | awk '{print int($2/$3)}'`
  endif
else if ( $?GRD_TOTAL_MPI_TASKS ) then
# GRD batch system
  if      ($NMPI == 0) then
    echo "error - NMPI=0, but running in a MPI batch queue"
    exit
  else
    setenv NMPI $GRD_TOTAL_MPI_TASKS
  endif
else if ( $?NSLOTS ) then
# codine or GRD batch system
  if      ($NMPI == 0) then
    setenv NOMP $NSLOTS
  else if ($NOMP == 0) then
    setenv NMPI $NSLOTS
  else
    setenv NMPI `echo $NSLOTS $NOMP | awk '{print int($1/$2)}'`
  endif
endif
echo "NOMP is " $NOMP " and NMPI is " $NMPI
C
C --- R is region name.
C --- V is source code version number.
C --- T is topography number.
C --- K is number of layers.
C --- E is expt number.
C --- P is primary path.
C --- D is permanent directory.
C --- S is scratch   directory, must not be the permanent directory.
C
setenv R ATLb2.00
setenv V 2.1.00
setenv T 01
setenv K 22
setenv E 010
setenv P hycom/${R}/expt_01.0/data
setenv D ~/$P
C
switch ($OS)
case 'SunOS':
    if (-e /net/hermes/scrb) then
#                  NRLSSC
      setenv S     /net/hermes/scrb/${user}/$P
    else
#                  NAVO MSRC
      mkdir        /scr/${user}
      chmod a+rx   /scr/${user}
      setenv S     /scr/${user}/$P
    endif
    breaksw
case 'Linux':
    if (-e /external/fast) then
#                  NRLSSC
      setenv S /external/fast/${user}/$P
    else
#                  Single Disk
      setenv S ~/$P/SCRATCH
    endif
    breaksw
case 'OSF1':
    mkdir        ~/scratch
    chmod a+rx   ~/scratch
    setenv S     ~/scratch/$P
    breaksw
case 'IRIX64':
    mkdir        /scr/${user}
    chmod a+rx   /scr/${user}
    setenv S     /scr/${user}/$P
    breaksw
case 'AIX':
    if      (-e /gpfs/work) then
#                  ERDC MSRC, under PBS
      mkdir        /gpfs/work/${user}
      chmod a+rx   /gpfs/work/${user}
      setenv S     /gpfs/work/${user}/$P
      setenv POE  pbspoe
    else if (-e /scr) then
#                  NAVO MSRC, under LoadLeveler
      mkdir        /scr/${user}
      chmod a+rx   /scr/${user}
      setenv S     /scr/${user}/$P
      setenv POE  poe
    else
#                  ARL MSRC, under GRD
      mkdir        /usr/var/tmp/${user}
      chmod a+rx   /usr/var/tmp/${user}
      setenv S     /usr/var/tmp/${user}/$P
      setenv POE  grd_poe
    endif
    breaksw
case 'unicos':
case 'unicosmk':
    mkdir        /tmp/${user}
    chmod a+rx   /tmp/${user}
    mkdir        /tmp/${user}/$R
    chgrp $ACCT  /tmp/${user}/$R
    setenv S     /tmp/${user}/$P
    breaksw
endsw
C
mkdir -p $S
cd       $S
C
C --- For whole year runs.
C ---   ymx number of years per model run.
C ---   Y01 initial model year of this run.
C ---   YXX is the last model year of this run, and the first of the next run.
C ---   A and B are identical, typically blank.
C --- For part year runs.
C ---   A is this part of the year, B is next part of the year.
C ---   Y01 initial model year of this run.
C ---   YXX year at end of this part year run.
C ---   ymx is 1.
C --- Note that these variables and the .awk generating script must
C ---  be consistant to get a good run script.
C
C --- One year spin-up run.
C
@ ymx =  1
C
setenv A ""
setenv B ""
setenv Y01 "001"
C
switch ("${B}")
case "${A}":
    setenv YXX `echo $Y01 $ymx | awk '{printf("%03d", $1+$2)}'`
    breaksw
case "a":
    setenv YXX `echo $Y01 | awk '{printf("%03d", $1+1)}'`
    breaksw
default:
    setenv YXX $Y01
endsw
C
echo "Y01 =" $Y01 "YXX = " $YXX  "A =" ${A} "B =" ${B}
#C
#C --- turn on detailed debugging.
#C
#touch PIPE_DEBUG
C
C --- local input files.
C
if (-e ${D}/../${E}y${Y01}${A}.limits) then
  /bin/cp ${D}/../${E}y${Y01}${A}.limits limits
else
#  use "LIMITI"  when starting a run after day zero.
#  use "LIMITS9" (say) for a 9-day run
  echo "LIMITS" | awk -f ${D}/../${E}.awk y01=${Y01} ab=${A} >! limits
endif
cat limits
C
if (-e ${D}/../ports.input_${Y01}${A}) then
  /bin/cp ${D}/../ports.input_${Y01}${A} ports.input
else
  /bin/cp ${D}/../ports.input ports.input
endif
C
if (-e ${D}/../blkdat.input_${Y01}${A}) then
  /bin/cp ${D}/../blkdat.input_${Y01}${A} blkdat.input
else
  /bin/cp ${D}/../blkdat.input blkdat.input
endif
C
C --- check that iexpt from blkdat.input agrees with E from this script.
C
setenv EB `grep "'iexpt ' =" blk* | awk '{printf("%03d", $1)}'`
if ($EB != $E) then
  cd $D/..
  /bin/mv LIST LIST_BADRUN
  echo "BADRUN" > LIST
  exit
endif
C
C --- pget, pput "copy" files between scratch and permanent storage.
C --- Can both be cp if the permanent filesystem is mounted locally.
C
switch ($OS)
case 'SunOS':
#case 'Linux':
case 'OSF1':
case 'IRIX64':
case 'AIX':
case 'unicos':
case 'unicosmk':
    if (-e ~wallcraf/bin/pget) then
      setenv pget ~wallcraf/bin/pget
      setenv pput ~wallcraf/bin/pput
    else
      setenv pget /bin/cp
      setenv pput /bin/cp
    endif
    breaksw
default:
    setenv pget /bin/cp
    setenv pput /bin/cp
endsw
C
C --- input files from file server.
C
touch regional.depth.a regional.depth.b
if (-z regional.depth.a) then
   ${pget} ${D}/../../topo/depth_${R}_${T}.a regional.depth.a &
endif
if (-z regional.depth.b) then
   ${pget} ${D}/../../topo/depth_${R}_${T}.b regional.depth.b &
endif
C
touch regional.grid.a regional.grid.b
if (-z regional.grid.a) then
   ${pget} ${D}/../../topo/regional.grid.a regional.grid.a &
endif
if (-z regional.grid.b) then
   ${pget} ${D}/../../topo/regional.grid.b regional.grid.b &
endif
C
if (! -e ./wind) then
C
C --- Climatological atmospheric forcing.
C
  touch forcing.tauewd.a forcing.taunwd.a forcing.wndspd.a forcing.airtmp.a forcing.precip.a forcing.radflx.a forcing.shwflx.a forcing.vapmix.a
  touch forcing.tauewd.b forcing.taunwd.b forcing.wndspd.b forcing.airtmp.b forcing.precip.b forcing.radflx.b forcing.shwflx.b forcing.vapmix.b
  if (-z forcing.tauewd.a) then
     ${pget} ${D}/../../force/coads/tauewd.a      forcing.tauewd.a &
  endif
  if (-z forcing.tauewd.b) then
     ${pget} ${D}/../../force/coads/tauewd.b      forcing.tauewd.b &
  endif
  if (-z forcing.taunwd.a) then
     ${pget} ${D}/../../force/coads/taunwd.a      forcing.taunwd.a &
  endif
  if (-z forcing.taunwd.b) then
     ${pget} ${D}/../../force/coads/taunwd.b      forcing.taunwd.b &
  endif
  if (-z forcing.wndspd.a) then
     ${pget} ${D}/../../force/coads/wndspd.a      forcing.wndspd.a &
  endif
  if (-z forcing.wndspd.b) then
     ${pget} ${D}/../../force/coads/wndspd.b      forcing.wndspd.b &
  endif
  if (-z forcing.airtmp.a) then
     ${pget} ${D}/../../force/coads/airtmp.a      forcing.airtmp.a &
  endif
  if (-z forcing.airtmp.b) then
     ${pget} ${D}/../../force/coads/airtmp.b      forcing.airtmp.b &
  endif
  if (-z forcing.vapmix.a) then
     ${pget} ${D}/../../force/coads/vapmix.a      forcing.vapmix.a &
  endif
  if (-z forcing.vapmix.b) then
     ${pget} ${D}/../../force/coads/vapmix.b      forcing.vapmix.b &
  endif
  if (-z forcing.precip.a) then
     ${pget} ${D}/../../force/coads/precip_zero.a forcing.precip.a &
  endif
  if (-z forcing.precip.b) then
     ${pget} ${D}/../../force/coads/precip_zero.b forcing.precip.b &
  endif
  if (-z forcing.radflx.a) then
     ${pget} ${D}/../../force/coads/radflx.a      forcing.radflx.a &
  endif
  if (-z forcing.radflx.b) then
     ${pget} ${D}/../../force/coads/radflx.b      forcing.radflx.b &
  endif
  if (-z forcing.shwflx.a) then
     ${pget} ${D}/../../force/coads/shwflx.a      forcing.shwflx.a &
  endif
  if (-z forcing.shwflx.b) then
     ${pget} ${D}/../../force/coads/shwflx.b      forcing.shwflx.b &
  endif
endif
C
touch relax.rmu.a relax.saln.a relax.temp.a relax.intf.a
touch relax.rmu.b relax.saln.b relax.temp.b relax.intf.b
if (-z relax.rmu.a) then
   ${pget} ${D}/../../relax/${E}/relax_rmu.a relax.rmu.a  &
endif
if (-z relax.rmu.b) then
   ${pget} ${D}/../../relax/${E}/relax_rmu.b relax.rmu.b  &
endif
if (-z relax.saln.a) then
   ${pget} ${D}/../../relax/${E}/relax_sal.a relax.saln.a &
endif
if (-z relax.saln.b) then
   ${pget} ${D}/../../relax/${E}/relax_sal.b relax.saln.b &
endif
if (-z relax.temp.a) then
   ${pget} ${D}/../../relax/${E}/relax_tem.a relax.temp.a &
endif
if (-z relax.temp.b) then
   ${pget} ${D}/../../relax/${E}/relax_tem.b relax.temp.b &
endif
if (-z relax.intf.a) then
   ${pget} ${D}/../../relax/${E}/relax_int.a relax.intf.a &
endif
if (-z relax.intf.b) then
   ${pget} ${D}/../../relax/${E}/relax_int.b relax.intf.b &
endif
C
C --- restart input
C
touch   restart_out.a restart_out.b
/bin/mv restart_out.a restart_in.a
/bin/mv restart_out.b restart_in.b
if (-z restart_in.a) then
  ${pget} ${D}/restart_${Y01}${A}.a restart_in.a &
endif
if (-z restart_in.b) then
  ${pget} ${D}/restart_${Y01}${A}.b restart_in.b &
endif
C
C --- model executable
C
if      ($NMPI == 0 && $NOMP == 0) then
  setenv TYPE one
else if ($NMPI == 0) then
  setenv TYPE omp
else if ($NOMP == 0) then
  if ( ! $?TYPE ) then
    setenv TYPE mpi
  endif
else
  setenv TYPE ompi
endif
/bin/cp ${D}/../../src_${V}_${K}_${TYPE}/hycom . &
C
C --- summary printout
C
touch   summary_out
/bin/mv summary_out summary_old
C
C --- heat transport output
C
touch   flxdp_out.a flxdp_out.b
/bin/mv flxdp_out.a flxdp_old.a
/bin/mv flxdp_out.b flxdp_old.b
C
touch   ovrtn_out
/bin/mv ovrtn_out ovrtn_old
C
C --- let all file copies complete.
C
wait
C
C --- Just in time atmospheric forcing.
C
if (-e ./wind) then
  if (! -e ./flux) then
    echo './flux must exist if ./wind does'
    exit
  endif
C
C --- Check to see if wind and flux files exist, if not make them and wait.
C
  /bin/rm -f forcing.tauewd.a forcing.taunwd.a forcing.wndspd.a
  /bin/rm -f forcing.tauewd.b forcing.taunwd.b forcing.wndspd.b
  /bin/rm -f ./wind/${E}w${Y01}${A}
  if (-e ./wind/tauewd_${Y01}${A}.a && \
      -e ./wind/tauewd_${Y01}${A}.b && \
      -e ./wind/taunwd_${Y01}${A}.a && \
      -e ./wind/taunwd_${Y01}${A}.b && \
      -e ./wind/wndspd_${Y01}${A}.a && \
      -e ./wind/wndspd_${Y01}${A}.b    ) then
    /bin/ln ./wind/tauewd_${Y01}${A}.a forcing.tauewd.a
    /bin/ln ./wind/taunwd_${Y01}${A}.a forcing.taunwd.a
    /bin/ln ./wind/wndspd_${Y01}${A}.a forcing.wndspd.a
    /bin/ln ./wind/tauewd_${Y01}${A}.b forcing.tauewd.b
    /bin/ln ./wind/taunwd_${Y01}${A}.b forcing.taunwd.b
    /bin/ln ./wind/wndspd_${Y01}${A}.b forcing.wndspd.b
  else
    cd ./wind
    touch ${E}w${Y01}${A}
    /bin/rm -f ${E}w${Y01}${A}.com ${E}w${Y01}${A}.log
    awk -f $D/../${E}.awk y01=${Y01} ab=${A} $D/../${E}W.com > ${E}w${Y01}${A}.com
    csh ${E}w${Y01}${A}.com >& ${E}w${Y01}${A}.log &
    cd ..
  endif
  /bin/rm -f forcing.airtmp.a forcing.precip.a forcing.radflx.a forcing.shwflx.a forcing.vapmix.a
  /bin/rm -f forcing.airtmp.b forcing.precip.b forcing.radflx.b forcing.shwflx.b forcing.vapmix.b
  /bin/rm -f ./flux/${E}f${Y01}${A}
  if (-e ./flux/airtmp_${Y01}${A}.a && \
      -e ./flux/airtmp_${Y01}${A}.b && \
      -e ./flux/precip_${Y01}${A}.a && \
      -e ./flux/precip_${Y01}${A}.b && \
      -e ./flux/radflx_${Y01}${A}.a && \
      -e ./flux/radflx_${Y01}${A}.b && \
      -e ./flux/shwflx_${Y01}${A}.a && \
      -e ./flux/shwflx_${Y01}${A}.b && \
      -e ./flux/vapmix_${Y01}${A}.a && \
      -e ./flux/vapmix_${Y01}${A}.b    ) then
    /bin/ln ./flux/airtmp_${Y01}${A}.a forcing.airtmp.a
    /bin/ln ./flux/precip_${Y01}${A}.a forcing.precip.a
    /bin/ln ./flux/radflx_${Y01}${A}.a forcing.radflx.a
    /bin/ln ./flux/shwflx_${Y01}${A}.a forcing.shwflx.a
    /bin/ln ./flux/vapmix_${Y01}${A}.a forcing.vapmix.a
    /bin/ln ./flux/airtmp_${Y01}${A}.b forcing.airtmp.b
    /bin/ln ./flux/precip_${Y01}${A}.b forcing.precip.b
    /bin/ln ./flux/radflx_${Y01}${A}.b forcing.radflx.b
    /bin/ln ./flux/shwflx_${Y01}${A}.b forcing.shwflx.b
    /bin/ln ./flux/vapmix_${Y01}${A}.b forcing.vapmix.b
  else
    cd ./flux
    touch ${E}f${Y01}${A}
    /bin/rm -f ${E}f${Y01}${A}.com ${E}f${Y01}${A}.log
    awk -f $D/../${E}.awk y01=${Y01} ab=${A} $D/../${E}F.com > ${E}f${Y01}${A}.com
    csh ${E}f${Y01}${A}.com >& ${E}f${Y01}${A}.log &
    cd ..
  endif
  wait
  if (-e ./wind/${E}w${Y01}${A}) then
    /bin/ln ./wind/tauewd_${Y01}${A}.a forcing.tauewd.a
    /bin/ln ./wind/taunwd_${Y01}${A}.a forcing.taunwd.a
    /bin/ln ./wind/wndspd_${Y01}${A}.a forcing.wndspd.a
    /bin/ln ./wind/tauewd_${Y01}${A}.b forcing.tauewd.b
    /bin/ln ./wind/taunwd_${Y01}${A}.b forcing.taunwd.b
    /bin/ln ./wind/wndspd_${Y01}${A}.b forcing.wndspd.b
  endif
  if (-e ./flux/${E}f${Y01}${A}) then
    /bin/ln ./flux/airtmp_${Y01}${A}.a forcing.airtmp.a
    /bin/ln ./flux/precip_${Y01}${A}.a forcing.precip.a
    /bin/ln ./flux/radflx_${Y01}${A}.a forcing.radflx.a
    /bin/ln ./flux/shwflx_${Y01}${A}.a forcing.shwflx.a
    /bin/ln ./flux/vapmix_${Y01}${A}.a forcing.vapmix.a
    /bin/ln ./flux/airtmp_${Y01}${A}.b forcing.airtmp.b
    /bin/ln ./flux/precip_${Y01}${A}.b forcing.precip.b
    /bin/ln ./flux/radflx_${Y01}${A}.b forcing.radflx.b
    /bin/ln ./flux/shwflx_${Y01}${A}.b forcing.shwflx.b
    /bin/ln ./flux/vapmix_${Y01}${A}.b forcing.vapmix.b
  endif
C
C --- If the winds or fluxes for the next segment dont exist, 
C --- interpolate them to model grid while current segment is running.
C
  if (-e ./wind/tauewd_${YXX}${B}.a && \
      -e ./wind/tauewd_${YXX}${B}.b && \
      -e ./wind/taunwd_${YXX}${B}.a && \
      -e ./wind/taunwd_${YXX}${B}.b && \
      -e ./wind/wndspd_${YXX}${B}.a && \
      -e ./wind/wndspd_${YXX}${B}.b    ) then
C
C --- next segments winds already exist.
C
  else
    cd ./wind
    /bin/rm -f ${E}w${YXX}${B}.com ${E}w${YXX}${B}.log
    awk -f $D/../${E}.awk y01=${YXX} ab=${B} $D/../${E}W.com > ${E}w${YXX}${B}.com
    csh ${E}w${YXX}${B}.com >& ${E}w${YXX}${B}.log &
    cd ..
  endif
  if (-e ./flux/airtmp_${YXX}${B}.a && \
      -e ./flux/airtmp_${YXX}${B}.b && \
      -e ./flux/precip_${YXX}${B}.a && \
      -e ./flux/precip_${YXX}${B}.b && \
      -e ./flux/radflx_${YXX}${B}.a && \
      -e ./flux/radflx_${YXX}${B}.b && \
      -e ./flux/shwflx_${YXX}${B}.a && \
      -e ./flux/shwflx_${YXX}${B}.b && \
      -e ./flux/vapmix_${YXX}${B}.a && \
      -e ./flux/vapmix_${YXX}${B}.b    ) then
C
C --- next segments fluxes already exist.
C
  else
    cd ./flux
    /bin/rm -f ${E}f${YXX}${B}.com ${E}f${YXX}${B}.log
    awk -f $D/../${E}.awk y01=${YXX} ab=${B} $D/../${E}F.com > ${E}f${YXX}${B}.com
    csh ${E}f${YXX}${B}.com >& ${E}f${YXX}${B}.log &
    cd ..
  endif
endif
C
chmod ug+x hycom
/bin/ls -laFq
C
C ---  Check to make sure restart file is there
C
if ($Y01 != "001" && -z restart_in.a) then
  cd $D/..
  /bin/mv LIST LIST_BADRUN
  echo "BADRUN" > LIST
  exit
endif
C
if ($NMPI == 0) then
C
C --- run the model, without MPI or SHMEM
C
switch ($OS)
case 'SunOS':
C
C   --- $NOMP CPUs/THREADs, if compiled for OpenMP.
C
    /bin/rm -f core
    touch core
    env OMP_NUM_THREADS=$NOMP ./hycom
    breaksw
case 'Linux':
C
C   --- $NOMP CPUs/THREADs, if compiled for OpenMP.
C
    /bin/rm -f core
    touch core
    env OMP_NUM_THREADS=$NOMP MPSTKZ=8M ./hycom
    breaksw
case 'OSF1':
C
C   --- $NOMP CPUs/THREADs, if compiled for OpenMP.
C
    /bin/rm -f core
    touch core
    env OMP_NUM_THREADS=$NOMP ./hycom
    breaksw
case 'IRIX64':
C
C   --- $NOMP CPUs/THREADs, if compiled for OpenMP.
C
    /bin/rm -f core
    touch core
    setenv FILENV .assign
    assign -R
    assign -s sbin u:18
    assign -V
    env OMP_NUM_THREADS=$NOMP ./hycom
    assign -V
    assign -R
    breaksw
case 'AIX':
C
C   --- $NOMP CPUs/THREADs, if compiled for IBM OpenMP.
C
    /bin/rm -f core
    touch core
    setenv SPINLOOPTIME     500
    setenv YIELDLOOPTIME    500
    setenv XLSMPOPTS       "parthds=${NOMP}:spins=0:yields=0"
    ./hycom
    breaksw
#case 'AIX':
#C
#C   --- $NOMP CPUs/THREADs, if compiled for KAI OpenMP.
#C
#    /bin/rm -f core
#    touch core
#    env OMP_NUM_THREADS=$NOMP ./hycom
#    breaksw
case 'unicosmk':
C
C   --- ONE CPU ONLY.
C
    /bin/rm -f core
    touch core
    ulimit
    assign -V
    ./hycom
    if (! -z core)  debugview hycom core
    ulimit
    assign -V
    assign -R
    limit -v
    breaksw
case 'unicos':
C
C   --- $NOMP CPUs/THREADs, if compiled for OpenMP.
C
    /bin/rm -f core
    touch core
    assign -V
    env OMP_NUM_THREADS=$NOMP ./hycom
    if (! -z core)  debug -s hycom core
    assign -V
    assign -R
    breaksw
endsw
else
C
C --- run the model, with MPI or SHMEM and perhaps also with OpenMP.
C
setenv NPATCH `echo $NMPI | awk '{printf("%03d", $1)}'`
/bin/cp ${D}/../../topo/partit/depth_${R}_${T}.${NPATCH} patch.input
C
switch ($OS)
case 'SunOS':
C
C   --- $NMPI MPI tasks and $NOMP THREADs, if compiled for OpenMP.
C
    /bin/rm -f core
    touch core
    setenv OMP_NUM_THREADS $NOMP
#   mpirun -np $NMPI ./hycom
    pam ./hycom
    breaksw
case 'Linux':
C
C   --- $NMPI MPI tasks and $NOMP THREADs, if compiled for OpenMP.
C
    /bin/rm -f core
    touch core
    setenv OMP_NUM_THREADS $NOMP
    mpirun -np $NMPI ./hycom
    breaksw
case 'OSF1':
C
C   --- $NMPI MPI tasks and $NOMP THREADs, if compiled for OpenMP.
C
    /bin/rm -f core
    touch core
    setenv OMP_NUM_THREADS $NOMP
    mpirun -np $NMPI ./hycom
    breaksw
case 'IRIX64':
if ($TYPE == "shmem") then
C
C   --- $NMPI SHMEM tasks
C
    /bin/rm -f core
    touch core
    setenv FILENV .assign
    assign -R
    assign -s sbin u:18
    assign -V
    setenv OMP_NUM_THREADS	1
    setenv SMA_DSM_TOPOLOGY	free
    setenv SMA_DSM_VERBOSE	1
    setenv SMA_VERSION		1
    env NPES=$NMPI ./hycom
    assign -V
    assign -R
    breaksw
else
C
C   --- $NMPI MPI tasks and $NOMP THREADs, if compiled for OpenMP.
C
    /bin/rm -f core
    touch core
    setenv FILENV .assign
    assign -R
    assign -s sbin u:18
    assign -V
    setenv OMP_NUM_THREADS	$NOMP
    setenv MPI_DSM_VERBOSE	1
    setenv MPI_REQUEST_MAX	8192
    mpirun -np $NMPI ./hycom < /dev/zero
    assign -V
    assign -R
    breaksw
endif
case 'AIX':
C
C   --- $NMPI MPI tasks and $NOMP THREADs, if compiled for IBM OpenMP.
C
    /bin/rm -f core
    touch core
    setenv SPINLOOPTIME		500
    setenv YIELDLOOPTIME	500
    setenv XLSMPOPTS		"parthds=${NOMP}:spins=0:yields=0"
    setenv MP_SHARED_MEMORY	yes
    setenv MP_SINGLE_THREAD	yes
    setenv MP_EAGER_LIMIT	65536
    setenv MP_EUILIB		us
    setenv MP_EUIDEVICE		css0
#   list where the MPI job will run
    env MP_LABELIO=YES $POE hostname
    $POE ./hycom
    breaksw
#case 'AIX':
#C
#C   --- $NMPI MPI tasks and $NOMP THREADs, if compiled for KAI OpenMP.
#C
#    /bin/rm -f core
#    touch core
#    setenv OMP_NUM_THREADS	$NOMP
#    setenv MP_SHARED_MEMORY	yes
#    setenv MP_SINGLE_THREAD	yes
#    setenv MP_EAGER_LIMIT	65536
#    setenv MP_EUILIB		us
#    setenv MP_EUIDEVICE		css0
##   list where the MPI job will run
#    env MP_LABELIO=YES $POE hostname
#    $POE ./hycom
#    breaksw
case 'unicosmk':
C
C   --- $NMPI MPI/SHMEM tasks
C
    /bin/rm -f core
    touch core
    ulimit
    assign -V
    mpprun -n $NMPI ./hycom
    if (! -z core)  debugview hycom core
    ulimit
    assign -V
    assign -R
    limit -v
    breaksw
default:
    echo "This O/S," $OS ", is not configured for MPI/SHMEM"
    exit (1)
endsw
endif
C
touch   PIPE_DEBUG
/bin/rm PIPE_DEBUG
C
C --- archive output
C
touch arch.dummy.b
foreach f (arch*)
  ${pput} $f $D/$f
  /bin/rm $f
end
C
C --- heat transport statistics output
C
if (-e flxdp_out.a) then
  ${pput} flxdp_out.a ${D}/flxdp_${Y01}${A}.a
endif
if (-e flxdp_out.b) then
  ${pput} flxdp_out.b ${D}/flxdp_${Y01}${A}.b
endif
if (-e ovrtn_out) then
  ${pput} ovrtn_out ${D}/ovrtn_${Y01}${A}
endif
C
C --- restart output
C
if (-e restart_out.a) then
  ${pput} restart_out.a ${D}/restart_${YXX}${B}.a
endif
if (-e restart_out.b) then
  ${pput} restart_out.b ${D}/restart_${YXX}${B}.b
endif
C
if (-e ./wind) then
C
C --- Delete just in time wind and flux files.
C
  touch summary_out
  if ( `tail -1 summary_out | grep -c "^normal stop"` == 1 ) then
    /bin/rm -f ./wind/*_${Y01}${A}*
    /bin/rm -f ./flux/*_${Y01}${A}*
  endif
C
  if (-e ./wind/${E}w${Y01}${A}.com) then
    /bin/mv ./wind/${E}w${Y01}${A}.{com,log} $D/..
  endif
  if (-e ./flux/${E}f${Y01}${A}.com) then
    /bin/mv ./flux/${E}f${Y01}${A}.{com,log} $D/..
  endif
C
C --- Wait for wind and flux interpolation of next segment.
C
  wait
C
  if (-e ./wind/${E}w${YXX}${B}.com) then
    /bin/mv ./wind/${E}w${YXX}${B}.{com,log} $D/..
  endif
  if (-e ./flux/${E}f${YXX}${B}.com) then
    /bin/mv ./flux/${E}f${YXX}${B}.{com,log} $D/..
  endif
endif
C
C --- HYCOM error stop is implied by the absence of a normal stop.
C
touch summary_out
if ( `tail -1 summary_out | grep -c "^normal stop"` == 0 ) then
  cd $D/..
  /bin/mv LIST LIST_BADRUN
  echo "BADRUN" > LIST
endif
C
C  --- END OF MODEL RUN SCRIPT
C
