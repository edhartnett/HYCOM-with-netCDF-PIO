#! /bin/csh
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
#
set echo
set time = 1
set timestamp
C
C --- Experiment ATLb2.00 - 01.5
C --- 22 layer HYCOM on 2.00 degree Atlantic region
C
C --- 01.2 - KPP: COADS forcing, Levitus bndry and sur. sal. relax.
C --- 01.3 - twin of 1.2 with 10x thkdf4.
C --- 01.5 - repat 1.5 with src_2.1.03.
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
    setenv NOMP 0
    setenv NMPI 0
    breaksw
case 'AIX':
    setenv NOMP 4
    setenv NMPI 0
    breaksw
case 'sn6705':
    setenv OS unicosmk
case 'unicosmk':
    setenv ACCT `newacct -l | awk '{print $4}'`
#   setenv ACCT `groups | awk '{print $1}'`
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
  if ( $?LSB_INITIAL_NUM_PROCESSORS) then
    setenv NCPU $LSB_INITIAL_NUM_PROCESSORS
  else
    setenv NCPU `echo $LSB_MCPU_HOSTS | awk '{print $2+$4+$6+$8+$10+$12}'`
  endif
  if      ($NMPI == 0) then
    setenv NOMP $NCPU
  else if ($NOMP == 0) then
    setenv NMPI $NCPU
  else
    setenv NMPI `echo $NCPU $NOMP | awk '{print int($1/$2)}'`
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
setenv V 2.1.03
setenv T 01
setenv K 22
setenv E 015
setenv P hycom/${R}/expt_01.5/data
setenv D ~/$P
C
switch ($OS)
case 'SunOS':
    if (-e /net/hermes/scrb) then
#                  NRLSSC
      setenv S     /net/hermes/scrb/${user}/$P
    else if (-e /scr) then
#                  NAVO MSRC
      mkdir        /scr/${user}
      chmod a+rx   /scr/${user}
      setenv S     /scr/${user}/$P
    else
#                  Single Disk
      setenv S     ~/$P/SCRATCH
    endif
    breaksw
case 'Linux':
    if (-e /export/a/$user) then
#              NRLSSC
      setenv S /export/a/${user}/$P
    else
#              Single Disk
      setenv S ~/$P/SCRATCH
    endif
    breaksw
case 'OSF1':
    if      (-e /work) then
#                  ERDC MSRC
      mkdir        /work/${user}
      chmod a+rx   /work/${user}
      setenv S     /work/${user}/$P
    else if (-e /workspace) then
#                  ASC MSRC
      mkdir        /workspace/${user}
      chmod a+rx   /workspace/${user}
      setenv S     /workspace/${user}/$P
    else
#                  Single Disk
      setenv S     ~/$P/SCRATCH
    endif
    breaksw
case 'IRIX64':
    if      (-e /scr) then
#                  NAVO MSRC
      mkdir        /scr/${user}
      chmod a+rx   /scr/${user}
      setenv S     /scr/${user}/$P
    else
#                  Single Disk
      setenv S     ~/$P/SCRATCH
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
    if      (-e /work) then
#                  ERDC MSRC
      mkdir        /work/${user}
      chmod a+rx   /work/${user}
      mkdir        /work/${user}/$R
      chgrp $ACCT  /work/${user}/$R
      setenv S     /work/${user}/$P
    else
      mkdir        /tmp/${user}
      chmod a+rx   /tmp/${user}
      mkdir        /tmp/${user}/$R
      chgrp $ACCT  /tmp/${user}/$R
      setenv S     /tmp/${user}/$P
    endif
    breaksw
endsw
C
mkdir -p $S
cd       $S
C
C --- For whole year runs.
C ---   Y01 initial model year of this run.
C ---   YXX is the last model year of this run, and the first of the next run.
C ---   A and B are identical, typically blank.
C --- For part year runs.
C ---   A is this part of the year, B is next part of the year.
C ---   Y01 is the start model year of this run.
C ---   YXX is the end   model year of this run, usually Y01.
C --- For a few hour/day run
C ---   A   is the start day and hour, of form "dDDDhHH".
C ---   B   is the end   day and hour, of form "dXXXhYY".
C ---   Y01 is the start model year of this run.
C ---   YXX is the end   model year of this run, usually Y01.
C --- Note that these variables are set by the .awk generating script.
C
setenv A ""
setenv B ""
setenv Y01 "003"
setenv YXX "004"
C
echo "Y01 =" $Y01 "YXX = " $YXX  "A =" ${A} "B =" ${B}
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
if ($NMPI != 0) then
  setenv NPATCH `echo $NMPI | awk '{printf("%03d", $1)}'`
  /bin/rm -f patch.input
  /bin/cp ${D}/../../topo/partit/depth_${R}_${T}.${NPATCH} patch.input
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
#C
#C --- turn on detailed debugging.
#C
#touch PIPE_DEBUG
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
    if      (-e ~wallcraf/bin/pget_navo) then
      setenv pget ~wallcraf/bin/pget_navo
      setenv pput ~wallcraf/bin/pput_navo
    else if (-e ~wallcraf/bin/pget) then
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
  setenv FN coads
  touch forcing.tauewd.a forcing.taunwd.a forcing.wndspd.a forcing.airtmp.a forcing.precip.a forcing.radflx.a forcing.shwflx.a forcing.vapmix.a
  touch forcing.tauewd.b forcing.taunwd.b forcing.wndspd.b forcing.airtmp.b forcing.precip.b forcing.radflx.b forcing.shwflx.b forcing.vapmix.b
  if (-z forcing.tauewd.a) then
     ${pget} ${D}/../../force/${FN}/tauewd.a      forcing.tauewd.a &
  endif
  if (-z forcing.tauewd.b) then
     ${pget} ${D}/../../force/${FN}/tauewd.b      forcing.tauewd.b &
  endif
  if (-z forcing.taunwd.a) then
     ${pget} ${D}/../../force/${FN}/taunwd.a      forcing.taunwd.a &
  endif
  if (-z forcing.taunwd.b) then
     ${pget} ${D}/../../force/${FN}/taunwd.b      forcing.taunwd.b &
  endif
  if (-z forcing.wndspd.a) then
     ${pget} ${D}/../../force/${FN}/wndspd.a      forcing.wndspd.a &
  endif
  if (-z forcing.wndspd.b) then
     ${pget} ${D}/../../force/${FN}/wndspd.b      forcing.wndspd.b &
  endif
  if (-z forcing.airtmp.a) then
     ${pget} ${D}/../../force/${FN}/airtmp.a      forcing.airtmp.a &
  endif
  if (-z forcing.airtmp.b) then
     ${pget} ${D}/../../force/${FN}/airtmp.b      forcing.airtmp.b &
  endif
  if (-z forcing.vapmix.a) then
     ${pget} ${D}/../../force/${FN}/vapmix.a      forcing.vapmix.a &
  endif
  if (-z forcing.vapmix.b) then
     ${pget} ${D}/../../force/${FN}/vapmix.b      forcing.vapmix.b &
  endif
# setenv PO ""
  setenv PO "_zero"
  if (-z forcing.precip.a) then
     ${pget} ${D}/../../force/${FN}/precip${PO}.a forcing.precip.a &
  endif
  if (-z forcing.precip.b) then
     ${pget} ${D}/../../force/${FN}/precip${PO}.b forcing.precip.b &
  endif
  setenv FO ""
# setenv FO "-25w"
  if (-z forcing.radflx.a) then
     ${pget} ${D}/../../force/${FN}/radflx${FO}.a forcing.radflx.a &
  endif
  if (-z forcing.radflx.b) then
     ${pget} ${D}/../../force/${FN}/radflx${FO}.b forcing.radflx.b &
  endif
  if (-z forcing.shwflx.a) then
     ${pget} ${D}/../../force/${FN}/shwflx${FO}.a forcing.shwflx.a &
  endif
  if (-z forcing.shwflx.b) then
     ${pget} ${D}/../../force/${FN}/shwflx${FO}.b forcing.shwflx.b &
  endif
endif
C
touch forcing.rivers.a
touch forcing.rivers.b
if (-z forcing.rivers.a) then
   ${pget} ${D}/../../force/rivers/rivers_${T}.a forcing.rivers.a &
endif
if (-z forcing.rivers.b) then
   ${pget} ${D}/../../force/rivers/rivers_${T}.b forcing.rivers.b &
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
touch   restart_in.a restart_in.b restart_out.a restart_out.b
if (-z restart_in.b) then
  setenv RI "       0.00"
else
  setenv RI `head -2 restart_in.b | tail -1 | awk  '{printf("%11.2f\n", $5)}'`
endif
if (-z restart_out.b) then
  setenv RO "       0.00"
else
  setenv RO `head -2 restart_out.b | tail -1 | awk  '{printf("%11.2f\n", $5)}'`
endif
setenv LI `awk  '{printf("%11.2f\n", $1)}' limits`
C
if (`echo $LI | awk '{print ($1 <= 0.0)}'`) then
C --- no restart needed
  /bin/rm restart_in.a restart_in.b restart_out.a restart_out.b
else
  if (`echo $LI $RI | awk '{print ($1-0.1 < $2 && $1+0.1 > $2)}'`) then
C --- restart is already in restart_in
    /bin/rm restart_out.a restart_out.b
  else
    if (`echo $LI $RO | awk '{print ($1-0.1 < $2 && $1+0.1 > $2)}'`) then
C --- restart is in restart_out
      /bin/mv restart_out.a restart_in.a
      /bin/mv restart_out.b restart_in.b
    else
C --- get restart from permenant storage
      /bin/rm restart_out.a restart_out.b
      ${pget} ${D}/restart_${Y01}${A}.a restart_in.a &
      ${pget} ${D}/restart_${Y01}${A}.b restart_in.b &
    endif
  endif
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
C --- clean up old archive files, typically from batch system rerun.
C
touch archv.dummy.b
foreach f (arch*)
  /bin/rm $f
end
C
C --- Nesting input archive files.
C
if (-e ./nest) then
  cd ./nest
  touch rmu.a rmu.b
  if (-z rmu.a) then
     ${pget} ${D}/../../relax/${E}/nest_rmu.a rmu.a &
  endif
  if (-z rmu.b) then
     ${pget} ${D}/../../relax/${E}/nest_rmu.b rmu.b &
  endif
  touch   arch.dummy.b
  /bin/rm arch*.[ab]
  touch archv_${Y01}${A}.tar
  if (-z archv_${Y01}${A}.tar) then
    ${pget} ${D}/nest/archv_${Y01}${A}.tar archv_${Y01}${A}.tar
  endif
  tar xvf archv_${Y01}${A}.tar
  cd ..
endif
C
C --- let all file copies complete.
C
wait
C
C --- zero file length means no rivers.
C
if (-z forcing.rivers.a) then
   /bin/rm forcing.rivers.[ab]
endif
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
  if (-e ./pcip) then
    /bin/rm -f forcing.precip.a
    /bin/rm -f forcing.precip.b
    /bin/rm -f ./pcip/${E}p${Y01}${A}
    if (-e ./pcip/precip_${Y01}${A}.a && \
        -e ./pcip/precip_${Y01}${A}.b    ) then
      /bin/ln ./pcip/precip_${Y01}${A}.a forcing.precip.a
      /bin/ln ./pcip/precip_${Y01}${A}.b forcing.precip.b
    else
      cd ./pcip
      touch ${E}p${Y01}${A}
      /bin/rm -f ${E}p${Y01}${A}.com ${E}p${Y01}${A}.log
      awk -f $D/../${E}.awk y01=${Y01} ab=${A} $D/../${E}P.com > ${E}p${Y01}${A}.com
      csh ${E}p${Y01}${A}.com >& ${E}p${Y01}${A}.log &
      cd ..
    endif
    /bin/rm -f forcing.airtmp.a forcing.radflx.a forcing.shwflx.a forcing.vapmix.a
    /bin/rm -f forcing.airtmp.b forcing.radflx.b forcing.shwflx.b forcing.vapmix.b
    /bin/rm -f ./flux/${E}f${Y01}${A}
    if (-e ./flux/airtmp_${Y01}${A}.a && \
        -e ./flux/airtmp_${Y01}${A}.b && \
        -e ./flux/radflx_${Y01}${A}.a && \
        -e ./flux/radflx_${Y01}${A}.b && \
        -e ./flux/shwflx_${Y01}${A}.a && \
        -e ./flux/shwflx_${Y01}${A}.b && \
        -e ./flux/vapmix_${Y01}${A}.a && \
        -e ./flux/vapmix_${Y01}${A}.b    ) then
      /bin/ln ./flux/airtmp_${Y01}${A}.a forcing.airtmp.a
      /bin/ln ./flux/radflx_${Y01}${A}.a forcing.radflx.a
      /bin/ln ./flux/shwflx_${Y01}${A}.a forcing.shwflx.a
      /bin/ln ./flux/vapmix_${Y01}${A}.a forcing.vapmix.a
      /bin/ln ./flux/airtmp_${Y01}${A}.b forcing.airtmp.b
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
  else
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
  if (-e ./pcip) then
    if (-e ./pcip/${E}p${Y01}${A}) then
      /bin/ln ./pcip/precip_${Y01}${A}.a forcing.precip.a
      /bin/ln ./pcip/precip_${Y01}${A}.b forcing.precip.b
    endif
    if (-e ./flux/${E}f${Y01}${A}) then
      /bin/ln ./flux/airtmp_${Y01}${A}.a forcing.airtmp.a
      /bin/ln ./flux/radflx_${Y01}${A}.a forcing.radflx.a
      /bin/ln ./flux/shwflx_${Y01}${A}.a forcing.shwflx.a
      /bin/ln ./flux/vapmix_${Y01}${A}.a forcing.vapmix.a
      /bin/ln ./flux/airtmp_${Y01}${A}.b forcing.airtmp.b
      /bin/ln ./flux/radflx_${Y01}${A}.b forcing.radflx.b
      /bin/ln ./flux/shwflx_${Y01}${A}.b forcing.shwflx.b
      /bin/ln ./flux/vapmix_${Y01}${A}.b forcing.vapmix.b
    endif
  else
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
  if (-e ./pcip) then
    if (-e ./pcip/precip_${YXX}${B}.a && \
        -e ./pcip/precip_${YXX}${B}.b    ) then
C
C ---   next segments pcip already exist.
C
    else
      cd ./pcip
      /bin/rm -f ${E}p${YXX}${B}.com ${E}p${YXX}${B}.log
      awk -f $D/../${E}.awk y01=${YXX} ab=${B} $D/../${E}P.com > ${E}p${YXX}${B}.com
      csh ${E}p${YXX}${B}.com >& ${E}p${YXX}${B}.log &
      cd ..
    endif
    if (-e ./flux/airtmp_${YXX}${B}.a && \
        -e ./flux/airtmp_${YXX}${B}.b && \
        -e ./flux/radflx_${YXX}${B}.a && \
        -e ./flux/radflx_${YXX}${B}.b && \
        -e ./flux/shwflx_${YXX}${B}.a && \
        -e ./flux/shwflx_${YXX}${B}.b && \
        -e ./flux/vapmix_${YXX}${B}.a && \
        -e ./flux/vapmix_${YXX}${B}.b    ) then
C
C ---   next segments fluxes already exist.
C
    else
      cd ./flux
      /bin/rm -f ${E}f${YXX}${B}.com ${E}f${YXX}${B}.log
      awk -f $D/../${E}.awk y01=${YXX} ab=${B} $D/../${E}F.com > ${E}f${YXX}${B}.com
      csh ${E}f${YXX}${B}.com >& ${E}f${YXX}${B}.log &
      cd ..
    endif
  else
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
C ---   next segments fluxes already exist.
C
    else
      cd ./flux
      /bin/rm -f ${E}f${YXX}${B}.com ${E}f${YXX}${B}.log
      awk -f $D/../${E}.awk y01=${YXX} ab=${B} $D/../${E}F.com > ${E}f${YXX}${B}.com
      csh ${E}f${YXX}${B}.com >& ${E}f${YXX}${B}.log &
      cd ..
    endif
  endif
endif
C
C --- Nesting input archive files for next segment.
C
if (-e ./nest) then
  cd ./nest
  touch archv_${YXX}${B}.tar
  if (-z archv_${YXX}${B}.tar) then
    ${pget} ${D}/nest/archv_${YXX}${B}.tar archv_${YXX}${B}.tar &
  endif
  cd ..
endif
C
chmod ug+x hycom
/bin/ls -laFq
C
if (-e ./nest) then
  ls -laFq nest
endif
C
C ---  Check to make sure restart file is there
C
if (`echo $LI | awk '{print ($1 > 0.0)}'` && -z restart_in.a) then
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
if ($NOMP == 0) then
  setenv NOMP 1
endif
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
touch patch.input
if (-z patch.input) then
C
C --- patch.input is always required for MPI or SHMEM.
C
  cd $D/..
  /bin/mv LIST LIST_BADRUN
  echo "BADRUN" > LIST
  exit
endif
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
#   mpirun -np $NMPI ./hycom
    time prun -n $NMPI ./hycom
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
touch archv.dummy.b
foreach f (arch*)
  ${pput} $f $D/$f
  /bin/rm $f
end
#foreach f (archv.*)
#  /bin/mv ${f} ${E}_${f}
#end
#tar cvf ${E}_archv_${Y01}${A}.tar ${E}_archv.*
#/bin/rm ${E}_archv.*
#${pput} ${E}_archv_${Y01}${A}.tar $D/${E}_archv_${Y01}${A}.tar
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
  if (-e ./pcip/${E}p${YXX}${B}.com) then
    /bin/mv ./pcip/${E}p${YXX}${B}.{com,log} $D/..
  endif
endif
C
C --- wait for nesting .tar file.
C
if (-e ./nest) then
  wait
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
