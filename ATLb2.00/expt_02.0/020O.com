#!/bin/csh
#
#@ job_name         = 011o
#@ output           = $(job_name).log
#@ error            = $(job_name).log
#@ restart          = yes
#@ job_type         = serial
#@ node_usage       = shared
#@ wall_clock_limit = 12:00:00
#@ account_no       = NRLSS018
#@ notification     = never
#@ class            = share
#@ resources        = ConsumableCpus(1) ConsumableMemory(1 gb)
#@ queue
#
set echo
set time = 1
set timestamp
C
C --- Create model interpolated seatmp for HYCOM.
C --- MODAS version, first daily then 6-hrly or 3-hrly.
C ---
C --- prebuild this script similar to a model script
# --- awk -f 284.awk y01=098 ab=a 284O.com > 284o098a.com
C
C --- Preamble, script keys on O/S name.
C
setenv OS `uname`
switch ($OS)
case 'SunOS':
C   assumes /usr/5bin is before /bin and /usr/bin in PATH.
    breaksw
case 'Linux':
    breaksw
case 'OSF1':
    breaksw
case 'IRIX64':
    breaksw
case 'AIX':
    breaksw
case 'unicos':
    setenv ACCT `newacct -l | awk '{print $4}'`
    breaksw
default:
    echo 'Unknown Operating System: ' $OS
    exit (1)
endsw
C
C --- pget, pput "copy" files between scratch and permanent storage.
C --- Can both be cp if the permanent filesystem is mounted locally.
C
switch ($OS)
case 'SunOS':
case 'Linux':
case 'OSF1':
case 'IRIX64':
case 'AIX':
case 'unicos':
    if (-e ~wallcraf/bin/pget) then
      setenv pget ~wallcraf/bin/pget
      setenv pput ~wallcraf/bin/pput
    else
      setenv pget cp
      setenv pput cp
    endif
    breaksw
default:
    setenv pget cp
    setenv pput cp
endsw
C
C --- E is expt number.
C --- P is primary path.
C --- D is permanent directory.
C --- S is scratch   directory, must not be the permanent directory.
C --- N is data-set name, e.g. ecmwf-reanal_ds111.6
C --- W is permanent native seatmp directory
C
setenv E 020
setenv P hycom/ATLb2.00/expt_02.0/data
setenv D ~/$P
setenv N MODAS
C
switch ($OS)
case 'SunOS':
    if (-e /net/hermes/scrb) then
#                  NRLSSC
      setenv S     /net/hermes/scrb/${user}/$P
      setenv W     /net/hermes/scrb/metzger/temp_ieee/$N
    else
#                  NAVO MSRC
      mkdir        /scr/${user}
      chmod a+rx   /scr/${user}
      setenv S     /scr/${user}/$P
      setenv W     /u/home/metzger/temp_ieee/$N
    endif
    breaksw
case 'Linux':
    if (-e /external/fast) then
#                  NRLSSC
      setenv S /external/fast/${user}/$P
      setenv W ~/temp_ieee/$N
    else if (-e /work) then
#                  ERDC
      setenv S /work/${user}/$P
      setenv W ~/temp_ieee/$N
    else
#                  Single Disk
      setenv S ~/$P/SCRATCH
      setenv W ~/temp_ieee/$N
    endif
    breaksw
case 'OSF1':
#                 ERDC MSRC
    mkdir        /work/${user}
    chmod a+rx   /work/${user}
    setenv S     /work/${user}/$P
    setenv W     /u/home/metzger/temp_ieee/$N
    breaksw
case 'IRIX64':
    mkdir        /scr/${user}
    chmod a+rx   /scr/${user}
    setenv S     /scr/${user}/$P
    setenv W     /u/home/metzger/temp_ieee/$N
    breaksw
case 'AIX':
    if (-e /gpfs/work) then
#                  ERDC MSRC
      mkdir        /gpfs/work/${user}
      chmod a+rx   /gpfs/work/${user}
      setenv S     /gpfs/work/${user}/$P
      setenv W     /u/home/metzger/temp_ieee/$N
    else if (-e /scr) then
#                  NAVO MSRC
      mkdir        /scr/${user}
      chmod a+rx   /scr/${user}
      setenv S     /scr/${user}/$P
      setenv W     /u/home/metzger/temp_ieee/$N
    else
#                  ARL MSRC, under GRD
      mkdir        /usr/var/tmp/${user}
      chmod a+rx   /usr/var/tmp/${user}
      setenv S     /usr/var/tmp/${user}/$P
      setenv W     /u/home/metzger/temp_ieee/$N
    endif
    breaksw
case 'unicos':
    mkdir        /tmp/${user}
    chmod a+rx   /tmp/${user}
    mkdir        /tmp/${user}/MEda0.04
    chgrp $ACCT  /tmp/${user}/MEda0.04
    setenv S     /tmp/${user}/$P
    setenv W     /u/home/metzger/temp_ieee/$N
    breaksw
endsw
C
mkdir -p $S/ssto
cd       $S/ssto
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
C --- For seatmp, only Y01 and A are used.
C
setenv A "a"
setenv B "b"
setenv Y01 "099"
setenv YXX "099"
C
echo "Y01 =" $Y01 "YXX = " $YXX  "A =" ${A} "B =" ${B}
setenv YY `echo ${Y01}  | awk '{printf("%04d", $1+1900)}'`
C
C --- time limits.
C
if (-e ${D}/../${E}y${Y01}${A}.limits) then
  setenv TS `sed -e "s/-/ /g" ${D}/../${E}y${Y01}${A}.limits | awk '{print $1}'`
  setenv TM `cat              ${D}/../${E}y${Y01}${A}.limits | awk '{print $2}'`
else
# use "LIMITI" when starting a run after day zero.
# use "LIMITS9" (say) for a 9-day run.
  setenv TS `echo "LIMITS" | awk -f ${D}/../${E}.awk y01=${Y01} ab=${A} | awk '{print $1}'`
  setenv TM `echo "LIMITS" | awk -f ${D}/../${E}.awk y01=${Y01} ab=${A} | awk '{print $2}'`
endif
C
echo "TS =" $TS "TM =" $TM
C
C --- input files from file server.
C
touch regional.grid.a regional.grid.b
if (-z regional.grid.b) then
  ${pget} ${D}/../../topo/regional.grid.b regional.grid.b &
endif
if (-z regional.grid.a) then
  ${pget} ${D}/../../topo/regional.grid.a regional.grid.a &
endif
C
touch  fort.30a
if (-z fort.30a) then
  ${pget} $D/../../relax/PF_SST/sst_all.a fort.30a &
endif
C
@ i = 70
foreach y ( 2003 2004 2005 2006 2007 )
  @ i = $i + 1
  setenv N `echo $i | awk '{printf("%02d\n", $1)}'`
  touch  fort.${N}
  if (-z fort.${N}) then
    ${pget} ${W}/modas_${y}_Ts_daily.d fort.${N} &
  endif
end
C
C --- executable, part 1
C
###${pget} ${D}/../../../ALL/force/src/kp . &
/bin/cp /u/home/wallcraf/hycom/ALL/force/src/kp . &
wait
chmod ug+rx kp
ls -laFq
C
C --- NAMELIST input, part 1.
C
cat <<E-o-D  >! fort.05i
 &AFTITL
  CTITLE = '1234567890123456789012345678901234567890',
  CTITLE = 'MODAS SST, daily, degC',
  CNAME  = 'seatmp',
 /
 &AFTIME
  PARMIN = -999.0,  ! disable parmin
  PARMAX =  999.0,  ! disable parmax
  PAROFF =    0.0,  ! degC to degC
  FSTART = ${TS},
  TSTART = ${TS},
  TMAX   = ${TM},
 /
 &AFFLAG
  IFFILE =   5,  !3:monthly-climo; 5:actual-day;
  NGLOBE =   1,  !1:near-global
  INTERP =   0,  !0:piecewise-linear; 1:cubic-spline;
  INTMSK =   0,  !0:no mask; 1:land/sea=0/1; 2:land/sea=1/0;
 /
E-o-D
C
C --- run the seatmp interpolation.
C
touch      fort.10 fort.10a
/bin/rm -f fort.10 fort.10a
C
setenv FOR010A fort.10a
C
switch ($OS)
case 'SunOS':
case 'Linux':
case 'OSF1':
case 'AIX':
    /bin/rm -f core
    touch core
    ./kp < fort.05i
    breaksw
case 'IRIX64':
    /bin/rm -f core
    touch core
    assign -V
    ./kp < fort.05i
    assign -V
    assign -R
    breaksw
case 'unicos':
    /bin/rm -f core
    touch core
    assign -V
    ./kp < fort.05i
    if (! -z core)  debug -s kp core
    assign -V
    assign -R
    breaksw
endsw
C
C --- executable, part 2
C
###${pget} ${D}/../../../ALL/force/src/aphf_extend . &
/bin/cp /u/home/wallcraf/hycom/ALL/force/src/aphf_extend . &
wait
chmod ug+rx aphf_extend
ls -laFq
C
C --- NAMELIST input, part 2.
C
cat <<E-o-D  >! fort.05i
 &AFTIME
  FINC   = 0.125,  !3-hrly
  FINC   = 0.25,   !6-hrly
  FSTART = ${TS},
  TSTART = ${TS},
  TMAX   = ${TM},
 /
E-o-D
C
C --- run the seatmp interpolation.
C
/bin/mv fort.10  fort.20
/bin/mv fort.10a fort.20a
C
setenv FOR010A fort.10a
setenv FOR020A fort.20a
setenv FOR030A fort.30a
C
switch ($OS)
case 'SunOS':
case 'Linux':
case 'OSF1':
case 'AIX':
    /bin/rm -f core
    touch core
    ./aphf_extend < fort.05i
    breaksw
case 'IRIX64':
    /bin/rm -f core
    touch core
    assign -V
    ./aphf_extend < fort.05i
    assign -V
    assign -R
    breaksw
case 'unicos':
    /bin/rm -f core
    touch core
    assign -V
    ./aphf_extend < fort.05i
    if (! -z core)  debug -s aphf_extend core
    assign -V
    assign -R
    breaksw
endsw
C
C --- Output.
C
/bin/mv fort.10  seatmp_${Y01}${A}.b
/bin/mv fort.10a seatmp_${Y01}${A}.a
C
if (-e ./SAVE) then
  ln seatmp_${Y01}${A}.a ./SAVE/seatmp_${Y01}${A}.a
  ln seatmp_${Y01}${A}.b ./SAVE/seatmp_${Y01}${A}.b
endif
C
C --- build monthly mean SSTs
C
touch seatmp*.[ab]
cd /scr/metzger/$P
C
hycom_mean ssto/seatmp_${Y01}${A}.a 4500 3298 1 1 1 0 ssto/seatmp_${Y01}${A}_mn.a | head -1 >! ssto/seatmp_${Y01}${A}_mn.b
C
C --- combine all monthly files if last month
C
if (${A} == l) then
  cat ssto/seatmp_${Y01}?_mn.a >! ssto/seatmp_${Y01}_mon.a
  cat ssto/seatmp_${Y01}?_mn.b >! ssto/seatmp_${Y01}_mon.b
  /site/unsupported/bin/gtar --format=posix -cvf modas_sst_${YY}_mon.tar ssto/seatmp_${Y01}?_mn.? >&! modas_sst_${YY}_mon.tar.lis
  rcp modas_sst_${YY}_mon.tar     vincent:${D}/../../relax/MODAS_SST
  rcp modas_sst_${YY}_mon.tar.lis vincent:${D}/../../relax/MODAS_SST
endif
C
C  --- END OF JUST IN TIME SEATMP GENERATION SCRIPT.
C
